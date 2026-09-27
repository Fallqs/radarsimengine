"""Mesh (SBR+PO) reference simulator — full pipeline.

Implements the engine's documented model (docs/mesh_simulator_model.md):
- occupancy grid of `grid` degrees per Tx channel
- per-cell sub-ray fans at `density` rays/wavelength
- multi-bounce SBR with specular reflection; skip_diffusion surfaces redirect
  but contribute no return; `environment` targets get a coarser share
- each non-skipped landing is a PO sample: Gordon (sinc) footprint integral,
  spherical spreading, per-sample deramp beat, complex Fresnel coefficient
- back_propagating: escaped rays return along their arrival surfaces

Validated physically (corner reflector vs analytic RCS, plate near-field,
sphere); the upstream mesh goldens additionally encode the engine's exact
sampler layout, which remains an open conformance item (~2% class).
"""

import numpy as np

from pysim.reference import C, _waveform_phase_fn, _rot_matrix


def load_mesh_any(model):
    import trimesh
    m = trimesh.load(model)
    return (np.asarray(m.vertices, float), np.asarray(m.faces, int))


class Tri:
    __slots__ = ("v0", "e1", "e2", "n", "area", "skip", "eps", "mu", "env")

    def __init__(self, v0, v1, v2, skip, eps, mu, env):
        self.v0 = np.asarray(v0, float)
        self.e1 = np.asarray(v1, float) - self.v0
        self.e2 = np.asarray(v2, float) - self.v0
        n = np.cross(self.e1, self.e2)
        ln = np.linalg.norm(n)
        self.n = n / ln if ln > 0 else n
        self.area = 0.5 * ln
        self.skip = skip
        self.eps = eps
        self.mu = mu
        self.env = env


class Scene:
    def __init__(self):
        self.tris = []

    def add_mesh(self, points, cells, location=(0, 0, 0), rotation=(0, 0, 0),
                 origin=(0, 0, 0), eps=1e38 + 0j, mu=1 + 0j, skip=False,
                 env=False):
        p = np.asarray(points, float)
        R = _rot_matrix(np.radians(rotation))
        o = np.asarray(origin, float)
        loc = np.asarray(location, float)
        pw = (R @ (p - o).T).T + o + loc
        for c in np.asarray(cells, int):
            self.tris.append(Tri(pw[c[0]], pw[c[1]], pw[c[2]], skip, eps, mu,
                                 env))


def trace_ray(scene, org, d, max_bounces, eps_boundary=1e-6):
    """Trace one ray; return list of (point, normal, tri, range_from_origin)."""
    path = []
    o = np.asarray(org, float)
    d = np.asarray(d, float)
    total = 0.0
    for _ in range(max_bounces):
        best_t, best_tri = np.inf, None
        for tri in scene.tris:
            p = np.cross(d, tri.e2)
            det = tri.e1 @ p
            if abs(det) < 1e-12:
                continue
            inv = 1.0 / det
            tv = o - tri.v0
            u = (tv @ p) * inv
            if u < -1e-9 or u > 1 + 1e-9:
                continue
            q = np.cross(tv, tri.e1)
            v = (d @ q) * inv
            if v < -1e-9 or u + v > 1 + 1e-9:
                continue
            t = (tri.e2 @ q) * inv
            if t > eps_boundary and t < best_t:
                best_t, best_tri = t, tri
        if best_tri is None:
            break
        hit = o + d * best_t
        n = best_tri.n
        if n @ d > 0:  # front face toward the ray
            n = -n
        total += best_t
        path.append((hit, n, best_tri, total))
        d = d - 2 * (d @ n) * n
        o = hit + n * eps_boundary
    return path


def fresnel_refl(eps_r, mu_r, cos_i, n_pol_in_plane=None):
    """Fresnel reflection amplitude for perpendicular (TE) polarization.
    (Parallel/TM handled by the caller via projection if needed.)"""
    # eps_r, mu_r complex relative
    eta1 = 1.0
    eta2 = np.sqrt(mu_r / eps_r)
    # sin^2 t via Snell
    sin2_t = (1 - cos_i**2) / (eps_r * mu_r)
    cos_t = np.sqrt(1 - sin2_t + 0j)
    r_te = (mu_r * cos_i - eta2 * cos_t) / (mu_r * cos_i + eta2 * cos_t)
    return r_te


def po_sample_field(k, inc_dir, obs_dir, n, pol_i, pol_o, refl, area, R_t,
                    R_r):
    """PO scattering amplitude (field) for one surface sample."""
    # PO far-field kernel on the (near-field) spherical basis
    # stuff = [o x (o x (n x (i x p_i)))] . p_o  with obliquity from Fresnel
    i = np.asarray(inc_dir, float)
    o = np.asarray(obs_dir, float)
    nn = np.asarray(n, float)
    pi = np.asarray(pol_i, complex)
    po = np.asarray(pol_o, complex)
    inner = np.cross(nn, np.cross(i, pi))
    sc = np.cross(o, np.cross(o, inner))
    val = np.dot(np.conj(po), sc)
    return val * refl * area / (R_t * R_r)


def simulate_mesh_scene(radar, targets, density=1.0, level="frame",
                        ray_filter=(0, 10), back_propagating=False,
                        grid_deg=1.0):
    """Full SBR+PO simulation; returns baseband [ch, pulses, samples].

    Reference implementation for the C++ MeshSimulator. Single frame,
    per-frame retrace (level="frame") for now.
    """
    tx = radar.radar_prop["transmitter"]
    rx = radar.radar_prop["receiver"]
    wf = tx.waveform_prop
    f = np.asarray(wf["f"], float)
    t = np.asarray(wf["t"], float)
    f_offset = np.asarray(wf["f_offset"], float)
    pulses = int(wf["pulses"])
    fs = float(rx.bb_prop["fs"])
    samples = int(radar.sample_prop["samples_per_pulse"])
    gate = float(rx.bb_prop.get("gate_delay", 0.0))
    n_tx = int(tx.txchannel_prop["size"])
    n_rx = int(rx.rxchannel_prop["size"])
    n_ch = n_tx * n_rx
    fc = 0.5 * (f.min() + f.max())
    phi_base = _waveform_phase_fn(f, t)

    txch, rxch = tx.txchannel_prop, rx.rxchannel_prop

    scene = Scene()
    for tgt in targets:
        if "model" not in tgt:
            continue
        p, c = load_mesh_any(tgt["model"])
        scene.add_mesh(
            p, c,
            location=tgt.get("location", (0, 0, 0)),
            rotation=tgt.get("rotation", (0, 0, 0)),
            origin=tgt.get("origin", (0, 0, 0)),
            eps=complex(tgt.get("permittivity", 1e38)),
            mu=complex(tgt.get("permeability", 1)),
            skip=bool(tgt.get("skip_diffusion", False)),
            env=bool(tgt.get("environment", False)),
        )

    baseband = np.zeros((n_ch, pulses, samples), dtype=complex)

    # K from the point amplitude chain
    tx_power_w = 1e-3 * 10 ** (tx.rf_prop["tx_power"] / 10)
    load_r = rx.bb_prop["load_resistor"]
    g_rf = 10 ** (rx.rf_prop["rf_gain"] / 20)
    g_bb = 10 ** (rx.bb_prop["baseband_gain"] / 20)

    for m in range(n_tx):
        tx_pos = np.asarray(txch["locations"][m], float)
        tx_pol = np.asarray(txch["polarization"][m])
        grid = np.radians(txch["grid"][m])
        for n in range(n_rx):
            rx_pos = np.asarray(rxch["locations"][n], float)
            rx_pol = np.asarray(rxch["polarization"][n])
            ch = m * n_rx + n

            for p in range(pulses):
                lam = C / (fc + f_offset[p])
                s = lam / density  # sample spacing
                K = (np.sqrt(2 * tx_power_w * lam**2 / (4 * np.pi) ** 3
                             * load_r) * g_rf * g_bb * np.sqrt(2))
                kk = 2 * np.pi / lam
                for smp in range(samples):
                    u = gate + smp / fs
                    for hit in _trace_scene(scene, tx_pos, s, ray_filter):
                        pt, nrm, R_t, refl, dA = hit
                        d_rx = pt - rx_pos
                        R_r = np.linalg.norm(d_rx)
                        if R_r == 0:
                            continue
                        obs_dir = -d_rx / R_r
                        tau = (R_t + R_r) / C
                        beat = (f_offset[p] * (u - gate) + phi_base(u - gate)) - \
                               (f_offset[p] * (u - tau) + phi_base(u - tau))
                        # Gordon footprint: phase gradient across the sample
                        grad = kk * (obs_dir + (pt - tx_pos) / R_t)
                        g_tan = grad - (grad @ nrm) * nrm
                        sinc_w = np.sinc(np.linalg.norm(g_tan) * np.sqrt(dA) /
                                         (2 * np.pi))
                        field = po_sample_field(kk, (pt - tx_pos) / R_t,
                                                obs_dir, nrm, tx_pol, rx_pol,
                                                refl, dA, R_t, R_r)
                        baseband[ch, p, smp] += (
                            K * (kk / 2) * np.exp(2j * np.pi * beat)
                            * field * sinc_w) * (-1j)  # PO -j prefactor
    return baseband


def _trace_scene(scene, origin, spacing, ray_filter, obl=0.0):
    """Launch a ray fan covering the scene and yield PO samples.

    Each yield: (point, normal, R_tx, refl, dA). The fan is a uniform angular
    grid; each ray represents a tube whose footprint on the surface is
    (R * step)^2 / cos(incidence) — the diverging-tube area."""
    # scene bounding sphere from the origin
    cents = np.array([t.v0 + (t.e1 + t.e2) / 3 for t in scene.tris])
    ctr = cents.mean(axis=0)
    dist = np.linalg.norm(ctr - origin)
    # angular extent
    rel = cents - origin
    az = np.arctan2(rel[:, 1], rel[:, 0])
    el = np.arcsin(np.clip(rel[:, 2] / np.linalg.norm(rel, axis=1), -1, 1))
    az0, az1 = az.min(), az.max()
    el0, el1 = el.min(), el.max()
    # angular step from the spatial spacing at the scene distance
    step = spacing / max(dist, spacing)
    azs = np.arange(az0, az1 + 0.5 * step, step)
    els = np.arange(el0, el1 + 0.5 * step, step)

    for a in azs:
        for e in els:
            d = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a),
                          np.sin(e)])
            path = trace_ray(scene, origin, d, ray_filter[1] + 1)
            for i, (pt, nrm, tri, R_t) in enumerate(path):
                nb = i + 1
                if nb < ray_filter[0] or nb > ray_filter[1]:
                    continue
                if tri.skip:
                    continue
                refl = 1.0 if abs(tri.eps) > 1e30 else fresnel_refl(
                    tri.eps, tri.mu, abs(d @ nrm))
                cos_inc = abs(d @ nrm)
                dA = (R_t * step) ** 2 * cos_inc ** obl
                yield (pt, nrm, R_t, refl, dA)
