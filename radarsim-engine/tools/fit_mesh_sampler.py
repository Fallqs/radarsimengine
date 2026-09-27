"""Fit the mesh sampler against the plate golden (test_scene_single_target)."""
import itertools
import numpy as np

C = 299792458.0
KS = 1.25e12
FS = 6e4
F0 = 24.075e9
FC = 24.125e9
GOLD = np.array([-0.035401657154601196 - 0.03127654608528041j,
                  0.04677351677154601 + 0.000666519439074309j,
                 -0.035729766119737105 + 0.02989363903760489j,
                  0.007594792597190564 - 0.04604849701138625j])
PEAK = np.abs(GOLD).max()


def run(sampler, anchor, lam_src, obliquity, gordon):
    lam = C / (F0 if lam_src == "f0" else FC)
    k = 2 * np.pi / lam
    s = lam / 0.4
    tot = np.zeros(4, dtype=complex)
    if sampler == "tube":
        # Cartesian tube grid over the plate
        n = int(np.ceil(5.0 / s)) + 1
        ys = -2.5 + (np.arange(n) + anchor) * s
        zs = -2.5 + (np.arange(n) + anchor) * s
        Y, Z = np.meshgrid(ys, zs, indexing="ij")
        mask = (np.abs(Y) <= 2.5) & (np.abs(Z) <= 2.5)
        Y, Z = Y[mask], Z[mask]
        R = np.sqrt(100 + Y**2 + Z**2)
        w = np.full_like(R, s * s) * (10.0 / R) ** obliquity
        pts = np.stack([Y, Z], -1)
    else:  # fan: angular grid from origin
        h = np.radians(1.0) / int(round(np.radians(1.0) * 10 / s))  # cells per deg rule
        lim = np.arctan(2.5 / 10)
        n = int(np.ceil(lim / h)) + 1
        gs = (np.arange(-n, n + 1) + anchor - 0.5) * h
        A, E = np.meshgrid(gs, gs, indexing="ij")
        ca, ce = np.cos(A), np.cos(E)
        tt = 10.0 / (ce * ca)
        Y, Z = ce * np.sin(A) * tt, np.sin(E) * tt
        mask = (np.abs(Y) <= 2.5) & (np.abs(Z) <= 2.5)
        Y, Z, tt = Y[mask], Z[mask], tt[mask]
        R = tt
        cosa = ce[mask] * ca[mask]  # dot(d, normal)
        dOm = h * h
        w = dOm * R * R / cosa ** obliquity
        pts = np.stack([Y, Z], -1)
    tau = 2 * R / C
    if gordon:
        gy = 2 * k * pts[:, 0] / R
        gz = 2 * k * pts[:, 1] / R
        w = w * np.sinc(gy * s / (2 * np.pi)) * np.sinc(gz * s / (2 * np.pi))
    for smp in range(4):
        u = smp / FS
        beat = F0 * tau + KS * u * tau - 0.5 * KS * tau * tau
        tot[smp] = np.sum(w * np.exp(2j * np.pi * beat) / (R * R))
    return tot


best = None
for sampler, anchor, lam_src, obl, gordon in itertools.product(
        ("tube", "fan"), (0.0, 0.25, 0.5), ("f0", "fc"), (0, 1, 3), (True, False)):
    tot = run(sampler, anchor, lam_src, obl, gordon)
    # free complex scale: isolate sampler structure from normalization
    alpha = np.sum(GOLD * np.conj(tot)) / np.sum(np.abs(tot) ** 2)
    pred = tot * alpha
    err = np.abs(pred - GOLD).max() / PEAK
    tag = f"{sampler:4s} anch={anchor:.2f} lam={lam_src} obl={obl} gordon={int(gordon)}"
    if best is None or err < best[0]:
        best = (err, tag)
    if err < 0.02:
        print(f"{tag}: max err {err*100:.2f}% of peak (free scale)")
print("BEST:", best)
