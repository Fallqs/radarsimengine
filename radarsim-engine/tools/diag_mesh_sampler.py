"""Search the SBR ray-sampler parameter space against the plate mesh golden.

Model (validated structure):
  bb[ch,p,s] += stuff * K * (k/2) * exp(j*2*pi*beat) / (R_t*R_r) * dA
with K the point-chain constant, beat the per-sample deramp phase, dA the
per-ray area. The sampler: 1-degree occupancy grid, per-cell N x N sub-rays.
Unknowns searched: N rounding rule, sub-ray layout, occupancy rule, weights.
"""
import numpy as np

C = 299792458.0
f0, f1, T = 24.075e9, 24.175e9, 80e-6
KS = (f1 - f0) / T
FS = 6e4
LAM = C / f0  # or fc? test both

GOLD = np.array([
    -0.035401657154601196 - 0.03127654608528041j,
    0.04677351677154601 + 0.000666519439074309j,
    -0.035729766119737105 + 0.02989363903760489j,
    0.007594792597190564 - 0.04604849701138625j,
])

# plate: x=10, y/z in [-2.5, 2.5]; tx=rx=(0,0,0)
def plate_hits(n_sub, grid_deg=1.0, density=0.4, layout="center"):
    """Ray hit points on the plate + per-ray area.
    Occupancy: 1-deg az/el cells the plate touches; N=n_sub sub-rays per cell
    per axis; weight = the ray's area on the plate (projected)."""
    g = np.radians(grid_deg)
    hits = []
    # plate angular extent from origin: az/el in [-atan(2.5/10), +..]
    lim = np.arctan(2.5 / 10.0)
    n_cells = int(np.ceil(2 * lim / g))
    i0 = int(np.floor(-lim / g))
    for iaz in range(i0, i0 + n_cells + 1):
        for iel in range(i0, i0 + n_cells + 1):
            for a in range(n_sub):
                for b in range(n_sub):
                    if layout == "center":
                        az = (iaz + (a + 0.5) / n_sub) * g
                        el = (iel + (b + 0.5) / n_sub) * g
                    else:
                        az = (iaz + a / n_sub) * g
                        el = (iel + b / n_sub) * g
                    d = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
                    # hit plane x=10
                    if d[0] <= 0:
                        continue
                    tt = 10.0 / d[0]
                    p = d * tt
                    if abs(p[1]) <= 2.5 and abs(p[2]) <= 2.5:
                        # area: angular subcell solid angle / cos^3 (projected)
                        dOm = (g / n_sub) ** 2
                        cosa = d[0]  # angle to plate normal
                        dA = dOm * tt * tt / cosa / cosa / cosa * 0  # placeholder
                        hits.append((p, tt, cosa, dOm))
    return hits


def evaluate(hits, weight_mode, density=0.4):
    tot = np.zeros(4, dtype=complex)
    for (p, tt, cosa, dOm) in hits:
        R = tt
        tau = 2 * R / C
        lam = C / 24.125e9  # fc
        # per-ray area
        if weight_mode == "solid":
            dA = dOm * R * R / cosa / cosa / cosa  # area on plate
        else:
            dA = dOm * R * R / cosa
        for s in range(4):
            u = s / FS
            beat = (24.075e9) * tau + KS * u * tau - 0.5 * KS * tau ** 2
            tot[s] += np.exp(2j * np.pi * beat) / (R * R) * dA
    return tot


for n_sub in (5, 6):
    for layout in ("center", "corner"):
        hits = plate_hits(n_sub, layout=layout)
        for wm in ("solid", "proj"):
            tot = evaluate(hits, wm)
            print(f"n={n_sub} {layout:6s} {wm:5s}: rays={len(hits):6d} "
                  f"|E0|={abs(tot[0]):.4e} ph={np.angle(tot[0])/(2*np.pi)%1:.4f}")
print("gold |s1|:", abs(GOLD[0]), "phase:", np.angle(GOLD[0])/(2*np.pi)%1)
EOF_MARKER = None
