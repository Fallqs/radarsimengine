#!/usr/bin/env python
"""Demo 4 — UHF EWR on hilly coastal upland, with terrain clutter (30 s).

Scenario
--------
- Same UHF EWR as demo 3 (fc 435 MHz, 200 kHz chirp, PRF 100 Hz,
  16-pulse CPI, 1 MW tx, 35+35 dB gains), now sited on a ~430 m
  coastal hilltop overlooking the sea.  1 frame/s over 30 s.
- Procedural terrain: analytic height field (Gaussian hills + coastal
  ridge), sea level at 0 m east of a wiggly coastline.
- Terrain-sourced interference: the surface is tiled into 2 km x 2 km
  clutter patches inside the +/-8 deg surveillance wedge.  Each patch
  gets an RCS from a reflectivity (sigma0) model:
      land:  sigma0 = gamma * sin(psi),        gamma = -15 dB
      sea:   sigma0_dB = -45 + 15*log10(psi_deg)   (sea state ~4)
  times patch area, times a lognormal speckle factor (sigma 2 dB), with
  a random phase; sea patches drift 0.3-1 m/s, land patches <=0.15 m/s
  (wind).  Patches shadowed by intervening terrain (earth-curvature LOS
  check, 4/3 earth radius) are excluded.
- Antenna elevation pattern peaks at +1..+2 deg and rolls off below the
  horizon, so low-elevation clutter is partially rejected while the
  raid axis keeps near-full gain.
- Raid: B-52 high altitude (always visible), B-1B at 200 m and B-2 at
  150 m pop over the radar horizon mid-scenario (terrain masking);
  two airliners as background traffic.  Receiver noise is on.

Outputs (saved to out/demo4_ewr_terrain.png):
- Terrain map with coastline, visible/shadowed clutter patches, tracks
- Range-Doppler map at mid-scenario (clutter ridge + aliased targets)
- Zero-Doppler clutter profile vs range (CNR)
- SNR per contact vs time with horizon pop-up markers

Run from the repo root:
    python radarsim-engine/demos/demo4_ewr_terrain.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from radarsimpy import Radar, Transmitter, Receiver
from radarsimpy.simulator import sim_radar
import radarsimpy.processing as proc

C = 299792458.0
RE = 4.0 / 3.0 * 6371e3          # effective earth radius (4/3 model)
OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(7)

# --- radar ------------------------------------------------------------------
FC = 435e6
BW = 200e3               # 750 m range resolution
PULSE = 3e-3
PRP = 10e-3              # PRF 100 Hz
PULSES = 16
FS = 500e3               # 1500 samples
FRAMES = 30              # 1 frame/s, 30 s

LAM = C / (FC + BW / 2)
SLOPE = BW / PULSE

# elevation pattern: beam peaked +1..+2 deg, rolled off below horizon
EL_ANG = [-90, -3, -1.5, -0.5, 1.0, 2.0, 3.0, 90]
EL_PAT = [-40, -40, -12, -3, 0, 0, -3, -40]
# azimuth pattern: 10-deg main lobe at 35 dB, ~10 dB sidelobes elsewhere
# (ring patches therefore model sidelobe clutter, wedge patches main-beam)
AZ_ANG = [-90, -7, -5, 5, 7, 90]
AZ_PAT = [10, 10, 35, 35, 10, 10]

# --- terrain ----------------------------------------------------------------
def terrain_h(x, y):
    """Analytic height field: hills inland (x<0), sea (h=0) past the coast."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    hills = (
        180 * np.exp(-((x + 3e3) ** 2 + y ** 2) / (2 * 2.0e3 ** 2)) +
        250 * np.exp(-((x + 35e3) ** 2 + (y - 20e3) ** 2) / (2 * 14e3 ** 2)) +
        460 * np.exp(-((x + 50e3) ** 2 + (y + 25e3) ** 2) / (2 * 17e3 ** 2)) +
        330 * np.exp(-((x + 15e3) ** 2 + (y - 30e3) ** 2) / (2 * 10e3 ** 2)) +
        200 * np.exp(-((x + 22e3) ** 2 + (y + 8e3) ** 2) / (2 * 9e3 ** 2))
    )
    base = 40.0 + hills
    xc = 2e3 * np.sin(y / 10e3)          # coastline near x=0
    # 1 inland, 0 at sea (clip keeps exp finite far out to sea)
    w = 1.0 / (1.0 + np.exp(np.clip((x - xc) / 600.0, -60, 60)))
    return base * w


RADAR_XY = (-3e3, 0.0)
RADAR_Z = float(terrain_h(*RADAR_XY)) + 15.0     # hilltop + 15 m mast


def los_blocked(p_radar, p_target, n=400):
    """Terrain/earth-curvature LOS check via effective-earth flat mapping."""
    dx, dy = p_target[0] - p_radar[0], p_target[1] - p_radar[1]
    r = np.hypot(dx, dy)
    s = np.linspace(0, r, n)
    xs = p_radar[0] + dx * s / r
    ys = p_radar[1] + dy * s / r
    h_terr = np.maximum(terrain_h(xs, ys), 0.0) - s ** 2 / (2 * RE)
    h_tgt_flat = p_target[2] - r ** 2 / (2 * RE)
    h_line = p_radar[2] + (h_tgt_flat - p_radar[2]) * s / r
    return bool(np.any(h_terr[1:-1] > h_line[1:-1]))


# --- clutter patches ----------------------------------------------------------
PATCH = 2e3
WEDGE = np.deg2rad(8.0)
R_MIN, R_MAX = 10e3, 85e3

patches = []   # (x, y, z, rcs_m2, phase_deg, vx, vy)

# dense main-beam wedge over the sea + coarse all-azimuth ring for
# sidelobe clutter from the surrounding hills
cells = []
gx = np.arange(RADAR_XY[0] + R_MIN, RADAR_XY[0] + R_MAX, PATCH)
gy = np.arange(-13e3, 13e3 + 1, PATCH)
cells += [(px, py) for px in gx for py in gy]
ring_r = np.arange(R_MIN, 40e3, 3.5e3)
for r in ring_r:
    n = max(6, int(2 * np.pi * r / 3.5e3))
    for a in np.linspace(-np.pi, np.pi, n, endpoint=False):
        cells.append((RADAR_XY[0] + r * np.cos(a), RADAR_XY[1] + r * np.sin(a)))

for px, py in cells:
    rx_, ry_ = px - RADAR_XY[0], py - RADAR_XY[1]
    r = np.hypot(rx_, ry_)
    if r < R_MIN or r > R_MAX:
        continue
    h = float(terrain_h(px, py))
    p3 = (px, py, max(h, 0.0) + 2.0)
    if los_blocked((RADAR_XY[0], RADAR_XY[1], RADAR_Z), p3):
        continue
    psi = max(np.arctan2(RADAR_Z - max(h, 0.0), r), np.deg2rad(0.3))
    if h < 1.0:    # sea
        s0_db = -45 + 15 * np.log10(np.rad2deg(psi))
        spd = rng.uniform(0.3, 1.0)
    else:          # land
        s0_db = -15 + 10 * np.log10(np.sin(psi))
        spd = rng.uniform(0.0, 0.15)
    s0 = 10 ** (s0_db / 10)
    speckle = 10 ** (rng.normal(0, 2.0) / 10)
    rcs = s0 * PATCH ** 2 * speckle
    vdir = rng.uniform(0, 2 * np.pi)
    patches.append((px, py, max(h, 0.0), rcs, rng.uniform(0, 360),
                    spd * np.cos(vdir), spd * np.sin(vdir)))

print(f"clutter patches: {len(patches)} "
      f"(sea {sum(1 for p in patches if p[2] < 1)}, "
      f"land {sum(1 for p in patches if p[2] >= 1)})")

# --- contacts -----------------------------------------------------------------
CONTACTS = [
    # name, rcs m2, start (x,y,z) m, velocity m/s
    ("B-52",      100.0, (430e3,  20e3, 10e3), (-250.0, 0.0, 0.0)),
    ("B-1B",       10.0, (124e3,  -8e3,  200), (-300.0, 0.0, 0.0)),
    ("B-2",         0.1, (116e3,   4e3,  150), (-260.0, 0.0, 0.0)),
    ("airliner 1", 30.0, (330e3,  25e3, 10.5e3), (0.0, -240.0, 0.0)),
    ("airliner 2", 30.0, (360e3, -25e3, 11e3), (0.0,  240.0, 0.0)),
]

patch_targets = [
    {"location": (p[0], p[1], p[2]), "speed": (p[5], p[6], 0.0),
     "rcs": 10 * np.log10(p[3]), "phase": p[4]}
    for p in patches
]


def contact_pos(ci, t):
    _, _, p0, v = CONTACTS[ci]
    return np.array(p0) + np.array(v) * t


def visible(ci, t):
    return not los_blocked((RADAR_XY[0], RADAR_XY[1], RADAR_Z),
                           tuple(contact_pos(ci, t)))


# predicted pop-up time from the LOS model itself
popup_pred = {}
for ci, (name, *_rest) in enumerate(CONTACTS):
    ts = np.arange(0, 30.001, 0.5)
    vis = [visible(ci, t) for t in ts]
    popup_pred[name] = ts[np.argmax(vis)] if any(vis) else None
print("predicted visibility: " +
      ", ".join(f"{n} @ {popup_pred[n]:.0f} s" if popup_pred[n] is not None
                else f"{n} never" for n, *_ in CONTACTS))

RADAR_LOC = (RADAR_XY[0], RADAR_XY[1], RADAR_Z)

# --- simulation ---------------------------------------------------------------
samples = int(FS * PULSE)
bb = np.zeros((FRAMES, PULSES, samples), dtype=complex)
print(f"running sim_radar: {FRAMES} frames x {PULSES} pulses x {samples} "
      f"samples x ({len(patch_targets)} patches + visible contacts) ...")
for fi in range(FRAMES):
    t = float(fi)
    tgts = list(patch_targets)
    for ci, (name, rcs, p0, v) in enumerate(CONTACTS):
        if visible(ci, t):
            tgts.append({"location": tuple(contact_pos(ci, t)),
                         "speed": tuple(v),
                         "rcs": 10 * np.log10(rcs), "phase": 0})
    tx = Transmitter(
        f=[FC, FC + BW], t=PULSE, tx_power=90, pulses=PULSES, prp=PRP,
        channels=[{"location": RADAR_LOC,
                   "azimuth_angle": AZ_ANG, "azimuth_pattern": AZ_PAT,
                   "elevation_angle": EL_ANG, "elevation_pattern": EL_PAT}])
    rx = Receiver(
        fs=FS, noise_figure=3, rf_gain=20, baseband_gain=40,
        load_resistor=1000,
        channels=[{"location": RADAR_LOC,
                   "azimuth_angle": AZ_ANG, "azimuth_pattern": AZ_PAT,
                   "elevation_angle": EL_ANG, "elevation_pattern": EL_PAT}])
    radar = Radar(transmitter=tx, receiver=rx, frame_time=0.0, seed=100 + fi)
    data = sim_radar(radar, tgts)
    bb[fi] = data["baseband"][0] + data["noise"][0]
print("baseband shape:", bb.shape)

# --- processing ---------------------------------------------------------------
# Blackman windows: clutter peaks exceed 100 dB, so the 2D window
# sidelobes must be deep or their residue raises the far-range noise floor
rwin = np.blackman(samples)
dwin = np.blackman(PULSES)
range_axis = np.arange(samples) * FS / samples * C / (2 * SLOPE)
vel_axis = np.fft.fftshift(np.fft.fftfreq(PULSES, d=PRP)) * LAM / 2
VA = LAM / (4 * PRP)

NOISE_DBM = -174 + 3 + 10 * np.log10(FS)
GAIN_DB = 35 + 35
PG_DB = 10 * np.log10(samples * PULSES)


def predict_snr(rcs_m2, rng_m):
    pr_dbm = (90 + 10 * np.log10(LAM ** 2 / (4 * np.pi) ** 3)
              + 10 * np.log10(rcs_m2) - 40 * np.log10(rng_m) + GAIN_DB)
    return pr_dbm - NOISE_DBM + PG_DB


rds = []
for fi in range(FRAMES):
    rd = np.fft.fftshift(
        proc.range_doppler_fft(bb[fi:fi + 1], rwin=rwin, dwin=dwin)[0],
        axes=0)
    rds.append(np.abs(rd))

snr_t = np.full((len(CONTACTS), FRAMES), np.nan)
for ci, (name, rcs, p0, vel) in enumerate(CONTACTS):
    for fi in range(FRAMES):
        if not visible(ci, float(fi)):
            continue
        mag = rds[fi]
        pos = contact_pos(ci, float(fi)) - np.array(
            [RADAR_XY[0], RADAR_XY[1], RADAR_Z])
        r_true = np.linalg.norm(pos)
        vr_true = pos @ np.array(vel) / r_true
        vr_wrapped = (vr_true + VA) % (2 * VA) - VA
        ri = np.argmin(np.abs(range_axis - r_true))
        vi = np.argmin(np.abs(vel_axis - vr_wrapped))
        lo_r, hi_r = max(0, ri - 6), min(samples, ri + 7)
        ri_det = lo_r + int(np.argmax(mag[vi, lo_r:hi_r]))
        ring = np.delete(mag[:, lo_r:hi_r].ravel(),
                         vi * (hi_r - lo_r) + ri_det - lo_r)
        noise = np.median(ring)
        snr_t[ci, fi] = 20 * np.log10(mag[vi, ri_det] / noise)

# --- summary table ------------------------------------------------------------
print(f"\n{'contact':11s} {'RCS(m2)':>8s} {'pop_pred':>9s} {'det_frames':>10s}"
      f" {'SNR_meas':>9s} {'SNR_pred':>9s}")
for ci, (name, rcs, p0, vel) in enumerate(CONTACTS):
    det = np.where(snr_t[ci] > 15)[0]
    s_m = np.nanmean(snr_t[ci])
    s_p = predict_snr(rcs, np.linalg.norm(
        contact_pos(ci, 15.0) - [RADAR_XY[0], RADAR_XY[1], RADAR_Z]))
    first = f"{det[0]:d} s" if len(det) else "none"
    pop = (f"{popup_pred[name]:.0f} s" if popup_pred[name] is not None
           else "never")
    print(f"{name:11s} {rcs:8.3f} {pop:>9s} {len(det):4d}/{FRAMES} "
          f"(1st {first:>4s}) {s_m:8.1f}  {s_p:8.1f}")

# clutter profile: zero-Doppler rows vs range, mid-scenario
fi_show = FRAMES // 2
mag0 = rds[fi_show]
zero_rows = np.abs(vel_axis) < 1.0
clutter_prof = 20 * np.log10(mag0[zero_rows].max(axis=0) /
                             np.median(mag0) + 1e-12)

# --- plots --------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(15, 11))

# (1) terrain map
xg = np.linspace(-70e3, 90e3, 220)
yg = np.linspace(-40e3, 40e3, 160)
XG, YG = np.meshgrid(xg, yg)
ZG = terrain_h(XG, YG)
cs = ax[0, 0].contourf(XG / 1e3, YG / 1e3, ZG, levels=18, cmap="terrain")
ax[0, 0].contour(XG / 1e3, YG / 1e3, ZG, levels=[1.0], colors="navy",
                 linewidths=1.5)
fig.colorbar(cs, ax=ax[0, 0], label="terrain height (m)")
ax[0, 0].plot(RADAR_XY[0] / 1e3, RADAR_XY[1] / 1e3, "r^", ms=12,
              label="EWR site")
pv = np.array([(p[0], p[1]) for p in patches])
ax[0, 0].plot(pv[:, 0] / 1e3, pv[:, 1] / 1e3, ".", color="orange", ms=2,
              alpha=0.6, label="clutter patches (visible)")
for ci, (name, *_r) in enumerate(CONTACTS):
    tr = np.array([contact_pos(ci, t) for t in np.arange(0, 30.1, 1.0)])
    ax[0, 0].plot(tr[:, 0] / 1e3, tr[:, 1] / 1e3, lw=1.5, label=name)
ax[0, 0].set(xlabel="x (km)", ylabel="y (km)",
             title="Terrain, clutter patches, contact tracks", xlim=(-70, 90),
             ylim=(-40, 40), aspect="equal")
ax[0, 0].legend(fontsize=7, loc="upper left")

# (2) range-Doppler map
im = ax[0, 1].pcolormesh(range_axis / 1e3, vel_axis,
                         20 * np.log10(mag0 / mag0.max() + 1e-12),
                         cmap="jet", vmin=-60, vmax=0, shading="auto")
ax[0, 1].set(xlabel="range (km)", ylabel="radial velocity (m/s, aliased)",
             title=f"Range-Doppler map, t={fi_show} s", xlim=(0, 450))
fig.colorbar(im, ax=ax[0, 1], label="dB")

# (3) clutter profile
ax[1, 0].plot(range_axis / 1e3, clutter_prof)
ax[1, 0].axvline(R_MAX / 1e3, color="g", ls="--",
                 label=f"patch limit {R_MAX/1e3:.0f} km")
ax[1, 0].set(xlabel="range (km)", ylabel="level over median (dB)",
             title=f"Zero-Doppler clutter profile, t={fi_show} s",
             xlim=(0, 200))
ax[1, 0].legend(); ax[1, 0].grid(alpha=0.3)

# (4) SNR per contact
t_ax = np.arange(FRAMES)
for ci, (name, *_r) in enumerate(CONTACTS):
    ax[1, 1].plot(t_ax, snr_t[ci], ".-", ms=4, label=name)
ax[1, 1].axhline(15, color="r", ls="--", label="detection ~15 dB")
ax[1, 1].set(xlabel="time (s)", ylabel="SNR (dB)",
             title="SNR per contact (gaps = terrain-masked)")
ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.3)

fig.suptitle("Demo 4 — UHF EWR on hilly coastal upland: terrain clutter, "
             "sea clutter, horizon masking")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "demo4_ewr_terrain.png"), dpi=130)
print(f"figure: {os.path.join(OUT, 'demo4_ewr_terrain.png')}")
