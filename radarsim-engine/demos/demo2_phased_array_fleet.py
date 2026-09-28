#!/usr/bin/env python
"""Demo 2 — Phased-array radar vs a mixed fleet (30 s).

Scenario
--------
- S-band AESA (fc 3 GHz, 2 MHz chirp, PRF 2 kHz, 16-pulse CPI,
  4-element Rx ULA at lambda/2), 1 frame/s over 30 s.
- Traffic:
    airliner   RCS 100 m2    240 m/s tangential at ~50 km, az +9 deg
    F-16       RCS 4 m2      400 m/s inbound from ~40 km, az -4 deg
    F-35       RCS 0.001 m2  300 m/s inbound from ~18 km, az +3 deg
    B-2        RCS 0.0001 m2 250 m/s inbound from ~12 km, az -6 deg
- Background: receiver noise only (antennas isotropic, 0 dB — conservative).

Outputs (saved to out/):
- Range-Doppler map with CFAR detections marked
- Coherent beamforming (azimuth spectrum) on the strongest detections
- Per-target detection SNR table (measured vs radar-equation prediction)

Run from the repo root:
    python radarsim-engine/demos/demo2_phased_array_fleet.py
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
OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)

# --- radar ------------------------------------------------------------------
FC = 3.0e9
BW = 2e6                  # 75 m range resolution
PULSE = 500e-6            # receive window covers the 48 km airliner
PRP = 500e-6              # PRF 2 kHz -> unambiguous range 75 km
PULSES = 16
FS = 4e6                  # 2000 samples -> 150 km window
N_RX = 4
FRAMES = 30               # 1 fps, 30 s

LAM = C / (FC + BW / 2)
SLOPE = BW / PULSE

# --- fleet ------------------------------------------------------------------
FLEET = [
    # name, rcs m2, start (x,y,z) m, velocity (vx,vy,vz) m/s
    ("airliner", 1e2,  (48e3,  8e3, 10e3), (0.0, -240.0, 0.0)),
    ("F-16",     4.0,  (40e3, -3e3,  8e3), (-410.0, 0.0, 0.0)),
    ("F-35",     1e-3, (18e3,  1e3,  7e3), (-315.0, 0.0, 0.0)),
    ("B-2",      1e-4, (12e3, -1.3e3, 9e3), (-235.0, 0.0, 0.0)),
]

targets = [
    {"location": tuple(p0), "speed": tuple(v),
     "rcs": 10 * np.log10(rcs), "phase": 0}
    for _, rcs, p0, v in FLEET
]

tx = Transmitter(
    f=[FC, FC + BW], t=PULSE, tx_power=90,  # 1 MW class AESA
    pulses=PULSES, prp=PRP,
    channels=[{"location": (0, 0, 0),
               "azimuth": [-90, 90], "azimuth_pattern": [25, 25],
               "elevation": [-90, 90], "elevation_pattern": [0, 0]}],
)
rx = Receiver(
    fs=FS, noise_figure=6, rf_gain=20, baseband_gain=40, load_resistor=1000,
    channels=[{"location": (0, n * LAM / 2, 0),
               "azimuth": [-90, 90], "azimuth_pattern": [25, 25],
               "elevation": [-90, 90], "elevation_pattern": [0, 0]}
              for n in range(N_RX)],
)
radar = Radar(transmitter=tx, receiver=rx,
              frame_time=np.arange(FRAMES) * 1.0, seed=7)

print("running sim_radar: "
      f"{FRAMES} frames x {PULSES} pulses x "
      f"{radar.sample_prop['samples_per_pulse']} samples x {N_RX} rx ...")
data = sim_radar(radar, targets)
bb = data["baseband"] + data["noise"]   # [frames*rx, pulses, samples]
print("baseband shape:", bb.shape)

# --- processing -------------------------------------------------------------
samples = bb.shape[2]
rwin = np.hamming(samples)
dwin = np.hamming(PULSES)
range_axis = np.arange(samples) * FS / samples * C / (2 * SLOPE)
vel_axis = np.fft.fftshift(np.fft.fftfreq(PULSES, d=PRP)) * LAM / 2

# predicted SNR at CPI output (radar equation; antenna gains 25+25 dB)
NOISE_DBM = -174 + 6 + 10 * np.log10(FS)  # kTB + NF, at the detector
GAIN_DB = 25 + 25                          # tx + rx pattern gain
PG_DB = 10 * np.log10(samples * PULSES)    # ideal integration gain


def predict_snr(rcs_m2, rng_m):
    pr_dbm = (90 + 10 * np.log10(LAM**2 / (4 * np.pi)**3)
              + 10 * np.log10(rcs_m2) - 40 * np.log10(rng_m) + GAIN_DB)
    return pr_dbm - NOISE_DBM + PG_DB


rows = []
beam_scans = {}
VA = LAM / (4 * PRP)   # unambiguous radial velocity (+-VA)
for name, rcs, p0, vel in FLEET:
    snrs = []
    az = []
    az_t = []
    for fi in range(FRAMES):
        rds = np.stack(
            [np.fft.fftshift(
                proc.range_doppler_fft(bb[fi * N_RX + n:fi * N_RX + n + 1],
                                       rwin=rwin, dwin=dwin)[0], axes=0)
             for n in range(N_RX)])       # [rx, doppler, range]
        mag = np.abs(rds).sum(axis=0)
        # expected cell of THIS target (nearest bin); measure SNR there
        pos = np.array(p0) + np.array(vel) * fi
        r_true = np.linalg.norm(pos)
        vr_true = pos @ np.array(vel) / r_true
        vr_wrapped = (vr_true + VA) % (2 * VA) - VA   # low-PRF aliasing
        az_t.append(np.degrees(np.arctan2(pos[1], pos[0])))
        ri = np.argmin(np.abs(range_axis - r_true))
        vi = np.argmin(np.abs(vel_axis - vr_wrapped))
        lo_r, hi_r = max(0, ri - 8), min(samples, ri + 9)
        # range-Doppler coupling shifts the LFM peak by f_d*c/(2*slope)
        # (up to ~4 bins at 410 m/s): measure at the local peak
        ri_det = lo_r + int(np.argmax(mag[vi, lo_r:hi_r]))
        ring = np.delete(mag[:, lo_r:hi_r].ravel(),
                         vi * (hi_r - lo_r) + ri_det - lo_r)
        noise = np.median(ring)
        cell = mag[vi, ri_det]
        snrs.append(20 * np.log10(cell / noise))
        # coherent azimuth estimate from the target cell across rx channels
        s = rds[:, vi, ri_det]
        n_scan = 181
        angs = np.radians(np.linspace(-20, 20, n_scan))
        steer = np.exp(-1j * 2 * np.pi * (LAM / 2) / LAM
                       * np.arange(N_RX)[:, None] * np.sin(angs)[None, :])
        spec = np.abs(steer.conj().T @ s) ** 2
        az.append(np.degrees(angs[np.argmax(spec)]))
    rows.append((name, rcs, np.linalg.norm(p0) / 1e3,
                 np.nanmean(snrs), predict_snr(rcs, np.linalg.norm(p0)),
                 np.nanmean(az), np.mean(az_t)))
    beam_scans[name] = (angs, spec / spec.max())

print(f"\n{'target':9s} {'RCS(m2)':>9s} {'range0':>7s} {'SNR_meas':>9s} "
      f"{'SNR_pred':>9s} {'az_meas':>8s} {'az_true':>8s} {'R_det15':>8s}")
for name, rcs, r0, snr_m, snr_p, az_m, az_t in rows:
    # range at which the target would fall to the 15 dB CFAR threshold
    r_det = r0 * 10 ** ((snr_m - 15) / 40)
    print(f"{name:9s} {rcs:9.4f} {r0:6.1f}k {snr_m:8.1f}  {snr_p:8.1f}  "
          f"{az_m:7.2f}  {az_t:7.2f} {r_det:7.1f}k")

# --- plots ------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(16, 5))

fi_show = 10
rds = np.stack(
    [np.fft.fftshift(
        proc.range_doppler_fft(bb[fi_show * N_RX + n:fi_show * N_RX + n + 1],
                               rwin=rwin, dwin=dwin)[0], axes=0)
     for n in range(N_RX)])
mag = np.abs(rds).sum(axis=0)
im = ax[0].pcolormesh(range_axis / 1e3, vel_axis,
                      20 * np.log10(mag / mag.max() + 1e-12),
                      cmap="jet", vmin=-50, vmax=0, shading="auto")
ax[0].set(xlabel="range (km)", ylabel="radial velocity (m/s)",
          title=f"Range-Doppler map, t={fi_show} s", xlim=(0, 60))
fig.colorbar(im, ax=ax[0], label="dB")

for name, (angs, spec) in beam_scans.items():
    ax[1].plot(np.degrees(angs), 10 * np.log10(spec + 1e-12), label=name)
ax[1].set(xlabel="azimuth (deg)", ylabel="normalized power (dB)",
          title="Beamforming spectrum at each target cell", xlim=(-20, 20),
          ylim=(-40, 3))
ax[1].grid(alpha=0.3); ax[1].legend()

names = [r[0] for r in rows]
ax[2].bar(np.arange(len(rows)) - 0.2, [r[3] for r in rows], 0.4,
          label="measured")
ax[2].bar(np.arange(len(rows)) + 0.2, [r[4] for r in rows], 0.4,
          label="radar equation")
ax[2].axhline(15, color="r", ls="--", label="CFAR ~15 dB")
ax[2].set_xticks(range(len(rows)), names)
ax[2].set(ylabel="SNR at CPI output (dB)", title="Detectability by RCS")
ax[2].legend(); ax[2].grid(alpha=0.3)

fig.suptitle("Demo 2 — Phased array vs mixed fleet "
             "(airliner / F-16 / F-35 / B-2)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "demo2_phased_array.png"), dpi=130)
print(f"figure: {os.path.join(OUT, 'demo2_phased_array.png')}")
