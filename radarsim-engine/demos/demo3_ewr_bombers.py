#!/usr/bin/env python
"""Demo 3 — UHF early-warning radar (EWR) vs bombers (30 s).

Scenario
--------
- UHF EWR (fc 435 MHz, 200 kHz chirp, PRF 100 Hz, 16-pulse CPI,
  1 MW tx, 35+35 dB antenna gains), 0.5 frames/s over 30 s.
- Inbound raid at altitude:
    B-52   RCS 100 m2   250 m/s from 400 km
    B-1B   RCS 10 m2    300 m/s from 330 km
    B-2    RCS 0.1 m2   250 m/s from 260 km (UHF resonance region:
           shaping optimized for microwave bands is ineffective)
- Background traffic: two civil airliners (30 m2) on distant airways,
  plus receiver noise.

Outputs (saved to out/):
- Range-Doppler map at mid-scenario
- Range tracks of all contacts vs time
- SNR vs time per contact
- Detection table: measured vs radar-equation SNR, 15 dB detection range

Run from the repo root:
    python radarsim-engine/demos/demo3_ewr_bombers.py
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
FC = 435e6
BW = 200e3               # 750 m range resolution (typical for EWR)
PULSE = 3e-3             # receive window covers 420 km
PRP = 10e-3              # PRF 100 Hz -> unambiguous range 1500 km
PULSES = 16
FS = 500e3               # 1500 samples
FRAMES = 15              # 2 s per frame, 30 s

LAM = C / (FC + BW / 2)
SLOPE = BW / PULSE

# --- raid + background traffic ----------------------------------------------
CONTACTS = [
    # name, rcs m2, start (x,y,z) m, velocity m/s
    ("B-52",      100.0, (400e3,  20e3, 12e3), (-250.0, 0.0, 0.0)),
    ("B-1B",       10.0, (330e3, -15e3, 11e3), (-300.0, 0.0, 0.0)),
    ("B-2",         0.1, (260e3,  10e3, 12e3), (-250.0, 0.0, 0.0)),
    ("airliner 1", 30.0, (180e3,  40e3, 10e5 * 0 + 10e3), (0.0, -240.0, 0.0)),
    ("airliner 2", 30.0, (220e3, -50e3,  9e3), (0.0,  240.0, 0.0)),
]

targets = [
    {"location": tuple(p0), "speed": tuple(v),
     "rcs": 10 * np.log10(rcs), "phase": 0}
    for _, rcs, p0, v in CONTACTS
]

tx = Transmitter(
    f=[FC, FC + BW], t=PULSE, tx_power=90,  # 1 MW
    pulses=PULSES, prp=PRP,
    channels=[{"location": (0, 0, 0),
               "azimuth": [-90, 90], "azimuth_pattern": [35, 35],
               "elevation": [-90, 90], "elevation_pattern": [0, 0]}],
)
rx = Receiver(
    fs=FS, noise_figure=3, rf_gain=20, baseband_gain=40, load_resistor=1000,
    channels=[{"location": (0, 0, 0),
               "azimuth": [-90, 90], "azimuth_pattern": [35, 35],
               "elevation": [-90, 90], "elevation_pattern": [0, 0]}],
)
radar = Radar(transmitter=tx, receiver=rx,
              frame_time=np.arange(FRAMES) * 2.0, seed=11)

print("running sim_radar: "
      f"{FRAMES} frames x {PULSES} pulses x "
      f"{radar.sample_prop['samples_per_pulse']} samples ...")
data = sim_radar(radar, targets)
bb = data["baseband"] + data["noise"]
print("baseband shape:", bb.shape)

# --- processing -------------------------------------------------------------
samples = bb.shape[2]
rwin = np.hamming(samples)
dwin = np.hamming(PULSES)
range_axis = np.arange(samples) * FS / samples * C / (2 * SLOPE)
vel_axis = np.fft.fftshift(np.fft.fftfreq(PULSES, d=PRP)) * LAM / 2
VA = LAM / (4 * PRP)

NOISE_DBM = -174 + 3 + 10 * np.log10(FS)
GAIN_DB = 35 + 35
PG_DB = 10 * np.log10(samples * PULSES)


def predict_snr(rcs_m2, rng_m):
    pr_dbm = (90 + 10 * np.log10(LAM**2 / (4 * np.pi)**3)
              + 10 * np.log10(rcs_m2) - 40 * np.log10(rng_m) + GAIN_DB)
    return pr_dbm - NOISE_DBM + PG_DB


rds_all = []
for fi in range(FRAMES):
    rd = np.fft.fftshift(
        proc.range_doppler_fft(bb[fi:fi + 1], rwin=rwin, dwin=dwin)[0],
        axes=0)
    rds_all.append(np.abs(rd))

snr_t = np.zeros((len(CONTACTS), FRAMES))
rng_meas = np.zeros((len(CONTACTS), FRAMES))
rng_true = np.zeros((len(CONTACTS), FRAMES))
for ci, (name, rcs, p0, vel) in enumerate(CONTACTS):
    for fi in range(FRAMES):
        mag = rds_all[fi]
        pos = np.array(p0) + np.array(vel) * fi * 2.0
        r_true = np.linalg.norm(pos)
        vr_true = pos @ np.array(vel) / r_true
        vr_wrapped = (vr_true + VA) % (2 * VA) - VA
        rng_true[ci, fi] = r_true
        ri = np.argmin(np.abs(range_axis - r_true))
        vi = np.argmin(np.abs(vel_axis - vr_wrapped))
        lo_r, hi_r = max(0, ri - 6), min(samples, ri + 7)
        # range-Doppler coupling shifts the LFM peak by f_d*c/(2*slope)
        ri_det = lo_r + int(np.argmax(mag[vi, lo_r:hi_r]))
        ring = np.delete(mag[:, lo_r:hi_r].ravel(),
                         vi * (hi_r - lo_r) + ri_det - lo_r)
        noise = np.median(ring)
        snr_t[ci, fi] = 20 * np.log10(mag[vi, ri_det] / noise)
        rng_meas[ci, fi] = range_axis[ri_det]

# --- summary table ----------------------------------------------------------
print(f"\n{'contact':11s} {'RCS(m2)':>8s} {'range0':>7s} {'SNR_meas':>9s} "
      f"{'SNR_pred':>9s} {'R_det15':>8s}")
for ci, (name, rcs, p0, vel) in enumerate(CONTACTS):
    r0 = np.linalg.norm(p0)
    s_m, s_p = snr_t[ci].mean(), predict_snr(rcs, r0)
    print(f"{name:11s} {rcs:8.3f} {r0/1e3:6.0f}k {s_m:8.1f}  {s_p:8.1f}  "
          f"{r0/1e3 * 10**((s_m-15)/40):7.0f}k")

# --- plots ------------------------------------------------------------------
t = np.arange(FRAMES) * 2.0
fig, ax = plt.subplots(1, 3, figsize=(16, 5))

fi_show = FRAMES // 2
mag = rds_all[fi_show]
im = ax[0].pcolormesh(range_axis / 1e3, vel_axis,
                      20 * np.log10(mag / mag.max() + 1e-12),
                      cmap="jet", vmin=-50, vmax=0, shading="auto")
ax[0].set(xlabel="range (km)", ylabel="radial velocity (m/s, aliased)",
          title=f"Range-Doppler map, t={fi_show*2} s", xlim=(100, 450))
fig.colorbar(im, ax=ax[0], label="dB")

for ci, (name, *_ ) in enumerate(CONTACTS):
    ax[1].plot(t, rng_true[ci] / 1e3, ls="--", alpha=0.5)
    ax[1].plot(t, rng_meas[ci] / 1e3, ".", ms=4, label=name)
ax[1].set(xlabel="time (s)", ylabel="range (km)",
          title="Range tracks (dots = measured, dashed = truth)")
ax[1].legend(); ax[1].grid(alpha=0.3)

for ci, (name, *_ ) in enumerate(CONTACTS):
    ax[2].plot(t, snr_t[ci], ".-", ms=4, label=name)
ax[2].axhline(15, color="r", ls="--", label="CFAR ~15 dB")
ax[2].set(xlabel="time (s)", ylabel="SNR (dB)", title="SNR per contact")
ax[2].legend(); ax[2].grid(alpha=0.3)

fig.suptitle("Demo 3 — UHF EWR vs B-52 / B-1B / B-2 "
             "(+ civil background traffic)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "demo3_ewr.png"), dpi=130)
print(f"figure: {os.path.join(OUT, 'demo3_ewr.png')}")
