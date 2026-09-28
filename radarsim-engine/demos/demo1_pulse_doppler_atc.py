#!/usr/bin/env python
"""Demo 1 — Pulse-Doppler surveillance radar vs a civil airliner (30 s).

Scenario
--------
- S-band pulse-Doppler radar (fc 3.3 GHz, 5 MHz chirp bandwidth,
  PRF 12 kHz, 64-pulse CPI, 10 frames/s over 30 s).
- A civil airliner (RCS 100 m^2) flies along an airway at 240 m/s and
  1 km height, crossing the radar's vicinity with a 3 km offset.
- Background: 25 static ground-clutter patches plus receiver noise.

Outputs (saved to out/):
- Range-Doppler map at the midpoint frame (target vs zero-Doppler clutter)
- Measured vs true range track
- Measured vs true radial-velocity track
- Estimated SNR vs time

Run from the repo root (after tools/build_and_integrate.sh):
    python radarsim-engine/demos/demo1_pulse_doppler_atc.py
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

# --- radar parameters -------------------------------------------------------
FC = 3.3e9
BW = 5e6                 # 30 m range resolution
PULSE = 60e-6
PRP = 1 / 12e3           # PRF 12 kHz -> unambiguous velocity +-272 m/s
PULSES = 64              # CPI
FS = 10e6                # 600 samples -> 9 km instrumented window
FRAMES = 300
FRAME_TIME = 0.1         # 30 s at 10 Hz

LAM = C / (FC + BW / 2)

# --- airway: airliner crossing with a 3 km offset at 1 km height ------------
V_AIR = 240.0            # m/s (cruise, low level for the demo window)
air0 = np.array([-5000.0, 3000.0, 1000.0])
air_v = np.array([V_AIR, 0.0, 0.0])

rng = np.random.default_rng(7)
clutter = [
    {"location": (float(x), float(y), 0.0), "rcs": float(r), "phase": 0}
    for x, y, r in zip(
        rng.uniform(500, 6000, 25),
        rng.uniform(-4000, 4000, 25),
        rng.uniform(-20, -10, 25),
    )
]

targets = [
    {"location": tuple(air0), "speed": tuple(air_v), "rcs": 20.0, "phase": 0},
    *clutter,
]

tx = Transmitter(
    f=[FC, FC + BW],
    t=PULSE,
    tx_power=60,             # 1 kW
    pulses=PULSES,
    prp=PRP,
    channels=[{"location": (0, 0, 0)}],
)
rx = Receiver(
    fs=FS,
    noise_figure=5,
    rf_gain=20,
    baseband_gain=40,
    load_resistor=1000,
    channels=[{"location": (0, 0, 0)}],
)
radar = Radar(
    transmitter=tx,
    receiver=rx,
    frame_time=np.arange(FRAMES) * FRAME_TIME,
    seed=42,
)

print("running sim_radar: "
      f"{FRAMES} frames x {PULSES} pulses x "
      f"{radar.sample_prop['samples_per_pulse']} samples ...")
data = sim_radar(radar, targets)
bb = data["baseband"] + data["noise"]  # [frames*channels, pulses, samples]
ts = data["timestamp"]
print("baseband shape:", bb.shape)

# --- processing -------------------------------------------------------------
samples = bb.shape[2]
rwin = np.hamming(samples)
dwin = np.hamming(PULSES)

# deramp (stretch) processing: beat bin k -> f_b = k*fs/N -> R = f_b*c/(2*slope)
SLOPE = BW / PULSE
range_axis = np.arange(samples) * FS / samples * C / (2 * SLOPE)
vel_axis = np.fft.fftshift(
    np.fft.fftfreq(PULSES, d=PRP)) * LAM / 2

det_range = np.full(FRAMES, np.nan)
det_vel = np.full(FRAMES, np.nan)
det_snr = np.full(FRAMES, np.nan)

# per-frame noise floor from the first frame's empty far range cells
for fi in range(FRAMES):
    rd = proc.range_doppler_fft(bb[fi:fi + 1], rwin=rwin, dwin=dwin)[0]
    rd = np.fft.fftshift(rd, axes=0)
    mag = np.abs(rd)
    # restrict to plausible airway ranges; the PD clutter notch excludes
    # zero Doppler (mainlobe clutter rejection)
    rmask = (range_axis > 2000) & (range_axis < 8000)
    vmask = np.abs(vel_axis) > 20       # exclude zero-Doppler clutter
    sub = mag[np.ix_(vmask, rmask)]
    pk = np.unravel_index(np.argmax(sub), sub.shape)
    noise = np.median(sub)
    snr = 20 * np.log10(sub[pk] / noise)
    det_snr[fi] = snr
    if snr >= 15.0:  # detection threshold; below it, no report (PD fade)
        det_range[fi] = range_axis[rmask][pk[1]]
        det_vel[fi] = vel_axis[vmask][pk[0]]
    if fi == FRAMES // 2:
        rd_show = 20 * np.log10(mag / mag.max() + 1e-12)

# truth
t = np.arange(FRAMES) * FRAME_TIME
true_pos = air0[None, :] + t[:, None] * air_v[None, :]
true_range = np.linalg.norm(true_pos, axis=1)
true_vr = (true_pos @ air_v) / true_range

# --- plots ------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(13, 9))

im = ax[0, 0].pcolormesh(range_axis / 1e3, vel_axis, rd_show,
                         cmap="jet", vmin=-60, vmax=0, shading="auto")
ax[0, 0].set(xlabel="range (km)", ylabel="radial velocity (m/s)",
             title=f"Range-Doppler map, t={t[FRAMES//2]:.1f} s")
fig.colorbar(im, ax=ax[0, 0], label="dB (peak-normalized)")

ax[0, 1].plot(t, true_range / 1e3, label="truth")
ax[0, 1].plot(t, det_range / 1e3, ".", ms=3, label="CFAR peak")
ax[0, 1].set(xlabel="time (s)", ylabel="range (km)",
             title="Range track"); ax[0, 1].legend(); ax[0, 1].grid(alpha=0.3)

ax[1, 0].plot(t, true_vr, label="truth")
ax[1, 0].plot(t, det_vel, ".", ms=3, label="CFAR peak")
ax[1, 0].set(xlabel="time (s)", ylabel="radial velocity (m/s)",
             title="Doppler track"); ax[1, 0].legend(); ax[1, 0].grid(alpha=0.3)

ax[1, 1].plot(t, det_snr)
ax[1, 1].set(xlabel="time (s)", ylabel="SNR (dB)",
             title="Estimated SNR at detection"); ax[1, 1].grid(alpha=0.3)

fig.suptitle("Demo 1 — Pulse-Doppler radar, civil airliner on airway "
             "(RCS 100 m$^2$, 240 m/s, clutter + noise)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "demo1_pulse_doppler.png"), dpi=130)

# --- summary ----------------------------------------------------------------
ok = ~np.isnan(det_range)
rng_err = det_range[ok] - true_range[ok]
vel_err = det_vel[ok] - true_vr[ok]
print(f"detected in {ok.sum()}/{FRAMES} frames "
      f"(fade near tangent crossing at ~21 s is the PD clutter notch)")
print(f"SNR (detected frames): min/med/max = {np.nanmin(det_snr):.1f}/"
      f"{np.nanmedian(det_snr):.1f}/{np.nanmax(det_snr):.1f} dB")
print(f"range track rms error:    {np.sqrt(np.mean(rng_err**2)):.2f} m "
      f"(resolution {C/(2*BW):.0f} m)")
print(f"velocity track rms error: {np.sqrt(np.mean(vel_err**2)):.2f} m/s "
      f"(resolution {LAM/(2*PULSES*PRP):.1f} m/s)")
print(f"figure: {os.path.join(OUT, 'demo1_pulse_doppler.png')}")
