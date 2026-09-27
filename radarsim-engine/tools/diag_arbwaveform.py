"""Compare candidate waveform-phase integration schemes against two golden
datasets: the ideal arbitrary-waveform test (4 samples) and the system
arbitrary-waveform test (16 samples, nonlinear 100-point sweep)."""
import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
REPO = TOOLS.parents[1]
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(REPO))

from pysim import shim
shim.install()
import numpy as np
from radarsimpy import Radar, Transmitter, Receiver

C = 299792458.0

# ---------------- dataset 1: ideal test ----------------
f1 = np.array([24.075e9, 24.175e9, 26e9, 28e9, 26e9])
t1 = np.array([0, 20e-6, 40e-6, 60e-6, 80e-6])
fs1 = 6e4
tau1 = 20 / C
gold1 = np.array([0.020911267027258873 + 0.015191296115517616j,
                  0.0165250264108181 - 0.019874105229973793j,
                  -0.022241275757551193 + 0.013167467899620533j,
                  0.025563737377524376 + 0.0038147002924233675j])

# ---------------- dataset 2: system test ----------------
# extract the literal frequency array from the test source
import re
txt = (REPO / "tests" / "test_system_arbitrary_waveform.py").read_text(
    encoding="utf-8"
)
m = re.search(r"freq_nonlinear = np\.array\(\s*\[([\d.eE\s,+-]+)\]", txt)
f2 = np.array([float(x) for x in m.group(1).split(",") if x.strip()])
t2 = np.linspace(0, 80e-6, 100)
fs2 = 2e5
# moving target: pos(T) = 200 - 5 T  -> tau varies per sample
def tau2(s):
    T = s / fs2
    return 2 * (200 - 5 * T) / C
gold2 = np.array([
    2.18508231e-03 + 0.00029834j, 1.94347688e-03 + 0.00104235j,
    2.15462260e-03 - 0.00047032j, 2.74110357e-04 - 0.00218826j,
    -2.20141465e-03 + 0.00013182j, 1.53208282e-03 + 0.00158629j,
    -8.78069136e-04 - 0.00202302j, 1.23647175e-03 + 0.00182613j,
    -2.12812091e-03 - 0.00057854j, 1.25044390e-03 - 0.00181659j,
    2.01748455e-03 + 0.00089071j, 6.40810057e-04 + 0.00211021j,
    -9.06653710e-05 + 0.0022035j, 3.22449600e-04 + 0.00218166j,
    1.67007645e-03 + 0.0014403j, 1.90978212e-03 - 0.00110289j,
])

# ---------------- candidate phase models ----------------
def make_models(f, t, fs):
    kseg = np.diff(f) / np.diff(t)
    dt = np.diff(t)
    phi_exact_knots = np.concatenate(([0.0], np.cumsum(f[:-1] * dt + 0.5 * kseg * dt**2)))

    def phi_exact(u):
        u = np.asarray(u, float)
        idx = np.clip(np.searchsorted(t, u, side="right") - 1, 0, len(t) - 2)
        du = u - t[idx]
        return phi_exact_knots[idx] + f[idx] * du + 0.5 * kseg[idx] * du**2

    def phi_trapz_linear(u):  # knots exact, linear between
        u = np.asarray(u, float)
        # beyond grid: linear extrapolation of f (needed for chirp tests)
        out = np.interp(u, t, phi_exact_knots)
        below = u < t[0]
        if np.any(below):
            out = np.asarray(out)
            out[below] = phi_exact_knots[0] + f[0] * (u[below] - t[0]) + \
                0.5 * kseg[0] * (u[below] - t[0])**2
        above = u > t[-1]
        if np.any(above):
            out[above] = phi_exact_knots[-1] + f[-1] * (u[above] - t[-1]) + \
                0.5 * kseg[-1] * (u[above] - t[-1])**2
        return out

    phi_left_knots = np.concatenate(([0.0], np.cumsum(f[:-1] * dt)))
    def phi_left(u):
        return np.interp(u, t, phi_left_knots)

    phi_right_knots = np.concatenate(([0.0], np.cumsum(f[1:] * dt)))
    def phi_right(u):
        return np.interp(u, t, phi_right_knots)

    # ADC grid: exact phase at ADC sample times, linear interp between
    t_adc = np.arange(0, t[-1] + 0.5 / fs, 1 / fs)
    phi_adc = phi_exact(t_adc)
    def phi_adcgrid(u):
        return np.interp(u, t_adc, phi_adc)

    def beat_midpoint(u, tau):
        um = np.asarray(u, float) - tau / 2
        idx = np.clip(np.searchsorted(t, um, side="right") - 1, 0, len(t) - 2)
        fm = f[idx] + kseg[idx] * (um - t[idx])
        return fm * tau

    return {
        "exact": lambda u, tau: phi_exact(u) - phi_exact(u - tau),
        "trapz_linear": lambda u, tau: phi_trapz_linear(u) - phi_trapz_linear(u - tau),
        "left": lambda u, tau: phi_left(u) - phi_left(u - tau),
        "right": lambda u, tau: phi_right(u) - phi_right(u - tau),
        "adcgrid": lambda u, tau: phi_adcgrid(u) - phi_adcgrid(u - tau),
        "midpoint": beat_midpoint,
    }

print("=== dataset 1 (ideal, 4 samples) ===")
g1 = np.angle(gold1) / (2 * np.pi)
for name, beat in make_models(f1, t1, fs1).items():
    ph = np.array([beat(s / fs1, tau1) for s in range(4)]) % 1
    err = np.abs(((ph - g1 + 0.5) % 1) - 0.5).max()
    print(f"  {name:14s} max phase err (cycles): {err:.6f}   {ph}")

print("=== dataset 2 (system, 16 samples) ===")
g2 = np.unwrap(np.angle(gold2)) / (2 * np.pi)
for name, beat in make_models(f2, t2, fs2).items():
    ph = np.array([beat(s / fs2, tau2(s)) for s in range(16)])
    err = np.abs(((ph - g2 + 0.5) % 1) - 0.5)
    print(f"  {name:14s} max phase err (cycles): {err.max():.6f}")
