"""Diag round 2: verify modulation rule, freq_offset failure, arbitrary waveform."""
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
from pysim import reference
from pysim.reference import _waveform_phase_fn, C

# ---------- modulation rule check ----------
mod_t = np.array([0, 10e-6, 20e-6, 30e-6, 40e-6])
amp = np.array([0, 1, 0, 3, 4]); phs = np.array([0, 90, 180, -90, -180])
var = amp * np.exp(1j * np.radians(phs))
step = mod_t[1] - mod_t[0]
tau = 20 / C
print("mod rule idx = floor((u-tau)/step)+1, out-of-range -> 0:")
for s in range(4):
    u = s / 6e4
    idx = int(np.floor((u - tau - mod_t[0]) / step)) + 1
    v = var[idx] if 0 <= idx < len(var) else 0.0
    print(f"  s{s}: u={u*1e6:7.2f}us idx={idx} var={v}")

# ---------- freq_offset ----------
tx = Transmitter(f=[24.075e9, 24.175e9], t=80e-6, tx_power=10, prp=100e-6,
                 pulses=3, f_offset=[0, 1e6, 2e6], channels=[{'location': (0,0,0)}])
rx = Receiver(fs=6e4, noise_figure=12, rf_gain=20, load_resistor=500,
              baseband_gain=30, channels=[{'location': (0,0,0)}])
radar = Radar(transmitter=tx, receiver=rx)
res = reference.sim_radar_reference(radar, [{'location': np.array([10,0,0]), 'rcs': 20}])
gold = np.array([
    [0.02167871966958046 + 0.017555851489305496j, -0.027893975377082825 + 0.00031773850787431j,
     0.021273192018270493 - 0.018045110628008842j, -0.004863050766289234 + 0.02746862918138504j],
    [0.012656153179705143 + 0.02485823817551136j, -0.025607381016016006 - 0.011062202043831348j,
     0.026774795725941658 - 0.007824353873729706j, -0.015621167607605457 + 0.023110372945666313j],
    [0.0014430786250159144 + 0.027856115251779556j, -0.018887367099523544 - 0.020525911822915077j,
     0.027640212327241898 + 0.0037502485793083906j, -0.023673780262470245 + 0.014751197770237923j],
])
err = np.abs(res['baseband'][0] - gold)
print("\nfreq_offset max abs err:", err.max())
print("mine[1,0]:", res['baseband'][0,1,0], " gold:", gold[1,0])
print("wf f_offset:", radar.radar_prop['transmitter'].waveform_prop['f_offset'])

# ---------- arbitrary waveform ----------
tx = Transmitter(f=[24.075e9, 24.175e9, 26e9, 28e9, 26e9],
                 t=[0, 20e-6, 40e-6, 60e-6, 80e-6],
                 tx_power=10, prp=100e-6, pulses=3, channels=[{'location': (0,0,0)}])
radar2 = Radar(transmitter=tx, receiver=rx)
res2 = reference.sim_radar_reference(radar2, [{'location': np.array([10,0,0]), 'rcs': 20}])
import re
txt = (REPO / 'tests' / 'test_module_sim_radar_ideal.py').read_text(encoding='utf-8')
m = re.search(r'def test_simc_arbitrary_waveform\(\):.*?result\["baseband"\],\s*np\.array\((\[.*?\])\s*,?\s*\)', txt, re.S)
gold2 = np.array(eval(m.group(1)))[0]
print("\narbitrary waveform max abs err:", np.abs(res2['baseband'][0] - gold2).max())
print("mine  p0:", res2['baseband'][0,0,:])
print("gold  p0:", gold2[0])
f = np.array([24.075e9, 24.175e9, 26e9, 28e9, 26e9]); t = np.array([0, 20e-6, 40e-6, 60e-6, 80e-6])
phi = _waveform_phase_fn(f, t)
print("wf props f:", radar2.radar_prop['transmitter'].waveform_prop['f'])
print("wf props t:", radar2.radar_prop['transmitter'].waveform_prop['t'])
for s in range(4):
    u = s / 6e4
    print(f"s{s}: my beat mod1 = {(phi(u)-phi(u-tau)) % 1:.8f}  gold phase mod1 = {(np.angle(gold2[0,s])/(2*np.pi)) % 1:.8f}")
