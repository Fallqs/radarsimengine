# Point Simulator — Recovered Model and Known Differences

Status of `tools/pysim/reference.py` against `tests/test_module_sim_radar_ideal.py`:
**20 / 21 pass.** The one failure is `test_simc_arbitrary_waveform` (phase only;
amplitudes match). See "Open issue" below.

## Confirmed model (validated against goldens to ≲1e-6 of peak)

Per (frame f, tx m, rx n, pulse p, sample s), channel index
`ch = f·(M·N) + m·N + n`:

- **Timestamps** come from `radar.time_prop["timestamp"]` (absolute, include
  tx channel delay, pulse starts, gate delay, frame starts).
- **Target position** is evaluated at the *reception* timestamp:
  `pos = loc + speed·T` or from the time-varying `[ch, p, s]` arrays.
- **Delay** `τ = (R_tx + R_rx)/c`, per-leg free-space spreading.
- **Amplitude** (matches `src/simpy.py` @4fe7b70):
  `pr_dbm = tx_power + g_tx + g_rx − 10·log10(4π·R_tx²) + rcs_dBsm
            − 10·log10(4π·R_rx²) + 10·log10(λ²/(4π)) + rf_gain`
  `amp = sqrt(2 · 1e-3·10^(pr_dbm/10) · R_load) · 10^(bb_gain/20)`
  with `λ = c / (fc + f_offset[p])`, `fc = (min(f)+max(f))/2`.
  (The per-pulse `f_offset` rides the wavelength — pinned by
  `test_simc_freq_offset`, err < 2.4e-6 → 0.)
- **Phase**: `beat = φ_p(u) − φ_p(u−τ)`, `u` = pulse-local time
  (`T − frame_start − pulse_start[p] − tx_delay[m]`; the tx delay cancels in
  the beat), `φ_p` integrates `f_offset[p] + f(t)` piecewise-quadratically
  (linear f between breakpoints), with **linear-f extrapolation beyond the
  grid ends** (clamping is ruled out by sample 0 of every chirp test).
- **Antenna patterns**: *nearest* table entry, not interpolation
  (az/el independently), plus the channel gain:
  `gain_db = az_pat[argmin|az_angles−az|] + el_pat[argmin|el_angles−el|] + gain`.
  Pinned by `test_simc_tx_az_pattern` (45° query on a ±46° table applies the
  ±46° value exactly).
- **Waveform modulation** (fast-time `amp`/`phs` table): zero-order hold at
  the *echo transmit time*: `idx = floor((u−τ − mod_t[0])/step) + 1`,
  `step = mod_t[1]−mod_t[0]`, out-of-range → 0. Pinned by
  `test_simc_waveform_modulation` (only sample 2 survives, value = var[4])
  and `test_system_pulsed_radar` (21-sample echo window at samples 40–60).
- **Pulse modulation**: per-pulse complex multiplier `pulse_mod[p]`.
- **Polarization**: `|vdot(rx_pol, tx_pol)|`.
- **Platform motion**: `loc(T) = loc0 + speed·T`, channel world positions
  `loc + R(rot0 + rate·T)·ch_loc` with `R = Rz(yaw)·Ry(−pitch)·Rx(roll)`
  (degrees).
- **Interference** (`sim_interference_reference`): one-way Friis with an
  extra `1/(4π)`: `P_r = P_t·G·λ²/((4π)³·R²)`; phase
  `2π·(φ_victim(u_v) − φ_interf(u_i))` (local minus emission, same structure
  as the target beat); hard passband gate `|f_v − f_i| ≤ noise_bandwidth`;
  interferer's modulations apply, victim's do not gate; interferer pulses are
  scanned and rejected outside their `[t0, t_end]` window.

## Open issue: multi-segment (arbitrary) waveform phase

`test_simc_arbitrary_waveform` (5-point non-monotonic `f`) and
`test_system_arbitrary_waveform` (100-point concave sweep, moving target)
fail on **phase only** at interior samples (s0 and the last sample match).
Findings so far (`tools/diag_arbwaveform.py`, `tools/diag_phase1.py`):

- No variant of {left/trapz/right knots} × {linear/Hermite phase eval} ×
  {upsampling 1–500×} × {cubic/PCHIP f} reproduces the goldens.
- DS2 implied effective-frequency error vs the exact piecewise-linear
  integral: −78 kHz at u=5µs rising linearly at ≈1.56e9 Hz/s (zero crossing
  ≈ 55µs), plus a −150 kHz outlier at u=0 (which involves extrapolation).
  I.e. the engine's phase carries a smooth quadratic-ish deviation from the
  exact integral of the given samples.
- DS1 (static target) shows the same *kind* of failure, so it is not a
  motion/Doppler handling difference.

Hypotheses not yet tried: phase table built on a grid tied to `pulse_length`
rather than the breakpoints; a specific float32 staging inside the phase
integrator; engine re-deriving segment slopes by finite differences with a
boundary rule that leaks into interior segments.

The C++ port implements the physically exact model; conformance on these two
tests is tracked as a known divergence until the quirk is identified.
