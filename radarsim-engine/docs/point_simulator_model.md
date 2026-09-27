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

## Additional confirmed semantics (round 2)

- **Range gate = stretch processing**: the deramp reference is the transmit
  chirp delayed by `gate_delay`, so `beat = φ(u − gate_delay) − φ(u − τ)`;
  a target at the gate range beats at DC (`test_target_at_gate_produces_dc`).
  Gate delay also applies to the victim LO in the interference path.
  Range-gate suite: **17/17**.
- **Modulation table is PERIODIC**: `idx = floor((u−τ − mod_t[0])/step) + 1`
  taken with floor-modulo over the table length — folded (beyond-PRP) echoes
  wrap into the table instead of being zeroed or reading OOB
  (`test_pulsed_radar_ambiguous_range_no_ghost`; upstream had a C++ `%`
  negative-index bug here). Pulsed suite: **4/4**.
- **Real baseband** (`bb_type="real"`): baseband/interference/noise outputs
  are the real part of the complex computation.
- **float32 geometry staging**: platform rotation/location and channel/target
  positions are narrowed to float32 at the marshalling boundary (they are
  `float_t` in the engine). The rotation matrix is built from the float32
  angles — `sin(f32(π)) = −8.74e-8`, which shifts a rotated channel by
  ~6e-7 m and shows up as a +7.8e-4 rad phase offset at 60 GHz in
  `test_system_interference`. With this staging the interference amplitude
  and phase match to 6.7e-7 / ~1e-5 rad respectively.
- **Noise model**: per-sample Gaussian keyed by (seed=0, rx index, float64
  timestamp bits) — channels sharing an Rx and identical timestamps get
  identical noise; complex mode = level·(n_re + j·n_im)/√2; real mode =
  level·n_re. Noise suite: **16/16**.
- **Phase noise**: mask shaping + IFFT generator is bit-exact against
  `radar.cal_phase_noise` (validation mode is deterministic: AWGN ≡
  (1+1j)/√2). Application: `bb ×= exp(j(φ_pn(u) − φ_pn(u−τ)))` with linear
  interpolation between samples.

## Open issues

1. **Multi-segment waveform phase** (`test_simc_arbitrary_waveform`,
   `test_system_arbitrary_waveform`): interior-sample phases of
   non-two-point waveforms carry a smooth deviation from the exact
   piecewise-linear integral (see below).
2. **Phase-noise golden** (`test_fmcw_phase_noise`): the range-profile-delta
   golden is hypersensitive at spectral nulls; my generator is bit-exact to
   the Python reference algorithm, so the C++ generator differs in some
   unobservable detail (candidate: normalization `M` vs `2M−2`, sign, or a
   different interpolation). Structure and magnitude are right; exact
   sequence match open.
3. **`test_system_interference`**: structure, window, amplitude all match;
   residual ~1e-5 rad/sample vs the engine's exact float32 rounding profile
   in the 60 GHz real-mode scenario (complex 24 GHz ideal interference test
   passes exactly).

The C++ port implements the physically exact model; conformance on these
tests is tracked as known divergence until the quirks are identified.

