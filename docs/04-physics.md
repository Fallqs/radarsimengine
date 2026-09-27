# 04 — Physics & Algorithms

Per-simulator algorithm specifications. Citations refer to papers in
`references/` and to the upstream documentation that fixes the semantics.

## 4.1 Common: signal model

Baseband sample for Tx channel `m`, Rx channel `n`, pulse `p`, sample `s`:

```
bb[m,n,p,s] = Σ_paths A_path · M_tx(t) · M_pulse(p) · exp(j·φ_path)
φ_path = 2π · [ f(t − τ)·(t − τ) − f(t)·t ]  terms from the sampled waveform
τ      = R_path / c
```

- Waveform `f(t)` is user-sampled (arbitrary: FMCW, PMCW, pulsed, CW) and
  interpolated in double precision (`transmitter.hpp` stores `freq`,
  `freq_time`, per-pulse `freq_offset`, `pulse_start_time`).
- Doppler sign convention follows `gen_docs/user_guide/doppler_convention.rst`
  — implement to the convention doc, not to intuition, and verify against
  `test_module_sim_radar_ideal.py`.
- Antenna gains: separable azimuth/elevation pattern tables, interpolated and
  converted dB→linear; complex polarization 3-vectors project via dot product
  at both Tx and Rx. Pattern normalization (max of azimuth pattern) is done in
  Python (`transmitter.py:510-528`) — the engine receives normalized tables.
- Range gate: only delays within `[gate_delay, gate_delay + pulse_length)`
  contribute; the gate also shifts phase reference (`receiver.hpp` gate_delay,
  `test_system_range_gate.py`).

## 4.2 Point simulator

For each point target with time-varying `location(t)`, `rcs(t)`, `phase(t)`:

1. Interpolate point state onto each timestamp of the `[channels, pulses,
   samples]` grid (linear).
2. Path: Tx_m → point → Rx_n; range = |Tx−P| + |P−Rx|.
3. Amplitude from the radar equation with RCS in dBsm, both antenna patterns,
   polarization projection, and the explicit phase offset.
4. Accumulate into baseband.

This simulator is the numerical ground truth: it has no geometric
approximation, so its golden references carry the tightest tolerance (1e-6 of
peak) and it must be validated first.

## 4.3 Noise simulator

- Complex baseband: independent Gaussian I/Q at the given RMS `noise_level`;
  real baseband: real Gaussian. (`bb_type` shapes `noise_bandwidth` upstream;
  the engine gets `is_complex` and the level.)
- Seeded PRNG; the binding always passes seed 0 → output must be a pure
  function of (shape, level, seed) (`simulator_radar.pyx:403,418`).
- Correlation structure across channels/pulses must match the reference
  captures in `test_noise_simulation.py`.

## 4.4 Phase noise

- Input: SSB phase-noise mask `(pn_f, pn_power)` in dBc/Hz, target sample rate,
  sample count, seed.
- Method: shape white Gaussian noise in the frequency domain with the mask,
  IFFT to a time-domain phase deviation; reference implementation readable in
  `radar.py:121` (`cal_phase_noise`) and `references/add_phase_noise.m`.
- The second `Transmitter` constructor (`radarsimc.pxd:176-187`) implies
  per-frame deferred generation — generate once per frame and reuse across
  channels with the correct correlation behavior; pin semantics to
  `test_config_transmitter.py` / phase-noise system tests.

## 4.5 Mesh simulator — SBR + Physical Optics

The core of the engine. Semantics fixed by
`gen_docs/user_guide/ray_tracing_simulation.rst`; physics by Gordon 1975 and
Hwang 2015 (`references/`).

### 4.5.1 Pipeline (per ray-tracing pass)

```
1. Move all targets to the pass instant (kinematics; §4.5.5)
2. Build/rebuild BVH over the moved scene
3. Pyramid ray generation from each Tx (§4.5.2)
4. SBR trace: bounce rays through the scene (§4.5.3)
5. At each contributing bounce: PO scattering back to every Rx (§4.5.4)
6. Accumulate baseband over passes
```

### 4.5.2 Ray generation

- Rays launched from the Tx through a grid over the scene's subtended solid
  angle ("pyramid"). `density` = rays per wavelength: at 77 GHz (λ ≈ 3.9 mm),
  density 1.0 is one ray per 3.9 mm of grid spacing.
- Directions that see only `skip_diffusion` surfaces are culled *before*
  launch — the largest single saving on ground-plane scenes
  (`ray_tracing_simulation.rst:224-231`).
- `environment=True` targets claim a coarser share of the ray budget; the
  target behind them sets the sampling (`ray_tracing_simulation.rst:272-299`).
- Grid capacity overflow → `TOO_MANY_RAYS_PER_GRID`.

### 4.5.3 Shooting and bouncing rays

- Specular reflection at each triangle hit; direction via the reflection law,
  amplitude via Fresnel coefficients for complex εr/μr with polarization
  decomposition (`references/Reflection and Refraction of Plane EM Waves.pdf`).
- Bounce count capped by `ray_filter[1]`; rays outside `[min, max]` are
  excluded from baseband. A ray stopped by the cap is treated as still
  travelling (no outgoing contributions) (`ray_tracing_simulation.rst:146-149`).
- `skip_diffusion` surfaces redirect rays but contribute no return; a ray
  touching only skipped surfaces contributes nothing
  (`ray_tracing_simulation.rst:204-218`).

### 4.5.4 Physical Optics scattering

- Each ray landing = one surface-current sample of the PO integral:
  `J = 2(n̂ × H_i)` on the lit side; the scattered field at the Rx follows
  the far-field Kirchhoff–Helmholtz / Gordon formulation
  (`references/Gordon_1975_*.pdf`, `references/Hwang_2015_*.pdf`).
- Facet size is set by landing density (rays per wavelength) — ray density is
  surface sampling density (`ray_tracing_simulation.rst:302-321`).
- Scattering points are quantized into a bounded PO lookup table; overflow →
  `PO_LUT_OVERFLOW`. Size the LUT from (facet count, density); make the bound
  a named constant and document it.
- Field → baseband: same accumulation kernel as §4.1, with the PO complex
  amplitude in place of the point-target radar-equation gain.

### 4.5.5 Fidelity levels and motion

- `level`: 0 = one pass per frame, 1 = per pulse, 2 = per ADC sample; cost
  tracks pass count (frame of P pulses × S samples = 1, P, or P·S passes).
- Between passes each target is carried forward at its **range rate** —
  straight-line extrapolation, exact only for constant-velocity translation
  (`ray_tracing_simulation.rst:78-91`). Rotating/vibrating targets need
  `level=sample` (micro-Doppler).

### 4.5.6 Back-propagation

- When enabled, a ray that **escapes** the scene also contributes returns that
  reflect their way back out over the surfaces it arrived over.
- Each return reflection must land on the surface (not past its edge) and be
  unoccluded. Chains are never built for rays stopped by the depth cap
  (`ray_tracing_simulation.rst:151-178`).
- Acceptance: `test_system_back_propagation.py`; invariant: a single convex
  target returns identical baseband with the flag on or off.

## 4.6 RCS simulator

- PO far-field integral over the mesh as seen from `inc_dir`; observation
  toward `obs_dir`; complex polarization vectors on both ends
  (`radarsimc.pxd:397-409`).
- Monostatic = obs_dir −inc_dir; bistatic supported natively.
- Analytical validation targets: flat plate (PO closed form), trihedral corner
  reflector (`models/plate.stl`, `models/cr.stl`, plus RCS tests).

## 4.7 LiDAR simulator

- Ray fan over `(phi, theta)` from the sensor position; closest hit against
  the BVH; per-hit record position, direction, normal, range; intensity ∝
  cos(incidence)/r² (Lambertian) (`simulator_lidar.pyx:230-239`).
- Kinematics evaluated once at `frame_time` (no time-varying motion).
- Trivial once §4.5's BVH exists — schedule it immediately after.

## 4.8 Interference simulator

- Victim radar and interfering radar share no scene; the interferer's
  waveform (its own `f(t)`, offsets, PRP, modulation, Tx power, patterns) is
  evaluated inside the victim's receive window.
- Writes into the victim's baseband buffers, which the binding re-points to
  separate interference buffers before the call (`simulator_radar.pyx:427`).
- Path loss interferer→victim by free-space range; both radars may be moving
  (platform kinematics apply).

## 4.9 Numerical validation anchors

| Check | Reference | Tolerance |
|-------|-----------|-----------|
| Point-target baseband vs golden | `test_module_sim_radar_ideal.py` | 1e-6 of peak |
| Mesh baseband vs golden (CPU) | `test_module_sim_radar_mesh.py` | 3e-5 of peak |
| Plate/CR RCS vs analytic | RCS test files + PO plate formula | per-test |
| CPU↔GPU consistency | `tests/conftest.py:66-76` | ~8.6e-6 of peak (mesh) |
| Regression scenes | `benchmarks/capture_reference.py` (6 scenes) | diff vs `baseline/*.npz` |
