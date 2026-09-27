# 02 — Requirements

Derived from `src/radarsimpy/includes/radarsimc.pxd` (the binding contract),
`src/radarsimpy/simulator_*.pyx` and `src/radarsimpy/lib/cp_radarsimc_*.pyx`
(usage semantics), and `gen_docs/user_guide/` (documented behavior).

## 2.1 Template and type requirements

| ID | Requirement | Source |
|----|-------------|--------|
| T1 | All scene/simulator classes are templates over `<H, L>` (high/low precision) or `<T>`; the Python build instantiates `H=double`, `L=float` | `radarsimc.pxd:14-17`, `type_def.pxd:62` |
| T2 | `L` is simultaneously the geometry/BVH type and the mesh baseband kernel type | `type_def.pxd:52-62` |
| T3 | Range, delay, range gate, and waveform phase are evaluated in `H` even inside `L`-precision kernels | `type_def.pxd:53-57` |
| T4 | `Receiver` takes `gate_delay` as `double` regardless of `T` | `radarsimc.pxd:211-219` |
| T5 | Custom `Vec2[T]`/`Vec3[T]` value types with raw-pointer and per-element constructors | `rsvector.pxd` |

## 2.2 Core objects

| ID | Requirement | Source |
|----|-------------|--------|
| C1 | `Transmitter[H,L]`: waveform as `(freq, freq_time)` samples, per-pulse `freq_offset` and `pulse_start_time`; optional SSB phase-noise constructor (mask + fs + sample count + seed + validation flag) | `radarsimc.pxd:165-187` |
| C2 | `Transmitter.AddChannel`: location, complex polarization vector, az/el pattern tables, antenna gain, fast-time complex modulation `(mod_t, mod_var)`, per-pulse complex `pulse_mod`, channel delay, time-grid resolution | `radarsimc.pxd:189-200` |
| C3 | `Receiver[T]`: fs, RF gain, load resistor, baseband gain, baseband bandwidth, gate delay; `AddChannel` with pattern tables and polarization | `radarsimc.pxd:206-227` |
| C4 | `Radar[H,L]`: Tx+Rx, frame start times, platform location/rotation arrays (scalar or per-timestamp), speed and rotation rate; `InitBaseband(H*, H*)` hands caller-owned buffers to the engine; `SyncBaseband()` for device memory; exposes `sample_size_` | `radarsimc.pxd:233-252` |
| C5 | `PointsManager[T]`: `AddPoint` (time-varying location/rcs/phase arrays) and `AddPointSimple` (static) | `radarsimc.pxd:258-273` |
| C6 | `TargetsManager[T]`: `AddTarget`/`AddTargetSimple` from raw `(points, cells)` arrays with origin, kinematics, complex permittivity/permeability, `skip_diffusion`, per-target ray density, `environment` flag | `radarsimc.pxd:320-350` |
| C7 | `Target[T]`: `Move(index, time)` kinematics evaluation; exposes `vect_mesh_`, `array_size_` for scene-state queries | `radarsimc.pxd:282-313` |
| C8 | `Ray[T,L]`: direction/location/normal/range arrays + reflection count | `radarsimc.pxd:356-363` |

## 2.3 Simulators (functional requirements)

| ID | Simulator | Requirement | Source |
|----|-----------|-------------|--------|
| S1 | `PointSimulator[H,L,P]` | `Run(radar, points)` — ideal point-target baseband accumulation | `radarsimc.pxd:372-377` |
| S2 | `MeshSimulator[H,L,P]` | `Run(radar, targets, level, density, ray_filter, back_propagating, log_path, dry_run)`; `level` 0/1/2 = re-trace per frame/pulse/sample; `density` = rays per wavelength; `ray_filter` = bounce-count window (also caps trace depth); `back_propagating` = returns reflected back out along arrival surfaces; `log_path` = HDF5 debug dump; `dry_run` = validate without tracing | `radarsimc.pxd:382-393`, `ray_tracing_simulation.rst`, `simulator_radar.pyx:154` |
| S3 | `RcsSimulator[T,P,L]` | `Run(targets, inc_dirs, obs_dirs, inc_pol, obs_pol, frequency, density)` + `GetRcs()`; monostatic and bistatic; complex polarization vectors | `radarsimc.pxd:397-409` |
| S4 | `LidarSimulator[T,P]` | `Run(targets, phi, theta, position)`; fills `cloud_` with `Ray` entries; Lambertian intensity | `radarsimc.pxd:413-421`, `simulator_lidar.pyx:230-239` |
| S5 | `InterferenceSimulator[H,L,P]` | `Run(victim_radar, interfering_radar)` — writes interference into the victim's (re-pointed) baseband buffers | `radarsimc.pxd:426-431`, `simulator_radar.pyx:427` |
| S6 | `NoiseSimulator[H,L,P]` | `Run(radar, noise_level, is_complex, timestamps, sizes, out_real, out_imag, seed)`; channel/pulse-correlated thermal noise | `radarsimc.pxd:435-448` |

Execution order inside one `sim_radar` call is fixed: point → mesh → (GPU sync)
→ interference (optional) → noise, all accumulating into the same buffers
(`simulator_radar.pyx:392-435`).

## 2.4 Error model

E1. All `Run` methods return `RadarSimErrorCode` — never throw across the
binding boundary. The enum includes generic codes plus per-simulator CUDA
sync/kernel codes in 100-blocks (`radarsimc.pxd:101-134`).

E2. Two codes have bespoke Python messages and must be producible:
`TOO_MANY_RAYS_PER_GRID` and `PO_LUT_OVERFLOW` (`simulator_radar.pyx:120-148`).

## 2.5 Execution and platform requirements

| ID | Requirement |
|----|-------------|
| X1 | Template-based execution policies `cpu_policy`/`gpu_policy` with static `is_gpu/is_cpu/name()/device_id()`; global `cpu`/`gpu` instances |
| X2 | Runtime `gpu_available()` probe, cached for process lifetime, honoring `CUDA_VISIBLE_DEVICES`; always false in CPU-only builds |
| X3 | A CUDA build must compile both policies and fall back to CPU at runtime when no device is usable (CI explicitly tests this: `unit_test_ubuntu.yml:22-29`) |
| X4 | CPU parallelism via OpenMP (all cores by default, `OMP_NUM_THREADS` to limit) |
| X5 | Platforms: Windows (MSVC), Linux (GCC), macOS (Clang, x86-64 + ARM64) |
| X6 | HDF5 output for `log_path` (debug ray dumps) — consumed as a prebuilt static lib upstream |

## 2.6 Non-functional requirements

| ID | Requirement | Measure |
|----|-------------|---------|
| N1 | CPU output bit-reproducible run-to-run with fixed seeds | rerun diff |
| N2 | Mesh baseband within 3e-5 of golden-reference peak | `tests/conftest.py:80-89` |
| N3 | Point/interference baseband within 1e-6 of peak | `tests/conftest.py` |
| N4 | Noise deterministic: `NoiseSimulator` is always invoked with seed 0 by the binding | `simulator_radar.pyx:403,418` |
| N5 | Stateless across calls: no cross-call engine state; every `Run` is independent | code review + interleaving test |
| N6 | Performance within 2× of `benchmarks/baseline/` on the bench_sbr scenes | `bench_sbr.py --compare` |
| N7 | Zero-copy interop: engine writes into caller-owned baseband buffers (`InitBaseband`) | `radarsimc.pxd:247-249` |

## 2.7 Out of scope (from perspective doc)

License management, bit-exact GPU parity, full-wave solvers, features beyond
the exposed interface.
