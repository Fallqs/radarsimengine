# 03 — Architecture & Decomposition

## 3.1 Layer view

```
┌─────────────────────────────────────────────────────────────┐
│ radarsimpy (unchanged, upstream)                            │
│   radar.py / transmitter.py / receiver.py  →  dicts of numpy │
│   simulator_*.pyx, lib/cp_radarsimc_*.pyx  →  Cython bridge  │
└──────────────────────────┬──────────────────────────────────┘
                           │ C++ ABI declared in includes/*.pxd
┌──────────────────────────▼──────────────────────────────────┐
│ radarsim-engine (this project)                              │
│                                                             │
│  api/        classes named in radarsimc.pxd — the only     │
│              headers the bridge includes                    │
│  core/       types, Vec2/Vec3, error codes, exec policies   │
│  geom/       mesh, kinematics, BVH, intersection            │
│  rf/         waveform eval, antenna patterns, propagation   │
│  sim/        point / mesh / rcs / lidar / interf / noise    │
│  libs/       mem copy, hdf5 logging, (license stub)         │
└─────────────────────────────────────────────────────────────┘
```

The binding includes headers by these exact relative paths
(`radarsimc.pxd:54-453`), so the public header layout is itself a contract:

```
core/execution_policy.hpp   core/enums.hpp
libs/mem_lib.hpp            libs/license_manager.hpp   libs/motion_lib.hpp
transmitter.hpp             receiver.hpp               radar.hpp
points_manager.hpp          target.hpp                 targets_manager.hpp
triangle.hpp                ray.hpp
simulator_point.hpp         simulator_mesh.hpp         simulator_rcs.hpp
simulator_lidar.hpp         simulator_interference.hpp simulator_noise.hpp
rsvector.hpp (include dir: includes/rsvector/)
```

## 3.2 Module decomposition

### 3.2.1 `core/` — foundations

- `types.hpp`: `int_t`, `uint_t`, and the `H`/`L` template discipline
  (req T1–T4).
- `execution_policy.hpp`: `cpu_policy` / `gpu_policy` tags, global instances,
  `gpu_available()` probe (cached, honors `CUDA_VISIBLE_DEVICES`). In a CPU-only
  build `gpu_policy` aliases the CPU implementation so templated code still
  compiles — mirroring the upstream trick (`radarsimc.pxd:85-95`).
- `enums.hpp`: `RadarSimErrorCode`, including the per-simulator CUDA 100-blocks.
- No third-party dependencies.

### 3.2.2 `geom/` — scene geometry

- `Triangle[T]`, `Mesh` (points + cells, flat arrays).
- `Kinematics`: location/rotation either constant or time-indexed arrays;
  `Move(index, time)` evaluates position/orientation; straight-line
  range-rate extrapolation between ray-tracing passes
  (`ray_tracing_simulation.rst:78-91`).
- `Rotation`: `Rz(yaw)·Ry(−pitch)·Rx(roll)` composition, degrees; the pitch
  sign convention is load-bearing (README Coordinate Systems, note at line 180);
  `libs/motion_lib.hpp::Rotate` is part of the public API.
- `BVH`: binary BVH over triangles; CPU build first, LBVH (Karras 2012,
  `references/`) for GPU. Möller–Trumbore or the watertight ray–triangle
  variant; ray–box per the Williams et al. paper in `references/`.
- `Target[T]`/`TargetsManager[T]`: material (complex εr, μr), `skip_diffusion`,
  `environment`, per-target density; target-level mesh caching so repeated
  passes within one `Run` do not rebuild.

### 3.2.3 `rf/` — radar front-end model

- `Waveform`: arbitrary `f(t)` sampled waveform; interpolation in `H`;
  per-pulse frequency offset and start time; bandwidth/pulse-length derived
  upstream in Python.
- `PhaseNoise`: SSB mask → time-domain phase noise. Reference algorithm:
  `radar.py:121` (`cal_phase_noise`), with per-frame deferred generation
  implied by the second Transmitter constructor.
- `AntennaPattern`: separable az/el tables, 1-D interpolation, dB → linear;
  polarization as complex 3-vector with dot-product projection.
- `Modulation`: fast-time complex modulation interpolated on `(mod_t, mod_var)`
  with time-grid resolution `grid`; per-pulse complex `pulse_mod`.
- `Receiver`: noise budget is computed in Python; engine needs fs, gains,
  baseband bandwidth, gate delay, and the output buffer convention.
- `Radar`: owns baseband pointers (`InitBaseband`), timestamp semantics:
  channel index = `frame·(M·N) + tx·N + rx` (`simulator_radar.pyx:240-252`).

### 3.2.4 `sim/` — the six simulators

Shared kernel: **baseband accumulation**. For a scatterer at range path `R`
seen by Tx `m`, Rx `n` at time `t`:

```
delay  = R / c                       (evaluated in H)
phase  = 2π · f(delay, t−delay)      (waveform in H)
gain   = radar-equation(RCS or PO E-field, patterns, polarization)
bb[n, p, s] += gain · mod(t) · exp(j·phase)
```

- `point`: apply per point target; interpolation of time-varying
  location/rcs/phase onto the timestamp grid.
- `mesh`: SBR + PO (see 04-physics). Cost driver; the only simulator whose
  inner loop justifies mixed precision.
- `rcs`: PO far-field integral over illuminated facets; polarization vectors;
  monostatic shortcut when `obs == inc`.
- `lidar`: ray fan `(phi, theta)` from `position`, closest-hit into `cloud_`.
- `interference`: evaluate the interferer's waveform inside the victim's time
  gate, including frequency-offset overlap; writes into the victim's buffers.
- `noise`: Gaussian complex noise at `noise_level`, shaped by baseband type,
  seeded; deterministic across calls (binding always passes seed 0).

### 3.2.5 `libs/`

- `mem_lib.hpp`: `Mem_Copy`, `Mem_Copy_Complex`, `Mem_Copy_Vec3` — raw pointer
  → `std::vector` helpers. Trivial but part of the contract.
- `hdf5_log.hpp`: `log_path` ray dump (HDF5, static-linked like upstream).
- `license_manager.hpp`: stub implementing `GetInstance/SetLicense/IsLicensed/
  IsFreeTier/GetLicenseInfo` so `license.pyx` compiles unmodified;
  `IsLicensed()=true`, `IsFreeTier()=false` (non-goal N3).

## 3.3 Precision policy (requirement T1–T4)

| Quantity | Precision | Rationale |
|----------|-----------|-----------|
| Range, delay, gate, waveform phase, timestamps | `H` (double) | float32 mantissa resolves ~5.8e-11 s at a 741 µs gate ≈ 0.5 carrier cycles at 9 GHz (`radarsimc.pxd:211-213`) |
| Geometry, BVH, patterns, PO kernel accumulation | `L` (float) | ~4× GPU throughput; ≤ −113 dBc error vs FP64 on range-Doppler (`type_def.pxd:52-62`) |

Changing `L` changes which rays are launched (grid quantization), so it is a
compile-time typedef, not a runtime switch. Provide `RSIM_LOW_PRECISION=float|
double` CMake option, default `float`, mirroring `type_def.pxd`.

## 3.4 Execution policy & parallelism

```
simulator<H, L, Policy>  — one implementation, two instantiations
  Policy = cpu_policy  → OpenMP over channels × pulses (or rays)
  Policy = gpu_policy  → CUDA kernels; async launch + sync, error mapped to
                         the per-simulator 100-block codes
```

- Ray tracing parallelizes over launched rays; baseband accumulation uses
  per-thread private buffers + reduction (CPU) or atomics (GPU, the known
  source of GPU/CPU divergence — document bounds, do not chase bit parity).
- `OMP_NUM_THREADS` honored at process start (OpenMP runtime reads it once).
- GPU build without device → CPU fallback at `Run` time (requirement X3).

## 3.5 Memory & data flow

```
Python (numpy, float64)                    Engine
─────────────────────────────────────────────────────────────
bb_real, bb_imag  ──InitBaseband(H*,H*)──► accumulates in place
targets dicts     ──AddTarget(points,    ► copies into Target
                     cells, ...)             (engine-owned)
timestamp arrays  ──AddPoint / Radar     ► read-only views during Run
noise buffers     ──NoiseSimulator.Run   ► writes in place
```

- The engine copies mesh data on `AddTarget` (Python may free after).
- Baseband is caller-owned; the engine never allocates the output.
- One `Run` call allocates its scratch (rays, BVH instance per pass, PO LUT)
  and frees at return — statelessness (N5).
- PO LUT: scattering points quantized into a bounded table; overflow is a
  first-class error (`PO_LUT_OVERFLOW`), sized by a CMake/runtime parameter.

## 3.6 Build system

- CMake ≥ 3.18, C++20; options: `GPU_BUILD`, `GTEST`, `ENABLE_LICENSE` (stub
  only), `RSIM_LOW_PRECISION`, `RADARSIMX_CUDA_ARCHITECTURES` passthrough.
- Output: `libradarsimcpp.{so,dylib,dll}` — **same library name** so the
  upstream `setup.py` links unmodified (`setup.py:319` links `radarsimcpp`).
- Header include dirs must match `setup.py:274-276`
  (`includes/`, `includes/rsvector/`).
- Layout the repo as a drop-in: checking this project out *as*
  `src/radarsimcpp` must make the upstream `build.sh`/`build.bat` work.

## 3.7 Extension seams (post-parity)

Designed-in, but only built after G2/G3 are met:

1. `nogil` regions so Python threads overlap (upstream holds the GIL for the
   whole call).
2. In-memory mesh targets (skip file re-read per call — the service-caching
   hook; upstream always reloads).
3. Device index + memory budget on the GPU policy.
4. Batch `Run` over target sets sharing one BVH.
