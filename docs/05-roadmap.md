# 05 — Roadmap

## 5.1 Strategy

Bottom-up, continuously verifiable. Every phase ends with a runnable artifact
that passes a defined slice of the existing RadarSimPy test suite — there is
no "big bang" integration at the end, because the interface contract and the
conformance suite exist from day one.

Guiding rules:

- Bind early (Phase 0). The Cython bridge is the cheapest place to discover
  contract misunderstandings.
- Never let a phase end red: the suite subset for each phase is its exit gate.
- CPU first, GPU second, always behind the execution-policy seam (X3).

## 5.2 Phases

### Phase 0 — Skeleton and binding (2–4 weeks)

- CMake project (C++20, options per architecture doc §3.6), repo layout that
  drops in as `src/radarsimcpp`.
- `rsvector`, type aliases, `RadarSimErrorCode`, execution-policy stubs,
  `mem_lib`, license stub.
- All `api/` headers with the exact class/method signatures from
  `radarsimc.pxd`; empty `Run` implementations returning `SUCCESS`.
- **Exit gate**: upstream `build.sh` compiles and links; `import radarsimpy`
  works; simulators run and return zeros; `test_package_api.py` passes.

### Phase 1 — Point, noise, interference (4–6 weeks)

- `rf/` waveform interpolation, antenna patterns, polarization, modulation.
- `PointSimulator`, `NoiseSimulator`, `InterferenceSimulator` (CPU/OpenMP).
- Phase noise generator (SSB mask → time series), seeded determinism.
- **Exit gate**: `test_module_sim_radar_ideal.py`, `test_noise_simulation.py`,
  `test_system_interference.py`, phase-noise and modulation tests pass
  (1e-6-of-peak tolerance).

### Phase 2 — Geometry core (4–6 weeks)

- Mesh container, kinematics (`Move`, rotation conventions), `TargetsManager`.
- BVH (CPU build + traversal), watertight ray–triangle intersection.
- `cp_GetTargetMesh` / scene-state support (`Target.Move`, `vect_mesh_`).
- **Exit gate**: `test_module_mesh_kit.py`, `test_module_scene.py`,
  `test_module_animation_kit.py` pass; BVH validated against brute-force
  intersection on all `models/*.stl`.

### Phase 3 — RCS simulator (6–8 weeks)

- PO far-field integral, polarization, monostatic/bistatic.
- Ray-density-based facet sampling reusing the Phase-2 ray casting.
- **Exit gate**: RCS tests pass; plate and corner-reflector RCS match analytic
  PO solutions across the azimuth sweeps used in the tests.

### Phase 4 — Mesh simulator, SBR (8–12 weeks)

- Pyramid ray generation (density, environment, skip_diffusion culling).
- Multi-bounce SBR with Fresnel/material handling; `ray_filter`.
- PO scattering into baseband; PO LUT with overflow error.
- Fidelity levels (frame/pulse/sample) with range-rate extrapolation.
- `back_propagating`; `log_path` (HDF5); `dry_run`.
- **Exit gate**: `test_module_sim_radar_mesh.py`,
  `test_system_back_propagation.py`, `test_system_range_gate.py`,
  `test_system_cpu_fallback.py` (CPU part) pass at 3e-5-of-peak;
  `capture_reference.py` regression diff clean.

### Phase 5 — LiDAR and polish (2 weeks)

- `LidarSimulator` on the Phase-2 BVH; structured cloud output.
- Free-tier error paths (if keeping the enum contract), error-message parity
  for `TOO_MANY_RAYS_PER_GRID` / `PO_LUT_OVERFLOW`.
- **Exit gate**: `test_module_sim_lidar.py` passes; full `pytest tests/` green
  on CPU.

### Phase 6 — CUDA backend (8–12 weeks)

- LBVH GPU construction (Karras 2012), GPU ray tracing and PO kernels,
  mixed-precision baseband kernel (H-accumulation of phase/delay in float
  kernel, per §3.3 of the architecture doc).
- `gpu_available()` probe, runtime CPU fallback, per-simulator CUDA error
  codes.
- **Exit gate**: GPU-run test suite passes within documented GPU tolerances;
  `bench_sbr.py --compare` within 2× of `benchmarks/baseline/` GPU numbers;
  CPU numbers within 2× of baseline CPU.

### Phase 7 — Extensions (unscheduled, post-parity)

Architecture doc §3.7 items, in priority order: `nogil` regions, in-memory
mesh targets, device index control, batch API. Each is an independent,
opt-in enhancement that must not regress the conformance suite.

## 5.3 Timeline summary

| Phase | Deliverable | Effort |
|-------|-------------|--------|
| 0 | Compilable drop-in skeleton | 2–4 wk |
| 1 | Point/noise/interference | 4–6 wk |
| 2 | Geometry + BVH | 4–6 wk |
| 3 | RCS | 6–8 wk |
| 4 | SBR + PO mesh simulator | 8–12 wk |
| 5 | LiDAR + polish | 2 wk |
| 6 | CUDA | 8–12 wk |
| **Total** | CPU parity ≈ 6–8 eng-months; with GPU ≈ 10–14 eng-months | |

- **M3a** (end Phase 3, as executed): the RCS golden values embed the engine's
  exact PO surface-sampling grid, which is the SBR ray generator's grid.
  Reverse-engineering it from 5 RCS numbers proved underdetermined, so Phase 3
  was reordered after Phase 4: the mesh suite's dense goldens pin the ray
  generator, and RCS then reuses it.

## 5.4 Risk register

| # | Risk | Impact | Mitigation |
|---|------|--------|------------|
| R1 | PO quantization/LUT semantics under-specified by the interface | Golden mesh tests fail at 3e-5 | Treat `capture_reference.py` + mesh tests as the spec; iterate against `cr.stl`/plate scenes first; use `log_path` outputs as oracles |
| R2 | Ray-generation grid details (pyramid layout, culling) not fully documented | Same as R1 | `RADARSIMCPP_ERROR_TOO_MANY_RAYS_PER_GRID` and docs pin the structure; reverse-engineer grid size from error thresholds in tests |
| R3 | Fresnel/polarization sign conventions differ from upstream | RCS/PO phase errors | Coordinate-system doc (pitch sign!) + polarized RCS tests are discriminators |
| R4 | GPU/CPU divergence larger than suite tolerance | Phase 6 gate fails | Match upstream's mixed-precision split exactly; use deterministic reduction where the suite requires it (their `RADARSIMX_DEVICE_PARITY` flag proves it is achievable) |
| R5 | Performance gap on mesh scenes | N6 missed | Profile against `bench_sbr.py` from Phase 4 on; BVH quality and PO-kernel vectorization are the levers |
| R6 | Upstream interface evolves (new `.pxd` symbols) | Drop-in breaks | Pin to a released radarsimpy version (current: 15.4.0, `__init__.py:80`); track upstream releases |
| R7 | License-stub `IsLicensed()=true` diverges from tiered behavior | Free-tier tests fail | Those tests enforce limits in Python (`simulator_radar.pyx:76-118`), so the stub is safe — verify in Phase 0 |

## 5.5 Milestones

- **M1** (end Phase 1): first physically meaningful output — ideal point-target
  baseband matching golden references.
- **M2** (end Phase 3): full far-field EM capability (RCS) — publishable
  validation against analytic solutions.
- **M3** (end Phase 5): **CPU feature parity** — entire suite green; the
  project is usable.
- **M4** (end Phase 6): **performance parity** — GPU build competitive with the
  closed engine.
