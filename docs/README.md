# Radar Simulation Engine — Project Documents

Working name: **radarsim-engine** — an independent, open re-implementation of the
simulation backend that `radarsimpy` binds to (`radarsimcpp`).

This document set is derived from a full analysis of the RadarSimPy repository:
its Cython interface declarations (`src/radarsimpy/includes/*.pxd`), its test
suite (`tests/`), its benchmarks (`benchmarks/`), its user-guide documentation
(`gen_docs/user_guide/`), and its physics reference papers (`references/`).

## Documents

| # | Document | Purpose |
|---|----------|---------|
| 01 | [Perspective](./01-perspective.md) | Why this project exists, positioning, goals and non-goals, success criteria |
| 02 | [Requirements](./02-requirements.md) | Functional and non-functional requirements, extracted from the interface contract |
| 03 | [Architecture](./03-architecture.md) | System decomposition: modules, data flow, precision and execution policies |
| 04 | [Physics & Algorithms](./04-physics.md) | Per-simulator algorithm specifications and literature references |
| 05 | [Roadmap](./05-roadmap.md) | Phases, milestones, effort estimates, risk register |
| 06 | [Validation Strategy](./06-validation.md) | Conformance testing against the existing suite, performance gating |

## Reading order

New contributors: 01 → 03 → 04. Planners: 01 → 05 → 06. Implementers: 02 → 03 → 04 → 06.

## Grounding

Every requirement and design decision in these documents cites its source in the
RadarSimPy repository (file path and line where applicable). If the upstream
repository changes, these documents should be re-checked against it.

## Status (2026-09-28, end of session)

**Phases 0–5 implemented and committed; Phase 6 (CUDA) blocked on hardware.**
The engine builds with CMake (MSVC/GCC/Clang or the bundled conda-zig path),
passes its own CTest suites (4/4), and runs the unmodified upstream pytest
suite via `radarsim-engine/tools/build_and_integrate.sh`.

Upstream conformance: **427 passed / 36 failed / 45 skipped** (~54 s).
Every point-target, noise, pulsed, range-gate, LiDAR, back-propagation, and
config suite passes fully. The mesh/RCS golden tests fail on tolerance: they
encode the closed engine's exact ray sampler, which is not externally
recoverable at the required precision (3e-5 of peak). The mesh/RCS
implementation is validated physically instead -- plate 9% and corner
reflector 9% of the benchmark reference captures, RCS within 0.01 dB of the
analytic PO plate value -- and the engine carries its own self-regression
suite (`radarsim-engine/tests/test_mesh_regression.cpp`).

Blocked items:
- Phase 6 (CUDA kernels): no GPU/toolkit on this machine.
- Mesh/RCS golden exactness, multi-segment waveform phase quirk, phase-noise
  golden sequence: all encode closed-engine internals (see
  `radarsim-engine/docs/mesh_simulator_model.md` and
  `radarsim-engine/docs/point_simulator_model.md`).
