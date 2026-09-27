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

## Status (2026-09-28)

**Phases 0–4 implemented (point/noise/interference/phase-noise/geometry/LiDAR/
mesh-SBR/RCS); GPU (Phase 6) not started.** The engine builds with CMake and
passes its CTest suites; via `radarsim-engine/tools/build_and_integrate.sh` it
drops into the upstream package and runs the unmodified upstream pytest suite.

Conformance against the full upstream suite: **425 passed / 38 failed /
45 skipped** (45 skips = glTF animation, pygltflib not installed here).
All point-target, noise, pulsed, range-gate, and LiDAR suites pass fully.
The mesh/RCS goldens fail on tolerance: they encode the engine's exact ray
sampler, which is not externally documented at the precision needed — the
implementation is physically validated instead (plate 9%, corner reflector 9%
of benchmark captures; RCS within 0.01 dB of the analytic PO plate value).
See `radarsim-engine/docs/mesh_simulator_model.md`.

Open items: multi-segment waveform phase quirk; phase-noise golden sequence;
mesh/RCS exact-sampler conformance; grazing-rim treatment on curved surfaces;
`back_propagating` return legs; CUDA backend.
