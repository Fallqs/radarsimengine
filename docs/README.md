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

**Phases 0–2 complete, plus LiDAR.** The engine builds with CMake and passes its
CTest suites; via `radarsim-engine/tools/build_and_integrate.sh` it drops into
the upstream package and runs the unmodified upstream pytest suite.

Conformance against the full upstream suite: **420 passed / 43 failed /
45 skipped** (45 skips = glTF animation, pygltflib not installed here).
Failures by cause:

| Cause | Tests | Status |
|-------|-------|--------|
| Mesh (SBR+PO) simulator not yet implemented | 34 | Phase 4 |
| RCS simulator not yet implemented | 2 | Phase 3 |
| Multi-segment waveform phase quirk | 2 | open (see `radarsim-engine/docs/point_simulator_model.md`) |
| Phase-noise golden (exact-sequence) | 1 | open |
| Interference f32 rounding profile (60 GHz real mode) | 1 | open — matches to ~1e-5 rad |
| `test_scene_interference` (mesh + interference) | 1 | Phase 4 |
| `test_sim_cw_raytracing` (mesh) | 1 | Phase 4 |
| range-gate mesh tests | 2 | Phase 4 |

All point-target, noise, pulsed, range-gate, and LiDAR suites pass fully.
