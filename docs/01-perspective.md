# 01 — Perspective

## 1.1 The situation

RadarSimPy is a capable radar simulator whose Python layer, tests, benchmarks,
and documentation are public, but whose compute backend — `radarsimcpp`, a
C++20/CUDA engine — is a private, closed-source submodule distributed only as
prebuilt binaries under a commercial license. The build instructions state this
plainly: without access to `radarsimcpp`, the repository cannot be compiled
(`build_instructions.md:7`).

Anyone who wants to audit the physics, modify the numerics, port the simulator
to a new platform (ARM Linux, ROCm, WebGPU), embed it in a larger system, or
simply guarantee long-term availability of their simulation toolchain currently
cannot do so.

## 1.2 The opportunity

Unusually for a closed engine, its complete specification is public:

- **Interface contract** — `src/radarsimpy/includes/radarsimc.pxd` declares every
  class, method, template parameter, and error code the Python layer uses.
- **Executable acceptance tests** — `tests/` contains 155 tests with golden
  baseband references and documented tolerances (`tests/conftest.py:58-113`).
- **Performance targets** — `benchmarks/baseline/` holds recorded CPU/GPU timing
  sweeps and `.npz` reference captures; the README publishes cross-release
  timings.
- **Physics provenance** — `references/` contains the exact papers the engine
  implements (Gordon 1975, Hwang 2015, Karras 2012, and others).
- **Algorithm semantics** — `gen_docs/user_guide/ray_tracing_simulation.rst`
  documents the ray-tracing model in detail.

Building an independent engine that satisfies this specification is therefore a
well-bounded engineering problem, not a research problem.

## 1.3 Vision

A radar simulation backend that is:

1. **Correct** — passes the RadarSimPy test suite as a drop-in `radarsimcpp`
   replacement, within the tolerances the suite itself defines.
2. **Open** — fully source-available; the physics can be read, cited, audited,
   and taught.
3. **Portable** — CPU everywhere; GPU through an abstraction that is not
   welded to one vendor.
4. **Embeddable** — stateless, re-entrant, and callable from Python without
   holding the GIL, so it can serve as a backend for services, not just
   notebooks.

## 1.4 Positioning

This project does **not** clone RadarSimPy. The Python configuration layer
(`transmitter.py`, `receiver.py`, `radar.py`), the signal-processing toolbox
(`processing.py`), and the mesh/animation utilities are already open and are
reused as-is. The project replaces only the compiled core below the Cython
boundary. Two compatibility postures are possible, in order of ambition:

| Posture | Meaning | Consequence |
|---------|---------|-------------|
| **API-compatible** | Implements the classes in `radarsimc.pxd` with the same semantics; binds through the existing, unmodified Cython layer | Can run the unmodified `radarsimpy` Python package and test suite |
| **Behaviorally comparable** | Same algorithms and outputs through its own API | Free to improve the interface (GIL release, batch API, in-memory meshes, device selection) |

The plan targets API-compatible first (it buys the conformance suite for free),
then layers the improvements of the second posture on top as extensions.

## 1.5 Goals

- G1. Implement all six simulators exposed by the interface: point, mesh
  (SBR + Physical Optics), RCS, LiDAR, interference, noise.
- G2. Pass the existing pytest suite unmodified, within its documented
  tolerances (mesh baseband within 3e-5 of peak; point/interference within
  1e-6 of peak).
- G3. Reach performance within 2× of the published `radarsimcpp` baselines on
  the `benchmarks/bench_sbr.py` scenes, CPU first, GPU second.
- G4. CPU parity on Windows, Linux, and macOS (x86-64 and ARM64); GPU via CUDA
  first, with the execution-policy seam kept thin enough for a second backend.
- G5. Numerical reproducibility: bit-reproducible CPU output run-to-run;
  documented GPU divergence bounds.

## 1.6 Non-goals

- N1. Re-implementing the Python configuration or processing layers.
- N2. Bit-exact match with `radarsimcpp` GPU output (its own GPU/CPU results
  differ by up to ~8.6e-6 of peak; see `tests/conftest.py:66-76`).
- N3. License-manager compatibility. The `LicenseManager` and free-tier
  enforcement are commercial mechanisms; the replacement either omits them or
  provides a trivial always-licensed stub so the binding compiles.
- N4. Full-wave EM solvers (MoM/FDTD/FEM). The engine is a high-frequency
  asymptotic simulator (ray tracing + PO), matching the original's regime.
- N5. Competing on features the interface does not expose (e.g. atmospheric
  propagation, clutter statistics) until after parity.

## 1.7 Success criteria

The project is successful when, on a clean machine:

1. `git clone` → `./build.sh --arch=cpu` produces a `libradarsim*` binary from
   open source only.
2. Dropping it into the unmodified RadarSimPy build makes `pytest tests/` pass
   with the same skip patterns as the genuine engine.
3. `python benchmarks/bench_sbr.py --compare` shows results within the G3
   envelope on the baseline scenes.
