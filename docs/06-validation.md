# 06 — Validation Strategy

## 6.1 Principle

The RadarSimPy repository already contains a complete, executable definition
of "correct" and "fast enough". Validation is therefore not a phase at the end
of the project — it is the steering mechanism for every phase. Three
independent oracles are used:

1. **Conformance** — the upstream pytest suite (behavioral correctness).
2. **Analytic solutions** — closed-form physics where they exist (absolute
   correctness, independent of upstream).
3. **Benchmarks** — `benchmarks/` timing and reference captures (performance
   and regression).

## 6.2 Conformance harness

### 6.2.1 Setup

The engine repository is laid out to be checked out *as* `src/radarsimcpp`
inside a RadarSimPy clone. CI then runs the unmodified upstream flow:

```
git clone <radarsimpy> && cd radarsimpy
git clone <radarsim-engine> src/radarsimcpp
./build.sh --arch=cpu --test=on --license=off --deps=release
pytest tests/
```

Because `setup.py` links `-lradarsimcpp` by name (`setup.py:319`) and includes
`src/radarsimcpp/includes` (`setup.py:274-276`), keeping the library name and
header layout identical makes the swap transparent.

### 6.2.2 Tolerance policy (from `tests/conftest.py`)

| Suite class | Metric | Tolerance | Origin of reference |
|-------------|--------|-----------|---------------------|
| Point / interference | peak-relative abs error | 1e-6 | recorded from upstream GPU run |
| Mesh | peak-relative abs error | 3e-5 | recorded from upstream CPU run |
| CPU↔GPU (mesh) | peak-relative | ~8.6e-6 | documented device-parity bound |

References are golden arrays baked into the tests; `assert_baseband_close`
normalizes by the peak of the expected array. We do not re-record references
to fit our engine — the recorded values are the target.

### 6.2.3 Phase-gate mapping

Each roadmap phase has an explicit suite subset as its exit gate
(`docs/05-roadmap.md` §5.2). Full-suite green is the M3 milestone.

### 6.2.4 Marks and optional deps

`pytest.ini` defines `mesh`, `gltf`, `slow` marks; `conftest.py:145` skips
marked tests when the optional library is absent. CI must run the full matrix
(trimesh + pygltflib installed) so nothing is silently skipped.

## 6.3 Analytic validation (independent of upstream)

| Target | Closed form | Test |
|--------|-------------|------|
| Flat plate RCS (broadside sweep) | PO plate: σ = 4πA²/λ²·sinc² terms | `models/plate.stl`, `plate5x5.stl` |
| Trihedral corner reflector | σ_max = 12πa⁴/λ² (per-face a) | `models/cr.stl` |
| Point target range/Doppler | radar equation + Doppler convention doc | synthetic, cross-checked with `processing.py` FFT peaks |
| Sphere (far field, PO limit) | σ → πr² | `models/ball_1m.stl` |
| Free-space path loss | Friis | interference tests |

These guard against the failure mode of matching upstream's numbers by
matching upstream's bugs: agreement with both the golden references *and* the
analytics is required where both exist.

## 6.4 Determinism and device parity

- **CPU reproducibility**: two runs of the same configuration must be
  bitwise identical (fixed seeds; no atomics on the CPU path). Test: rerun
  diff = 0.
- **Seeded noise**: `NoiseSimulator` output is a pure function of
  (shape, level, seed); the binding passes seed 0 always.
- **GPU divergence**: bounded and documented, not eliminated. The suite's
  tolerance encodes the acceptable bound; where a stricter check is wanted,
  provide a deterministic-reduction build option (the upstream
  `RADARSIMX_DEVICE_PARITY` env flag proves this is feasible).
- **Fallback**: `test_system_cpu_fallback.py` runs a child process with
  `CUDA_VISIBLE_DEVICES=-1` and requires `device="gpu"` to produce exactly
  the CPU result with a warning, and `"auto"`/`"cpu"` to stay silent. This
  pins the `gpu_available()` caching semantics.

## 6.5 Performance validation

- Harness: `benchmarks/bench_sbr.py` — sweeps `density, level, tx, rx, model,
  pulses, rcs`; records `wall_s`, `cpu_s`, `threads_busy`; `--compare` renders
  a speedup table against saved JSON.
- Baseline: `benchmarks/baseline/` (pre/post CPU and GPU sweeps).
- Gate (N6): median ratio ours/upstream ≤ 2.0 on every baseline scene at M4;
  track per-scene ratios in CI history to catch regressions.
- Thread scaling: `bench_sbr.py --omp-scaling` (spawns subprocesses per
  `OMP_NUM_THREADS`); expect near-linear scaling to physical cores on the
  mesh scenes.
- Note the upstream observation that cost follows subtended solid angle, not
  triangle count (`benchmarks/scenes.py`) — perf investigations should start
  from ray counts, not mesh size.

## 6.6 Regression infrastructure

- `benchmarks/capture_reference.py`: six scenes chosen for code-path coverage;
  diff a fresh capture against `baseline/*.npz`. Runs in CI on every merge to
  main (CPU) and nightly (GPU, where hardware exists).
- C++ unit tests (GoogleTest) live in the engine repo for geometry kernels
  (intersection, BVH vs brute force), waveform interpolation, and PRNG
  statistics — the layers below the reach of the Python suite.
- Platforms in CI: Windows (MSVC), Ubuntu 22.04/24.04 (GCC 11/13), macOS
  (Clang, ARM64 + x86-64); one CUDA-compile job that verifies the CPU
  fallback path, mirroring upstream's GPU-less GPU CI.

## 6.7 What we cannot validate (and say so)

- Bit-exact equality with upstream GPU output (upstream itself does not
  achieve this; N2 in the perspective doc).
- Behavior under valid licenses vs free tier — the engine ships a permissive
  stub; tier enforcement lives in the Python layer and is unaffected.
- `log_path` HDF5 internal schema — validate by round-tripping our own dumps;
  byte-compatibility with upstream's debug files is not a goal.
