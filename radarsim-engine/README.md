# radarsim-engine

An independent, open re-implementation of the radar simulation backend that
[RadarSimPy](https://github.com/radarsimx/radarsimpy) binds to as
`radarsimcpp`. Project documents (perspective, requirements, architecture,
physics, roadmap, validation) live in `../docs/`.

**Status: Phase 0 — compilable drop-in skeleton.** All public types and
simulators from the binding contract exist with exact signatures; simulator
bodies validate inputs and return `SUCCESS` without computing physics. See
`../docs/05-roadmap.md` for the phase plan.

## Layout (contractual)

The Cython bridge includes headers by relative path and links the library by
name, so this layout is part of the interface:

```
includes/
  core/    execution_policy.hpp  enums.hpp  types.hpp  export.hpp
  libs/    mem_lib.hpp  license_manager.hpp  motion_lib.hpp
  rsvector/rsvector.hpp          (include dir of its own: includes/rsvector)
  transmitter.hpp  receiver.hpp  radar.hpp
  points_manager.hpp  target.hpp  targets_manager.hpp
  triangle.hpp  ray.hpp
  simulator_{point,mesh,rcs,lidar,interference,noise}.hpp
src/       non-template implementations (license stub, GPU probe)
tests/     CTest smoke tests
```

## Build standalone

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGTEST=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Options: `GPU_BUILD`, `GTEST`, `ENABLE_LICENSE` (stub), `RSIM_LOW_PRECISION`
(`float` default, `double` for the all-FP64 reference build).

### Windows without Visual Studio (verified path)

The conda zig toolchain is a self-contained C++20 compiler (clang + mingw-w64):

```bash
conda create -y -n rsimdev -c conda-forge zig cmake ninja
export ZIGCXX_ENV=<conda>/envs/rsimdev   # or set CONDA_PREFIX adjacency
cmake -S . -B build -G Ninja \
    -DCMAKE_CXX_COMPILER=$PWD/tools/zigcxx.bat \
    -DCMAKE_RC_COMPILER=$PWD/tools/zigwindres.bat -DGTEST=ON
cmake --build build && ctest --test-dir build --output-on-failure
```

### Integrated build (engine + Python package, verified)

`tools/build_and_integrate.sh` builds the engine, cythonizes and compiles the
three extension modules against it, and assembles `./radarsimpy/` so the
**unmodified upstream pytest suite** runs against this engine:

```bash
bash radarsim-engine/tools/build_and_integrate.sh
python -m pytest tests/
```

## Drop in as the RadarSimPy backend

```bash
git clone https://github.com/radarsimx/radarsimpy.git
cd radarsimpy
rm -rf src/radarsimcpp          # empty submodule mount point
cp -r /path/to/radarsim-engine src/radarsimcpp
./build.sh --arch=cpu --test=on --license=off
```

The library name (`radarsimcpp`), the include directories (`includes/`,
`includes/rsvector/`), and the CMake option names match what the upstream
`build.sh` / `build.bat` / `setup.py` expect.

## Compatibility notes

- `radarsimx::gpu_policy::is_gpu` / `is_cpu` are static constexpr data
  members, not functions — the binding aliases `gpu_policy::is_gpu` as the
  compile-time `CUDA_BUILD` value.
- The license manager is a permissive stub (`IsLicensed() == true`); usage
  limits are enforced in the Python layer, not here.
- Precision policy: `H = double` (range/delay/phase), `L = float` (geometry/
  BVH/PO kernel) — see `../docs/03-architecture.md` §3.3.
