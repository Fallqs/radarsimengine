#!/usr/bin/env bash
# Build radarsim-engine and integrate it into a working radarsimpy package on
# machines without MSVC, using the conda zig toolchain.
#
#   conda install -n rsimdev -c conda-forge zig cmake ninja
#   conda install -n base -c conda-forge cython trimesh
#   bash radarsim-engine/tools/build_and_integrate.sh
#
# Produces ./radarsimpy/ (gitignored upstream) with the compiled extensions,
# then `python -m pytest tests/` runs the upstream suite against the engine.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENGINE_SRC="$REPO/radarsim-engine"
BUILD="$ENGINE_SRC/build"
ENVDIR="${ZIGCXX_ENV:-/c/Users/27382/miniconda3/envs/rsimdev}"
export PATH="$(cygpath -u "$ENVDIR")/Library/bin:$PATH"
export ZIGCXX_ENV="$ENVDIR"
export PREFIX="$ENVDIR"

PYINC="$(python -c 'import sysconfig; print(sysconfig.get_paths()["include"])')"
PYLIBS="$(python -c 'import sysconfig; print(sysconfig.get_config_var("prefix"))')/libs"
NUMPYINC="$(python -c 'import numpy; print(numpy.get_include())')"
PYVER="$(python -c 'import sys; print(f"python{sys.version_info.major}{sys.version_info.minor}")')"

CXX="$ENVDIR/Library/bin/x86_64-w64-mingw32-zig-cxx.exe"
TARGET="-target x86_64-windows-gnu"
CXXFLAGS="-std=c++20 -O2 -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION"
INCLUDES="-I$ENGINE_SRC/includes -I$ENGINE_SRC/includes/rsvector -I$PYINC -I$NUMPYINC"

echo "== 1. engine (cmake + zig) =="
cmake -S "$ENGINE_SRC" -B "$BUILD" -G Ninja \
    -DCMAKE_CXX_COMPILER="$ENGINE_SRC/tools/zigcxx.bat" \
    -DCMAKE_RC_COMPILER="$ENGINE_SRC/tools/zigwindres.bat" \
    -DGTEST=ON > /dev/null
cmake --build "$BUILD" > /dev/null
ctest --test-dir "$BUILD" --output-on-failure

echo "== 2. cythonize =="
CYBUILD="$REPO/build-cython"
mkdir -p "$CYBUILD"
for mod in simulator license; do
    cython --cplus --directive language_level=3 -I "$REPO/src" \
        "$REPO/src/radarsimpy/$mod.pyx" -o "$CYBUILD/$mod.cpp"
done
mkdir -p "$CYBUILD/lib"
cython --cplus --directive language_level=3 -I "$REPO/src" \
    "$REPO/src/radarsimpy/lib/cp_radarsimc.pyx" -o "$CYBUILD/lib/cp_radarsimc.cpp"

echo "== 3. compile extensions =="
mkdir -p "$REPO/radarsimpy/lib"
for mod in simulator license; do
    "$CXX" $TARGET $CXXFLAGS $INCLUDES -shared "$CYBUILD/$mod.cpp" \
        -L"$BUILD" -lradarsimcpp -L"$PYLIBS" -l"$PYVER" \
        -o "$REPO/radarsimpy/$mod.pyd"
done
"$CXX" $TARGET $CXXFLAGS $INCLUDES -shared "$CYBUILD/lib/cp_radarsimc.cpp" \
    -L"$BUILD" -lradarsimcpp -L"$PYLIBS" -l"$PYVER" \
    -o "$REPO/radarsimpy/lib/cp_radarsimc.pyd"

echo "== 4. assemble package =="
cp "$REPO"/src/radarsimpy/*.py "$REPO/radarsimpy/"
cp "$REPO"/src/radarsimpy/lib/__init__.py "$REPO/radarsimpy/lib/"
cp "$BUILD/radarsimcpp.dll" "$REPO/radarsimpy/"
echo "done: $REPO/radarsimpy"
