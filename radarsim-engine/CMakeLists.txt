# ==============================================================================
# radarsim-engine — Phase 0 skeleton
# Drop-in replacement source tree for RadarSimPy's src/radarsimcpp.
#
# Consumed by the upstream build:
#   - build.sh / build.bat run CMake on this directory with:
#       -DCMAKE_BUILD_TYPE=Release
#       -DGPU_BUILD=ON|OFF
#       -DENABLE_LICENSE=ON|OFF
#       -DGTEST=ON|OFF
#       -DRADARSIMCPP_DEPS_PREFER_DOWNLOAD=ON|OFF   (accepted, unused here)
#   - setup.py then links -lradarsimcpp from ./build (./build/Release on Windows)
#     and includes ./includes and ./includes/rsvector.
# ==============================================================================
cmake_minimum_required(VERSION 3.18)

project(radarsimcpp VERSION 15.4.0 LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------
option(GPU_BUILD "Build with CUDA support" OFF)
option(GTEST "Build unit tests" OFF)
option(ENABLE_LICENSE "Enable license verification (stub in Phase 0)" OFF)
set(RSIM_LOW_PRECISION "float" CACHE STRING "Low-precision type L: float or double")
set_property(CACHE RSIM_LOW_PRECISION PROPERTY STRINGS float double)

# Accepted for compatibility with the upstream build scripts; the Phase-0
# skeleton has no prebuilt dependencies to download.
set(RADARSIMCPP_DEPS_PREFER_DOWNLOAD OFF CACHE BOOL "Unused in radarsim-engine Phase 0")

# ------------------------------------------------------------------------------
# Language settings
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ------------------------------------------------------------------------------
# Precision policy (docs/03-architecture.md §3.3)
# H (time/frequency/delay/phase) is always double in the Python build.
# L (geometry/BVH/PO kernel) defaults to float; switching it changes which
# rays are launched, so it is a compile-time choice, not a runtime one.
# ------------------------------------------------------------------------------
if(RSIM_LOW_PRECISION STREQUAL "double")
    set(RSIM_PRECISION_DEF "RSIM_LOW_PRECISION_DOUBLE=1")
else()
    set(RSIM_PRECISION_DEF "RSIM_LOW_PRECISION_DOUBLE=0")
endif()

# ------------------------------------------------------------------------------
# OpenMP (optional; macOS without libomp falls back to serial)
# ------------------------------------------------------------------------------
find_package(OpenMP QUIET)

# ------------------------------------------------------------------------------
# CUDA
# ------------------------------------------------------------------------------
if(GPU_BUILD)
    enable_language(CUDA)
    if(DEFINED ENV{RADARSIMX_CUDA_ARCHITECTURES})
        set(CMAKE_CUDA_ARCHITECTURES "$ENV{RADARSIMX_CUDA_ARCHITECTURES}")
    elseif(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "75")  # Turing, matching upstream minimum
    endif()
    set(RSIM_CUDA_SOURCES src/execution_policy_gpu.cu)
else()
    set(RSIM_CUDA_SOURCES "")
endif()

# ------------------------------------------------------------------------------
# Shared library
# The library name must stay "radarsimcpp": upstream setup.py links it by name.
# ------------------------------------------------------------------------------
add_library(radarsimcpp SHARED
    src/execution_policy.cpp
    src/license_manager.cpp
    ${RSIM_CUDA_SOURCES}
)

target_include_directories(radarsimcpp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/includes>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/includes/rsvector>
)

target_compile_definitions(radarsimcpp
    PRIVATE radarsimcpp_EXPORTS
    PUBLIC ${RSIM_PRECISION_DEF}
)

if(GPU_BUILD)
    target_compile_definitions(radarsimcpp PUBLIC _CUDA_)
    set_target_properties(radarsimcpp PROPERTIES
        CUDA_STANDARD 17
        CUDA_SEPARABLE_COMPILATION ON
    )
endif()

if(ENABLE_LICENSE)
    target_compile_definitions(radarsimcpp PUBLIC ENABLE_LICENSE)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(radarsimcpp PUBLIC OpenMP::OpenMP_CXX)
endif()

# Windows: export symbols from the DLL
set_target_properties(radarsimcpp PROPERTIES
    WINDOWS_EXPORT_ALL_SYMBOLS OFF
)

# ------------------------------------------------------------------------------
# Tests (ctest is driven by upstream build scripts when --test=on)
# ------------------------------------------------------------------------------
if(GTEST)
    enable_testing()
    add_subdirectory(tests)
endif()
