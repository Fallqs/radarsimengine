// ==============================================================================
// radarsim-engine — core/execution_policy.hpp
// Template-based CPU/GPU execution policy, namespace radarsimx.
//
// Contract notes (src/radarsimpy/includes/radarsimc.pxd:54-95):
//  - is_gpu / is_cpu must be static constexpr DATA members, not functions:
//    the binding aliases `radarsimx::gpu_policy::is_gpu` as a compile-time
//    value (CUDA_BUILD) and uses it as `if CUDA_BUILD:`.
//  - In a CPU-only build, gpu_policy is a genuine type whose is_gpu is false;
//    templated simulators still instantiate it, and Run() must fall back to
//    the CPU implementation.
//  - gpu_available() is the runtime probe: cached for process lifetime and
//    honors CUDA_VISIBLE_DEVICES. Always false without _CUDA_.
// ==============================================================================
#pragma once

#include "core/export.hpp"

namespace radarsimx {

struct execution_policy_base {};

struct cpu_policy : execution_policy_base {
    static constexpr bool is_gpu = false;
    static constexpr bool is_cpu = true;
    static const char *name() { return "cpu"; }
    static int device_id() { return -1; }
};

struct gpu_policy : execution_policy_base {
#ifdef _CUDA_
    static constexpr bool is_gpu = true;
#else
    static constexpr bool is_gpu = false;
#endif
    static constexpr bool is_cpu = false;
    static const char *name() { return "gpu"; }
    static int device_id() { return 0; }
};

// Global policy instances used by the binding.
inline const cpu_policy cpu{};
inline const gpu_policy gpu{};

// Runtime probe: is a usable CUDA device present? Result is cached.
RADARSIMCPP_API bool gpu_available();

}  // namespace radarsimx
