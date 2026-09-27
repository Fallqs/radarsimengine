// ==============================================================================
// radarsim-engine — src/execution_policy_gpu.cu
// CUDA-build runtime probe: a usable device is one that cudaGetDeviceCount
// reports (which honors CUDA_VISIBLE_DEVICES). The result is cached for the
// process lifetime, matching the semantics pinned by
// test_system_cpu_fallback.py (device selection must not change mid-process).
// ==============================================================================
#include <cuda_runtime_api.h>

#include "core/execution_policy.hpp"

namespace radarsimx {

bool gpu_available() {
    static const bool available = [] {
        int count = 0;
        const cudaError_t err = cudaGetDeviceCount(&count);
        return err == cudaSuccess && count > 0;
    }();
    return available;
}

}  // namespace radarsimx
