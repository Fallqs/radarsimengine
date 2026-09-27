// ==============================================================================
// radarsim-engine — src/execution_policy.cpp
// CPU-build implementation of the runtime GPU probe: always unavailable.
// The CUDA build replaces this translation unit's probe with
// execution_policy_gpu.cu.
// ==============================================================================
#include "core/execution_policy.hpp"

#ifndef _CUDA_

namespace radarsimx {

bool gpu_available() {
    static const bool available = false;
    return available;
}

}  // namespace radarsimx

#endif
