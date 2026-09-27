// ==============================================================================
// radarsim-engine — core/types.hpp
// Primitive type aliases matching src/radarsimpy/includes/type_def.pxd.
// Global namespace: the Cython binding declares these headers without one.
// ==============================================================================
#pragma once

using int_t = int;
using uint_t = unsigned int;

// Low-precision type L for the Python build (geometry/BVH/PO kernel).
// See docs/03-architecture.md §3.3 — this is a compile-time policy.
#if RSIM_LOW_PRECISION_DOUBLE
using float_t = double;
#else
using float_t = float;
#endif
