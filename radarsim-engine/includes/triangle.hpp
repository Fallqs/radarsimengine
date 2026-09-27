// ==============================================================================
// radarsim-engine — triangle.hpp
// A triangle as three vertices in one contiguous allocation; vertex_ points at
// the first of the three (radarsimc.pxd:277-280).
// ==============================================================================
#pragma once

#include "rsvector.hpp"

template <typename T>
struct Triangle {
    Triangle() = default;

    rsv::Vec3<T> *vertex_ = nullptr;
};
