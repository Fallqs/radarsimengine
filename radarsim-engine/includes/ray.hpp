// ==============================================================================
// radarsim-engine — ray.hpp
// Ray record for tracing and LiDAR output (radarsimc.pxd:356-363). The array
// members index the bounce history: location_/normal_/range_ per bounce.
// ==============================================================================
#pragma once

#include "rsvector.hpp"

template <typename T, typename L>
struct Ray {
    Ray() = default;

    rsv::Vec3<T> *direction_ = nullptr;  // Ray direction (unit)
    rsv::Vec3<T> *location_ = nullptr;   // Origin / hit locations per bounce
    rsv::Vec3<T> *normal_ = nullptr;     // Surface normal at each hit
    T *range_ = nullptr;                 // Cumulative range at each bounce
    int reflections_ = 0;                // Number of reflections encountered
};
