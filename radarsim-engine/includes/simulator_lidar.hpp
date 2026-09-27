// ==============================================================================
// radarsim-engine — simulator_lidar.hpp
// LiDAR point-cloud simulator (radarsimc.pxd:413-421). Ray fan over
// (phi, theta) from the sensor position; closest hit fills cloud_.
// Phase 0: clears cloud_; ray casting lands with Phase 5 on the Phase-2 BVH
// (docs/04-physics.md §4.7).
// ==============================================================================
#pragma once

#include <memory>
#include <vector>

#include "core/enums.hpp"
#include "ray.hpp"
#include "rsvector.hpp"
#include "targets_manager.hpp"

template <typename T, typename ExecutionPolicy>
class LidarSimulator {
public:
    LidarSimulator() = default;

    RadarSimErrorCode Run(
        const std::shared_ptr<TargetsManager<T>> &targets_manager,
        const std::vector<T> &phi, const std::vector<T> &theta,
        const rsv::Vec3<T> &position) {
        if (!targets_manager) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (phi.empty() || theta.empty()) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        (void)position;
        cloud_.clear();
        return SUCCESS;
    }

    std::vector<Ray<T, T>> cloud_;
};
