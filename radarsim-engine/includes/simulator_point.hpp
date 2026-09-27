// ==============================================================================
// radarsim-engine — simulator_point.hpp
// Ideal point-target simulator (radarsimc.pxd:372-377).
// Phase 0: validates inputs and returns SUCCESS; physics lands with Phase 1
// (docs/04-physics.md §4.2).
// ==============================================================================
#pragma once

#include <memory>

#include "core/enums.hpp"
#include "points_manager.hpp"
#include "radar.hpp"

template <typename H, typename L, typename ExecutionPolicy>
class PointSimulator {
public:
    PointSimulator() = default;

    RadarSimErrorCode Run(
        const std::shared_ptr<Radar<H, L>> &radar,
        const std::shared_ptr<PointsManager<L>> &points_manager) {
        if (!radar || !points_manager) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (radar->bb_real_ == nullptr || radar->bb_imag_ == nullptr) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        return SUCCESS;
    }
};
