// ==============================================================================
// radarsim-engine — simulator_interference.hpp
// Radar-to-radar interference simulator (radarsimc.pxd:426-431). Evaluates
// the interferer's waveform inside the victim's receive window and writes
// into the victim's (re-pointed) baseband buffers.
// Phase 0: validates inputs; physics lands with Phase 1
// (docs/04-physics.md §4.8).
// ==============================================================================
#pragma once

#include <memory>

#include "core/enums.hpp"
#include "radar.hpp"

template <typename H, typename L, typename ExecutionPolicy>
class InterferenceSimulator {
public:
    InterferenceSimulator() = default;

    RadarSimErrorCode Run(const std::shared_ptr<Radar<H, L>> &radar,
                          const std::shared_ptr<Radar<H, L>> &interf_radar) {
        if (!radar || !interf_radar) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (radar->bb_real_ == nullptr || radar->bb_imag_ == nullptr) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        return SUCCESS;
    }
};
