// ==============================================================================
// radarsim-engine — simulator_rcs.hpp
// Physical-Optics RCS simulator (radarsimc.pxd:397-409): monostatic and
// bistatic, complex polarization vectors on incidence and observation.
// Phase 0: fills GetRcs() with zeros; physics lands with Phase 3
// (docs/04-physics.md §4.6).
// ==============================================================================
#pragma once

#include <complex>
#include <memory>
#include <vector>

#include "core/enums.hpp"
#include "rsvector.hpp"
#include "targets_manager.hpp"

template <typename T, typename ExecutionPolicy, typename L>
class RcsSimulator {
public:
    RcsSimulator() = default;

    RadarSimErrorCode Run(
        const std::shared_ptr<TargetsManager<L>> &targets_manager,
        std::vector<rsv::Vec3<T>> inc_dir_array,
        std::vector<rsv::Vec3<T>> obs_dir_array,
        rsv::Vec3<std::complex<T>> inc_polarization,
        rsv::Vec3<std::complex<T>> obs_polarization, T frequency, T density) {
        if (!targets_manager) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (inc_dir_array.empty() ||
            inc_dir_array.size() != obs_dir_array.size()) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        if (frequency <= T(0) || density <= T(0)) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        (void)inc_polarization;
        (void)obs_polarization;
        rcs_.assign(inc_dir_array.size(), T(0));
        return SUCCESS;
    }

    const std::vector<T> &GetRcs() { return rcs_; }

private:
    std::vector<T> rcs_;
};
