// ==============================================================================
// radarsim-engine — targets_manager.hpp
// Scene container for mesh targets (radarsimc.pxd:320-350). Targets are
// owned by the manager; AddTarget copies the mesh data, so the caller may
// free its arrays immediately.
// ==============================================================================
#pragma once

#include <complex>
#include <memory>
#include <vector>

#include "core/types.hpp"
#include "rsvector.hpp"
#include "target.hpp"

template <typename T>
class TargetsManager {
public:
    TargetsManager() = default;

    void AddTarget(const T *points, const int_t *cells, const int_t &cell_size,
                   const rsv::Vec3<T> &origin,
                   const std::vector<rsv::Vec3<T>> &location_array,
                   const std::vector<rsv::Vec3<T>> &speed_array,
                   const std::vector<rsv::Vec3<T>> &rotation_array,
                   const std::vector<rsv::Vec3<T>> &rotrate_array,
                   const std::complex<T> &ep, const std::complex<T> &mu,
                   const bool &skip_diffusion, const T &density,
                   const bool &environment) {
        targets_.push_back(std::make_shared<Target<T>>(
            points, cells, cell_size, origin, location_array, speed_array,
            rotation_array, rotrate_array, ep, mu, skip_diffusion, density,
            environment));
    }

    void AddTargetSimple(const T *points, const int_t *cells,
                         const int_t &cell_size, const rsv::Vec3<T> &origin,
                         const rsv::Vec3<T> &location,
                         const rsv::Vec3<T> &speed,
                         const rsv::Vec3<T> &rotation,
                         const rsv::Vec3<T> &rotation_rate,
                         const bool &skip_diffusion, const T &density,
                         const bool &environment) {
        targets_.push_back(std::make_shared<Target<T>>(
            points, cells, cell_size, origin, location, speed, rotation,
            rotation_rate, skip_diffusion, density, environment));
    }

    const std::vector<std::shared_ptr<Target<T>>> &targets() const {
        return targets_;
    }

private:
    std::vector<std::shared_ptr<Target<T>>> targets_;
};
