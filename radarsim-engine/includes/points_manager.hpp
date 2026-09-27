// ==============================================================================
// radarsim-engine — points_manager.hpp
// Container for ideal point targets (radarsimc.pxd:258-273). Time-varying
// location/rcs/phase arrays are indexed by the flattened
// [channels, pulses, samples] timestamp grid.
// ==============================================================================
#pragma once

#include <vector>

#include "rsvector.hpp"

template <typename T>
class PointsManager {
public:
    PointsManager() = default;

    struct Point {
        std::vector<rsv::Vec3<T>> location_array;
        rsv::Vec3<T> speed;
        std::vector<T> rcs_array;    // dBsm
        std::vector<T> phase_array;  // rad
    };

    void AddPoint(const std::vector<rsv::Vec3<T>> &location_array,
                  const rsv::Vec3<T> &speed_array,
                  const std::vector<T> &rcs_array,
                  const std::vector<T> &phase_array) {
        points_.push_back(
            Point{location_array, speed_array, rcs_array, phase_array});
    }

    void AddPointSimple(const rsv::Vec3<T> &location, const rsv::Vec3<T> &speed,
                        const T &rcs, const T &phase) {
        points_.push_back(Point{{location}, speed, {rcs}, {phase}});
    }

    const std::vector<Point> &points() const { return points_; }

private:
    std::vector<Point> points_;
};
