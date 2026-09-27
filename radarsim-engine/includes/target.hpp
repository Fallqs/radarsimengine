// ==============================================================================
// radarsim-engine — target.hpp
// One mesh target: geometry (copied from caller), material, kinematics, and
// the sampling/return flags (radarsimc.pxd:282-313).
//
// Kinematics: location/speed/rotation/rotation-rate are either single values
// (static constructor) or per-timestamp arrays (time-varying constructor).
// Move(index, time) evaluates the pose at one timeline point.
// ==============================================================================
#pragma once

#include <complex>
#include <vector>

#include "core/types.hpp"
#include "rsvector.hpp"
#include "triangle.hpp"

template <typename T>
struct Target {
    Target() = default;

    // Time-varying kinematics, with material properties.
    Target(const T *points, const int_t *cells, const int_t &cell_size,
           const rsv::Vec3<T> &origin,
           const std::vector<rsv::Vec3<T>> &location_array,
           const std::vector<rsv::Vec3<T>> &speed_array,
           const std::vector<rsv::Vec3<T>> &rotation_array,
           const std::vector<rsv::Vec3<T>> &rotrate_array,
           const std::complex<T> &ep, const std::complex<T> &mu,
           const bool &skip_diffusion, const T &density,
           const bool &environment)
        : origin_(origin),
          location_array_(location_array),
          speed_array_(speed_array),
          rotation_array_(rotation_array),
          rotrate_array_(rotrate_array),
          permittivity_(ep),
          permeability_(mu),
          skip_diffusion_(skip_diffusion),
          density_(density),
          environment_(environment) {
        InitMesh(points, cells, cell_size);
    }

    // Static kinematics, PEC material (no ep/mu in this overload).
    Target(const T *points, const int_t *cells, const int_t &cell_size,
           const rsv::Vec3<T> &origin, const rsv::Vec3<T> &location,
           const rsv::Vec3<T> &speed, const rsv::Vec3<T> &rotation,
           const rsv::Vec3<T> &rotation_rate, const bool &skip_diffusion,
           const T &density, const bool &environment)
        : origin_(origin),
          location_array_{location},
          speed_array_{speed},
          rotation_array_{rotation},
          rotrate_array_{rotation_rate},
          permittivity_(T(1e38), T(0)),  // PEC, matches cp_radarsimc default
          permeability_(T(1), T(0)),
          skip_diffusion_(skip_diffusion),
          density_(density),
          environment_(environment) {
        InitMesh(points, cells, cell_size);
    }

    // Evaluate pose at timeline entry `index` / absolute `time`.
    // Phase 0: stores the query; kinematics evaluation lands with Phase 2.
    void Move(int index, double time) {
        last_move_index_ = index;
        last_move_time_ = time;
    }

    std::vector<Triangle<T>> vect_mesh_;  // Triangle soup (3 verts per entry)
    int array_size_ = 0;                  // Number of vertex Vec3s (3 x cells)

    rsv::Vec3<T> origin_;
    std::vector<rsv::Vec3<T>> location_array_;
    std::vector<rsv::Vec3<T>> speed_array_;
    std::vector<rsv::Vec3<T>> rotation_array_;
    std::vector<rsv::Vec3<T>> rotrate_array_;
    std::complex<T> permittivity_{T(1e38), T(0)};
    std::complex<T> permeability_{T(1), T(0)};
    bool skip_diffusion_ = false;
    T density_ = T(0);  // 0 = use the global sim_radar density
    bool environment_ = false;

private:
    void InitMesh(const T *points, const int_t *cells, const int_t &cell_size) {
        vertices_.resize(static_cast<size_t>(cell_size) * 3);
        for (int_t c = 0; c < cell_size; ++c) {
            for (int_t v = 0; v < 3; ++v) {
                const int_t vi = cells[c * 3 + v];
                vertices_[c * 3 + v] = rsv::Vec3<T>(
                    points[vi * 3], points[vi * 3 + 1], points[vi * 3 + 2]);
            }
        }
        vect_mesh_.resize(cell_size);
        for (int_t c = 0; c < cell_size; ++c) {
            vect_mesh_[c].vertex_ = &vertices_[c * 3];
        }
        array_size_ = cell_size * 3;
    }

    std::vector<rsv::Vec3<T>> vertices_;  // backing store for vect_mesh_
    int last_move_index_ = -1;
    double last_move_time_ = 0.0;
};
