// ==============================================================================
// radarsim-engine — target.hpp
// One mesh target: geometry (copied from caller), material, kinematics, and
// the sampling/return flags (radarsimc.pxd:282-313).
//
// Kinematics contract (cp_radarsimc_mesh.pyx):
//   - kinematics arrays are either size 1 (static: value + rate * time) or
//     one entry per timeline point (time-varying: use array entry `index`,
//     rates ignored)
//   - array_size_ is the kinematics timeline length; the binding branches on
//     array_size_ > 1 to pick the query mode of Move()
//   - Move(index, time) rewrites vect_mesh_ vertices:
//       world = Rotate(local - origin, rot) + origin + loc
//     with Rotate the Rz(yaw) Ry(-pitch) Rx(roll) convention, angles in
//     radians (the marshalling layer converts).
// ==============================================================================
#pragma once

#include <complex>
#include <vector>

#include "core/types.hpp"
#include "libs/motion_lib.hpp"
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
        : Target(points, cells, cell_size, origin,
                 std::vector<rsv::Vec3<T>>{location},
                 std::vector<rsv::Vec3<T>>{speed},
                 std::vector<rsv::Vec3<T>>{rotation},
                 std::vector<rsv::Vec3<T>>{rotation_rate},
                 std::complex<T>(T(1e38), T(0)),  // PEC
                 std::complex<T>(T(1), T(0)), skip_diffusion, density,
                 environment) {}

    // Move the mesh to its pose at timeline `index` / absolute `time`.
    // Time-varying (array_size_ > 1): use entry `index`, ignore `time`.
    // Static: extrapolate base + rate * time.
    void Move(int index, double time) {
        rsv::Vec3<T> loc, rot;
        if (array_size_ > 1) {
            const size_t i =
                static_cast<size_t>(index) < location_array_.size()
                    ? static_cast<size_t>(index)
                    : location_array_.size() - 1;
            loc = location_array_[i];
            rot = rotation_array_[i];
        } else {
            const T t = static_cast<T>(time);
            const rsv::Vec3<T> &spd = speed_array_[0];
            const rsv::Vec3<T> &rrt = rotrate_array_[0];
            loc = location_array_[0] + spd * t;
            rot = rotation_array_[0] + rrt * t;
        }
        // rotation matrix (Rz.Ry(-pitch).Rx(roll)) in double for accuracy
        double R[3][3];
        {
            const double yaw = rot[0], pitch = rot[1], roll = rot[2];
            const double cy = std::cos(yaw), sy = std::sin(yaw);
            const double cp = std::cos(pitch), sp = std::sin(pitch);
            const double cr = std::cos(roll), sr = std::sin(roll);
            R[0][0] = cy * cp;
            R[0][1] = cy * sp * sr - sy * cr;
            R[0][2] = cy * sp * cr + sy * sr;
            R[1][0] = sy * cp;
            R[1][1] = sy * sp * sr + cy * cr;
            R[1][2] = sy * sp * cr - cy * sr;
            R[2][0] = -sp;
            R[2][1] = cp * sr;
            R[2][2] = cp * cr;
        }
        const size_t n = base_vertices_.size();
        for (size_t i = 0; i < n; ++i) {
            const rsv::Vec3<T> &v = base_vertices_[i];
            // world = R . (v - origin) + origin + location
            const double x = static_cast<double>(v[0] - origin_[0]);
            const double y = static_cast<double>(v[1] - origin_[1]);
            const double z = static_cast<double>(v[2] - origin_[2]);
            vertices_[i] = rsv::Vec3<T>(
                static_cast<T>(R[0][0] * x + R[0][1] * y + R[0][2] * z +
                               origin_[0] + loc[0]),
                static_cast<T>(R[1][0] * x + R[1][1] * y + R[1][2] * z +
                               origin_[1] + loc[1]),
                static_cast<T>(R[2][0] * x + R[2][1] * y + R[2][2] * z +
                               origin_[2] + loc[2]));
        }
    }

    std::vector<Triangle<T>> vect_mesh_;  // moved triangle soup
    int array_size_ = 0;                  // kinematics timeline length

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

    // Original (unmoved) vertex soup; vertices_ is the moved copy that
    // vect_mesh_ points into.
    std::vector<rsv::Vec3<T>> base_vertices_;
    std::vector<rsv::Vec3<T>> vertices_;

private:
    void InitMesh(const T *points, const int_t *cells, const int_t &cell_size) {
        base_vertices_.resize(static_cast<size_t>(cell_size) * 3);
        for (int_t c = 0; c < cell_size; ++c) {
            for (int_t v = 0; v < 3; ++v) {
                const int_t vi = cells[c * 3 + v];
                base_vertices_[c * 3 + v] =
                    rsv::Vec3<T>(points[vi * 3], points[vi * 3 + 1],
                                 points[vi * 3 + 2]);
            }
        }
        vertices_ = base_vertices_;
        vect_mesh_.resize(cell_size);
        for (int_t c = 0; c < cell_size; ++c) {
            vect_mesh_[c].vertex_ = &vertices_[c * 3];
        }
        array_size_ = static_cast<int>(location_array_.size());
    }
};
