// ==============================================================================
// radarsim-engine — simulator_lidar.hpp
// LiDAR point-cloud simulator (radarsimc.pxd:413-421).
//
// Ray fan over (phi, theta) — phi outer, theta inner — with directions
//   d = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta))
// (z-up; theta from zenith, per the project coordinate convention). Per ray
// that hits, cloud_ gets a Ray with two timeline entries:
//   [0] = sensor (origin, incident direction)   [1] = hit (point, reflected
//   direction, face normal flipped toward the sensor, range)
// Kinematics are evaluated once at frame_time (the marshalling layer already
// folds frame_time into the target locations/rotations), i.e. Move(0, 0).
// ==============================================================================
#pragma once

#include <array>
#include <cmath>
#include <deque>
#include <memory>
#include <vector>

#include "core/enums.hpp"
#include "geom/bvh.hpp"
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

        // gather all moved targets into one scene soup (two passes so the
        // triangle pointers never dangle: vertex_store_ is reserved up front)
        for (const auto &tgt : targets_manager->targets()) {
            tgt->Move(0, 0.0);
        }
        size_t total_verts = 0;
        for (const auto &tgt : targets_manager->targets()) {
            total_verts += tgt->vertices_.size();
        }
        vertex_store_.clear();
        vertex_store_.reserve(total_verts);
        scene_.clear();
        for (const auto &tgt : targets_manager->targets()) {
            const size_t base = vertex_store_.size();
            vertex_store_.insert(vertex_store_.end(), tgt->vertices_.begin(),
                                 tgt->vertices_.end());
            for (size_t c = 0; c < tgt->vect_mesh_.size(); ++c) {
                Triangle<T> tri;
                tri.vertex_ = vertex_store_.data() + base + c * 3;
                scene_.push_back(tri);
            }
        }

        Bvh scene_bvh(scene_);

        cloud_.clear();
        loc_store_.clear();
        dir_store_.clear();
        nrm_store_.clear();
        rng_store_.clear();

        for (const T &ph : phi) {
            for (const T &th : theta) {
                const rsv::Vec3<T> dir(std::sin(th) * std::cos(ph),
                                       std::sin(th) * std::sin(ph),
                                       std::cos(th));
                typename rsim::Bvh<T>::Hit hit;
                if (!scene_bvh.ClosestHit(position, dir, hit)) {
                    continue;
                }
                const Triangle<T> &tri = scene_[hit.tri];
                const rsv::Vec3<T> p0 = tri.vertex_[0];
                const rsv::Vec3<T> e1 = tri.vertex_[1] - p0;
                const rsv::Vec3<T> e2 = tri.vertex_[2] - p0;
                rsv::Vec3<T> n = e1.Cross(e2);
                const T nl = std::sqrt(n.Dot(n));
                if (nl > T(0)) {
                    n = n * (T(1) / nl);
                }
                if (n.Dot(dir) > T(0)) {  // face the sensor
                    n = n * T(-1);
                }
                // reflected direction: d - 2 (d.n) n
                const T dn = dir.Dot(n);
                const rsv::Vec3<T> refl = dir - n * (T(2) * dn);
                const rsv::Vec3<T> hitp = position + dir * hit.t;

                loc_store_.emplace_back(std::array<rsv::Vec3<T>, 2>{
                    rsv::Vec3<T>(position[0], position[1], position[2]), hitp});
                dir_store_.emplace_back(
                    std::array<rsv::Vec3<T>, 2>{dir, refl});
                nrm_store_.emplace_back(std::array<rsv::Vec3<T>, 2>{
                    rsv::Vec3<T>(T(0), T(0), T(0)), n});
                rng_store_.emplace_back(std::array<T, 2>{T(0), hit.t});

                Ray<T, T> ray;
                ray.location_ = loc_store_.back().data();
                ray.direction_ = dir_store_.back().data();
                ray.normal_ = nrm_store_.back().data();
                ray.range_ = rng_store_.back().data();
                ray.reflections_ = 1;
                cloud_.push_back(ray);
            }
        }
        return SUCCESS;
    }

    std::vector<Ray<T, T>> cloud_;

private:
    using Bvh = rsim::Bvh<T>;
    std::vector<Triangle<T>> scene_;
    std::vector<rsv::Vec3<T>> vertex_store_;
    // stable backing for the Ray pointer fields (deque never moves elements)
    std::deque<std::array<rsv::Vec3<T>, 2>> loc_store_;
    std::deque<std::array<rsv::Vec3<T>, 2>> dir_store_;
    std::deque<std::array<rsv::Vec3<T>, 2>> nrm_store_;
    std::deque<std::array<T, 2>> rng_store_;
};
