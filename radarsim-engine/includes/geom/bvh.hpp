// ==============================================================================
// radarsim-engine — geom/bvh.hpp
// BVH over triangles with Möller–Trumbore intersection.
//
// CPU build: recursive median split on the longest centroid axis (simple and
// correct; the GPU LBVH lands with Phase 6, Karras 2012). Both queries needed
// by the simulators are provided:
//   - ClosestHit: nearest triangle along a ray (LiDAR, ray tracing)
//   - Occluded:   any hit in (t0, t1) (back-propagation occlusion checks)
// Geometry type is the low-precision type L (float), matching the engine's
// BVH precision policy.
// ==============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "rsvector.hpp"
#include "triangle.hpp"

namespace rsim {

template <typename T>
struct Aabb {
    rsv::Vec3<T> lo{std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
                    std::numeric_limits<T>::max()};
    rsv::Vec3<T> hi{std::numeric_limits<T>::lowest(),
                    std::numeric_limits<T>::lowest(),
                    std::numeric_limits<T>::lowest()};

    void Expand(const rsv::Vec3<T> &p) {
        for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], p[a]);
            hi[a] = std::max(hi[a], p[a]);
        }
    }
    rsv::Vec3<T> Centroid() const {
        return rsv::Vec3<T>((lo[0] + hi[0]) * T(0.5), (lo[1] + hi[1]) * T(0.5),
                            (lo[2] + hi[2]) * T(0.5));
    }
};

// Möller–Trumbore. Returns t > 0 on hit.
template <typename T>
inline bool RayTriangle(const rsv::Vec3<T> &org, const rsv::Vec3<T> &dir,
                        const rsv::Vec3<T> &v0, const rsv::Vec3<T> &v1,
                        const rsv::Vec3<T> &v2, T &t_out) {
    const T eps = T(1e-9);
    const rsv::Vec3<T> e1 = v1 - v0;
    const rsv::Vec3<T> e2 = v2 - v0;
    const rsv::Vec3<T> p = dir.Cross(e2);
    const T det = e1.Dot(p);
    if (std::fabs(det) < eps) {
        return false;
    }
    const T inv = T(1) / det;
    const rsv::Vec3<T> tv = org - v0;
    const T u = tv.Dot(p) * inv;
    if (u < T(0) || u > T(1)) {
        return false;
    }
    const rsv::Vec3<T> q = tv.Cross(e1);
    const T v = dir.Dot(q) * inv;
    if (v < T(0) || u + v > T(1)) {
        return false;
    }
    const T t = e2.Dot(q) * inv;
    if (t <= T(0)) {
        return false;
    }
    t_out = t;
    return true;
}

template <typename T>
class Bvh {
public:
    Bvh() = default;

    // triangles: soup of 3-vertex triangles (target.vect_mesh_ layout)
    explicit Bvh(const std::vector<Triangle<T>> &tris) : tris_(&tris) {
        const size_t n = tris.size();
        prim_idx_.resize(n);
        std::iota(prim_idx_.begin(), prim_idx_.end(), 0);
        nodes_.reserve(2 * n);
        if (n > 0) {
            BuildNode(0, n);
        }
    }

    struct Hit {
        T t;
        size_t tri;
    };

    bool ClosestHit(const rsv::Vec3<T> &org, const rsv::Vec3<T> &dir,
                    Hit &hit) const {
        hit.t = std::numeric_limits<T>::max();
        hit.tri = size_t(-1);
        if (nodes_.empty()) {
            return false;
        }
        size_t stack[64];
        int sp = 0;
        stack[sp++] = 0;
        while (sp > 0) {
            const Node &node = nodes_[stack[--sp]];
            T t_box;
            if (!RayAabb(org, dir, node.box, hit.t, t_box)) {
                continue;
            }
            if (node.count > 0) {  // leaf
                for (size_t i = node.first; i < node.first + node.count; ++i) {
                    const Triangle<T> &tri = (*tris_)[prim_idx_[i]];
                    T t;
                    if (RayTriangle(org, dir, tri.vertex_[0], tri.vertex_[1],
                                    tri.vertex_[2], t) &&
                        t < hit.t) {
                        hit.t = t;
                        hit.tri = prim_idx_[i];
                    }
                }
            } else {
                stack[sp++] = node.left;
                stack[sp++] = node.right;
            }
        }
        return hit.tri != size_t(-1);
    }

    // Any hit with t in (t_min, t_max) — occlusion query.
    bool Occluded(const rsv::Vec3<T> &org, const rsv::Vec3<T> &dir, T t_min,
                  T t_max) const {
        if (nodes_.empty()) {
            return false;
        }
        size_t stack[64];
        int sp = 0;
        stack[sp++] = 0;
        while (sp > 0) {
            const Node &node = nodes_[stack[--sp]];
            T t_box;
            if (!RayAabb(org, dir, node.box, t_max, t_box)) {
                continue;
            }
            if (node.count > 0) {
                for (size_t i = node.first; i < node.first + node.count; ++i) {
                    const Triangle<T> &tri = (*tris_)[prim_idx_[i]];
                    T t;
                    if (RayTriangle(org, dir, tri.vertex_[0], tri.vertex_[1],
                                    tri.vertex_[2], t) &&
                        t > t_min && t < t_max) {
                        return true;
                    }
                }
            } else {
                stack[sp++] = node.left;
                stack[sp++] = node.right;
            }
        }
        return false;
    }

private:
    struct Node {
        Aabb<T> box;
        size_t first = 0;  // leaf: first prim index; inner: unused
        size_t count = 0;  // leaf: prim count (0 => inner)
        size_t left = 0, right = 0;
    };

    size_t BuildNode(size_t begin, size_t end) {
        const size_t idx = nodes_.size();
        nodes_.emplace_back();
        Node &node = nodes_.back();

        Aabb<T> box, cbox;
        for (size_t i = begin; i < end; ++i) {
            const Triangle<T> &tri = (*tris_)[prim_idx_[i]];
            Aabb<T> tb;
            tb.Expand(tri.vertex_[0]);
            tb.Expand(tri.vertex_[1]);
            tb.Expand(tri.vertex_[2]);
            for (int a = 0; a < 3; ++a) {
                box.lo[a] = std::min(box.lo[a], tb.lo[a]);
                box.hi[a] = std::max(box.hi[a], tb.hi[a]);
            }
            cbox.Expand(tb.Centroid());
        }
        node.box = box;

        const size_t count = end - begin;
        if (count <= 4) {
            node.first = begin;
            node.count = count;
            return idx;
        }
        // split on the longest centroid axis at the median
        int axis = 0;
        T longest = cbox.hi[0] - cbox.lo[0];
        for (int a = 1; a < 3; ++a) {
            if (cbox.hi[a] - cbox.lo[a] > longest) {
                longest = cbox.hi[a] - cbox.lo[a];
                axis = a;
            }
        }
        const size_t mid = begin + count / 2;
        std::nth_element(prim_idx_.begin() + begin, prim_idx_.begin() + mid,
                         prim_idx_.begin() + end, [&](size_t a, size_t b) {
                             return CentroidOf(a)[axis] < CentroidOf(b)[axis];
                         });
        // nodes_ may reallocate during recursion; index, don't reference
        nodes_[idx].left = BuildNode(begin, mid);
        nodes_[idx].right = BuildNode(mid, end);
        return idx;
    }

    rsv::Vec3<T> CentroidOf(size_t prim) const {
        const Triangle<T> &tri = (*tris_)[prim];
        return rsv::Vec3<T>(
            (tri.vertex_[0][0] + tri.vertex_[1][0] + tri.vertex_[2][0]) / T(3),
            (tri.vertex_[0][1] + tri.vertex_[1][1] + tri.vertex_[2][1]) / T(3),
            (tri.vertex_[0][2] + tri.vertex_[1][2] + tri.vertex_[2][2]) /
                T(3));
    }

    // Ray vs AABB (slab); t_limit prunes boxes beyond the current best hit.
    static bool RayAabb(const rsv::Vec3<T> &org, const rsv::Vec3<T> &dir,
                        const Aabb<T> &box, T t_limit, T &t_hit) {
        T t_lo = T(0);
        T t_hi = t_limit;
        for (int a = 0; a < 3; ++a) {
            if (dir[a] == T(0)) {  // ray parallel to the slab
                if (org[a] < box.lo[a] || org[a] > box.hi[a]) {
                    return false;
                }
                continue;
            }
            const T inv = T(1) / dir[a];
            T t0 = (box.lo[a] - org[a]) * inv;
            T t1 = (box.hi[a] - org[a]) * inv;
            if (t0 > t1) {
                std::swap(t0, t1);
            }
            t_lo = std::max(t_lo, t0);
            t_hi = std::min(t_hi, t1);
            if (t_lo > t_hi) {
                return false;
            }
        }
        t_hit = t_lo;
        return true;
    }

    const std::vector<Triangle<T>> *tris_ = nullptr;
    std::vector<size_t> prim_idx_;
    std::vector<Node> nodes_;
};

}  // namespace rsim
