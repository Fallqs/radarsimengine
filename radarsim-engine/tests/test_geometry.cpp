// ==============================================================================
// radarsim-engine — tests/test_geometry.cpp
// Phase 2: Target::Move kinematics and BVH correctness (vs brute force).
// ==============================================================================
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "geom/bvh.hpp"
#include "target.hpp"

using T = float;
int g_failures = 0;

void Check(bool ok, const char *msg) {
    if (!ok) {
        std::printf("FAIL: %s\n", msg);
        ++g_failures;
    }
}

bool Vec3Close(const rsv::Vec3<T> &a, const rsv::Vec3<T> &b, double tol) {
    return std::fabs(a[0] - b[0]) < tol && std::fabs(a[1] - b[1]) < tol &&
           std::fabs(a[2] - b[2]) < tol;
}

// unit triangle in the z=0 plane, vertices (0,0,0) (1,0,0) (0,1,0)
void MakeTriangle(std::vector<T> &points, std::vector<int_t> &cells) {
    points = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    cells = {0, 1, 2};
}

void TestStaticMove() {
    std::vector<T> points;
    std::vector<int_t> cells;
    MakeTriangle(points, cells);
    // translate only
    Target<T> tgt(points.data(), cells.data(), 1, rsv::Vec3<T>(0, 0, 0),
                  rsv::Vec3<T>(10, 0, 0), rsv::Vec3<T>(0, 0, 0),
                  rsv::Vec3<T>(0, 0, 0), rsv::Vec3<T>(0, 0, 0), false, T(0),
                  false);
    tgt.Move(0, 0.0);
    Check(Vec3Close(tgt.vect_mesh_[0].vertex_[1], {11, 0, 0}, 1e-5),
          "static translate");
    Check(tgt.array_size_ == 1, "array_size_ == 1 for static");

    // translate + speed * time
    Target<T> tgt2(points.data(), cells.data(), 1, rsv::Vec3<T>(0, 0, 0),
                   rsv::Vec3<T>(10, 0, 0), rsv::Vec3<T>(2, 0, 0),
                   rsv::Vec3<T>(0, 0, 0), rsv::Vec3<T>(0, 0, 0), false, T(0),
                   false);
    tgt2.Move(0, 3.0);
    Check(Vec3Close(tgt2.vect_mesh_[0].vertex_[1], {17, 0, 0}, 1e-5),
          "static translate + speed * t");
}

void TestRotationOrigin() {
    std::vector<T> points;
    std::vector<int_t> cells;
    MakeTriangle(points, cells);
    const T half_pi = T(1.5707963267948966);
    // yaw 90 deg about origin=(1,0,0): vertex (1,0,0) stays, (0,0,0) -> (1,-1,0)
    Target<T> tgt(points.data(), cells.data(), 1, rsv::Vec3<T>(1, 0, 0),
                  rsv::Vec3<T>(0, 0, 0), rsv::Vec3<T>(0, 0, 0),
                  rsv::Vec3<T>(half_pi, 0, 0), rsv::Vec3<T>(0, 0, 0), false,
                  T(0), false);
    tgt.Move(0, 0.0);
    Check(Vec3Close(tgt.vect_mesh_[0].vertex_[1], {1, 0, 0}, 1e-5),
          "rotate about origin: pivot fixed");
    Check(Vec3Close(tgt.vect_mesh_[0].vertex_[0], {1, -1, 0}, 1e-5),
          "rotate about origin: (0,0,0) -> (1,-1,0)");
}

void TestTimeVaryingMove() {
    std::vector<T> points;
    std::vector<int_t> cells;
    MakeTriangle(points, cells);
    std::vector<rsv::Vec3<T>> locs = {rsv::Vec3<T>(0, 0, 0),
                                      rsv::Vec3<T>(10, 0, 0),
                                      rsv::Vec3<T>(20, 0, 0)};
    std::vector<rsv::Vec3<T>> zeros = {rsv::Vec3<T>(0, 0, 0),
                                       rsv::Vec3<T>(0, 0, 0),
                                       rsv::Vec3<T>(0, 0, 0)};
    Target<T> tgt(points.data(), cells.data(), 1, rsv::Vec3<T>(0, 0, 0), locs,
                  zeros, zeros, zeros, std::complex<T>(T(1e38), T(0)),
                  std::complex<T>(T(1), T(0)), false, T(0), false);
    Check(tgt.array_size_ == 3, "array_size_ == timeline length");
    tgt.Move(2, 0.0);
    Check(Vec3Close(tgt.vect_mesh_[0].vertex_[1], {21, 0, 0}, 1e-5),
          "time-varying: index 2 location");
    tgt.Move(1, 99.0);  // time ignored for time-varying
    Check(Vec3Close(tgt.vect_mesh_[0].vertex_[1], {11, 0, 0}, 1e-5),
          "time-varying: index 1 location, time ignored");
}

void TestBvhVsBruteForce() {
    // random triangle soup
    std::mt19937 rng(42);
    std::uniform_real_distribution<T> dist(-5.0f, 5.0f);
    const int n_tri = 500;
    std::vector<Triangle<T>> tris(n_tri);
    std::vector<rsv::Vec3<T>> soup(n_tri * 3);
    for (auto &v : soup) {
        v = rsv::Vec3<T>(dist(rng), dist(rng), dist(rng));
    }
    for (int i = 0; i < n_tri; ++i) {
        tris[i].vertex_ = &soup[i * 3];
    }
    rsim::Bvh<T> bvh(tris);

    std::uniform_real_distribution<T> rdir(-1.0f, 1.0f);
    int n_checked = 0;
    for (int r = 0; r < 2000; ++r) {
        rsv::Vec3<T> org(dist(rng), dist(rng), dist(rng));
        rsv::Vec3<T> dir(rdir(rng), rdir(rng), rdir(rng));
        const T len = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] +
                                dir[2] * dir[2]);
        dir = rsv::Vec3<T>(dir[0] / len, dir[1] / len, dir[2] / len);

        // brute force
        T best_t = std::numeric_limits<T>::max();
        size_t best_tri = size_t(-1);
        for (size_t i = 0; i < tris.size(); ++i) {
            T t;
            if (rsim::RayTriangle(org, dir, tris[i].vertex_[0],
                                  tris[i].vertex_[1], tris[i].vertex_[2], t) &&
                t < best_t) {
                best_t = t;
                best_tri = i;
            }
        }
        rsim::Bvh<T>::Hit hit;
        const bool bvh_hit = bvh.ClosestHit(org, dir, hit);
        Check(bvh_hit == (best_tri != size_t(-1)), "bvh hit flag matches");
        if (bvh_hit && best_tri != size_t(-1)) {
            Check(std::fabs(hit.t - best_t) < 1e-3, "bvh t matches brute force");
        }
        ++n_checked;
    }
    Check(n_checked == 2000, "bvh checks ran");
}

int main() {
    TestStaticMove();
    TestRotationOrigin();
    TestTimeVaryingMove();
    TestBvhVsBruteForce();
    if (g_failures == 0) {
        std::printf("test_geometry: all passed\n");
        return 0;
    }
    std::printf("test_geometry: %d failures\n", g_failures);
    return 1;
}
