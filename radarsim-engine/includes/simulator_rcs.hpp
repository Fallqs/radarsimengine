// ==============================================================================
// radarsim-engine — simulator_rcs.hpp
// Physical-Optics RCS simulator (radarsimc.pxd:397-409).
//
// Far-field PO sum over the lit, unoccluded surface:
//   sigma = (k^2/pi) |sum_samples val * exp(j k (inc+obs) . r) dA|^2
// with the vector-PO kernel val = conj(p_obs) . (J - (J.o) o),
// J = n x (i x p_inc), i = propagation direction = -inc_dir.
// Facets are sampled at `density` points per wavelength; occlusion against
// the incident direction uses the scene BVH.
//
// Divergence notice: the upstream goldens embed the engine's exact sampling
// grid; this implementation reproduces the PO model but not the exact
// recorded values (docs/mesh_simulator_model.md).
// ==============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "core/enums.hpp"
#include "core/types.hpp"
#include "geom/bvh.hpp"
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

        // gather all triangles (RCS targets are posed once; Move(0,0))
        size_t total_verts = 0;
        for (const auto &tgt : targets_manager->targets()) {
            total_verts += tgt->vertices_.size();
        }
        std::vector<rsv::Vec3<L>> verts;
        verts.reserve(total_verts);
        std::vector<Triangle<L>> tris;
        for (const auto &tgt : targets_manager->targets()) {
            tgt->Move(0, 0.0);
            const size_t base = verts.size();
            verts.insert(verts.end(), tgt->vertices_.begin(),
                         tgt->vertices_.end());
            for (size_t ci = 0; ci < tgt->vect_mesh_.size(); ++ci) {
                Triangle<L> t2;
                t2.vertex_ = verts.data() + base + ci * 3;
                tris.push_back(t2);
            }
        }
        if (tris.empty()) {
            rcs_.assign(inc_dir_array.size(), T(0));
            return SUCCESS;
        }
        rsim::Bvh<L> bvh(tris);

        const double lam = 299792458.0 / static_cast<double>(frequency);
        const double k = 2.0 * 3.14159265358979323846 / lam;
        const double ds = lam / static_cast<double>(density);

        rcs_.assign(inc_dir_array.size(), T(0));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int idx = 0; idx < static_cast<int>(inc_dir_array.size()); ++idx) {
            // i: propagation direction (wave travels from source to target)
            const double i[3] = {-static_cast<double>(inc_dir_array[idx][0]),
                                 -static_cast<double>(inc_dir_array[idx][1]),
                                 -static_cast<double>(inc_dir_array[idx][2])};
            const double o[3] = {static_cast<double>(obs_dir_array[idx][0]),
                                 static_cast<double>(obs_dir_array[idx][1]),
                                 static_cast<double>(obs_dir_array[idx][2])};
            const std::complex<double> pi[3] = {
                std::complex<double>(inc_polarization[0]),
                std::complex<double>(inc_polarization[1]),
                std::complex<double>(inc_polarization[2])};
            const std::complex<double> po[3] = {
                std::conj(std::complex<double>(obs_polarization[0])),
                std::conj(std::complex<double>(obs_polarization[1])),
                std::conj(std::complex<double>(obs_polarization[2]))};
            std::complex<double> sum(0.0, 0.0);

            for (size_t ti = 0; ti < tris.size(); ++ti) {
                const rsv::Vec3<L> &v0 = tris[ti].vertex_[0];
                const rsv::Vec3<L> &v1 = tris[ti].vertex_[1];
                const rsv::Vec3<L> &v2 = tris[ti].vertex_[2];
                const rsv::Vec3<L> e1 = v1 - v0;
                const rsv::Vec3<L> e2 = v2 - v0;
                rsv::Vec3<L> n = e1.Cross(e2);
                const double nl = std::sqrt(n.Dot(n));
                if (nl == 0.0) {
                    continue;
                }
                n = n * L(1.0 / nl);
                const double area = 0.5 * nl;
                // lit side: the outward normal faces the incident wave
                const double facing =
                    -(n[0] * i[0] + n[1] * i[1] + n[2] * i[2]);
                if (facing <= 0.0) {
                    continue;
                }
                // H_inc = i x p_inc; PO current J = n x H_inc
                const std::complex<double> hx[3] = {i[1] * pi[2] - i[2] * pi[1],
                                                    i[2] * pi[0] - i[0] * pi[2],
                                                    i[0] * pi[1] - i[1] * pi[0]};
                const double nd[3] = {static_cast<double>(n[0]),
                                      static_cast<double>(n[1]),
                                      static_cast<double>(n[2])};
                const std::complex<double> J[3] = {nd[1] * hx[2] - nd[2] * hx[1],
                                                   nd[2] * hx[0] - nd[0] * hx[2],
                                                   nd[0] * hx[1] - nd[1] * hx[0]};
                const std::complex<double> Jdo =
                    J[0] * o[0] + J[1] * o[1] + J[2] * o[2];
                const std::complex<double> Jp[3] = {J[0] - Jdo * o[0],
                                                    J[1] - Jdo * o[1],
                                                    J[2] - Jdo * o[2]};
                const std::complex<double> val =
                    Jp[0] * po[0] + Jp[1] * po[1] + Jp[2] * po[2];
                if (std::abs(val) == 0.0) {
                    continue;
                }

                // per-facet samples at density points per wavelength
                const int n_s = std::max(
                    1, static_cast<int>(std::ceil(std::sqrt(area) / ds)));
                const int n_samples = (n_s + 1) * (n_s + 2) / 2;
                const double da = area / n_samples;
                for (int a = 0; a <= n_s; ++a) {
                    for (int b = 0; b + a <= n_s; ++b) {
                        const double u = (a + 1.0 / 3.0) / (n_s + 1.0);
                        const double w = (b + 1.0 / 3.0) / (n_s + 1.0);
                        const rsv::Vec3<L> pt = v0 + e1 * L(u) + e2 * L(w);
                        // occlusion: visible from the source direction?
                        // step off the lit face toward the source
                        const rsv::Vec3<L> dir_i(L(-i[0]), L(-i[1]), L(-i[2]));
                        const rsv::Vec3<L> src = pt + dir_i * L(1e-4);
                        if (bvh.Occluded(src, dir_i, L(1e-6), L(1e9))) {
                            continue;
                        }
                        const double phase =
                            k * ((static_cast<double>(inc_dir_array[idx][0]) +
                                  o[0]) *
                                     pt[0] +
                                 (static_cast<double>(inc_dir_array[idx][1]) +
                                  o[1]) *
                                     pt[1] +
                                 (static_cast<double>(inc_dir_array[idx][2]) +
                                  o[2]) *
                                     pt[2]);
                        sum += val *
                               std::exp(std::complex<double>(0, phase)) * da;
                    }
                }
            }
            rcs_[idx] = static_cast<T>(k * k / 3.14159265358979323846 *
                                       std::norm(sum));
        }
        return SUCCESS;
    }

    const std::vector<T> &GetRcs() { return rcs_; }

private:
    std::vector<T> rcs_;
};
