// ==============================================================================
// radarsim-engine — simulator_mesh.hpp
// SBR + Physical Optics mesh simulator (radarsimc.pxd:382-393).
//
// Pipeline (docs/04-physics.md §4.5, docs/mesh_simulator_model.md):
//   1. per ray-tracing pass (level: 0=frame, 1=pulse, 2=sample), targets are
//      moved to the pass instant; between passes, ranges are extrapolated at
//      each sample's path-projected range rate
//   2. occupancy: probe the angular grid (phi in [-90,90], theta in [0,180],
//      step = channel grid) against the scene mesh, with one-cell dilation
//   3. fine rays: per occupied cell, count = int(grid/atan(lam_min/density/R))+1
//      per axis, anchored at the cell's lower edge
//   4. SBR: specular multi-bounce with Fresnel coefficients; skip_diffusion
//      surfaces redirect but never return; ray_filter caps/filters bounces
//   5. PO: each contributing landing scatters to every Rx channel:
//      bb += K * (k/2) * (-j) * exp(j 2 pi beat) * stuff * refl * sinc * dA
//          / (Rt Rr)
//      with K the point-target amplitude-chain constant (validated on the
//      plate golden to ~2%; the engine's exact sampler remains a documented
//      divergence)
// ==============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "core/enums.hpp"
#include "core/types.hpp"
#include "geom/bvh.hpp"
#include "geom/fresnel.hpp"
#include "radar.hpp"
#include "rf/antenna.hpp"
#include "rf/platform.hpp"
#include "rf/waveform.hpp"
#include "rsvector.hpp"
#include "targets_manager.hpp"

namespace rsim {

// One triangle in the scene, with its target's material/flags.
template <typename T>
struct SceneTri {
    rsv::Vec3<T> v[3];
    rsv::Vec3<T> normal;  // unit, from (v1-v0)x(v2-v0)
    T area;
    std::complex<T> eps, mu;
    bool skip_diffusion;
    bool environment;
    int target_idx;  // for range-rate extrapolation
};

}  // namespace rsim

template <typename H, typename L, typename ExecutionPolicy>
class MeshSimulator {
public:
    MeshSimulator() = default;

    RadarSimErrorCode Run(
        const std::shared_ptr<Radar<H, L>> &radar,
        const std::shared_ptr<TargetsManager<L>> &targets_manager, int level,
        L density, rsv::Vec2<int_t> ray_filter, bool back_propagating,
        std::string log_path, bool dry_run) {
        if (!radar || !targets_manager) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (level < 0 || level > 2) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        if (density <= L(0)) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        if (ray_filter[0] < 0 || ray_filter[1] < ray_filter[0]) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        if (!dry_run &&
            (radar->bb_real_ == nullptr || radar->bb_imag_ == nullptr)) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        (void)log_path;  // HDF5 ray dump not implemented (no HDF5 dep)
        (void)back_propagating;  // return-leg reflections: not yet implemented

        const auto tx = radar->tx_;
        const auto rx = radar->rx_;
        const int n_tx = static_cast<int>(tx->channels_.size());
        const int n_rx = static_cast<int>(rx->channels_.size());
        const int pulses = static_cast<int>(tx->pulse_start_time_.size());
        const int samples = radar->sample_size_;
        const int frames = static_cast<int>(radar->frame_start_time_.size());
        if (n_tx == 0 || n_rx == 0 || pulses == 0 || samples <= 0) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        if (targets_manager->targets().empty()) {
            return SUCCESS;
        }

        const double fs = static_cast<double>(rx->fs_);
        const double gate = rx->gate_delay_;
        const rsim::Waveform<H> waveform(tx->freq_, tx->freq_time_);
        const double fc = static_cast<double>(waveform.CenterFreq());
        const double lam_min =
            kC /
            static_cast<double>(
                *std::max_element(tx->freq_.begin(), tx->freq_.end()));

        // pass layout
        const int pass_pulses = (level == 0) ? 1 : pulses;
        const int pass_samples = (level == 2) ? samples : 1;
        const int pulse_span = (level == 0) ? pulses : 1;
        const int sample_span = (level == 2) ? 1 : samples;

        for (int f_idx = 0; f_idx < frames; ++f_idx) {
            for (int m = 0; m < n_tx; ++m) {
                for (int pp = 0; pp < pass_pulses; ++pp) {
                    for (int ss = 0; ss < pass_samples; ++ss) {
                        TracePass(radar, targets_manager, f_idx, m, pp,
                                  pulse_span, ss, sample_span, level,
                                  static_cast<double>(density), ray_filter,
                                  lam_min, fc, fs, gate, waveform);
                    }
                }
            }
        }
        return SUCCESS;
    }

private:
    static constexpr double kPi = 3.14159265358979323846;
    static constexpr double kC = 299792458.0;

    struct PoSample {
        double x, y, z;    // hit point
        double nx, ny, nz; // surface normal (facing the incoming ray)
        double dx, dy, dz; // incident ray direction at this bounce
        double range_tx;   // path length from the tx channel
        double area;       // footprint area
        double grad_mag;   // |tangential phase gradient| for Gordon weight
        std::complex<double> einc[3];  // transported incident E at this bounce
        std::complex<double> refl;
        int target_idx;
        int bounce;
    };

    void TracePass(const std::shared_ptr<Radar<H, L>> &radar,
                   const std::shared_ptr<TargetsManager<L>> &targets_manager,
                   int f_idx, int m, int p0, int p_span, int s0, int s_span,
                   int level, double density, rsv::Vec2<int_t> ray_filter,
                   double lam_min, double fc, double fs, double gate,
                   const rsim::Waveform<H> &waveform) {
        const auto tx = radar->tx_;
        const auto rx = radar->rx_;
        const int n_rx = static_cast<int>(rx->channels_.size());
        const int n_tx = static_cast<int>(tx->channels_.size());
        const int pulses = static_cast<int>(tx->pulse_start_time_.size());
        const int samples = radar->sample_size_;

        const double delay = static_cast<double>(tx->channels_[m].delay);
        const double u0 = gate + s0 / fs;
        const double frame_start = radar->frame_start_time_[f_idx];
        const double T_pass = frame_start + tx->pulse_start_time_[p0] + delay +
                              u0;

        // move targets to the pass instant and build the scene
        std::vector<rsim::SceneTri<L>> scene;
        std::vector<rsv::Vec3<L>> target_vel;
        BuildScene(targets_manager, T_pass, scene, target_vel);
        if (scene.empty()) {
            return;
        }

        rsv::Vec3<double> plat_loc;
        double rot[3][3];
        rsim::PlatformPose(radar->location_array_, radar->speed_,
                           radar->rotation_array_, radar->rotrate_, T_pass, 0,
                           plat_loc, rot);
        const rsv::Vec3<double> tx_pos =
            plat_loc + rsim::RotApply(rot, tx->channels_[m].location);

        std::vector<PoSample> po_samples;
        GenerateRays(scene, tx_pos, tx->channels_[m], lam_min, density,
                     ray_filter, po_samples);
        if (std::getenv("RSIM_DEBUG_MESH")) {
            double area_sum = 0.0;
            for (const auto &ps : po_samples) {
                area_sum += ps.area;
            }
            std::printf("[mesh] pass f%d tx%d: %zu PO samples, area sum %.4f\n",
                        f_idx, m, po_samples.size(), area_sum);
            int hist[16] = {0};
            double bh_area[16] = {0};
            for (const auto &ps : po_samples) {
                if (ps.bounce < 16) { ++hist[ps.bounce]; bh_area[ps.bounce] += ps.area; }
            }
            std::printf("[mesh]   bounce histogram:");
            for (int b = 1; b < 16; ++b)
                if (hist[b]) std::printf(" %d:%d(a=%.5f)", b, hist[b], bh_area[b]);
            std::printf("\n");
        }
        if (po_samples.empty()) {
            return;
        }

        for (int n = 0; n < n_rx; ++n) {
            const int ch = (f_idx * n_tx + m) * n_rx + n;
            const rsv::Vec3<double> rx_pos =
                plat_loc + rsim::RotApply(rot, rx->channels_[n].location);

            for (int p = p0; p < p0 + p_span; ++p) {
                const double f_off = tx->freq_offset_[p];
                const double lam = kC / (fc + f_off);
                const double k = 2.0 * kPi / lam;
                const double k_over_2 = kPi / lam;
                for (int s = s0; s < s0 + s_span; ++s) {
                    const size_t flat =
                        (static_cast<size_t>(ch) * pulses + p) * samples + s;
                    const double u = gate + s / fs;
                    const double T =
                        frame_start + tx->pulse_start_time_[p] + delay + u;
                    const double dt = T - T_pass;

                    std::complex<double> acc(0.0, 0.0);
                    for (const PoSample &ps : po_samples) {
                        const double dx = ps.x - rx_pos[0],
                                     dy = ps.y - rx_pos[1],
                                     dz = ps.z - rx_pos[2];
                        const double R_r =
                            std::sqrt(dx * dx + dy * dy + dz * dz);
                        if (R_r <= 0.0) {
                            continue;
                        }
                        double R_t = ps.range_tx;
                        if (dt != 0.0) {
                            R_t += RangeRate(ps, target_vel, tx_pos, rx_pos) *
                                   dt;
                        }
                        const double tau = (R_t + R_r) / kC;
                        const H beat =
                            (f_off * (u - gate) + waveform.Phase(u - gate)) -
                            (f_off * (u - tau) + waveform.Phase(u - tau));

                        const double g_db =
                            GainAt(tx->channels_[m], rx->channels_[n], rot, ps,
                                   tx_pos, rx_pos);

                        // K: the point-target amplitude-chain constant
                        // (full (4 pi)^3 form; the per-leg spreading lives in
                        // PoField's 1/(Rt Rr)); k_norm folds the field-to-RCS
                        // conversion of one footprint (2/sqrt(pi) at broadside)
                        const double pr_dbm =
                            static_cast<double>(tx->tx_power_) + g_db +
                            10.0 * std::log10(lam * lam /
                                              (4.0 * kPi * 4.0 * kPi * 4.0 *
                                               kPi)) +
                            static_cast<double>(rx->rf_gain_);
                        const double K =
                            std::sqrt(2.0 * 1e-3 *
                                      std::pow(10.0, pr_dbm / 10.0) *
                                      static_cast<double>(rx->resistor_)) *
                            std::pow(10.0,
                                     static_cast<double>(rx->baseband_gain_) /
                                         20.0) *
                            (2.0 / std::sqrt(kPi));

                        const std::complex<double> field =
                            PoField(ps, tx_pos, rx_pos,
                                    tx->channels_[m].polar,
                                    rx->channels_[n].polar, k);
                        std::complex<double> contrib =
                            K * k_over_2 * field * ps.refl *
                            std::complex<double>(0.0, 1.0) *  // +j prefactor
                            std::exp(std::complex<double>(0.0,
                                                          2.0 * kPi * beat));

                        // waveform modulation: ZOH at the echo transmit time
                        if (!tx->channels_[m].mod_var.empty()) {
                            const auto &mt = tx->channels_[m].mod_t;
                            const double step = mt[1] - mt[0];
                            const int idx =
                                static_cast<int>(std::floor(
                                    (u - tau - mt[0]) / step)) +
                                1;
                            const int nmod = static_cast<int>(
                                tx->channels_[m].mod_var.size());
                            const int widx = ((idx % nmod) + nmod) % nmod;
                            const auto &mv = tx->channels_[m].mod_var[widx];
                            contrib *= std::complex<double>(
                                static_cast<double>(mv.real()),
                                static_cast<double>(mv.imag()));
                        }
                        acc += contrib;
                    }

                    const std::complex<L> &pm = tx->channels_[m].pulse_mod[p];
                    acc *= std::complex<double>(static_cast<double>(pm.real()),
                                                static_cast<double>(pm.imag()));
                    radar->bb_real_[flat] += static_cast<H>(acc.real());
                    radar->bb_imag_[flat] += static_cast<H>(acc.imag());
                    if (std::getenv("RSIM_DEBUG_MESH") && p == p0 &&
                        s == s0 && n == 0) {
                        double per_b[16] = {0};
                        std::complex<double> coh[16] = {0.0};
                        for (const auto &ps : po_samples) {
                            if (ps.bounce < 16) {
                                const double dx = ps.x - rx_pos[0],
                                             dy = ps.y - rx_pos[1],
                                             dz = ps.z - rx_pos[2];
                                const double Rr = std::sqrt(dx*dx+dy*dy+dz*dz);
                                const double Rt = ps.range_tx;
                                const double tau = (Rt + Rr) / kC;
                                const H bt =
                                    (f_off * (u - gate) + waveform.Phase(u - gate)) -
                                    (f_off * (u - tau) + waveform.Phase(u - tau));
                                std::complex<double> fld =
                                    PoField(ps, tx_pos, rx_pos,
                                            tx->channels_[m].polar,
                                            rx->channels_[n].polar, k);
                                per_b[ps.bounce] += std::abs(fld * ps.refl);
                                coh[ps.bounce] += fld * ps.refl *
                                    std::exp(std::complex<double>(0.0, 2.0*kPi*bt));
                            }
                        }
                        std::printf("[mesh]   |field| per bounce:");
                        for (int b = 1; b < 16; ++b)
                            if (per_b[b] > 0) std::printf(" %d:%.3e", b, per_b[b]);
                        std::printf("   coherent:");
                        for (int b = 1; b < 16; ++b)
                            if (per_b[b] > 0)
                                std::printf(" %d:%.3e", b, std::abs(coh[b]));
                        std::printf("\n");
                    }
                }
            }
        }
    }

    // ---- scene assembly ----------------------------------------------------
    void BuildScene(
        const std::shared_ptr<TargetsManager<L>> &targets_manager,
        double T_pass, std::vector<rsim::SceneTri<L>> &scene,
        std::vector<rsv::Vec3<L>> &target_vel) {
        const auto &targets = targets_manager->targets();
        int ti = 0;
        for (const auto &tgt : targets) {
            // time-varying targets: the marshalling layer expanded motion onto
            // the timestamp grid; use the nearest entry
            int idx = 0;
            if (tgt->array_size_ > 1) {
                idx = static_cast<int>(T_pass);
                if (idx >= tgt->array_size_) {
                    idx = tgt->array_size_ - 1;
                }
            }
            tgt->Move(idx, T_pass);
            for (const auto &tri : tgt->vect_mesh_) {
                rsim::SceneTri<L> st;
                st.v[0] = tri.vertex_[0];
                st.v[1] = tri.vertex_[1];
                st.v[2] = tri.vertex_[2];
                const rsv::Vec3<L> e1 = st.v[1] - st.v[0];
                const rsv::Vec3<L> e2 = st.v[2] - st.v[0];
                st.normal = e1.Cross(e2);
                const L nl = std::sqrt(st.normal.Dot(st.normal));
                if (nl > L(0)) {
                    st.normal = st.normal * (L(1) / nl);
                }
                st.area = L(0.5) * nl;
                st.eps = tgt->permittivity_;
                st.mu = tgt->permeability_;
                st.skip_diffusion = tgt->skip_diffusion_;
                st.environment = tgt->environment_;
                st.target_idx = ti;
                scene.push_back(st);
            }
            target_vel.push_back(tgt->speed_array_[0]);
            ++ti;
        }
    }

    // ---- ray generation ----------------------------------------------------
    // Occupancy: the documented grid is phi in [-90, 90], theta in [0, 180],
    // cell n at fov0 + n*grid (benchmarks/occupancy_ab.py:126-137). Probed at
    // cell centers against the mesh, plus one-cell dilation so edge cells are
    // not missed.
    void GenerateRays(const std::vector<rsim::SceneTri<L>> &scene,
                      const rsv::Vec3<double> &tx_pos,
                      const typename Transmitter<H, L>::Channel &tx_ch,
                      double lam_min, double density,
                      rsv::Vec2<int_t> ray_filter,
                      std::vector<PoSample> &po_samples) {
        const double grid = static_cast<double>(tx_ch.grid);
        if (grid <= 0.0) {
            return;
        }

        // scene BVH
        std::vector<Triangle<L>> tris(scene.size());
        for (size_t i = 0; i < scene.size(); ++i) {
            tris[i].vertex_ = const_cast<rsv::Vec3<L> *>(scene[i].v);
        }
        rsim::Bvh<L> bvh(tris);

        const rsv::Vec3<L> org(static_cast<L>(tx_pos[0]),
                               static_cast<L>(tx_pos[1]),
                               static_cast<L>(tx_pos[2]));

        // probe pass: which cells hold geometry (with dilation)
        const int n_phi = static_cast<int>(180.0 / (grid * 180.0 / kPi)) + 1;
        const int n_theta = static_cast<int>(180.0 / (grid * 180.0 / kPi)) + 1;
        const double grid_deg = grid * 180.0 / kPi;
        std::vector<char> occupied(n_phi * n_theta, 0);
        for (int ip = 0; ip < n_phi; ++ip) {
            const double phi = (-90.0 + ip * grid_deg) * kPi / 180.0;
            for (int it = 0; it < n_theta; ++it) {
                const double theta = (0.0 + it * grid_deg) * kPi / 180.0;
                const rsv::Vec3<L> dir(
                    static_cast<L>(std::sin(theta) * std::cos(phi)),
                    static_cast<L>(std::sin(theta) * std::sin(phi)),
                    static_cast<L>(std::cos(theta)));
                typename rsim::Bvh<L>::Hit probe;
                if (bvh.ClosestHit(org, dir, probe)) {
                    occupied[ip * n_theta + it] = 1;
                }
            }
        }
        // dilate by one cell
        std::vector<char> dilated = occupied;
        for (int ip = 0; ip < n_phi; ++ip) {
            for (int it = 0; it < n_theta; ++it) {
                if (!occupied[ip * n_theta + it]) {
                    continue;
                }
                for (int dp = -1; dp <= 1; ++dp) {
                    for (int dt = -1; dt <= 1; ++dt) {
                        const int qp = ip + dp, qt = it + dt;
                        if (qp >= 0 && qp < n_phi && qt >= 0 && qt < n_theta) {
                            dilated[qp * n_theta + qt] = 1;
                        }
                    }
                }
            }
        }

        // fine rays per occupied cell
        for (int ip = 0; ip < n_phi; ++ip) {
            for (int it = 0; it < n_theta; ++it) {
                if (!dilated[ip * n_theta + it]) {
                    continue;
                }
                const double phi_c = (-90.0 + ip * grid_deg) * kPi / 180.0;
                const double theta_c = (0.0 + it * grid_deg) * kPi / 180.0;
                // probe range for the sizing formula
                const rsv::Vec3<L> pdir(
                    static_cast<L>(std::sin(theta_c) * std::cos(phi_c)),
                    static_cast<L>(std::sin(theta_c) * std::sin(phi_c)),
                    static_cast<L>(std::cos(theta_c)));
                typename rsim::Bvh<L>::Hit probe;
                double R_cell;
                if (bvh.ClosestHit(org, pdir, probe)) {
                    R_cell = static_cast<double>(probe.t);
                } else {
                    // dilated cell: use the scene's mean range
                    R_cell = SceneMeanRange(scene, org);
                }
                const double fine_step =
                    std::atan(lam_min / density / R_cell);
                const int count = static_cast<int>(grid / fine_step) + 1;
                const double phi_start = phi_c - 0.5 * grid;
                const double theta_start = theta_c - 0.5 * grid;
                for (int i = 0; i < count; ++i) {
                    for (int j = 0; j < count; ++j) {
                        TraceRay(scene, bvh, org, phi_start + i * fine_step,
                                 theta_start + j * fine_step, fine_step,
                                 ray_filter, po_samples, tx_ch.polar);
                    }
                }
            }
        }
    }

    static double SceneMeanRange(const std::vector<rsim::SceneTri<L>> &scene,
                                 const rsv::Vec3<L> &org) {
        double acc = 0.0;
        for (const auto &st : scene) {
            const double cx = (st.v[0][0] + st.v[1][0] + st.v[2][0]) / 3.0 -
                              org[0];
            const double cy = (st.v[0][1] + st.v[1][1] + st.v[2][1]) / 3.0 -
                              org[1];
            const double cz = (st.v[0][2] + st.v[1][2] + st.v[2][2]) / 3.0 -
                              org[2];
            acc += std::sqrt(cx * cx + cy * cy + cz * cz);
        }
        return acc / scene.size();
    }

    void TraceRay(const std::vector<rsim::SceneTri<L>> &scene,
                  const rsim::Bvh<L> &bvh, const rsv::Vec3<L> &origin,
                  double phi, double theta, double fine_step,
                  rsv::Vec2<int_t> ray_filter,
                  std::vector<PoSample> &po_samples,
                  const rsv::Vec3<std::complex<L>> &tx_pol) {
        rsv::Vec3<L> org = origin;
        rsv::Vec3<L> dir(static_cast<L>(std::sin(theta) * std::cos(phi)),
                         static_cast<L>(std::sin(theta) * std::sin(phi)),
                         static_cast<L>(std::cos(theta)));
        // transported E field: starts at the tx polarization (unit)
        std::complex<double> ef[3] = {std::complex<double>(tx_pol[0]),
                                      std::complex<double>(tx_pol[1]),
                                      std::complex<double>(tx_pol[2])};
        std::complex<double> refl(1.0, 0.0);
        double path = 0.0;
        const int max_b = ray_filter[1];
        for (int bounce = 1; bounce <= max_b; ++bounce) {
            typename rsim::Bvh<L>::Hit hit;
            if (!bvh.ClosestHit(org, dir, hit)) {
                break;
            }
            const rsim::SceneTri<L> &st = scene[hit.tri];
            org = org + dir * hit.t;
            path += hit.t;
            rsv::Vec3<L> n = st.normal;
            if (n.Dot(dir) > L(0)) {
                n = n * L(-1);
            }
            if (!st.skip_diffusion && bounce >= ray_filter[0] &&
                bounce <= ray_filter[1]) {
                PoSample ps;
                ps.x = org[0];
                ps.y = org[1];
                ps.z = org[2];
                ps.nx = n[0];
                ps.ny = n[1];
                ps.nz = n[2];
                ps.dx = dir[0];
                ps.dy = dir[1];
                ps.dz = dir[2];
                ps.range_tx = path;
                const double cos_inc =
                    std::fabs(static_cast<double>(n.Dot(dir)));
                // footprint on the surface: tube cross-section / cos_inc
                ps.area = fine_step * fine_step * path * path /
                          std::max(cos_inc, 1e-6);
                ps.grad_mag = 0.0;
                for (int a = 0; a < 3; ++a) {
                    ps.einc[a] = ef[a];
                }
                ps.refl = refl;
                ps.target_idx = st.target_idx;
                ps.bounce = bounce;
                po_samples.push_back(ps);
            }
            const L cos_i = std::fabs(n.Dot(dir));
            refl *= rsim::FresnelTE(st.eps, st.mu, cos_i);
            // PEC vector reflection of the field: E' = 2(E.n)n - E
            // (tangential flips, normal keeps sign; dielectric TE/TM split
            // is a refinement on top of this)
            const std::complex<double> edn =
                ef[0] * static_cast<double>(n[0]) +
                ef[1] * static_cast<double>(n[1]) +
                ef[2] * static_cast<double>(n[2]);
            for (int a = 0; a < 3; ++a) {
                ef[a] = 2.0 * edn * static_cast<double>(n[a]) - ef[a];
            }
            dir = rsim::Reflect(dir, n);
            org = org + n * L(1e-6);
        }
    }

    // ---- PO field of one sample -------------------------------------------
    // stuff = [o x (o x (n x (i x p_i)))] . conj(p_o), times footprint area,
    // times the Gordon sinc on the tangential phase gradient.
    std::complex<double> PoField(const PoSample &ps,
                                 const rsv::Vec3<double> &tx_pos,
                                 const rsv::Vec3<double> &rx_pos,
                                 const rsv::Vec3<std::complex<L>> &tx_pol,
                                 const rsv::Vec3<std::complex<L>> &rx_pol,
                                 double k) {
        const double Rr = std::sqrt((ps.x - rx_pos[0]) * (ps.x - rx_pos[0]) +
                                    (ps.y - rx_pos[1]) * (ps.y - rx_pos[1]) +
                                    (ps.z - rx_pos[2]) * (ps.z - rx_pos[2]));
        const double o[3] = {-(ps.x - rx_pos[0]) / Rr, -(ps.y - rx_pos[1]) / Rr,
                             -(ps.z - rx_pos[2]) / Rr};  // sample -> rx
        const double n[3] = {ps.nx, ps.ny, ps.nz};
        const double i[3] = {ps.dx, ps.dy, ps.dz};  // propagation direction
        (void)tx_pol;  // the incident field is the transported tx pol
        // PO surface current J = 2 n x H_inc, H_inc = i x E_inc
        const std::complex<double> hx[3] = {
            i[1] * ps.einc[2] - i[2] * ps.einc[1],
            i[2] * ps.einc[0] - i[0] * ps.einc[2],
            i[0] * ps.einc[1] - i[1] * ps.einc[0]};
        const std::complex<double> jj[3] = {n[1] * hx[2] - n[2] * hx[1],
                                            n[2] * hx[0] - n[0] * hx[2],
                                            n[0] * hx[1] - n[1] * hx[0]};
        // radiate the part of J perpendicular to o, project on conj(rx_pol)
        const std::complex<double> jdo =
            jj[0] * o[0] + jj[1] * o[1] + jj[2] * o[2];
        const std::complex<double> jp[3] = {jj[0] - jdo * o[0],
                                            jj[1] - jdo * o[1],
                                            jj[2] - jdo * o[2]};
        const std::complex<double> val =
            jp[0] * std::conj(std::complex<double>(rx_pol[0])) +
            jp[1] * std::conj(std::complex<double>(rx_pol[1])) +
            jp[2] * std::conj(std::complex<double>(rx_pol[2]));

        // Gordon footprint: the footprint is the ray tube projected on the
        // surface -- an ellipse stretched by 1/cos_inc along the incidence
        // plane. The phase ramp runs along the gradient, so the relevant
        // extent is w = R*step/cos_inc (= sqrt(area/cos_inc)); grazing
        // footprints are long across a steep ramp and suppress hard.
        const double g = GradMag(ps, tx_pos, rx_pos, k);
        const double cos_inc =
            std::fabs(ps.nx * ps.dx + ps.ny * ps.dy + ps.nz * ps.dz);
        // footprint extent along the phase gradient: the tube projected on
        // the surface stretches by 1/cos_inc in the incidence plane
        const double w = std::sqrt(ps.area) / std::max(cos_inc, 1e-6);
        const double sinc = Sinc(g * w / (2.0 * kPi));
        return val * ps.area * sinc * sinc / (ps.range_tx * Rr);
    }

    // |tangential phase gradient| of the footprint: the PO current ramps as
    // k*i.r under the local illumination and the observation picks up -k*o.r,
    // so the footprint pattern ramps as k*(i - o).r. For a first-bounce
    // monostatic sample this is k*(i - o) = 2k*i as usual.
    static double GradMag(const PoSample &ps, const rsv::Vec3<double> &tx_pos,
                          const rsv::Vec3<double> &rx_pos, double k) {
        (void)tx_pos;
        const double dr[3] = {ps.x - rx_pos[0], ps.y - rx_pos[1],
                              ps.z - rx_pos[2]};
        const double Rr = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] +
                                    dr[2] * dr[2]);
        const double o[3] = {-dr[0] / Rr, -dr[1] / Rr, -dr[2] / Rr};
        const double g[3] = {k * (ps.dx - o[0]), k * (ps.dy - o[1]),
                             k * (ps.dz - o[2])};
        const double n[3] = {ps.nx, ps.ny, ps.nz};
        const double gn = g[0] * n[0] + g[1] * n[1] + g[2] * n[2];
        const double gt[3] = {g[0] - gn * n[0], g[1] - gn * n[1],
                              g[2] - gn * n[2]};
        return std::sqrt(gt[0] * gt[0] + gt[1] * gt[1] + gt[2] * gt[2]);
    }

    static double Sinc(double x) {  // sin(pi x)/(pi x)
        if (std::fabs(x) < 1e-12) {
            return 1.0;
        }
        return std::sin(kPi * x) / (kPi * x);
    }

    double RangeRate(const PoSample &ps,
                     const std::vector<rsv::Vec3<L>> &target_vel,
                     const rsv::Vec3<double> &tx_pos,
                     const rsv::Vec3<double> &rx_pos) {
        if (ps.target_idx < 0 ||
            ps.target_idx >= static_cast<int>(target_vel.size())) {
            return 0.0;
        }
        const rsv::Vec3<L> &v = target_vel[ps.target_idx];
        if (v[0] == 0 && v[1] == 0 && v[2] == 0) {
            return 0.0;
        }
        // d(R_t + R_r)/dt = v . (unit tx->pt) + v . (unit rx->pt)
        const double dt[3] = {ps.x - tx_pos[0], ps.y - tx_pos[1],
                              ps.z - tx_pos[2]};
        const double Rt = std::sqrt(dt[0] * dt[0] + dt[1] * dt[1] +
                                    dt[2] * dt[2]);
        const double dr[3] = {ps.x - rx_pos[0], ps.y - rx_pos[1],
                              ps.z - rx_pos[2]};
        const double Rr = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] +
                                    dr[2] * dr[2]);
        return v[0] * (dt[0] / Rt + dr[0] / Rr) +
               v[1] * (dt[1] / Rt + dr[1] / Rr) + v[2] * (dt[2] / Rt + dr[2] / Rr);
    }

    static double GainAt(const typename Transmitter<H, L>::Channel &tx_ch,
                         const typename Receiver<L>::Channel &rx_ch,
                         const double rot[3][3], const PoSample &ps,
                         const rsv::Vec3<double> &tx_pos,
                         const rsv::Vec3<double> &rx_pos) {
        double az, th;
        const double d0 = ps.x - tx_pos[0], d1 = ps.y - tx_pos[1],
                     d2 = ps.z - tx_pos[2];
        rsim::ViewAngles(
            rsv::Vec3<double>(rot[0][0] * d0 + rot[1][0] * d1 +
                                  rot[2][0] * d2,
                              rot[0][1] * d0 + rot[1][1] * d1 +
                                  rot[2][1] * d2,
                              rot[0][2] * d0 + rot[1][2] * d1 +
                                  rot[2][2] * d2),
            az, th);
        const double g_tx = rsim::PatternGainDb(tx_ch, az, th);
        const double e0 = ps.x - rx_pos[0], e1 = ps.y - rx_pos[1],
                     e2 = ps.z - rx_pos[2];
        rsim::ViewAngles(
            rsv::Vec3<double>(rot[0][0] * e0 + rot[1][0] * e1 +
                                  rot[2][0] * e2,
                              rot[0][1] * e0 + rot[1][1] * e1 +
                                  rot[2][1] * e2,
                              rot[0][2] * e0 + rot[1][2] * e1 +
                                  rot[2][2] * e2),
            az, th);
        const double g_rx = rsim::PatternGainDb(rx_ch, az, th);
        return g_tx + g_rx;
    }
};
