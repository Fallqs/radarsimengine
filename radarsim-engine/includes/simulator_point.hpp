// ==============================================================================
// radarsim-engine — simulator_point.hpp
// Ideal point-target simulator (radarsimc.pxd:372-377).
//
// Model (validated against the upstream ideal suite, 20/21; full derivation
// in docs/point_simulator_model.md):
//   per (frame f, tx m, rx n, pulse p, sample s), ch = (f*M + m)*N + n:
//     T   = frame_start[f] + pulse_start[p] + tx_delay[m] + gate + s/fs
//     u   = gate + s/fs                      (pulse-local; tx delay cancels)
//     pos = target position at T             (float-staged geometry)
//     tau = (|pos - tx_pos| + |pos - rx_pos|) / c
//     beat = phi_p(u - gate) - phi_p(u - tau)  (cycles; phi integrates
//            f_offset[p] + f(t))
//     amp  = radar-equation voltage with per-leg 4 pi R^2 spreading,
//            lambda at (fc + f_offset[p]), nearest-entry pattern gains,
//            |vdot(rx_pol, tx_pol)|
//     bb  += amp * mod_zoh(u - tau) * pulse_mod[p] * exp(j(2 pi beat + phs))
//            [* pn_factor(u, tau) when phase noise is configured]
// The point simulator runs first and OVERWRITES the baseband buffers.
// ==============================================================================
#pragma once

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "core/enums.hpp"
#include "points_manager.hpp"
#include "radar.hpp"
#include "rf/antenna.hpp"
#include "rf/phase_noise.hpp"
#include "rf/platform.hpp"
#include "rf/waveform.hpp"

template <typename H, typename L, typename ExecutionPolicy>
class PointSimulator {
public:
    PointSimulator() = default;

    RadarSimErrorCode Run(
        const std::shared_ptr<Radar<H, L>> &radar,
        const std::shared_ptr<PointsManager<L>> &points_manager) {
        if (!radar || !points_manager) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (radar->bb_real_ == nullptr || radar->bb_imag_ == nullptr) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }

        const auto tx = radar->tx_;
        const auto rx = radar->rx_;
        const int n_tx = static_cast<int>(tx->channels_.size());
        const int n_rx = static_cast<int>(rx->channels_.size());
        const int n_ch = n_tx * n_rx;
        const int pulses = static_cast<int>(tx->pulse_start_time_.size());
        const int samples = radar->sample_size_;
        const int frames = static_cast<int>(radar->frame_start_time_.size());
        if (n_tx == 0 || n_rx == 0 || pulses == 0 || samples <= 0) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }

        const double fs = static_cast<double>(rx->fs_);
        const double gate = radar->rx_->gate_delay_;  // double by contract
        const double tx_power_dbm = static_cast<double>(tx->tx_power_);
        const double rf_gain_db = static_cast<double>(rx->rf_gain_);
        const double bb_gain_db = static_cast<double>(rx->baseband_gain_);
        const double resistor = static_cast<double>(rx->resistor_);

        const rsim::Waveform<H> waveform(tx->freq_, tx->freq_time_);
        const double fc = static_cast<double>(waveform.CenterFreq());

        // phase noise phase deviation per fast-time sample (empty if none)
        std::vector<double> phi_pn;
        if (tx->has_phase_noise_) {
            phi_pn = rsim::GeneratePhaseNoise(
                std::vector<double>(tx->pn_freq_.begin(), tx->pn_freq_.end()),
                std::vector<double>(tx->pn_power_.begin(), tx->pn_power_.end()),
                static_cast<double>(tx->pn_fs_), samples, tx->pn_seed_,
                tx->pn_validation_);
        }

        const int64_t total =
            static_cast<int64_t>(frames) * n_ch * pulses * samples;

        // point simulator runs first: overwrite, do not accumulate
        std::fill(radar->bb_real_, radar->bb_real_ + total, H(0));
        std::fill(radar->bb_imag_, radar->bb_imag_ + total, H(0));

        const auto &points = points_manager->points();
        const bool time_varying_platform = radar->location_array_.size() > 1;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ch = 0; ch < frames * n_ch; ++ch) {
            const int f_idx = ch / n_ch;
            const int m = (ch % n_ch) / n_rx;
            const int n = ch % n_rx;
            const auto &tx_ch = tx->channels_[m];
            const auto &rx_ch = rx->channels_[n];

            const double delay = static_cast<double>(tx_ch.delay);
            // polarization projection |vdot(rx_pol, tx_pol)|
            const double pol = std::abs(
                std::conj(rx_ch.polar[0]) * tx_ch.polar[0] +
                std::conj(rx_ch.polar[1]) * tx_ch.polar[1] +
                std::conj(rx_ch.polar[2]) * tx_ch.polar[2]);

            for (int p = 0; p < pulses; ++p) {
                const double f_off = tx->freq_offset_[p];
                const double lambda =
                    299792458.0 / (fc + f_off);
                const double lambda_term =
                    10.0 * std::log10(lambda * lambda / (4.0 * kPi));

                for (int s = 0; s < samples; ++s) {
                    const size_t flat =
                        (static_cast<size_t>(ch) * pulses + p) * samples + s;
                    const double u = gate + s / fs;  // pulse-local waveform time
                    const double T = radar->frame_start_time_[f_idx] +
                                     tx->pulse_start_time_[p] + delay + u;

                    rsv::Vec3<double> plat_loc;
                    double rot[3][3];
                    rsim::PlatformPose(radar->location_array_, radar->speed_,
                                       radar->rotation_array_, radar->rotrate_,
                                       T, time_varying_platform ? flat : 0,
                                       plat_loc, rot);

                    const rsv::Vec3<double> tx_pos =
                        plat_loc + rsim::RotApply(rot, tx_ch.location);
                    const rsv::Vec3<double> rx_pos =
                        plat_loc + rsim::RotApply(rot, rx_ch.location);

                    for (const auto &pt : points) {
                        const bool tv = pt.location_array.size() > 1;
                        const rsv::Vec3<L> &pl = tv ? pt.location_array[flat]
                                                    : pt.location_array[0];
                        const double pos_x = static_cast<double>(pl[0]) +
                                             pt.speed[0] * T;
                        const double pos_y = static_cast<double>(pl[1]) +
                                             pt.speed[1] * T;
                        const double pos_z = static_cast<double>(pl[2]) +
                                             pt.speed[2] * T;
                        const double rcs_db = static_cast<double>(
                            pt.rcs_array.size() > 1 ? pt.rcs_array[flat]
                                                    : pt.rcs_array[0]);
                        const double phs = static_cast<double>(
                            pt.phase_array.size() > 1 ? pt.phase_array[flat]
                                                      : pt.phase_array[0]);

                        const double dtx[3] = {pos_x - tx_pos[0],
                                               pos_y - tx_pos[1],
                                               pos_z - tx_pos[2]};
                        const double drx[3] = {pos_x - rx_pos[0],
                                               pos_y - rx_pos[1],
                                               pos_z - rx_pos[2]};
                        const double R_tx =
                            std::sqrt(dtx[0] * dtx[0] + dtx[1] * dtx[1] +
                                      dtx[2] * dtx[2]);
                        const double R_rx =
                            std::sqrt(drx[0] * drx[0] + drx[1] * drx[1] +
                                      drx[2] * drx[2]);
                        if (R_tx <= 0.0 || R_rx <= 0.0) {
                            continue;
                        }
                        const double tau = (R_tx + R_rx) / 299792458.0;

                        // antenna gains in each channel's body frame
                        double az, th;
                        rsim::ViewAngles(
                            rsv::Vec3<double>(
                                rot[0][0] * dtx[0] + rot[1][0] * dtx[1] +
                                    rot[2][0] * dtx[2],
                                rot[0][1] * dtx[0] + rot[1][1] * dtx[1] +
                                    rot[2][1] * dtx[2],
                                rot[0][2] * dtx[0] + rot[1][2] * dtx[1] +
                                    rot[2][2] * dtx[2]),
                            az, th);
                        const double g_tx =
                            rsim::PatternGainDb(tx_ch, az, th);
                        rsim::ViewAngles(
                            rsv::Vec3<double>(
                                rot[0][0] * drx[0] + rot[1][0] * drx[1] +
                                    rot[2][0] * drx[2],
                                rot[0][1] * drx[0] + rot[1][1] * drx[1] +
                                    rot[2][1] * drx[2],
                                rot[0][2] * drx[0] + rot[1][2] * drx[1] +
                                    rot[2][2] * drx[2]),
                            az, th);
                        const double g_rx =
                            rsim::PatternGainDb(rx_ch, az, th);

                        const double pr_dbm =
                            tx_power_dbm + g_tx + g_rx -
                            10.0 * std::log10(4.0 * kPi * R_tx * R_tx) +
                            rcs_db -
                            10.0 * std::log10(4.0 * kPi * R_rx * R_rx) +
                            lambda_term + rf_gain_db;
                        const double amp =
                            std::sqrt(1e-3 * std::pow(10.0, pr_dbm / 10.0) *
                                      resistor) *
                            std::pow(10.0, bb_gain_db / 20.0) * kSqrt2;

                        // beat phase: deramp reference delayed by the gate
                        const H beat =
                            (f_off * (u - gate) + waveform.Phase(u - gate)) -
                            (f_off * (u - tau) + waveform.Phase(u - tau));

                        double re = amp * pol *
                                    std::cos(2.0 * kPi * beat + phs);
                        double im = amp * pol *
                                    std::sin(2.0 * kPi * beat + phs);

                        // waveform modulation: ZOH at echo transmit time
                        if (!tx_ch.mod_var.empty()) {
                            const auto &mt = tx_ch.mod_t;
                            const double step = mt[1] - mt[0];
                            const int idx =
                                static_cast<int>(
                                    std::floor((u - tau - mt[0]) / step)) +
                                1;
                            const int nmod = static_cast<int>(tx_ch.mod_var.size());
                            const int widx = ((idx % nmod) + nmod) % nmod;
                            const std::complex<L> &mv = tx_ch.mod_var[widx];
                            const double mre = static_cast<double>(mv.real());
                            const double mim = static_cast<double>(mv.imag());
                            const double r2 = re * mre - im * mim;
                            im = re * mim + im * mre;
                            re = r2;
                        }
                        // pulse modulation
                        {
                            const std::complex<L> &pm = tx_ch.pulse_mod[p];
                            const double pre = static_cast<double>(pm.real());
                            const double pim = static_cast<double>(pm.imag());
                            const double r2 = re * pre - im * pim;
                            im = re * pim + im * pre;
                            re = r2;
                        }
                        // phase noise: echo carries oscillator phase at
                        // emission, LO carries it now
                        if (!phi_pn.empty()) {
                            const double f_hi = u * fs;
                            const double f_lo = (u - tau) * fs;
                            const double dphi =
                                PnAt(phi_pn, f_hi) - PnAt(phi_pn, f_lo);
                            const double c = std::cos(dphi), sn = std::sin(dphi);
                            const double r2 = re * c - im * sn;
                            im = re * sn + im * c;
                            re = r2;
                        }

                        radar->bb_real_[flat] += static_cast<H>(re);
                        radar->bb_imag_[flat] += static_cast<H>(im);
                    }
                }
            }
        }
        return SUCCESS;
    }

private:
    static constexpr double kPi = 3.14159265358979323846;
    static constexpr double kSqrt2 = 1.4142135623730950488;

    // Phase-noise phase deviation at fractional sample index (wraps).
    static double PnAt(const std::vector<double> &phi, double idx) {
        const int n = static_cast<int>(phi.size());
        const double fl = std::floor(idx);
        const int i0 = static_cast<int>(fl);
        const double frac = idx - fl;
        const int w0 = ((i0 % n) + n) % n;
        const int w1 = (w0 + 1) % n;
        return phi[w0] * (1.0 - frac) + phi[w1] * frac;
    }
};
