// ==============================================================================
// radarsim-engine — simulator_interference.hpp
// Radar-to-radar interference (radarsimc.pxd:426-431). Writes into the victim
// radar's currently-pointed baseband buffers (the binding re-points them to
// the interference buffers before the call), OVERWRITING them.
//
// Model (docs/point_simulator_model.md):
//   - one-way Friis path with an extra 1/(4 pi): P_r = P_t G lambda^2 /
//     ((4 pi)^3 R^2), lambda at the interferer's fc + f_offset
//   - phase = phi_v(u_v - gate) - phi_i(u_i): the victim's (gated) LO minus
//     the interferer's emission-time phase
//   - hard passband gate: |f_v(u_v) - f_i(u_i)| > noise_bandwidth -> nothing
//   - the interferer's waveform/pulse modulation applies; its pulses are
//     scanned and rejected when the emission falls outside [t0, t_end]
// ==============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "core/enums.hpp"
#include "radar.hpp"
#include "rf/antenna.hpp"
#include "rf/platform.hpp"
#include "rf/waveform.hpp"

template <typename H, typename L, typename ExecutionPolicy>
class InterferenceSimulator {
public:
    InterferenceSimulator() = default;

    RadarSimErrorCode Run(const std::shared_ptr<Radar<H, L>> &radar,
                          const std::shared_ptr<Radar<H, L>> &interf_radar) {
        if (!radar || !interf_radar) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (radar->bb_real_ == nullptr || radar->bb_imag_ == nullptr) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }

        const auto tx_v = radar->tx_;
        const auto rx_v = radar->rx_;
        const auto tx_i = interf_radar->tx_;

        const int n_tx_v = static_cast<int>(tx_v->channels_.size());
        const int n_rx = static_cast<int>(rx_v->channels_.size());
        const int n_tx_i = static_cast<int>(tx_i->channels_.size());
        const int n_ch = n_tx_v * n_rx;
        const int pulses_v = static_cast<int>(tx_v->pulse_start_time_.size());
        const int pulses_i = static_cast<int>(tx_i->pulse_start_time_.size());
        const int samples = radar->sample_size_;
        const int frames = static_cast<int>(radar->frame_start_time_.size());
        if (n_ch == 0 || pulses_v == 0 || pulses_i == 0 || samples <= 0) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }

        const double fs = static_cast<double>(rx_v->fs_);
        const double gate = rx_v->gate_delay_;
        const double noise_bw = static_cast<double>(rx_v->baseband_bw_);
        const double resistor = static_cast<double>(rx_v->resistor_);
        const double rf_gain_db = static_cast<double>(rx_v->rf_gain_);
        const double bb_gain_db = static_cast<double>(rx_v->baseband_gain_);

        const rsim::Waveform<H> wave_v(tx_v->freq_, tx_v->freq_time_);
        const rsim::Waveform<H> wave_i(tx_i->freq_, tx_i->freq_time_);
        const double t_i0 = static_cast<double>(wave_i.t_begin());
        const double t_i1 = static_cast<double>(wave_i.t_end());
        const double fc_i = static_cast<double>(wave_i.CenterFreq());

        const int64_t total =
            static_cast<int64_t>(frames) * n_ch * pulses_v * samples;
        // the interference buffers are caller-allocated but uninitialized:
        // overwrite, do not accumulate
        std::fill(radar->bb_real_, radar->bb_real_ + total, H(0));
        std::fill(radar->bb_imag_, radar->bb_imag_ + total, H(0));

        const auto &i_frames = interf_radar->frame_start_time_;
        const bool tv_plat_v = radar->location_array_.size() > 1;
        const bool tv_plat_i = interf_radar->location_array_.size() > 1;

        for (int ch = 0; ch < frames * n_ch; ++ch) {
            const int f_idx = ch / n_ch;
            const int m = (ch % n_ch) / n_rx;
            const int n = ch % n_rx;
            const auto &v_tx_ch = tx_v->channels_[m];
            const auto &v_rx_ch = rx_v->channels_[n];
            const double delay_v = static_cast<double>(v_tx_ch.delay);

            for (int mi = 0; mi < n_tx_i; ++mi) {
                const auto &i_ch = tx_i->channels_[mi];

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                for (int p = 0; p < pulses_v; ++p) {
                    for (int s = 0; s < samples; ++s) {
                        const size_t flat =
                            (static_cast<size_t>(ch) * pulses_v + p) * samples +
                            s;
                        const double u_v = gate + s / fs;
                        const double T = radar->frame_start_time_[f_idx] +
                                         tx_v->pulse_start_time_[p] + delay_v +
                                         u_v;

                        rsv::Vec3<double> plat_v;
                        double rot_v[3][3];
                        rsim::PlatformPose(radar->location_array_,
                                           radar->speed_,
                                           radar->rotation_array_,
                                           radar->rotrate_, T,
                                           tv_plat_v ? flat : 0, plat_v, rot_v);
                        const rsv::Vec3<double> rx_pos =
                            plat_v + rsim::RotApply(rot_v, v_rx_ch.location);

                        for (int q = 0; q < pulses_i; ++q) {
                            rsv::Vec3<double> plat_i;
                            double rot_i[3][3];
                            rsim::PlatformPose(
                                interf_radar->location_array_,
                                interf_radar->speed_,
                                interf_radar->rotation_array_,
                                interf_radar->rotrate_, T,
                                tv_plat_i ? flat : 0, plat_i, rot_i);
                            const rsv::Vec3<double> tx_pos_i =
                                plat_i +
                                rsim::RotApply(rot_i, i_ch.location);

                            const double dx = tx_pos_i[0] - rx_pos[0];
                            const double dy = tx_pos_i[1] - rx_pos[1];
                            const double dz = tx_pos_i[2] - rx_pos[2];
                            const double R =
                                std::sqrt(dx * dx + dy * dy + dz * dz);
                            if (R <= 0.0) {
                                continue;
                            }
                            const double tau = R / 299792458.0;
                            const double f_start_i =
                                i_frames[std::min(f_idx,
                                                  static_cast<int>(
                                                      i_frames.size()) - 1)];
                            const double u_i =
                                T - tau - f_start_i -
                                tx_i->pulse_start_time_[q] -
                                static_cast<double>(i_ch.delay);
                            if (u_i < t_i0 || u_i > t_i1) {
                                continue;
                            }
                            const double foff_v = tx_v->freq_offset_[p];
                            const double foff_i = tx_i->freq_offset_[q];
                            const double beat_f =
                                (static_cast<double>(wave_v.Freq(u_v - gate)) +
                                 foff_v) -
                                (static_cast<double>(wave_i.Freq(u_i)) +
                                 foff_i);
                            if (std::fabs(beat_f) > noise_bw) {
                                continue;
                            }

                            // antenna gains: interferer toward victim, victim
                            // toward interferer, each in its own body frame
                            double az, th;
                            rsim::ViewAngles(
                                rsv::Vec3<double>(
                                    -(rot_i[0][0] * dx + rot_i[1][0] * dy +
                                      rot_i[2][0] * dz),
                                    -(rot_i[0][1] * dx + rot_i[1][1] * dy +
                                      rot_i[2][1] * dz),
                                    -(rot_i[0][2] * dx + rot_i[1][2] * dy +
                                      rot_i[2][2] * dz)),
                                az, th);
                            const double g_i =
                                rsim::PatternGainDb(i_ch, az, th);
                            rsim::ViewAngles(
                                rsv::Vec3<double>(
                                    rot_v[0][0] * dx + rot_v[1][0] * dy +
                                        rot_v[2][0] * dz,
                                    rot_v[0][1] * dx + rot_v[1][1] * dy +
                                        rot_v[2][1] * dz,
                                    rot_v[0][2] * dx + rot_v[1][2] * dy +
                                        rot_v[2][2] * dz),
                                az, th);
                            const double g_r =
                                rsim::PatternGainDb(v_rx_ch, az, th);

                            const double lambda = 299792458.0 / (fc_i + foff_i);
                            const double p_w =
                                1e-3 *
                                std::pow(10.0,
                                         (static_cast<double>(
                                              tx_i->tx_power_) +
                                          g_i + g_r) /
                                             10.0) *
                                lambda * lambda /
                                (4.0 * kPi * 4.0 * kPi * 4.0 * kPi * R * R);
                            const double amp =
                                std::sqrt(2.0 * p_w * resistor) *
                                std::pow(10.0, rf_gain_db / 20.0) *
                                std::pow(10.0, bb_gain_db / 20.0);

                            const H ph =
                                (foff_v * (u_v - gate) +
                                 wave_v.Phase(u_v - gate)) -
                                (foff_i * u_i + wave_i.Phase(u_i));

                            const double pol = std::abs(
                                std::conj(v_rx_ch.polar[0]) * i_ch.polar[0] +
                                std::conj(v_rx_ch.polar[1]) * i_ch.polar[1] +
                                std::conj(v_rx_ch.polar[2]) * i_ch.polar[2]);

                            double re =
                                amp * pol * std::cos(2.0 * kPi * ph);
                            double im =
                                amp * pol * std::sin(2.0 * kPi * ph);

                            // interferer's waveform modulation at emission
                            if (!i_ch.mod_var.empty()) {
                                const auto &mt = i_ch.mod_t;
                                const double step = mt[1] - mt[0];
                                const int idx =
                                    static_cast<int>(std::floor(
                                        (u_i - mt[0]) / step)) +
                                    1;
                                const int nmod =
                                    static_cast<int>(i_ch.mod_var.size());
                                const int widx = ((idx % nmod) + nmod) % nmod;
                                const double mre = static_cast<double>(
                                    i_ch.mod_var[widx].real());
                                const double mim = static_cast<double>(
                                    i_ch.mod_var[widx].imag());
                                const double r2 = re * mre - im * mim;
                                im = re * mim + im * mre;
                                re = r2;
                            }
                            // interferer's pulse modulation
                            {
                                const std::complex<L> &pm = i_ch.pulse_mod[q];
                                const double pre =
                                    static_cast<double>(pm.real());
                                const double pim =
                                    static_cast<double>(pm.imag());
                                const double r2 = re * pre - im * pim;
                                im = re * pim + im * pre;
                                re = r2;
                            }

                            radar->bb_real_[flat] += static_cast<H>(re);
                            radar->bb_imag_[flat] += static_cast<H>(im);
                        }
                    }
                }
            }
        }
        return SUCCESS;
    }

private:
    static constexpr double kPi = 3.14159265358979323846;
};
