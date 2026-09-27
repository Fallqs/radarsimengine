// ==============================================================================
// radarsim-engine — simulator_noise.hpp
// Receiver thermal-noise simulator (radarsimc.pxd:435-448).
//
// Contract (tests/test_noise_simulation.py):
//   - noise is keyed by (seed, rx channel index, absolute timestamp): virtual
//     channels sharing an Rx and bit-identical timestamps get identical noise
//   - complex baseband: level * (n_re + j n_im) / sqrt(2); real: level * n_re
//   - deterministic: the binding always passes seed 0
//
// The timestamps argument is the ORIGIN timestamp grid [M*N, pulses, samples]
// (pre-frame); frame f adds radar->frame_start_time_[f]. The hash uses the
// float64 bit pattern of the absolute timestamp, so identical times produce
// identical noise bit-for-bit.
// ==============================================================================
#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>

#include "core/enums.hpp"
#include "radar.hpp"

template <typename H, typename L, typename ExecutionPolicy>
class NoiseSimulator {
public:
    NoiseSimulator() = default;

    RadarSimErrorCode Run(const std::shared_ptr<Radar<H, L>> &radar,
                          H noise_level, bool is_complex, const H *timestamps,
                          int ts_channel_size, int ts_pulse_size,
                          int ts_sample_size, H *noise_real, H *noise_imag,
                          unsigned long long seed) {
        if (!radar || !timestamps || !noise_real || !noise_imag) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (ts_channel_size <= 0 || ts_pulse_size <= 0 || ts_sample_size <= 0) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }

        const int frames =
            static_cast<int>(radar->frame_start_time_.size());
        const int n_rx =
            static_cast<int>(radar->rx_->channels_.size());
        if (n_rx <= 0) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        const double level = static_cast<double>(noise_level);
        const double inv_sqrt2 = 0.7071067811865475244;

        const int64_t grid =
            static_cast<int64_t>(ts_channel_size) * ts_pulse_size *
            ts_sample_size;

        for (int f = 0; f < frames; ++f) {
            const double frame_t =
                static_cast<double>(radar->frame_start_time_[f]);
            for (int ch = 0; ch < ts_channel_size; ++ch) {
                const int n = ch % n_rx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                for (int64_t idx = 0; idx <
                                       static_cast<int64_t>(ts_pulse_size) *
                                           ts_sample_size;
                     ++idx) {
                    const int64_t g =
                        (static_cast<int64_t>(f) * ts_channel_size + ch) *
                            (static_cast<int64_t>(ts_pulse_size) *
                             ts_sample_size) +
                        idx;
                    const double ts_abs =
                        static_cast<double>(timestamps[ch * ts_pulse_size *
                                                           ts_sample_size +
                                                       idx]) +
                        frame_t;
                    double nre, nim;
                    GaussianPair(seed, n, ts_abs, nre, nim);
                    if (is_complex) {
                        noise_real[g] = static_cast<H>(level * nre * inv_sqrt2);
                        noise_imag[g] = static_cast<H>(level * nim * inv_sqrt2);
                    } else {
                        noise_real[g] = static_cast<H>(level * nre);
                        noise_imag[g] = H(0);
                    }
                }
            }
        }
        (void)grid;
        return SUCCESS;
    }

private:
    static uint64_t Splitmix64(uint64_t x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

    // Deterministic Gaussian pair keyed by (seed, rx index, timestamp bits).
    // Must stay bit-compatible with tools/pysim/reference.py::_noise_normals.
    static void GaussianPair(unsigned long long seed, int rx_idx, double ts,
                             double &nre, double &nim) {
        uint64_t ts_bits;
        std::memcpy(&ts_bits, &ts, 8);
        uint64_t h = Splitmix64(seed ^ ts_bits);
        h = Splitmix64(h ^ (static_cast<uint64_t>(rx_idx + 1) *
                            0x9E3779B97F4A7C15ULL));
        const uint64_t h2 = Splitmix64(h);
        double u1 = static_cast<double>(h >> 11) /
                    static_cast<double>(1ULL << 53);
        if (u1 < 1e-300) {
            u1 = 1e-300;
        }
        const double u2 =
            static_cast<double>(h2 >> 11) / static_cast<double>(1ULL << 53);
        const double r = std::sqrt(-2.0 * std::log(u1));
        const double ang = 2.0 * 3.14159265358979323846 * u2;
        nre = r * std::cos(ang);
        nim = r * std::sin(ang);
    }
};
