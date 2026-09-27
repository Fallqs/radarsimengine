// ==============================================================================
// radarsim-engine — simulator_noise.hpp
// Receiver thermal-noise simulator (radarsimc.pxd:435-448).
//
// Determinism is contractual: the binding always passes seed 0
// (simulator_radar.pyx:403,418), so output must be a pure function of
// (shape, noise_level, is_complex, seed).
// Phase 0: fills the output buffers deterministically (zeros); the Gaussian
// generator lands with Phase 1 (docs/04-physics.md §4.3).
// ==============================================================================
#pragma once

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
        if (ts_channel_size <= 0 || ts_pulse_size <= 0 ||
            ts_sample_size <= 0) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        (void)noise_level;
        (void)is_complex;
        (void)seed;
        const long long total = static_cast<long long>(ts_channel_size) *
                                ts_pulse_size * ts_sample_size;
        for (long long i = 0; i < total; ++i) {
            noise_real[i] = H(0);
            noise_imag[i] = H(0);
        }
        return SUCCESS;
    }
};
