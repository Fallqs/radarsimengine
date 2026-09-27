// ==============================================================================
// radarsim-engine — rf/phase_noise.hpp
// SSB phase-noise generation, replicating radar.py's cal_phase_noise:
// log-linear mask interpolation -> shaped deterministic spectrum -> IFFT ->
// phi_pn = Re(x_t). With pn_validation the AWGN is exactly (1+1j)/sqrt(2);
// otherwise a splitmix64-driven Gaussian stream from pn_seed.
//
// The IFFT is a direct O(K^2) DFT: K = 2M-2 equals the per-pulse sample
// count, and generation happens once per Run. A mixed-radix FFT is a Phase-4+
// optimization.
// ==============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <utility>
#include <vector>

namespace rsim {

constexpr double kPi = 3.14159265358979323846;

namespace pn_detail {

inline uint64_t Splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Log-scale SSB mask interpolation (mirrors _interpolate_phase_noise_power,
// including the realmin offsets).
inline std::vector<double> InterpolatePnPower(const std::vector<double> &freq,
                                              const std::vector<double> &power,
                                              const std::vector<double> &f_grid) {
    constexpr double realmin = 2.2250738585072014e-308;  // nextafter(0,1)-scale
    const size_t m = f_grid.size();
    std::vector<double> log_p(m, 0.0);
    const size_t n = freq.size();
    for (size_t i = 0; i < n; ++i) {
        const double left = freq[i];
        const double t1 = power[i];
        double right, t2;
        bool last_incl;
        if (i + 1 == n) {
            right = f_grid.back() * 2.0;
            t2 = power.back();
            last_incl = true;
        } else {
            right = freq[i + 1];
            t2 = power[i + 1];
            last_incl = false;
        }
        const double denom =
            std::log10(right + 2.0 * realmin) - std::log10(left + realmin);
        for (size_t j = 0; j < m; ++j) {
            const bool inside =
                last_incl ? (f_grid[j] >= left && f_grid[j] <= right)
                          : (f_grid[j] >= left && f_grid[j] < right);
            if (inside) {
                log_p[j] = t1 +
                           (std::log10(f_grid[j] + realmin) -
                            std::log10(left + realmin)) /
                               denom * (t2 - t1);
            }
        }
    }
    std::vector<double> p(m);
    for (size_t j = 0; j < m; ++j) {
        p[j] = std::pow(10.0, log_p[j] / 10.0);
    }
    return p;
}

}  // namespace pn_detail

// Phase-noise phase deviation phi_pn[sample] (radians) for one pulse.
// Returns an empty vector when no phase noise is configured.
//   num_samples: samples per pulse; fs: sample rate (pn_fs).
inline std::vector<double> GeneratePhaseNoise(
    const std::vector<double> &pn_freq_in, const std::vector<double> &pn_pow_in,
    double fs, int num_samples, unsigned long long seed, bool validation) {
    if (pn_freq_in.empty()) {
        return {};
    }

    // sort by frequency, drop >= fs/2, prepend 0 dBc/Hz at DC
    std::vector<std::pair<double, double>> pairs;
    for (size_t i = 0; i < pn_freq_in.size(); ++i) {
        pairs.emplace_back(pn_freq_in[i], pn_pow_in[i]);
    }
    std::sort(pairs.begin(), pairs.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    std::vector<double> freq, power;
    for (const auto &pr : pairs) {
        if (pr.first < fs / 2) {
            freq.push_back(pr.first);
            power.push_back(pr.second);
        }
    }
    if (std::find(freq.begin(), freq.end(), 0.0) == freq.end()) {
        freq.insert(freq.begin(), 0.0);
        power.insert(power.begin(), 0.0);
    }

    const int num_f =
        (num_samples % 2) ? (num_samples + 1) / 2 + 1 : num_samples / 2 + 1;
    std::vector<double> f_grid(num_f);
    for (int i = 0; i < num_f; ++i) {
        f_grid[i] = fs / 2 * i / (num_f - 1);
    }
    std::vector<double> delta_f(num_f);
    for (int i = 0; i + 1 < num_f; ++i) {
        delta_f[i] = f_grid[i + 1] - f_grid[i];
    }
    delta_f[num_f - 1] = f_grid[num_f - 1] - f_grid[num_f - 2];

    const std::vector<double> p_interp =
        pn_detail::InterpolatePnPower(freq, power, f_grid);

    // shaped spectrum
    uint64_t state = seed;
    std::vector<std::complex<double>> spec(num_f);
    for (int i = 0; i < num_f; ++i) {
        std::complex<double> w;
        if (validation) {
            w = std::complex<double>(std::sqrt(0.5), std::sqrt(0.5));
        } else {
            state = pn_detail::Splitmix64(state);
            const double u1 =
                std::max((state >> 11) / double(1ULL << 53), 1e-300);
            state = pn_detail::Splitmix64(state);
            const double u2 = (state >> 11) / double(1ULL << 53);
            const double r = std::sqrt(-2.0 * std::log(u1));
            w = std::sqrt(0.5) * std::complex<double>(r * std::cos(2.0 * kPi * u2),
                                                      r * std::sin(2.0 * kPi * u2));
        }
        spec[i] = num_f * std::sqrt(delta_f[i] * p_interp[i]) * w;
    }

    // full spectrum [0, fs): symmetric negative part, DC removed
    const int k_len = 2 * num_f - 2;
    std::vector<std::complex<double>> full(k_len, {0.0, 0.0});
    for (int i = 0; i < num_f; ++i) {
        full[i] = spec[i];
    }
    for (int i = 0; i < num_f - 2; ++i) {
        full[num_f + i] = std::conj(spec[num_f - 2 - i]);
    }
    full[0] = {0.0, 0.0};

    // inverse DFT (direct; K = num_samples per pulse)
    std::vector<double> phi(num_samples);
    for (int n = 0; n < num_samples && n < k_len; ++n) {
        std::complex<double> acc(0.0, 0.0);
        for (int k = 0; k < k_len; ++k) {
            const double ang = 2.0 * kPi * k * n / k_len;
            acc += full[k] * std::complex<double>(std::cos(ang), std::sin(ang));
        }
        phi[n] = acc.real() / k_len;
    }
    return phi;
}

}  // namespace rsim
