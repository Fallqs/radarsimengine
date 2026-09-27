// ==============================================================================
// radarsim-engine — rf/waveform.hpp
// Sampled arbitrary waveform: frequency f(t) given at breakpoints, linearly
// interpolated between them (piecewise-quadratic phase) and extrapolated
// LINEARLY beyond the grid ends — the golden references rule out clamping
// (see docs/point_simulator_model.md).
//
// All phase work is in H (double): range, delay and waveform phase are the
// quantities the mixed-precision policy keeps in high precision.
// ==============================================================================
#pragma once

#include <algorithm>
#include <vector>

namespace rsim {

template <typename H>
class Waveform {
public:
    Waveform() = default;

    Waveform(const std::vector<H> &freq, const std::vector<H> &time)
        : f_(freq), t_(time) {
        const size_t n = t_.size();
        phi_.assign(n, H(0));
        k_.assign(n > 1 ? n - 1 : 0, H(0));
        for (size_t i = 0; i + 1 < n; ++i) {
            const H dt = t_[i + 1] - t_[i];
            k_[i] = (f_[i + 1] - f_[i]) / dt;
            phi_[i + 1] = phi_[i] + f_[i] * dt + H(0.5) * k_[i] * dt * dt;
        }
    }

    // Cumulative phase in cycles: phi(t[0]) = 0.
    H Phase(H u) const {
        if (t_.empty()) {
            return H(0);
        }
        if (t_.size() == 1) {
            return f_[0] * (u - t_[0]);
        }
        size_t i = SegmentIndex(u);
        const H du = u - t_[i];
        return phi_[i] + f_[i] * du + H(0.5) * k_[i] * du * du;
    }

    // Instantaneous frequency (linear between breakpoints, linearly extended).
    H Freq(H u) const {
        if (t_.empty()) {
            return H(0);
        }
        if (t_.size() == 1) {
            return f_[0];
        }
        size_t i = SegmentIndex(u);
        return f_[i] + k_[i] * (u - t_[i]);
    }

    // Center frequency used for the wavelength: (min(f) + max(f)) / 2.
    H CenterFreq() const {
        const auto mm = std::minmax_element(f_.begin(), f_.end());
        return H(0.5) * (*mm.first + *mm.second);
    }

    H t_end() const { return t_.empty() ? H(0) : t_.back(); }
    H t_begin() const { return t_.empty() ? H(0) : t_.front(); }

private:
    // Segment index: the last segment with t[i] <= u, clamped to the valid
    // range so out-of-grid u extrapolates from the end segments.
    size_t SegmentIndex(H u) const {
        if (u <= t_[0]) {
            return 0;
        }
        if (u >= t_.back()) {
            return t_.size() - 2;
        }
        // upper_bound: first element > u
        auto it = std::upper_bound(t_.begin(), t_.end(), u);
        size_t idx = static_cast<size_t>(it - t_.begin()) - 1;
        return std::min(idx, t_.size() - 2);
    }

    std::vector<H> f_;     // frequency at breakpoints (Hz)
    std::vector<H> t_;     // breakpoint times (s), zero-referenced
    std::vector<H> phi_;   // cumulative phase at breakpoints (cycles)
    std::vector<H> k_;     // per-segment frequency slope (Hz/s)
};

}  // namespace rsim
