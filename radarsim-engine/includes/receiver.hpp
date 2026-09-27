// ==============================================================================
// radarsim-engine — receiver.hpp
// Receiver[T] (radarsimc.pxd:206-227).
//
// gate_delay is deliberately double regardless of T: with T=float, a float32
// mantissa resolves only ~5.8e-11 s at a 741 us gate, i.e. ~0.5 cycles of
// carrier phase at 9 GHz. Keep it double everywhere it propagates.
// ==============================================================================
#pragma once

#include <complex>
#include <vector>

#include "rsvector.hpp"

template <typename T>
class Receiver {
public:
    Receiver() = default;

    Receiver(const T &fs, const T &rf_gain, const T &resistor,
             const T &baseband_gain, const T &baseband_bw,
             const double &gate_delay)
        : fs_(fs),
          rf_gain_(rf_gain),
          resistor_(resistor),
          baseband_gain_(baseband_gain),
          baseband_bw_(baseband_bw),
          gate_delay_(gate_delay) {}

    struct Channel {
        rsv::Vec3<T> location;
        rsv::Vec3<std::complex<T>> polar;
        std::vector<T> phi, phi_ptn;
        std::vector<T> theta, theta_ptn;
        T antenna_gain = 0;
    };

    void AddChannel(const rsv::Vec3<T> &location,
                    const rsv::Vec3<std::complex<T>> &polar,
                    const std::vector<T> &phi, const std::vector<T> &phi_ptn,
                    const std::vector<T> &theta,
                    const std::vector<T> &theta_ptn, const T &antenna_gain) {
        channels_.push_back(
            Channel{location, polar, phi, phi_ptn, theta, theta_ptn,
                    antenna_gain});
    }

    T fs_ = 0;
    T rf_gain_ = 0;
    T resistor_ = 0;
    T baseband_gain_ = 0;
    T baseband_bw_ = 0;
    double gate_delay_ = 0.0;
    std::vector<Channel> channels_;
};
