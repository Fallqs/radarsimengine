// ==============================================================================
// radarsim-engine — transmitter.hpp
// Transmitter[H, L] (radarsimc.pxd:165-200): sampled arbitrary waveform,
// optional SSB phase-noise spec, and per-channel antenna/modulation data.
// Phase 0 stores everything; waveform evaluation lands with Phase 1 (rf/).
// ==============================================================================
#pragma once

#include <complex>
#include <vector>

#include "rsvector.hpp"

template <typename H, typename L>
class Transmitter {
public:
    Transmitter() = default;

    Transmitter(const L &tx_power, const std::vector<H> &freq,
                const std::vector<H> &freq_time,
                const std::vector<H> &freq_offset,
                const std::vector<H> &pulse_start_time)
        : tx_power_(tx_power),
          freq_(freq),
          freq_time_(freq_time),
          freq_offset_(freq_offset),
          pulse_start_time_(pulse_start_time) {}

    // With SSB phase-noise specification (per-frame deferred generation).
    Transmitter(const L &tx_power, const std::vector<H> &freq,
                const std::vector<H> &freq_time,
                const std::vector<H> &freq_offset,
                const std::vector<H> &pulse_start_time,
                const std::vector<H> &pn_freq, const std::vector<H> &pn_power,
                const H &pn_fs, const int &pn_num_samples,
                const unsigned long long &pn_seed, const bool &pn_validation)
        : Transmitter(tx_power, freq, freq_time, freq_offset,
                      pulse_start_time) {
        has_phase_noise_ = true;
        pn_freq_ = pn_freq;
        pn_power_ = pn_power;
        pn_fs_ = pn_fs;
        pn_num_samples_ = pn_num_samples;
        pn_seed_ = pn_seed;
        pn_validation_ = pn_validation;
    }

    struct Channel {
        rsv::Vec3<L> location;
        rsv::Vec3<std::complex<L>> polar;
        std::vector<L> phi, phi_ptn;      // azimuth table (rad, dB)
        std::vector<L> theta, theta_ptn;  // elevation table (rad, dB)
        L antenna_gain = 0;
        std::vector<L> mod_t;
        std::vector<std::complex<L>> mod_var;
        std::vector<std::complex<L>> pulse_mod;
        L delay = 0;
        L grid = 0;
    };

    void AddChannel(const rsv::Vec3<L> &location,
                    const rsv::Vec3<std::complex<L>> &polar,
                    const std::vector<L> &phi, const std::vector<L> &phi_ptn,
                    const std::vector<L> &theta,
                    const std::vector<L> &theta_ptn, const L &antenna_gain,
                    const std::vector<L> &mod_t,
                    const std::vector<std::complex<L>> &mod_var,
                    const std::vector<std::complex<L>> &pulse_mod,
                    const L &delay, const L &grid) {
        channels_.push_back(Channel{location, polar, phi, phi_ptn, theta,
                                    theta_ptn, antenna_gain, mod_t, mod_var,
                                    pulse_mod, delay, grid});
    }

    L tx_power_ = 0;
    std::vector<H> freq_, freq_time_, freq_offset_, pulse_start_time_;
    bool has_phase_noise_ = false;
    std::vector<H> pn_freq_, pn_power_;
    H pn_fs_ = 0;
    int pn_num_samples_ = 0;
    unsigned long long pn_seed_ = 0;
    bool pn_validation_ = false;
    std::vector<Channel> channels_;
};
