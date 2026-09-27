// ==============================================================================
// radarsim-engine — radar.hpp
// Radar[H, L] (radarsimc.pxd:233-252): Tx + Rx + platform kinematics, and the
// baseband buffer contract.
//
// InitBaseband(H*, H*) hands caller-owned (numpy) buffers to the engine; the
// simulators accumulate into them in place. The engine never allocates the
// output and never frees it.
// ==============================================================================
#pragma once

#include <cmath>
#include <memory>
#include <vector>

#include "receiver.hpp"
#include "rsvector.hpp"
#include "transmitter.hpp"

template <typename H, typename L>
class Radar {
public:
    Radar() = default;

    Radar(const std::shared_ptr<Transmitter<H, L>> &tx,
          const std::shared_ptr<Receiver<L>> &rx,
          std::vector<H> &frame_start_time,
          std::vector<rsv::Vec3<L>> &location_array,
          rsv::Vec3<L> speed_array,
          std::vector<rsv::Vec3<L>> &rotation_array,
          rsv::Vec3<L> rotrate_array)
        : tx_(tx),
          rx_(rx),
          frame_start_time_(frame_start_time),
          location_array_(location_array),
          speed_(speed_array),
          rotation_array_(rotation_array),
          rotrate_(rotrate_array) {
        // Samples per pulse from the waveform extent and the sample rate.
        // Samples per pulse: pulse_length * fs, truncated, with a snap to the
        // nearest integer when within rounding noise of one. The product is
        // formed from the float-narrowed fs, mirroring the Python snap rule
        // (radar.py:432-448) — sim_radar rejects the radar if the two counts
        // disagree.
        if (tx_ && rx_ && !tx_->freq_time_.empty()) {
            const double pulse_length = static_cast<double>(
                tx_->freq_time_.back() - tx_->freq_time_.front());
            const double raw = pulse_length * static_cast<double>(rx_->fs_);
            const double nearest = std::round(raw);
            if (std::fabs(raw - nearest) <=
                1e-9 * std::max(1.0, std::fabs(raw))) {
                sample_size_ = static_cast<int>(nearest);
            } else {
                sample_size_ = static_cast<int>(raw);
            }
        }
    }

    // Caller-owned output buffers, each [frames*tx*rx, pulses, samples]
    // flattened. Must be called before any simulator Run().
    void InitBaseband(H *bb_real, H *bb_imag) {
        bb_real_ = bb_real;
        bb_imag_ = bb_imag;
    }

    // Device-to-host sync for the GPU policy; no-op on the CPU policy.
    void SyncBaseband() {}

    int sample_size_ = 0;

    std::shared_ptr<Transmitter<H, L>> tx_;
    std::shared_ptr<Receiver<L>> rx_;
    std::vector<H> frame_start_time_;
    std::vector<rsv::Vec3<L>> location_array_;
    rsv::Vec3<L> speed_;
    std::vector<rsv::Vec3<L>> rotation_array_;
    rsv::Vec3<L> rotrate_;

    H *bb_real_ = nullptr;
    H *bb_imag_ = nullptr;
};
