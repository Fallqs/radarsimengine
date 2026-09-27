// ==============================================================================
// radarsim-engine — rf/antenna.hpp
// Antenna view angles and pattern gain.
//
// Pattern semantics (pinned by test_simc_tx_az_pattern / test_simc_tx_el_
// pattern): the gain is the NEAREST table entry, never an interpolation:
//   gain_db = az_ptn[argmin|az_table - az|] + el_ptn[argmin|el_table - el|]
//             + antenna_gain
// Tables arrive from the marshalling layer as: azimuth in radians; elevation
// already transformed to the polar angle theta = 90deg - el, flipped and in
// radians (cp_radarsimc_radar.pyx:168). So the lookup key for elevation is
// theta = acos(dz), not the elevation angle itself.
// ==============================================================================
#pragma once

#include <cmath>

#include "rsvector.hpp"

namespace rsim {

// Nearest-entry table lookup. Tables are non-decreasing angle arrays (L),
// query in the same units.
template <typename L>
L NearestPattern(const L *angles, const L *values, int count, double query) {
    int best = 0;
    double best_d = std::fabs(static_cast<double>(angles[0]) - query);
    for (int i = 1; i < count; ++i) {
        const double d = std::fabs(static_cast<double>(angles[i]) - query);
        if (d < best_d) {
            best_d = d;
            best = i;
        }
    }
    return values[best];
}

// View angles of a unit (or unnormalized) direction in the channel frame:
// azimuth (rad, 0 at +x, + toward +y) and polar angle theta (rad from +z).
template <typename L>
inline void ViewAngles(const rsv::Vec3<L> &dir, double &az, double &theta) {
    const double x = dir[0], y = dir[1], z = dir[2];
    az = std::atan2(y, x);
    const double r = std::sqrt(x * x + y * y + z * z);
    theta = r > 0.0 ? std::acos(z / r) : 0.0;
}

// Total pattern gain (dB) for a channel at (az, theta). The channel struct is
// the Transmitter/Receiver AddChannel record; both have identical layout of
// the pattern fields, so this is templated on the channel type.
template <typename Channel>
double PatternGainDb(const Channel &ch, double az, double theta) {
    const double g_az = static_cast<double>(NearestPattern(
        ch.phi.data(), ch.phi_ptn.data(), static_cast<int>(ch.phi.size()), az));
    const double g_el = static_cast<double>(
        NearestPattern(ch.theta.data(), ch.theta_ptn.data(),
                       static_cast<int>(ch.theta.size()), theta));
    return g_az + g_el + static_cast<double>(ch.antenna_gain);
}

}  // namespace rsim
