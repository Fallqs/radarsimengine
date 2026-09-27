// ==============================================================================
// radarsim-engine — libs/motion_lib.hpp
// Orientation math. The rotation convention is load-bearing
// (README "Coordinate Systems"): [yaw, pitch, roll], applied as
// R = Rz(yaw) · Ry(-pitch) · Rx(roll) — the aerospace "nose up is positive"
// convention; pitch is NOT a right-handed rotation about +y.
//
// Angles are in RADIANS: the Python layer converts before marshalling
// (radar.py:844 stores platform rotation via np.radians; cp_radarsimc_mesh.pyx
// converts target rotation the same way).
// ==============================================================================
#pragma once

#include <cmath>

#include "rsvector.hpp"

template <typename T>
rsv::Vec3<T> Rotate(const rsv::Vec3<T> &vect, const rsv::Vec3<T> &rotation) {
    const double yaw = rotation[0];
    const double pitch = rotation[1];
    const double roll = rotation[2];

    const double cy = std::cos(yaw), sy = std::sin(yaw);
    const double cp = std::cos(pitch), sp = std::sin(pitch);
    const double cr = std::cos(roll), sr = std::sin(roll);

    // R = Rz(yaw) · Ry(-pitch) · Rx(roll)
    const double r00 = cy * cp;
    const double r01 = cy * sp * sr - sy * cr;
    const double r02 = cy * sp * cr + sy * sr;
    const double r10 = sy * cp;
    const double r11 = sy * sp * sr + cy * cr;
    const double r12 = sy * sp * cr - cy * sr;
    const double r20 = -sp;
    const double r21 = cp * sr;
    const double r22 = cp * cr;

    const double x = vect[0], y = vect[1], z = vect[2];
    return rsv::Vec3<T>(
        static_cast<T>(r00 * x + r01 * y + r02 * z),
        static_cast<T>(r10 * x + r11 * y + r12 * z),
        static_cast<T>(r20 * x + r21 * y + r22 * z));
}
