// ==============================================================================
// radarsim-engine — rf/platform.hpp
// Radar platform pose (location + rotation) at a given absolute time.
//
// Geometry is float (L) end to end — positions, rotation angles — matching
// the marshalling boundary, which narrows them to float_t. The rotation
// matrix is built from the float angles with double trig. This staging is
// observable in the goldens (sin(f32(pi)) != 0 shifts a rotated channel by
// ~6e-7 m); do not "fix" it.
// ==============================================================================
#pragma once

#include <vector>

#include "libs/motion_lib.hpp"
#include "rsvector.hpp"

namespace rsim {

// Rotation matrix from [yaw, pitch, roll] (radians): Rz(yaw) Ry(-pitch) Rx(roll).
// Double trig of the (float) angles.
template <typename L>
void RotMatrix(const rsv::Vec3<L> &rot, double R[3][3]) {
    const double yaw = rot[0], pitch = rot[1], roll = rot[2];
    const double cy = std::cos(yaw), sy = std::sin(yaw);
    const double cp = std::cos(pitch), sp = std::sin(pitch);
    const double cr = std::cos(roll), sr = std::sin(roll);
    R[0][0] = cy * cp;
    R[0][1] = cy * sp * sr - sy * cr;
    R[0][2] = cy * sp * cr + sy * sr;
    R[1][0] = sy * cp;
    R[1][1] = sy * sp * sr + cy * cr;
    R[1][2] = sy * sp * cr - cy * sr;
    R[2][0] = -sp;
    R[2][1] = cp * sr;
    R[2][2] = cp * cr;
}

// v_world = R . v
template <typename L>
rsv::Vec3<double> RotApply(const double R[3][3], const rsv::Vec3<L> &v) {
    return rsv::Vec3<double>(
        R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
        R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
        R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]);
}

// Platform location (double) and rotation matrix at absolute time T.
// location/rotation arrays are either size 1 (constant base + speed/rate) or
// one entry per flattened timestamp (time-varying, indexed by flat_idx).
template <typename L>
inline void PlatformPose(const std::vector<rsv::Vec3<L>> &location_array,
                         const rsv::Vec3<L> &speed,
                         const std::vector<rsv::Vec3<L>> &rotation_array,
                         const rsv::Vec3<L> &rotrate, double T, size_t flat_idx,
                         rsv::Vec3<double> &loc_out, double R_out[3][3]) {
    rsv::Vec3<L> base_loc, base_rot;
    if (location_array.size() == 1) {
        base_loc = location_array[0];
        base_rot = rotation_array[0];
        loc_out = rsv::Vec3<double>(base_loc[0] + speed[0] * T,
                                    base_loc[1] + speed[1] * T,
                                    base_loc[2] + speed[2] * T);
        rsv::Vec3<L> rot(base_rot[0] + rotrate[0] * static_cast<L>(T),
                         base_rot[1] + rotrate[1] * static_cast<L>(T),
                         base_rot[2] + rotrate[2] * static_cast<L>(T));
        RotMatrix(rot, R_out);
    } else {
        const rsv::Vec3<L> &l = location_array[flat_idx];
        const rsv::Vec3<L> &r = rotation_array[flat_idx];
        loc_out = rsv::Vec3<double>(l[0], l[1], l[2]);
        RotMatrix(r, R_out);
    }
}

}  // namespace rsim
