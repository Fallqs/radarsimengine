// ==============================================================================
// radarsim-engine — rsvector.hpp
// 2D/3D value vectors, namespace rsv. Matches
// src/radarsimpy/includes/rsvector.pxd: default / component / raw-pointer
// constructors, assignment, unchecked component access.
// T must also work as std::complex<L> (polarization vectors), so arithmetic
// operators are only instantiated on use.
// ==============================================================================
#pragma once

namespace rsv {

template <typename T>
struct Vec3 {
    T data_[3];

    Vec3() : data_{T{}, T{}, T{}} {}
    Vec3(const T &x, const T &y, const T &z) : data_{x, y, z} {}
    explicit Vec3(T *ptr) : data_{ptr[0], ptr[1], ptr[2]} {}

    Vec3 &operator=(const Vec3 &other) = default;

    T &operator[](const unsigned int &i) { return data_[i]; }
    const T &operator[](const unsigned int &i) const { return data_[i]; }

    Vec3 operator+(const Vec3 &rhs) const {
        return Vec3(data_[0] + rhs.data_[0], data_[1] + rhs.data_[1],
                    data_[2] + rhs.data_[2]);
    }
    Vec3 operator-(const Vec3 &rhs) const {
        return Vec3(data_[0] - rhs.data_[0], data_[1] - rhs.data_[1],
                    data_[2] - rhs.data_[2]);
    }
    Vec3 operator*(const T &s) const {
        return Vec3(data_[0] * s, data_[1] * s, data_[2] * s);
    }
    T Dot(const Vec3 &rhs) const {
        return data_[0] * rhs.data_[0] + data_[1] * rhs.data_[1] +
               data_[2] * rhs.data_[2];
    }
    Vec3 Cross(const Vec3 &rhs) const {
        return Vec3(data_[1] * rhs.data_[2] - data_[2] * rhs.data_[1],
                    data_[2] * rhs.data_[0] - data_[0] * rhs.data_[2],
                    data_[0] * rhs.data_[1] - data_[1] * rhs.data_[0]);
    }
};

template <typename T>
struct Vec2 {
    T data_[2];

    Vec2() : data_{T{}, T{}} {}
    Vec2(const T &x, const T &y) : data_{x, y} {}
    explicit Vec2(T *ptr) : data_{ptr[0], ptr[1]} {}

    Vec2 &operator=(const Vec2 &other) = default;

    T &operator[](const unsigned int &i) { return data_[i]; }
    const T &operator[](const unsigned int &i) const { return data_[i]; }
};

}  // namespace rsv
