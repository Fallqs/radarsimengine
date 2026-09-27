// ==============================================================================
// radarsim-engine — geom/fresnel.hpp
// Fresnel reflection for the SBR bounces. Complex relative permittivity /
// permeability; PEC is encoded as |eps| = 1e38 (the marshalling default).
// ==============================================================================
#pragma once

#include <cmath>
#include <complex>

#include "rsvector.hpp"

namespace rsim {

// Specular reflection direction: d - 2 (d.n) n.
template <typename T>
inline rsv::Vec3<T> Reflect(const rsv::Vec3<T> &d, const rsv::Vec3<T> &n) {
    const T dn = d.Dot(n);
    return d - n * (T(2) * dn);
}

// Fresnel amplitude for the perpendicular (TE) component. cos_i >= 0.
template <typename T>
inline std::complex<T> FresnelTE(const std::complex<T> &eps_r,
                                 const std::complex<T> &mu_r, T cos_i) {
    if (std::abs(eps_r) > T(1e30)) {  // PEC
        return std::complex<T>(T(-1), T(0));
    }
    const std::complex<T> eta2 = std::sqrt(mu_r / eps_r);
    const std::complex<T> sin2_t =
        (T(1) - cos_i * cos_i) / (eps_r * mu_r);
    const std::complex<T> cos_t = std::sqrt(std::complex<T>(T(1), T(0)) - sin2_t);
    return (mu_r * cos_i - eta2 * cos_t) / (mu_r * cos_i + eta2 * cos_t);
}

// Parallel (TM) component.
template <typename T>
inline std::complex<T> FresnelTM(const std::complex<T> &eps_r,
                                 const std::complex<T> &mu_r, T cos_i) {
    if (std::abs(eps_r) > T(1e30)) {  // PEC
        return std::complex<T>(T(1), T(0));
    }
    const std::complex<T> eta2 = std::sqrt(mu_r / eps_r);
    const std::complex<T> sin2_t =
        (T(1) - cos_i * cos_i) / (eps_r * mu_r);
    const std::complex<T> cos_t = std::sqrt(std::complex<T>(T(1), T(0)) - sin2_t);
    return (eps_r * cos_i - eta2 * cos_t) / (eps_r * cos_i + eta2 * cos_t);
}

}  // namespace rsim
