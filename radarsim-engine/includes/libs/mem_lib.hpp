// ==============================================================================
// radarsim-engine — libs/mem_lib.hpp
// Raw-pointer → std::vector copy helpers used by the Cython marshalling layer
// (src/radarsimpy/includes/radarsimc.pxd:140-145).
// ==============================================================================
#pragma once

#include <complex>
#include <vector>

#include "rsvector.hpp"
#include "core/types.hpp"

template <typename T>
void Mem_Copy(T *ptr, int_t size, std::vector<T> &vect) {
    vect.assign(ptr, ptr + size);
}

template <typename T>
void Mem_Copy_Complex(T *ptr_real, T *ptr_imag, int_t size,
                      std::vector<std::complex<T>> &vect) {
    vect.resize(size);
    for (int_t i = 0; i < size; ++i) {
        vect[i] = std::complex<T>(ptr_real[i], ptr_imag[i]);
    }
}

template <typename T>
void Mem_Copy_Vec3(T *ptr_x, T *ptr_y, T *ptr_z, int_t size,
                   std::vector<rsv::Vec3<T>> &vect) {
    vect.resize(size);
    for (int_t i = 0; i < size; ++i) {
        vect[i] = rsv::Vec3<T>(ptr_x[i], ptr_y[i], ptr_z[i]);
    }
}
