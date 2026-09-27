// ==============================================================================
// radarsim-engine — core/export.hpp
// Symbol visibility for the shared library. Only non-template symbols need it
// (LicenseManager, radarsimx::gpu_available); everything else is header-only
// template code instantiated in the consumer's translation unit.
// ==============================================================================
#pragma once

#if defined(_WIN32) || defined(__CYGWIN__)
    #if defined(radarsimcpp_EXPORTS)
        #define RADARSIMCPP_API __declspec(dllexport)
    #else
        #define RADARSIMCPP_API __declspec(dllimport)
    #endif
#else
    #define RADARSIMCPP_API __attribute__((visibility("default")))
#endif
