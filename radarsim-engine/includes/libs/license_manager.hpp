// ==============================================================================
// radarsim-engine — libs/license_manager.hpp
// Permissive stub honoring the LicenseManager contract
// (src/radarsimpy/includes/radarsimc.pxd:151-159) so license.pyx compiles and
// links unmodified. Always licensed, never free tier — tier enforcement lives
// in the Python layer (simulator_radar.pyx:76-118) and is unaffected.
// ==============================================================================
#pragma once

#include <string>
#include <vector>

#include "core/export.hpp"

class RADARSIMCPP_API LicenseManager {
public:
    static LicenseManager &GetInstance();

    void SetLicense(const std::string &license_file_path,
                    const std::string &product);
    void SetLicense(const std::vector<std::string> &license_file_paths,
                    const std::string &product);

    bool IsLicensed() const;
    bool IsFreeTier() const;
    std::string GetLicenseInfo() const;

private:
    LicenseManager() = default;
    std::vector<std::string> license_paths_;
    std::string product_;
};
