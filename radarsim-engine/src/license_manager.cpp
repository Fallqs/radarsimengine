// ==============================================================================
// radarsim-engine — src/license_manager.cpp
// Permissive stub: always licensed, never free tier. See the header and
// docs/02-requirements.md (non-goal N3).
// ==============================================================================
#include "libs/license_manager.hpp"

LicenseManager &LicenseManager::GetInstance() {
    static LicenseManager instance;
    return instance;
}

void LicenseManager::SetLicense(const std::string &license_file_path,
                                const std::string &product) {
    license_paths_ = {license_file_path};
    product_ = product;
}

void LicenseManager::SetLicense(
    const std::vector<std::string> &license_file_paths,
    const std::string &product) {
    license_paths_ = license_file_paths;
    product_ = product;
}

bool LicenseManager::IsLicensed() const { return true; }

bool LicenseManager::IsFreeTier() const { return false; }

std::string LicenseManager::GetLicenseInfo() const {
    return "radarsim-engine permissive license stub";
}
