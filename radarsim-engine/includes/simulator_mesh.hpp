// ==============================================================================
// radarsim-engine — simulator_mesh.hpp
// SBR + Physical Optics mesh simulator (radarsimc.pxd:382-393), the core of
// the engine. Phase 0: validates inputs, honors dry_run, returns SUCCESS;
// the pipeline lands with Phase 4 (docs/04-physics.md §4.5).
//
// Parameter semantics (gen_docs/user_guide/ray_tracing_simulation.rst):
//   level            0 = re-trace per frame, 1 = per pulse, 2 = per sample
//   density          rays per wavelength; per-target density overrides when
//                    the target's own density != 0
//   ray_filter       [min, max] bounce-count window; max also caps trace depth
//   back_propagating escaped rays also return via reflections on arrival
//   log_path         HDF5 ray dump for debugging (Phase 4/5)
//   dry_run          full setup and validation without ray tracing
// ==============================================================================
#pragma once

#include <memory>
#include <string>

#include "core/enums.hpp"
#include "core/types.hpp"
#include "radar.hpp"
#include "rsvector.hpp"
#include "targets_manager.hpp"

template <typename H, typename L, typename ExecutionPolicy>
class MeshSimulator {
public:
    MeshSimulator() = default;

    RadarSimErrorCode Run(
        const std::shared_ptr<Radar<H, L>> &radar,
        const std::shared_ptr<TargetsManager<L>> &targets_manager, int level,
        L density, rsv::Vec2<int_t> ray_filter, bool back_propagating,
        std::string log_path, bool dry_run) {
        if (!radar || !targets_manager) {
            return RADARSIMCPP_ERROR_NULL_POINTER;
        }
        if (level < 0 || level > 2) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        if (density <= L(0)) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        if (ray_filter[0] < 0 || ray_filter[1] < ray_filter[0]) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        if (!dry_run &&
            (radar->bb_real_ == nullptr || radar->bb_imag_ == nullptr)) {
            return RADARSIMCPP_ERROR_INVALID_PARAMETER;
        }
        (void)back_propagating;
        (void)log_path;
        return SUCCESS;
    }
};
