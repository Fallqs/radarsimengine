// ==============================================================================
// radarsim-engine — tests/test_smoke.cpp
// Phase-0 smoke test: instantiates every public type the Cython binding
// declares, with the same template arguments the Python build uses
// (H = double, L = float). Exits non-zero on failure so ctest catches it.
// ==============================================================================
#include <cassert>
#include <cmath>
#include <complex>
#include <cstring>
#include <memory>
#include <vector>

#include "core/execution_policy.hpp"
#include "libs/license_manager.hpp"
#include "libs/mem_lib.hpp"
#include "libs/motion_lib.hpp"
#include "points_manager.hpp"
#include "radar.hpp"
#include "simulator_interference.hpp"
#include "simulator_lidar.hpp"
#include "simulator_mesh.hpp"
#include "simulator_noise.hpp"
#include "simulator_point.hpp"
#include "simulator_rcs.hpp"
#include "targets_manager.hpp"

using H = double;
using L = float;

namespace {

void TestVecAndRotation() {
    rsv::Vec3<L> v(1.0f, 0.0f, 0.0f);
    const float pi_2 = 1.5707963267948966f;
    rsv::Vec3<L> rot(pi_2, 0.0f, 0.0f);  // yaw 90 deg (radians): +x -> +y
    rsv::Vec3<L> r = Rotate(v, rot);
    assert(std::fabs(r[0]) < 1e-6f);
    assert(std::fabs(r[1] - 1.0f) < 1e-6f);
    assert(std::fabs(r[2]) < 1e-6f);

    rsv::Vec2<int_t> f(0, 10);
    assert(f[0] == 0 && f[1] == 10);
}

void TestRadarConstruction() {
    std::vector<H> freq = {77e9, 78e9};
    std::vector<H> freq_time = {0.0, 40e-6};
    std::vector<H> freq_offset = {0.0};
    std::vector<H> pulse_start = {0.0};

    const rsv::Vec3<std::complex<L>> pol(std::complex<L>(0, 0),
                                         std::complex<L>(0, 0),
                                         std::complex<L>(1, 0));

    auto tx = std::make_shared<Transmitter<H, L>>(0.0f, freq, freq_time,
                                                  freq_offset, pulse_start);
    tx->AddChannel(rsv::Vec3<L>(0, 0, 0), pol, {0.0f}, {0.0f}, {0.0f}, {0.0f},
                   0.0f, {0.0f}, {std::complex<L>(1, 0)},
                   {std::complex<L>(1, 0)}, 0.0f, 0.0f);

    auto rx = std::make_shared<Receiver<L>>(10e6f, 0.0f, 500.0f, 0.0f, 10e6f,
                                            0.0);
    rx->AddChannel(rsv::Vec3<L>(0, 0, 0), pol, {0.0f}, {0.0f}, {0.0f}, {0.0f},
                   0.0f);

    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot_arr = {rsv::Vec3<L>(0, 0, 0)};
    auto radar = std::make_shared<Radar<H, L>>(tx, rx, frame_start, loc,
                                               rsv::Vec3<L>(0, 0, 0), rot_arr,
                                               rsv::Vec3<L>(0, 0, 0));

    // 40 us at 10 MHz -> 400 samples
    assert(radar->sample_size_ == 400);

    std::vector<H> bb_real(400, 0.0), bb_imag(400, 0.0);
    radar->InitBaseband(bb_real.data(), bb_imag.data());

    auto points = std::make_shared<PointsManager<L>>();
    points->AddPointSimple(rsv::Vec3<L>(10, 0, 0), rsv::Vec3<L>(0, 0, 0),
                           10.0f, 0.0f);

    PointSimulator<H, L, radarsimx::cpu_policy> point_sim;
    assert(point_sim.Run(radar, points) == SUCCESS);

    InterferenceSimulator<H, L, radarsimx::cpu_policy> interf_sim;
    assert(interf_sim.Run(radar, radar) == SUCCESS);

    std::vector<H> noise_real(400, 1.0), noise_imag(400, 1.0);
    std::vector<H> timestamps(400, 0.0);
    NoiseSimulator<H, L, radarsimx::cpu_policy> noise_sim;
    assert(noise_sim.Run(radar, 1.0, true, timestamps.data(), 1, 1, 400,
                         noise_real.data(), noise_imag.data(), 0) == SUCCESS);
    // real noise now: keyed by (rx, timestamp); all-zero timestamps share one
    // value, which must be nonzero
    assert(noise_real[0] != 0.0);
}

void TestMeshAndRcs() {
    // Single triangle
    std::vector<L> points = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    std::vector<int_t> cells = {0, 1, 2};

    auto targets = std::make_shared<TargetsManager<L>>();
    targets->AddTargetSimple(points.data(), cells.data(), 1,
                             rsv::Vec3<L>(0, 0, 0), rsv::Vec3<L>(10, 0, 0),
                             rsv::Vec3<L>(0, 0, 0), rsv::Vec3<L>(0, 0, 0),
                             rsv::Vec3<L>(0, 0, 0), false, 0.0f, false);
    assert(targets->targets().size() == 1);
    assert(targets->targets()[0]->array_size_ == 1);  // static kinematics
    targets->targets()[0]->Move(0, 0.0);

    std::vector<H> freq = {77e9};
    std::vector<H> freq_time = {0.0, 1e-6};
    std::vector<H> foff = {0.0};
    std::vector<H> pst = {0.0};
    auto tx = std::make_shared<Transmitter<H, L>>(0.0f, freq, freq_time, foff,
                                                  pst);
    tx->AddChannel(rsv::Vec3<L>(0, 0, 0),
                   rsv::Vec3<std::complex<L>>({0, 0}, {0, 0}, {1, 0}),
                   {0.0f}, {0.0f}, {0.0f}, {0.0f}, 0.0f, {0.0f},
                   {std::complex<L>(1, 0)}, {std::complex<L>(1, 0)}, 0.0f,
                   0.0f);
    auto rx = std::make_shared<Receiver<L>>(1e6f, 0.0f, 500.0f, 0.0f, 1e6f,
                                            0.0);
    rx->AddChannel(rsv::Vec3<L>(0, 0, 0),
                   rsv::Vec3<std::complex<L>>({0, 0}, {0, 0}, {1, 0}),
                   {0.0f}, {0.0f}, {0.0f}, {0.0f}, 0.0f);
    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot_arr = {rsv::Vec3<L>(0, 0, 0)};
    auto radar = std::make_shared<Radar<H, L>>(tx, rx, frame_start, loc,
                                               rsv::Vec3<L>(0, 0, 0), rot_arr,
                                               rsv::Vec3<L>(0, 0, 0));
    std::vector<H> bb_r(1, 0.0), bb_i(1, 0.0);
    radar->InitBaseband(bb_r.data(), bb_i.data());

    MeshSimulator<H, L, radarsimx::cpu_policy> mesh_sim;
    assert(mesh_sim.Run(radar, targets, 0, 1.0f, rsv::Vec2<int_t>(0, 10),
                        false, "", true) == SUCCESS);
    assert(mesh_sim.Run(radar, targets, 3, 1.0f, rsv::Vec2<int_t>(0, 10),
                        false, "", false) ==
           RADARSIMCPP_ERROR_INVALID_PARAMETER);

    RcsSimulator<L, radarsimx::cpu_policy, L> rcs_sim;
    assert(rcs_sim.Run(targets, {rsv::Vec3<L>(-1, 0, 0)},
                       {rsv::Vec3<L>(-1, 0, 0)},
                       rsv::Vec3<std::complex<L>>(std::complex<L>(0, 0),
                                                  std::complex<L>(0, 0),
                                                  std::complex<L>(1, 0)),
                       rsv::Vec3<std::complex<L>>(std::complex<L>(0, 0),
                                                  std::complex<L>(0, 0),
                                                  std::complex<L>(1, 0)),
                       77e9f, 1.0f) == SUCCESS);
    assert(rcs_sim.GetRcs().size() == 1);

    LidarSimulator<L, radarsimx::cpu_policy> lidar_sim;
    assert(lidar_sim.Run(targets, {0.0f}, {0.0f}, rsv::Vec3<L>(0, 0, 0)) ==
           SUCCESS);
}

void TestLicenseAndPolicy() {
    auto &lm = LicenseManager::GetInstance();
    lm.SetLicense("dummy.lic", "RadarSimPy");
    assert(lm.IsLicensed());
    assert(!lm.IsFreeTier());

    static_assert(radarsimx::cpu_policy::is_cpu);
    static_assert(!radarsimx::cpu_policy::is_gpu);
#ifndef _CUDA_
    static_assert(!radarsimx::gpu_policy::is_gpu);
    assert(!radarsimx::gpu_available());
#endif
    assert(std::strcmp(radarsimx::cpu.name(), "cpu") == 0);
    assert(std::strcmp(radarsimx::gpu.name(), "gpu") == 0);
}

}  // namespace

int main() {
    TestVecAndRotation();
    TestRadarConstruction();
    TestMeshAndRcs();
    TestLicenseAndPolicy();
    return 0;
}
