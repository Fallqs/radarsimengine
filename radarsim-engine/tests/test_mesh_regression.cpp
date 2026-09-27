// ==============================================================================
// radarsim-engine — tests/test_mesh_regression.cpp
// Self-captured regression references for the mesh and RCS paths.
//
// The upstream mesh/RCS goldens encode the closed engine's exact sampler,
// which this project does not reproduce (docs/mesh_simulator_model.md). These
// references are captured from THIS engine, so the suite guards the mesh/RCS
// pipeline against regressions without conflating that with the sampler
// conformance gap. Tolerance is 1e-3 relative -- a guard, not the contract.
// ==============================================================================
#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>
#include <vector>

#include "core/execution_policy.hpp"
#include "radar.hpp"
#include "simulator_mesh.hpp"
#include "simulator_rcs.hpp"
#include "targets_manager.hpp"

using H = double;
using L = float;

namespace {

int g_failures = 0;

void Check(bool ok, const char *msg) {
    if (!ok) {
        std::printf("FAIL: %s\n", msg);
        ++g_failures;
    }
}

std::shared_ptr<Radar<H, L>> MakeRadar24G() {
    std::vector<H> f = {24.075e9, 24.175e9};
    std::vector<H> t = {0.0, 80e-6};
    std::vector<H> foff = {0.0, 0.0, 0.0};
    std::vector<H> pstart = {0.0, 100e-6, 200e-6};
    auto tx = std::make_shared<Transmitter<H, L>>(L(10), f, t, foff, pstart);
    tx->AddChannel(rsv::Vec3<L>(0, 0, 0),
                   rsv::Vec3<std::complex<L>>({0, 0}, {0, 0}, {1, 0}),
                   {-1.5707963f, 1.5707963f}, {0.0f, 0.0f},
                   {0.0f, 3.1415927f}, {0.0f, 0.0f}, 0.0f, {0.0f},
                   {std::complex<L>(1, 0)},
                   {std::complex<L>(1, 0), std::complex<L>(1, 0),
                    std::complex<L>(1, 0)},
                   0.0f, L(1.0f * 3.1415927f / 180.0f));
    auto rx = std::make_shared<Receiver<L>>(L(6e4), L(20), L(500), L(30),
                                            L(6e4), 0.0);
    rx->AddChannel(rsv::Vec3<L>(0, 0, 0),
                   rsv::Vec3<std::complex<L>>({0, 0}, {0, 0}, {1, 0}),
                   {-1.5707963f, 1.5707963f}, {0.0f, 0.0f},
                   {0.0f, 3.1415927f}, {0.0f, 0.0f}, 0.0f);
    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot = {rsv::Vec3<L>(0, 0, 0)};
    return std::make_shared<Radar<H, L>>(tx, rx, frame_start, loc,
                                         rsv::Vec3<L>(0, 0, 0), rot,
                                         rsv::Vec3<L>(0, 0, 0));
}

std::shared_ptr<TargetsManager<L>> MakePlate() {
    std::vector<L> points = {0, -2.5f, -2.5f, 0, 2.5f, -2.5f,
                             0, 2.5f,  2.5f,  0, -2.5f, 2.5f};
    std::vector<int_t> cells = {0, 1, 2, 0, 2, 3};
    auto tgts = std::make_shared<TargetsManager<L>>();
    tgts->AddTargetSimple(points.data(), cells.data(), 2, rsv::Vec3<L>(0, 0, 0),
                          rsv::Vec3<L>(10, 0, 0), rsv::Vec3<L>(0, 0, 0),
                          rsv::Vec3<L>(0, 0, 0), rsv::Vec3<L>(0, 0, 0), false,
                          0.0f, false);
    return tgts;
}

void TestPlateMesh() {
    auto radar = MakeRadar24G();
    auto tgts = MakePlate();
    std::vector<H> re(12, 0), im(12, 0);
    radar->InitBaseband(re.data(), im.data());
    MeshSimulator<H, L, radarsimx::cpu_policy> sim;
    Check(sim.Run(radar, tgts, 0, 0.4f, rsv::Vec2<int_t>(0, 10), false, "",
                  false) == SUCCESS,
          "plate mesh Run");
    // self-captured reference (this engine), 1e-3 relative guard
    const std::complex<double> ref[4] = {
        {-0.033270, -0.028343},
        {0.043332, -0.001778},
        {-0.031157, 0.031146},
        {0.002454, -0.045549},
    };
    double peak = 0.0;
    for (const auto &v : ref) {
        peak = std::max(peak, std::abs(v));
    }
    for (int s = 0; s < 4; ++s) {
        const double err = std::abs(
            std::complex<double>(re[s], im[s]) - ref[s]);
        if (err > 1e-3 * peak) {
            std::printf("FAIL: plate mesh s%d: %.6f%+.6fj vs ref, err %.3e\n",
                        s, re[s], im[s], err);
            ++g_failures;
        }
    }
}

void TestPlateRcs() {
    auto tgts = MakePlate();  // at origin for the RCS case
    std::vector<L> points = {0, -2.5f, -2.5f, 0, 2.5f, -2.5f,
                             0, 2.5f,  2.5f,  0, -2.5f, 2.5f};
    std::vector<int_t> cells = {0, 1, 2, 0, 2, 3};
    auto tgts0 = std::make_shared<TargetsManager<L>>();
    tgts0->AddTargetSimple(points.data(), cells.data(), 2,
                           rsv::Vec3<L>(0, 0, 0), rsv::Vec3<L>(0, 0, 0),
                           rsv::Vec3<L>(0, 0, 0), rsv::Vec3<L>(0, 0, 0),
                           rsv::Vec3<L>(0, 0, 0), false, 0.0f, false);
    RcsSimulator<double, radarsimx::cpu_policy, L> sim;
    Check(sim.Run(tgts0, {rsv::Vec3<double>(1, 0, 0)},
                  {rsv::Vec3<double>(1, 0, 0)},
                  rsv::Vec3<std::complex<double>>({0, 0}, {0, 0}, {1, 0}),
                  rsv::Vec3<std::complex<double>>({0, 0}, {0, 0}, {1, 0}),
                  1e9, 1.0) == SUCCESS,
          "plate RCS Run");
    // analytic PO plate: 4 pi A^2 / lambda^2, A=25, lambda=0.299792458
    const double lam = 299792458.0 / 1e9;
    const double ideal = 4.0 * 3.14159265358979 * 625.0 / (lam * lam);
    const double got = sim.GetRcs()[0];
    if (std::fabs(got - ideal) / ideal > 1e-3) {
        std::printf("FAIL: plate RCS %.1f vs ideal %.1f\n", got, ideal);
        ++g_failures;
    }
}

void TestTimeVaryingTargetIndex() {
    // Guards the timestamp-grid index used to pose time-varying targets at
    // each pass instant (MeshSimulator::BuildScene). With the index
    // miscomputed (e.g. derived from seconds), every pass poses the target
    // from entry 0 and the per-sample outputs coincide.
    std::vector<H> f = {24.125e9, 24.125e9};
    std::vector<H> t = {0.0, 2e-6};
    std::vector<H> foff = {0.0};
    std::vector<H> pstart = {0.0};
    auto tx = std::make_shared<Transmitter<H, L>>(L(10), f, t, foff, pstart);
    tx->AddChannel(rsv::Vec3<L>(0, 0, 0),
                   rsv::Vec3<std::complex<L>>({0, 0}, {0, 0}, {1, 0}),
                   {-1.5707963f, 1.5707963f}, {0.0f, 0.0f},
                   {0.0f, 3.1415927f}, {0.0f, 0.0f}, 0.0f, {0.0f},
                   {std::complex<L>(1, 0)},
                   {std::complex<L>(1, 0)},
                   0.0f, L(1.0f * 3.1415927f / 180.0f));
    auto rx = std::make_shared<Receiver<L>>(L(1e6), L(20), L(500), L(30),
                                            L(1e6), 0.0);
    rx->AddChannel(rsv::Vec3<L>(0, 0, 0),
                   rsv::Vec3<std::complex<L>>({0, 0}, {0, 0}, {1, 0}),
                   {-1.5707963f, 1.5707963f}, {0.0f, 0.0f},
                   {0.0f, 3.1415927f}, {0.0f, 0.0f}, 0.0f);
    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot = {rsv::Vec3<L>(0, 0, 0)};
    auto radar = std::make_shared<Radar<H, L>>(tx, rx, frame_start, loc,
                                               rsv::Vec3<L>(0, 0, 0), rot,
                                               rsv::Vec3<L>(0, 0, 0));

    // 5x5 plate in the y-z plane; kinematics expanded on the flat timestamp
    // grid (1 frame x 1 channel x 1 pulse x 2 samples): entry 0 at x=10,
    // entry 1 at x=20
    std::vector<L> points = {0, -2.5f, -2.5f, 0, 2.5f, -2.5f,
                             0, 2.5f,  2.5f,  0, -2.5f, 2.5f};
    std::vector<int_t> cells = {0, 1, 2, 0, 2, 3};
    std::vector<rsv::Vec3<L>> loc2 = {rsv::Vec3<L>(10, 0, 0),
                                      rsv::Vec3<L>(20, 0, 0)};
    std::vector<rsv::Vec3<L>> zero2 = {rsv::Vec3<L>(0, 0, 0),
                                       rsv::Vec3<L>(0, 0, 0)};
    auto tgts = std::make_shared<TargetsManager<L>>();
    tgts->AddTarget(points.data(), cells.data(), 2, rsv::Vec3<L>(0, 0, 0),
                    loc2, zero2, zero2, zero2,
                    std::complex<L>(1e38f, 0.0f), std::complex<L>(1.0f, 0.0f),
                    false, 0.0f, false);

    std::vector<H> re(2, 0), im(2, 0);
    radar->InitBaseband(re.data(), im.data());
    MeshSimulator<H, L, radarsimx::cpu_policy> sim;
    Check(sim.Run(radar, tgts, 2, 1.0f, rsv::Vec2<int_t>(0, 10), false, "",
                  false) == SUCCESS,
          "time-varying mesh Run");

    const std::complex<double> s0(re[0], im[0]), s1(re[1], im[1]);
    // discriminating invariant: the two passes must see different kinematics
    // entries (a miscomputed index poses both from entry 0 and the outputs
    // coincide exactly for this CW waveform, whose beat is u-independent)
    if (std::abs(s1 - s0) < 0.1 * std::abs(s0)) {
        std::printf("FAIL: time-varying idx: samples coincide "
                    "(%.6e%+.6ej vs %.6e%+.6ej) -- both passes posed "
                    "entry 0\n", s0.real(), s0.imag(), s1.real(), s1.imag());
        ++g_failures;
    }
    // self-captured guard values (this engine), 1e-3 relative
    const std::complex<double> ref[2] = {
        {4.899248e-02, -1.751534e-02},
        {-1.926321e-02, 3.345133e-02},
    };
    const double peak = std::max(std::abs(ref[0]), std::abs(ref[1]));
    for (int s = 0; s < 2; ++s) {
        const std::complex<double> got(re[s], im[s]);
        if (std::abs(got - ref[s]) > 1e-3 * peak) {
            std::printf("FAIL: time-varying mesh s%d: %.6e%+.6ej vs ref, "
                        "err %.3e\n", s, got.real(), got.imag(),
                        std::abs(got - ref[s]));
            ++g_failures;
        }
    }
}

}  // namespace

int main() {
    TestPlateMesh();
    TestPlateRcs();
    TestTimeVaryingTargetIndex();
    if (g_failures == 0) {
        std::printf("test_mesh_regression: all passed\n");
        return 0;
    }
    std::printf("test_mesh_regression: %d failures\n", g_failures);
    return 1;
}
