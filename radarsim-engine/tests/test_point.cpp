// ==============================================================================
// radarsim-engine — tests/test_point.cpp
// Golden-value tests for the Phase-1 simulators (point, interference, noise),
// replicating scenarios and expected values from the upstream suites
// tests/test_module_sim_radar_ideal.py (goldens recorded from the reference
// engine; tolerance = 1e-6 of peak, IDEAL_BASEBAND_PEAK_ATOL).
// ==============================================================================
#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>
#include <vector>

#include "core/execution_policy.hpp"
#include "points_manager.hpp"
#include "radar.hpp"
#include "simulator_interference.hpp"
#include "simulator_noise.hpp"
#include "simulator_point.hpp"

using namespace std::literals::complex_literals;

using H = double;
using L = float;

namespace {

constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;

void Check(bool ok, const char *msg) {
    if (!ok) {
        std::printf("FAIL: %s\n", msg);
        ++g_failures;
    }
}

// Build a default (isotropic) channel config matching the marshalled form of
// {"location": loc}: az table [-90, 90] deg in radians with [0, 0] dB,
// theta table [0, 180] deg in radians with [0, 0] dB, pol (0, 0, 1).
struct DefaultChannel {
    rsv::Vec3<L> location;
    std::vector<L> az = {-1.5707963267948966f, 1.5707963267948966f};
    std::vector<L> az_ptn = {0.0f, 0.0f};
    std::vector<L> th = {0.0f, 3.141592653589793f};
    std::vector<L> th_ptn = {0.0f, 0.0f};
};

void AddDefaultTxChannel(Transmitter<H, L> &tx, rsv::Vec3<L> loc, L delay,
                         int pulses,
                         const std::vector<std::complex<L>> &pulse_mod,
                         const std::vector<L> &mod_t,
                         const std::vector<std::complex<L>> &mod_var) {
    DefaultChannel dc;
    tx.AddChannel(loc, rsv::Vec3<std::complex<L>>({0, 0}, {0, 0}, {1, 0}),
                  dc.az, dc.az_ptn, dc.th, dc.th_ptn, L(0), mod_t, mod_var,
                  pulse_mod.empty() ? std::vector<std::complex<L>>(pulses, {1, 0})
                                    : pulse_mod,
                  delay, L(0.017453292519943295));
}

void AddDefaultRxChannel(Receiver<L> &rx, rsv::Vec3<L> loc) {
    DefaultChannel dc;
    rx.AddChannel(loc, rsv::Vec3<std::complex<L>>({0, 0}, {0, 0}, {1, 0}),
                  dc.az, dc.az_ptn, dc.th, dc.th_ptn, L(0));
}

std::shared_ptr<Radar<H, L>> MakeRadar24G(L fs, int pulses, L prp,
                                          double gate = 0.0) {
    std::vector<H> f = {24.075e9, 24.175e9};
    std::vector<H> t = {0.0, 80e-6};
    std::vector<H> foff(pulses, 0.0);
    std::vector<H> pstart(pulses);
    for (int p = 0; p < pulses; ++p) {
        pstart[p] = p * static_cast<double>(prp);
    }
    auto tx = std::make_shared<Transmitter<H, L>>(L(10), f, t, foff, pstart);
    AddDefaultTxChannel(*tx, rsv::Vec3<L>(0, 0, 0), L(0), pulses, {}, {}, {});
    auto rx = std::make_shared<Receiver<L>>(fs, L(20), L(500), L(30), fs, gate);
    AddDefaultRxChannel(*rx, rsv::Vec3<L>(0, 0, 0));
    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot = {rsv::Vec3<L>(0, 0, 0)};
    return std::make_shared<Radar<H, L>>(tx, rx, frame_start, loc,
                                         rsv::Vec3<L>(0, 0, 0), rot,
                                         rsv::Vec3<L>(0, 0, 0));
}

// Peak-relative comparison against golden
void CheckBaseband(const std::vector<H> &re, const std::vector<H> &im,
                   const std::vector<std::complex<double>> &gold,
                   const char *name) {
    double peak = 0.0;
    for (const auto &g : gold) {
        peak = std::max(peak, std::abs(g));
    }
    const double atol = 1e-6 * peak;
    for (size_t i = 0; i < gold.size(); ++i) {
        const double err = std::abs(
            std::complex<double>(re[i], im[i]) - gold[i]);
        if (err > atol) {
            std::printf("FAIL: %s sample %zu: err %.3e > %.3e\n", name, i, err,
                        atol);
            ++g_failures;
            return;
        }
    }
}

// --- test_simc_single_target ---
void TestSingleTarget() {
    auto radar = MakeRadar24G(L(6e4), 3, L(100e-6));
    auto points = std::make_shared<PointsManager<L>>();
    points->AddPointSimple(rsv::Vec3<L>(10, 0, 0), rsv::Vec3<L>(0, 0, 0),
                           L(20), L(0));
    const int64_t total = 3 * 4;
    std::vector<H> re(total), im(total);
    radar->InitBaseband(re.data(), im.data());
    PointSimulator<H, L, radarsimx::cpu_policy> sim;
    Check(sim.Run(radar, points) == SUCCESS, "single target Run");
    const std::complex<double> g = 0.02167871966958046 + 0.017555851489305496i;
    std::vector<std::complex<double>> gold = {
        g,
        -0.027893975377082825 + 0.00031773850787431i,
        0.021273192018270493 - 0.018045110628008842i,
        -0.004863050766289234 + 0.02746862918138504i,
        g,
        -0.027893975377082825 + 0.00031773850787431i,
        0.021273192018270493 - 0.018045110628008842i,
        -0.004863050766289234 + 0.02746862918138504i,
        g,
        -0.027893975377082825 + 0.00031773850787431i,
        0.021273192018270493 - 0.018045110628008842i,
        -0.004863050766289234 + 0.02746862918138504i,
    };
    CheckBaseband(re, im, gold, "single_target");
}

// --- test_simc_varing_prp: moving target, nonuniform PRP ---
void TestVaryingPrp() {
    std::vector<H> f = {24.075e9, 24.175e9};
    std::vector<H> t = {0.0, 80e-6};
    std::vector<H> foff = {0.0, 0.0, 0.0};
    std::vector<H> pstart = {0.0, 110e-6, 240e-6};  // cumsum([100,110,130])-100
    auto tx = std::make_shared<Transmitter<H, L>>(L(10), f, t, foff, pstart);
    AddDefaultTxChannel(*tx, rsv::Vec3<L>(0, 0, 0), L(0), 3, {}, {}, {});
    auto rx = std::make_shared<Receiver<L>>(L(6e4), L(20), L(500), L(30),
                                            L(6e4), 0.0);
    AddDefaultRxChannel(*rx, rsv::Vec3<L>(0, 0, 0));
    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot = {rsv::Vec3<L>(0, 0, 0)};
    auto radar = std::make_shared<Radar<H, L>>(tx, rx, frame_start, loc,
                                               rsv::Vec3<L>(0, 0, 0), rot,
                                               rsv::Vec3<L>(0, 0, 0));
    auto points = std::make_shared<PointsManager<L>>();
    points->AddPointSimple(rsv::Vec3<L>(10, 0, 0), rsv::Vec3<L>(-10, 0, 0),
                           L(20), L(0));
    std::vector<H> re(12), im(12);
    radar->InitBaseband(re.data(), im.data());
    PointSimulator<H, L, radarsimx::cpu_policy> sim;
    Check(sim.Run(radar, points) == SUCCESS, "varying prp Run");
    std::vector<std::complex<double>> gold = {
        0.02167871966958046 + 0.017555851489305496i,
        -0.027447368949651718 + 0.004986839834600687i,
        0.014111651107668877 - 0.024065323173999786i,
        0.009057712741196156 + 0.026387274265289307i,
        0.025369323790073395 - 0.011615276336669922i,
        -0.007712406571954489 + 0.02681581676006317i,
        -0.015329970046877861 - 0.023315511643886566i,
        0.027679286897182465 + 0.0035397966857999563i,
        -0.0047342246398329735 - 0.02750471606850624i,
        0.023969972506165504 + 0.01429736614227295i,
        -0.02644996903836727 + 0.008912081830203533i,
        0.010422983206808567 - 0.025892846286296844i,
    };
    CheckBaseband(re, im, gold, "varying_prp");
}

// --- test_simc_tx_offset: tx channel at (5,0,0) ---
void TestTxOffset() {
    std::vector<H> f = {24.075e9, 24.175e9};
    std::vector<H> t = {0.0, 80e-6};
    std::vector<H> foff = {0.0, 0.0, 0.0};
    std::vector<H> pstart = {0.0, 100e-6, 200e-6};
    auto tx = std::make_shared<Transmitter<H, L>>(L(10), f, t, foff, pstart);
    AddDefaultTxChannel(*tx, rsv::Vec3<L>(5, 0, 0), L(0), 3, {}, {}, {});
    auto rx = std::make_shared<Receiver<L>>(L(6e4), L(20), L(500), L(30),
                                            L(6e4), 0.0);
    AddDefaultRxChannel(*rx, rsv::Vec3<L>(0, 0, 0));
    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot = {rsv::Vec3<L>(0, 0, 0)};
    auto radar = std::make_shared<Radar<H, L>>(tx, rx, frame_start, loc,
                                               rsv::Vec3<L>(0, 0, 0), rot,
                                               rsv::Vec3<L>(0, 0, 0));
    auto points = std::make_shared<PointsManager<L>>();
    points->AddPointSimple(rsv::Vec3<L>(10, 0, 0), rsv::Vec3<L>(0, 0, 0),
                           L(20), L(0));
    std::vector<H> re(12), im(12);
    radar->InitBaseband(re.data(), im.data());
    PointSimulator<H, L, radarsimx::cpu_policy> sim;
    Check(sim.Run(radar, points) == SUCCESS, "tx offset Run");
    const std::complex<double> g0 =
        -0.04858788475394249 - 0.027421100065112114i;
    std::vector<std::complex<double>> gold(12);
    const std::complex<double> row[4] = {
        g0,
        -0.03965778648853302 - 0.03924231231212616i,
        -0.027931272983551025 - 0.04829641059041023i,
        -0.014235264621675014 - 0.05394493788480759i,
    };
    for (int p = 0; p < 3; ++p) {
        for (int s = 0; s < 4; ++s) {
            gold[p * 4 + s] = row[s];
        }
    }
    CheckBaseband(re, im, gold, "tx_offset");
}

// --- test_simc_waveform_modulation: periodic ZOH table ---
void TestWaveformModulation() {
    std::vector<H> f = {24.075e9, 24.175e9};
    std::vector<H> t = {0.0, 80e-6};
    std::vector<H> foff = {0.0, 0.0, 0.0};
    std::vector<H> pstart = {0.0, 100e-6, 200e-6};
    auto tx = std::make_shared<Transmitter<H, L>>(L(10), f, t, foff, pstart);
    std::vector<L> mod_t = {0.0f, 10e-6f, 20e-6f, 30e-6f, 40e-6f};
    std::vector<std::complex<L>> mod_var = {
        {0, 0}, {0, 1}, {0, 0}, {0, -3}, {-4, 0}};
    AddDefaultTxChannel(*tx, rsv::Vec3<L>(0, 0, 0), L(0), 3, {}, mod_t,
                        mod_var);
    auto rx = std::make_shared<Receiver<L>>(L(6e4), L(20), L(500), L(30),
                                            L(6e4), 0.0);
    AddDefaultRxChannel(*rx, rsv::Vec3<L>(0, 0, 0));
    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot = {rsv::Vec3<L>(0, 0, 0)};
    auto radar = std::make_shared<Radar<H, L>>(tx, rx, frame_start, loc,
                                               rsv::Vec3<L>(0, 0, 0), rot,
                                               rsv::Vec3<L>(0, 0, 0));
    auto points = std::make_shared<PointsManager<L>>();
    points->AddPointSimple(rsv::Vec3<L>(10, 0, 0), rsv::Vec3<L>(0, 0, 0),
                           L(20), L(0));
    std::vector<H> re(12), im(12);
    radar->InitBaseband(re.data(), im.data());
    PointSimulator<H, L, radarsimx::cpu_policy> sim;
    Check(sim.Run(radar, points) == SUCCESS, "waveform mod Run");
    std::vector<std::complex<double>> gold(12, {0.0, 0.0});
    const std::complex<double> mv =
        -0.08509276807308197 + 0.07218044251203537i;
    for (int p = 0; p < 3; ++p) {
        gold[p * 4 + 2] = mv;
    }
    CheckBaseband(re, im, gold, "waveform_modulation");
}

// --- test_simc_interference (ideal): burst at sample 24 ---
void TestInterference() {
    // victim: chirp up, 1 pulse, fs 6e5 -> 48 samples
    std::vector<H> f = {24.075e9, 24.175e9};
    std::vector<H> t = {0.0, 80e-6};
    std::vector<H> foff = {0.0};
    std::vector<H> pstart = {0.0};
    auto tx_v = std::make_shared<Transmitter<H, L>>(L(10), f, t, foff, pstart);
    AddDefaultTxChannel(*tx_v, rsv::Vec3<L>(0, 0, 0), L(0), 1, {}, {}, {});
    auto rx_v = std::make_shared<Receiver<L>>(L(6e5), L(20), L(500), L(30),
                                              L(6e5), 0.0);
    AddDefaultRxChannel(*rx_v, rsv::Vec3<L>(0, 0, 0));
    std::vector<H> frame_start = {0.0};
    std::vector<rsv::Vec3<L>> loc0 = {rsv::Vec3<L>(0, 0, 0)};
    std::vector<rsv::Vec3<L>> rot0 = {rsv::Vec3<L>(0, 0, 0)};
    auto radar = std::make_shared<Radar<H, L>>(tx_v, rx_v, frame_start, loc0,
                                               rsv::Vec3<L>(0, 0, 0), rot0,
                                               rsv::Vec3<L>(0, 0, 0));

    // interferer: chirp down, 3 pulses, at (20,0,0) rotated 180 deg (radians)
    std::vector<H> fi = {24.175e9, 24.075e9};
    std::vector<H> foff3 = {0.0, 0.0, 0.0};
    std::vector<H> pstart3 = {0.0, 100e-6, 200e-6};
    auto tx_i = std::make_shared<Transmitter<H, L>>(L(10), fi, t, foff3,
                                                    pstart3);
    AddDefaultTxChannel(*tx_i, rsv::Vec3<L>(0, 0, 0), L(0), 3, {}, {}, {});
    std::vector<rsv::Vec3<L>> loc_i = {rsv::Vec3<L>(20, 0, 0)};
    std::vector<rsv::Vec3<L>> rot_i = {
        rsv::Vec3<L>(L(kPi), L(0), L(0))};
    auto iradar = std::make_shared<Radar<H, L>>(tx_i, rx_v, frame_start, loc_i,
                                                rsv::Vec3<L>(0, 0, 0), rot_i,
                                                rsv::Vec3<L>(0, 0, 0));

    std::vector<H> re(48, 0), im(48, 0);
    radar->InitBaseband(re.data(), im.data());
    InterferenceSimulator<H, L, radarsimx::cpu_policy> sim;
    Check(sim.Run(radar, iradar) == SUCCESS, "interference Run");
    std::vector<std::complex<double>> gold(48, {0.0, 0.0});
    gold[24] = -0.013252748176455498 + 0.004348371643573046i;
    CheckBaseband(re, im, gold, "interference");
}

// --- noise: correlation structure and determinism ---
void TestNoise() {
    auto radar = MakeRadar24G(L(6e4), 3, L(100e-6));
    const int pulses = 3, samples = 4, ch = 1, frames = 1;
    std::vector<H> ts(ch * pulses * samples);
    for (int p = 0; p < pulses; ++p) {
        for (int s = 0; s < samples; ++s) {
            ts[p * samples + s] = p * 100e-6 + s / 6e4;
        }
    }
    std::vector<H> nre(ch * pulses * samples), nim(ch * pulses * samples);
    NoiseSimulator<H, L, radarsimx::cpu_policy> sim;
    Check(sim.Run(radar, H(1e-3), true, ts.data(), ch, pulses, samples,
                  nre.data(), nim.data(), 0) == SUCCESS,
          "noise Run");
    // nonzero
    bool any = false;
    for (auto v : nre) {
        any = any || v != 0.0;
    }
    Check(any, "noise nonzero");
    // deterministic: second run identical
    std::vector<H> nre2(ch * pulses * samples), nim2(ch * pulses * samples);
    sim.Run(radar, H(1e-3), true, ts.data(), ch, pulses, samples, nre2.data(),
            nim2.data(), 0);
    Check(nre == nre2 && nim == nim2, "noise deterministic");
    // real mode: imag zero
    sim.Run(radar, H(1e-3), false, ts.data(), ch, pulses, samples, nre2.data(),
            nim2.data(), 0);
    bool imag_zero = true;
    for (auto v : nim2) {
        imag_zero = imag_zero && v == 0.0;
    }
    Check(imag_zero, "noise real mode imag zero");
}

}  // namespace

int main() {
    TestSingleTarget();
    TestVaryingPrp();
    TestTxOffset();
    TestWaveformModulation();
    TestInterference();
    TestNoise();
    if (g_failures == 0) {
        std::printf("test_point: all passed\n");
        return 0;
    }
    std::printf("test_point: %d failures\n", g_failures);
    return 1;
}
