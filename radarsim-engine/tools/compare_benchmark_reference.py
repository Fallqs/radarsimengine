"""Diff radarsim-engine against the recorded benchmark reference captures.

Runs the capture_reference.py scenes against the integrated package and
compares with benchmarks/baseline/cpu_reference.npz. Report per-scene
peak-relative max error.

Usage: python radarsim-engine/tools/compare_benchmark_reference.py
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks"))

import numpy as np  # noqa: E402
import scenes  # noqa: E402
from radarsimpy.simulator import sim_radar  # noqa: E402

# mirrors benchmarks/capture_reference.py CASES
CASES = {
    "plate_normal": {"radar": {"pulses": 2, "samples": 20},
                     "targets": {"model": "plate5x5", "distance": 20.0},
                     "sim": {"density": 0.3}},
    "corner_multibounce": {"radar": {"pulses": 2, "samples": 20},
                           "targets": {"model": "cr", "distance": 10.0},
                           "sim": {"density": 1.0}},
    "sphere_grazing": {"radar": {"pulses": 2, "samples": 20},
                       "targets": {"model": "ball_1m", "distance": 20.0},
                       "sim": {"density": 0.6}},
    "sphere_mimo": {"radar": {"pulses": 2, "samples": 20, "tx_channels": 2,
                              "rx_channels": 4},
                    "targets": {"model": "ball_1m", "distance": 20.0},
                    "sim": {"density": 0.4}},
    "turbine_moving": {"radar": {"pulses": 4, "samples": 20},
                       "targets": {"model": "turbine", "distance": 15.0,
                                   "speed": (-5, 0, 0)},
                       "sim": {"density": 0.5, "level": "pulse"}},
    "sphere_backprop": {"radar": {"pulses": 2, "samples": 20},
                        "targets": {"model": "ball_1m", "distance": 20.0},
                        "sim": {"density": 0.5, "back_propagating": True}},
}


def main():
    ref = np.load(REPO / "benchmarks" / "baseline" / "cpu_reference.npz")
    names = sys.argv[1:] if len(sys.argv) > 1 else list(CASES)
    for name in names:
        case = CASES[name]
        radar = scenes.make_radar(**case["radar"])
        targets = scenes.make_targets(**case["targets"])
        # model paths are absolute in scenes.py; rebase to cwd for the engine
        data = sim_radar(radar, targets, **case["sim"])
        mine = data["baseband"]
        theirs = ref[name]
        peak = np.abs(theirs).max()
        err = np.abs(mine - theirs).max() / peak
        print(f"{name:22s} peak_rel_max_err = {err:.4f}  "
              f"(|mine| peak {np.abs(mine).max():.4e} vs {peak:.4e})")


if __name__ == "__main__":
    main()
