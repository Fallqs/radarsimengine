"""Run the upstream ideal-point-target test suite against the Python reference
implementation of the engine model.

Each test in tests/test_module_sim_radar_ideal.py builds real
Radar/Transmitter/Receiver objects (pure Python) and calls sim_radar, which we
replace with tools.pysim.reference.sim_radar_reference. The golden assertions
inside the tests are the acceptance criteria.

Usage: python radarsim-engine/tools/run_ideal_tests.py [test_name_substring]
"""

import sys
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pysim import shim  # noqa: E402

shim.install()

from pysim import reference  # noqa: E402

# Point the simulator stub at the reference implementation before the test
# module imports it.
import radarsimpy  # noqa: E402

radarsimpy.simulator.sim_radar = reference.sim_radar_reference

import tests.test_module_sim_radar_ideal as suite  # noqa: E402


def main():
    name_filter = sys.argv[1] if len(sys.argv) > 1 else ""
    tests = [
        (name, fn)
        for name, fn in vars(suite).items()
        if name.startswith("test_") and callable(fn) and name_filter in name
    ]
    n_pass = n_fail = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS {name}")
            n_pass += 1
        except Exception as exc:  # noqa: BLE001
            print(f"FAIL {name}: {type(exc).__name__}: {exc}")
            if "-v" in sys.argv:
                traceback.print_exc()
            n_fail += 1
    print(f"\n{n_pass} passed, {n_fail} failed, {len(tests)} collected")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
