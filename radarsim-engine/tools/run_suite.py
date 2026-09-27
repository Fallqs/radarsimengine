"""Run an upstream test suite against the Python reference implementation.

Usage: python radarsim-engine/tools/run_suite.py <test_module_basename> [filter]
e.g.:  python radarsim-engine/tools/run_suite.py test_noise_simulation
"""

import sys
import importlib
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pysim import shim  # noqa: E402

shim.install()

from pysim import reference  # noqa: E402
import radarsimpy  # noqa: E402

radarsimpy.simulator.sim_radar = reference.sim_radar_reference


def main():
    module_name = sys.argv[1]
    name_filter = sys.argv[2] if len(sys.argv) > 2 else ""
    suite = importlib.import_module(f"tests.{module_name}")
    tests = [
        (name, fn)
        for name, fn in vars(suite).items()
        if name.startswith("test_") and callable(fn) and name_filter in name
    ]
    n_pass = n_fail = n_error = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS {name}")
            n_pass += 1
        except AssertionError as exc:
            print(f"FAIL {name}: {str(exc)[:200]}")
            n_fail += 1
        except Exception as exc:  # noqa: BLE001
            print(f"ERROR {name}: {type(exc).__name__}: {str(exc)[:200]}")
            if "-v" in sys.argv:
                traceback.print_exc()
            n_error += 1
    print(f"\n{n_pass} passed, {n_fail} failed, {n_error} errors, {len(tests)} collected")
    return 1 if (n_fail or n_error) else 0


if __name__ == "__main__":
    sys.exit(main())
