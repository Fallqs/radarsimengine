"""Shim that loads the pure-Python parts of radarsimpy without the compiled
extensions, so the reference implementation can build real Radar/Transmitter/
Receiver objects.

Usage:
    import shim
    shim.install()                      # registers the "radarsimpy" package
    from radarsimpy import Radar, Transmitter, Receiver
"""

import sys
import types
from pathlib import Path


def install(src_dir=None):
    """Create the ``radarsimpy`` package module from source without executing
    its ``__init__.py`` (which imports the compiled simulator/license)."""
    if src_dir is None:
        src_dir = Path(__file__).resolve().parents[3] / "src"
    pkg_dir = Path(src_dir) / "radarsimpy"
    if not pkg_dir.is_dir():
        raise RuntimeError(f"radarsimpy package not found at {pkg_dir}")

    pkg = types.ModuleType("radarsimpy")
    pkg.__path__ = [str(pkg_dir)]
    sys.modules["radarsimpy"] = pkg

    from radarsimpy.radar import Radar  # noqa: E402
    from radarsimpy.receiver import Receiver  # noqa: E402
    from radarsimpy.transmitter import Transmitter  # noqa: E402

    pkg.Radar = Radar
    pkg.Receiver = Receiver
    pkg.Transmitter = Transmitter

    # Stub the compiled simulator module; callers replace sim_radar.
    simulator = types.ModuleType("radarsimpy.simulator")

    def _unavailable(*args, **kwargs):
        raise NotImplementedError("compiled simulator not available in shim")

    simulator.sim_radar = _unavailable
    simulator.sim_lidar = _unavailable
    simulator.sim_rcs = _unavailable
    simulator.gpu_available = lambda: False
    sys.modules["radarsimpy.simulator"] = simulator
    pkg.simulator = simulator

    return pkg
