#!/usr/bin/env python
"""Analytic validation of the RCS simulator (radarsim-engine).

Unlike the upstream goldens (which record the closed engine's lossy discrete
sampler), these checks compare against analytic PO references:

- plate 5x5 at broadside: 4 pi A^2 / lambda^2 (exact)
- sphere d=1 m: pi r^2 in the optical region; tolerance +-1.5 dB absorbs the
  triangulation and hard facet-level silhouette of ball_1m.stl

Run after tools/build_and_integrate.sh, from the repo root:
    python radarsim-engine/tools/validate_rcs_analytic.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from radarsimpy.simulator import sim_rcs

POL = np.array([0, 0, 1])
failures = 0


def check(name, got_db, want_db, tol_db):
    global failures
    ok = abs(got_db - want_db) <= tol_db
    print(f"{'ok ' if ok else 'FAIL'} {name}: {got_db:8.3f} dBsm "
          f"(ref {want_db:8.3f}, tol {tol_db} dB)")
    if not ok:
        failures += 1


# plate broadside: exact
for f in [1e9, 3e9, 10e9]:
    rcs = sim_rcs([{"model": "./models/plate5x5.stl", "location": (0, 0, 0)}],
                  f, 0, 90, inc_pol=POL, density=1)
    lam = 299792458.0 / f
    check(f"plate broadside {f/1e9:.0f} GHz", 10 * np.log10(rcs),
          10 * np.log10(4 * np.pi * 25.0**2 / lam**2), 0.01)

# sphere: optical-limit pi r^2, sweep 2-30 GHz
sphere_db = 10 * np.log10(np.pi * 0.5**2)
worst = 0.0
for fghz in np.arange(2.0, 30.1, 1.0):
    rcs = sim_rcs([{"model": "./models/ball_1m.stl", "location": (0, 0, 0)}],
                  fghz * 1e9, 0, 90, inc_pol=POL, density=1)
    err = 10 * np.log10(rcs) - sphere_db
    worst = max(worst, abs(err))
print(f"sphere 2-30 GHz: worst |err| = {worst:.3f} dB (tol 1.5 dB)")
if worst > 1.5:
    failures += 1

print("FAILURES:", failures)
raise SystemExit(1 if failures else 0)
