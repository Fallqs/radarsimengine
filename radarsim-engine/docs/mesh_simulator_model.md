# Mesh (SBR + PO) Simulator — Model Findings

Status of the mesh/RCS reverse-engineering (see also `point_simulator_model.md`).

## Confirmed model structure (validated against goldens to ~2-3%)

The mesh baseband is a Physical-Optics surface sum with spherical (near-field)
propagation — NOT a far-field plate formula. For the 5x5 m plate at 10 m
(Fresnel number ~200, deep near field), the phase across the aperture spans
~100 pi rad, so contributions mostly cancel; the result is set by the
stationary-phase region.

Per PO surface sample s (one per ray landing):

    bb += stuff_s * K * (k/2) * exp(j * 2 pi * beat_s) / (R_tx_s * R_rx_s) * dA_s

- `K` = the point-target chain constant (tx power, patterns, rf/bb gains,
  load resistor — the exact amplitude chain validated by the point suite).
- `k/2 = pi / lambda` — pinned to <2% by test_scene_single_target
  (|bb| 0.04724 vs 0.04630 model).
- `beat_s = phi(u) - phi(u - tau_s)` — the same deramp beat as point targets,
  per-sample round-trip delay `tau_s`.
- The PO prefactor carries `j` (+0.25 cycles of phase) and NO Fresnel sign
  flip for PEC (golden phase 0.6152 vs integral+j 0.6170).
- density only matters through the sampling; the golden at density 0.4 sits
  ~2% off the converged integral, i.e. the engine's sampler is nearly
  converged there.

## The open problem: the exact ray sampler

The mesh test tolerance (3e-5 of peak) makes the golden a recording of the
engine's exact discrete sampling. Established facts about the sampler:

- rays launch from the Tx through an occupancy grid of `channel.grid` degrees
  (default 1°);
- `density` = rays per wavelength sets the per-cell sub-ray count;
- each landing = one PO sample with an area weight;
- `environment` targets get a coarser share; `skip_diffusion` cells are culled
  before launch.

What is NOT yet pinned: the per-cell sub-ray count rounding rule, the sub-ray
layout within a cell, the area weight's obliquity convention, and edge
truncation. The RCS goldens (48.3 dBsm at 1 GHz for the 5x5 plate vs 49.41
ideal) additionally show the RCS sampler covers 88% of the plate at density 1
in a frequency-independent way that no simple origin- or corner-anchored
lambda-grid reproduces.

## The plan

Implement the mesh simulator with a clean, documented sampler (angular
occupancy grid + per-cell sub-ray fans, Gordon/PO per sample), validate
against the *converged* integral and physical invariants (corner-reflector
multibounce, plate near-field, sphere ~pi r^2), and keep golden-exact
matching as a tracked divergence. The Python reference harness
(`tools/pysim/`) gets a mesh path first for fast iteration.

## Data assets for this work

- `tests/test_module_sim_radar_mesh.py` — 22 golden scenarios.
- `benchmarks/baseline/cpu_reference.npz` — six full baseband captures with
  documented scenes (`benchmarks/capture_reference.py`, `benchmarks/scenes.py`).
- `references/` — Gordon 1975 (PO far field), Hwang 2015 (PO RCS), Karras
  2012 (GPU BVH).
