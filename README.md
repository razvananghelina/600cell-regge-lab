# 600-cell Regge lab

Reproducible finite-geometry and Regge-calculus investigations built around
the regular 600-cell and explicitly labelled as **DERIVED**, **PATTERN**,
**FITTING**, or **OPEN**.

This repository is a research notebook, not a completed theory of nature.
Positive numerical patterns are kept separate from the claims actually proved
by the registered verifiers, and failed routes are retained as evidence.

## Current physical result

The best-supported dynamical result is deliberately narrow:

- the fixed regular 600-cell has an exact homogeneous closed-dust Regge tick
  law;
- its weak-lapse acceleration differs from closed Friedmann by about `7.90%`;
- on two projected red refinements, the direct non-averaged irregular Regge
  action reduces that error to about `1.91%` and `0.474%`;
- the reduction is close to a factor of four at each level, consistently
  across four preregistered resolutions of the octahedral-diagonal ambiguity.

This is **DERIVED NUMERICAL** finite evidence that the homogeneous result is a
Regge discretization of closed dust Friedmann.  The apparent second-order rate
is a **PATTERN**, not an infinite-refinement theorem.

It does not yet establish local general relativity, gravitational waves, a
limiting speed, a physical tick size, Planck units, or particle masses.  The
next discriminating calculation is the gauge-reduced quadratic action for
inhomogeneous edge perturbations.

## Start here

- [Latest refinement comparison](gravity_600cell_projected_refinement_acceleration_comparison_result.md)
- [Blind coefficients frozen before target comparison](gravity_600cell_projected_refinement_acceleration_blind_result.md)
- [Exact all-tick homogeneous theorem](gravity_600cell_cellular_weak_lapse_all_n_result.md)
- [Dimension reconciliation](dimension_reconciliation.md)
- [Consolidated status](consolidation_summary.md)
- [Step-by-step theory ledger](theory_step_by_step_master.md)
- [Reproducibility registry](reproducible/README.md)

The repository contains older speculative directions as well as later
refutations.  Read the dated status/result notes and executable verifier
outputs rather than inferring validity from a filename.

## Reproduce a targeted result

Use the project interpreter:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_projected_refinement_acceleration_comparison.py
```

That comparison consumes a coefficient artifact committed before its
continuum target was evaluated.  The relevant commit order is:

```text
98a1e1d  Freeze refined Regge coefficients before target comparison
33ee717  Register disclosed Regge refinement comparison
c10be7f  Compare refined Regge acceleration with Friedmann
```

The full registry is in `reproducible/run_all.py`, but individual missions
should normally run their targeted verifier while being developed.

## Repository layout

```text
commons/              shared exact 600-cell constructions
reproducible/         registered verifiers and frozen artifacts
legacy/experiments/   historical exploratory exp*.py scripts
*.md                   protocols, prior-art gates, results, and ledgers
```

The historical experiments were moved out of the root only to keep GitHub's
directory view below its 1,000-entry display limit.  Their Git history is
unchanged; they are not promoted to verified evidence by being archived.
