# 600-cell Regge lab

A reproducible research notebook for classical Regge dynamics on the regular
600-cell and its refinements.

This is not a completed theory of nature. Every public claim is labelled
**DERIVED**, **STRUCTURAL**, **PATTERN**, or **OPEN**. Numerical agreement is
not treated as derivation, and failed preregistered routes remain part of the
evidence.

## Current result

The strongest current statement is finite and deliberately narrow.

For the fixed homogeneous tetrahedral-frustum 600-cell action with zero
cosmological constant, conserved global dust, positive proper heights and
positive endpoint scales:

1. the exact one-slab canonical equations admit a nontrivial finite-height
   update on a derived interval of incoming states;
2. composing the update at the frozen state `v=3/2` initially produces two
   physical second slabs, so the action alone defines a set-valued relation;
3. causality, future orientation, real-action branch membership and local
   Legendre regularity do not select between those two branches;
4. one-step future extendibility removes one branch: branch A has no physical
   third slab, while branch B has exactly one;
5. the surviving history has exactly one physical fourth and fifth slab;
6. an exact scale-free compactification gives a continuous family of
   asymptotic fixed points, not one selected universal point.

The third-, fourth- and fifth-slab classifications were reproduced by
mechanically different eliminations and direct-action reconstructions.  The
asymptotic coefficients were also obtained both by a symbolic series and by
derivatives at the compactified boundary.

This is **DERIVED EXACT** for the scale-free limiting family and **DERIVED
COMPUTATIONAL / STRUCTURAL, five-slab scoped** for the finite history. It does
not prove an infinite history or a fundamental time quantum.

The latest physical continuation is

```text
q5    = 1006.53493784425818414892...
h5    = 0.002163977529932147004...
L5/L4 = 3.17811898858662493889...
```

The limiting map explains the large-slope self-similarity but also supplies a
negative: it fixes a continuum of boundary states and therefore selects no
universal ratio. The active gate is now an invariant-region theorem for the
exact normalized map, or a finite counterexample. Another isolated slab is
not enough.

## Independent geometric control

A separate preregistered refinement calculation compares the homogeneous
closed-dust acceleration with Friedmann:

```text
fixed 600-cell error     7.8979%
first refinement error  1.9134%
second refinement error 0.4744%
```

All four frozen central-octahedron diagonal conventions improve at both
levels. The near-factor-four reduction is a **PATTERN**, not an
infinite-refinement theorem. It supports the interpretation of the coarse
coefficient as discretization error, not new cosmological acceleration.

## What is not derived

The current repository does not derive:

- local general relativity or propagating gravitational waves;
- indefinite deterministic evolution;
- an absolute tick or limiting speed;
- `G`, Planck length, Planck time, or particle masses;
- a quantum theory of the Regge histories.

Global scale covariance already excludes an absolute classical tick from the
scale-free action without an additional dimensionful input.

## Start here

- [Authoritative current status](docs/gravity/CURRENT_STATUS.md)
- [Scale-free map and fifth-slab result](docs/gravity/gravity_600cell_finite_height_asymptotic_map_result.md)
- [Fourth-slab result](docs/gravity/gravity_600cell_finite_height_fourth_slab_result.md)
- [Third-slab branch selection](docs/gravity/gravity_600cell_finite_height_third_slab_result.md)
- [Two-slab nonuniqueness](docs/gravity/gravity_600cell_finite_height_composition_result.md)
- [Exact finite-height classification](docs/gravity/gravity_600cell_finite_height_classification_result.md)
- [Scale-covariance no-go](docs/gravity/gravity_600cell_tick_scale_covariance_result.md)
- [Projected-refinement comparison](docs/gravity/gravity_600cell_projected_refinement_acceleration_comparison_result.md)
- [Documentation index](docs/README.md)
- [Verifier index](reproducible/README.md)

## Reproduce the latest result

Use the project interpreter:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_asymptotic_map.py

/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_asymptotic_map_adversarial.py
```

Expected summaries are `11/11 PASS` and `12/12 PASS`. No full-suite result is
claimed for the current checkpoint.

## Repository layout

```text
commons/          shared exact 600-cell construction
docs/gravity/     prior-art gates, preregistrations, corrections and results
reproducible/     registered gravity verifiers and frozen artifacts
CLAUDE.md         binding scientific-method rules
```

Legacy fitted particle-physics papers, parameter scans and promotional
visualizations have been removed from the current tree. The deletion is
recoverable from Git history; the history was not rewritten.
