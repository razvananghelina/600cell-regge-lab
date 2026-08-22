# 600-cell Regge lab

A reproducible research notebook for classical Regge dynamics on the regular
600-cell and its refinements.

This is not a completed theory of nature. Every public claim is labelled
**DERIVED**, **STRUCTURAL**, **PATTERN**, or **OPEN**. Numerical agreement is
not treated as derivation, and failed preregistered routes remain part of the
evidence.

## Current result

The strongest current statement is global in discrete step number but still
deliberately narrow in physical scope.

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
   asymptotic fixed points, not one selected universal point;
7. a rigorously certified invariant half-strip contains the accepted branch,
   proving a unique physical successor at every later finite step.
8. a preregistered census over the complete one-slab incoming domain finds 36
   distinct depth-four tree signatures rather than one generic pattern; 1080
   of 3072 diagnostic inputs retain at least one branch outside the proved
   invariant region at the frozen depth.
9. two mechanically different rigorous root certificates prove that the
   `DEAD+ENTERED_D` tree at `v=3/2` persists on some nonzero open neighbourhood
   of the incoming curve, although no radius is selected.

The third-, fourth- and fifth-slab classifications were reproduced by
mechanically different eliminations and direct-action reconstructions.  The
asymptotic coefficients were also obtained both by a symbolic series and by
derivatives at the compactified boundary.

This is **DERIVED EXACT** for the scale-free limiting family and **DERIVED
COMPUTATIONAL / STRUCTURAL, representative-seed scoped** for the forward
history.  The invariant-region theorem proves, by induction, a unique
successor at every later finite step in the frozen homogeneous model.  It
does not prove convergence, infinite total proper duration, completeness, a
local evolution law or a fundamental time quantum.

The latest physical continuation is

```text
q5    = 1006.53493784425818414892...
h5    = 0.002163977529932147004...
L5/L4 = 3.17811898858662493889...
```

The limiting map explains the large-slope self-similarity but also supplies a
negative: it fixes a continuum of boundary states and therefore selects no
universal ratio. The invariant-region theorem removes the finite-horizon
weakness. The incoming census shows that the `v=3/2` pattern persists at
adjacent diagnostic inputs but is not representative of the full domain. The
continuum basin remains **OPEN**: the next mathematical gate is an interval
certificate, while the next physical gate is a geometry-selected
nonhomogeneous perturbation reconciled with the existing anisotropic no-go
results.
Another isolated slab is not enough.

## Latest nonhomogeneous canonicity gate

At the first positive-height slab from `v=3/2`, a preregistered primary
calculation and a mechanically different 180-digit orbit-kernel replication
find no resolved difference between the two staircase parities in the
one-sided quadratic form on the exact rank-240 scale-plus-strut tangent
carrier.  The adversarial verifier passes `18/18`; direct complete-action
second derivatives reproduce the pulled-back forms to `1e-85` or better.

The next rank gate has also been completed by two mechanically different
methods.  The `720 x 240` diagonal-only map has rank `119`; adding the 120
pole equations gives rank `239`.  The sole complete survivor is the same
homogeneous line in both parities, while every nonhomogeneous sector has full
column rank.  Direct nonlinear action-gradient secants corroborate the line.

An exact registered reconciliation now proves that the survivor is precisely
the homogeneous lapse-constraint tangent `dC=0`.  The fixed incoming
canonical momentum has nonzero derivative on it, so the complete fixed-input
kernel is zero.  This is a **DERIVED BOUNDED NEGATIVE** for interpreting the
line as a tick or free mode, not a claim of zero dynamics.

That next forced response has now been constructed.  A preregistered
group-reduced complex-ball calculation and a mechanically different complete
real-space calculation both find a regular action-generated 1440-dimensional
canonical boundary map at the first positive-height slab.  Its symplectic
identities pass, and no dependence on the two staircase schedules is
resolved.  The targeted verifiers pass `21/21` and `22/22`.

This is a **DERIVED COMPUTATIONAL** linearized map, not yet a physical mode
spectrum.  It connects tangent fibres at different background scales.  The
next gate is to derive a co-moving canonical trivialization, build the map on
the next accepted slab and compose the two.  No graviton, dispersion law,
physical tick, `c`, `G` or Planck scale has been derived.

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
- generic deterministic evolution over the original incoming states;
- an absolute tick or limiting speed;
- `G`, Planck length, Planck time, or particle masses;
- a quantum theory of the Regge histories.

Global scale covariance already excludes an absolute classical tick from the
scale-free action without an additional dimensionful input.

## Start here

- [Authoritative current status](docs/gravity/CURRENT_STATUS.md)
- [Invariant half-strip and finite-step forward continuation](docs/gravity/gravity_600cell_finite_height_invariant_region_result.md)
- [Incoming-state candidate skeleton](docs/gravity/gravity_600cell_finite_height_incoming_basin_discovery_result.md)
- [Adversarially corroborated local branch theorem](docs/gravity/gravity_600cell_finite_height_local_signature_result.md)
- [Scale-free map and fifth-slab result](docs/gravity/gravity_600cell_finite_height_asymptotic_map_result.md)
- [Fourth-slab result](docs/gravity/gravity_600cell_finite_height_fourth_slab_result.md)
- [Third-slab branch selection](docs/gravity/gravity_600cell_finite_height_third_slab_result.md)
- [Two-slab nonuniqueness](docs/gravity/gravity_600cell_finite_height_composition_result.md)
- [Exact finite-height classification](docs/gravity/gravity_600cell_finite_height_classification_result.md)
- [Scale-covariance no-go](docs/gravity/gravity_600cell_tick_scale_covariance_result.md)
- [Finite-height internal-carrier rank](docs/gravity/gravity_600cell_finite_height_internal_carrier_rank_result.md)
- [Canonical meaning of the internal kernel](docs/gravity/gravity_600cell_finite_height_internal_kernel_canonical_reconciliation_result.md)
- [Projected-refinement comparison](docs/gravity/gravity_600cell_projected_refinement_acceleration_comparison_result.md)
- [Documentation index](docs/README.md)
- [Verifier index](reproducible/README.md)

## Reproduce the latest result

Use the project interpreter:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_invariant_region.py

/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_invariant_region_adversarial_resolution.py
```

Expected summaries are `14/14 PASS` and `12/12 PASS`. No full-suite result is
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
