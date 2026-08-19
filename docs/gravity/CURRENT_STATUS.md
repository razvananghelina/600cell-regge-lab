# Current gravity status

Updated: 2026-08-19

This is the authoritative resume bookmark for the active gravity route.  Read
`CLAUDE.md`, this file, `git status --short`, and the commits listed below
before starting another calculation.  Do not infer acceptance from an
uncommitted artifact, a dimension match, or a modular rank alone.

## Last accepted checkpoint

Commit `e3a77fe` froze the target-blind canonical-data admissibility census for
the complete variable-face flat-frustum 600-cell slab.  At both
`(lambda,tau)=(2,5)` and `(3,11)`, and modulo both primes 1000003 and 1000033,

```text
fixed cell-flex rank       3600 / 3600
full augmented rank        4200 / 4440
full augmented nullity      240
```

**DERIVED (modular).** The fixed-data problem is injective and the augmented
system has a stable 240-dimensional modular kernel.

**OPEN.** The rational kernel dimension and its exact geometric carrier have
not yet been proved.  The modular 240 is an upper bound on the rational
dimension, not an accepted rational equality.

Primary files:

- `reproducible/verify_gravity_600cell_canonical_data_admissibility.py`
- `reproducible/gravity_600cell_canonical_data_admissibility.json`
- `docs/gravity/gravity_600cell_canonical_data_admissibility_protocol.md`

## Most recent falsification

The tempting interpretation

```text
240 = 120 vertex-scale modes + 120 vertex-lapse modes
```

was tested by exact substitution, not accepted from its dimension.  It fails
on all 3600 face rows for each of the two representatives, with the same exact
first residual under two different local right-inverse graphs.

**DERIVED NEGATIVE.** `Q^120 scale + Q^120 lapse` is not the complete
admissible carrier.

The first execution was preserved in commit `829bf29`.  Two verdict-only
classifier corrections were preregistered and committed separately; no
scientific matrix or residual was changed.  The final targeted run passed
10/10 with the refuted outcome and was frozen in commit `96b912d`.

Primary files:

- `docs/gravity/gravity_600cell_canonical_data_carrier_first_failure.md`
- `docs/gravity/gravity_600cell_canonical_data_carrier_classifier_correction.md`
- `reproducible/verify_gravity_600cell_canonical_data_carrier.py`
- `reproducible/gravity_600cell_canonical_data_carrier.json`

## What the gravity route has established

- **DERIVED.** The carrier is the boundary of the regular 600-cell, a
  triangulated 3-sphere with `f=(120,720,1200,600)`.
- **DERIVED.** The homogeneous cellular weak-lapse Regge equations give the
  all-step discrete Friedmann relation, but do not select an absolute tick.
- **PATTERN with exact finite controls.** Canonical projected refinements move
  the homogeneous acceleration toward the dust-FLRW value `-1/2`.
- **DERIVED.** A schedule-free local flat tetrahedral frustum has six shape
  flexes; variable face transitions are necessary for global closure.
- **DERIVED.** On unequal-scale global slabs the fixed-data system has full
  rank 3600.  On the static slab its 119-dimensional kernel is exactly the old
  vertex-gradient/prism-shift carrier, not a new physical sector.
- **DERIVED NEGATIVE.** The first natural 240-dimensional vertex
  scale/lapse carrier is not compatible globally.

## What is not yet derived

- **OPEN:** the exact rational meaning of the 240 modular modes;
- **OPEN:** the admissible symplectic/canonical phase space;
- **OPEN:** an action Hessian restricted to that phase space;
- **OPEN:** propagating tensor modes and a dispersion relation;
- **OPEN:** a dynamically selected lapse or physical tick duration;
- **OPEN:** an effective universal speed `c`;
- **OPEN:** Newton's constant `G`, an absolute length, and the Planck scales;
- **OPEN:** any particle mass formula or Standard-Model derivation from this
  gravity calculation.

## Active result and next calculation

The prior-art gate is commit `faf5b08`, the projection protocol is commit
`36ddebd`, and the registered source is commit `d259477`.  The first census
passed 11/11 and was frozen without a carrier comparison in commit `9b97775`.
All primes, representatives, graphs, and conventions give:

```text
(kernel, edge-only, strut-only, edge projection, strut projection)
= (240, 120, 120, 120, 120).
```

**DERIVED (modular).** The compatible data split as a direct sum of a
120-dimensional edge-only sector and the entire 120-dimensional strut sector.

**STRUCTURAL PROTOCOL DEVIATION.** Commit `a3fe2b9` records that a speculative
`119+121` target was disclosed after the complete source was committed but
before the artifact.  The result falsified that target, but strict blindness
is not claimed.

The next exact test is:

1. construct the derived unsigned 720-by-120 vertex-edge scale image;
2. test exact data-space inclusion/equality with the 120-dimensional edge
   sector rather than infer it from dimension;
3. retain arbitrary struts as the already established modular strut sector;
4. if data equality holds, construct the globally solved rational cell-flex
   lifts, which differ from the refuted local lift.

This requires a new target-disclosed protocol and a mechanically independent
implementation before any material interpretation is accepted.  The current
result does not select a tick: every strut variation is kinematically allowed
modulo both primes.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The last independently reported full baseline was 79/79,
but it is not evidence for later files.  Before stopping, update this file with
the last accepted commit, the current open question, and the exact next test.
