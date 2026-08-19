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

## Active calculation

The prior-art gate is commit `faf5b08`, the target-blind protocol is commit
`36ddebd`, and the registered source is commit `d259477`.  Before proposing
another carrier, compute modular projection dimensions for the 240-space:

1. the compatible subspace with all strut variations set to zero;
2. the compatible subspace with all upper-edge variations set to zero;
3. the ranks of the projections onto the 720 upper-edge and 120 strut data.

For `A=[F E S]`, with `rank(F)=3600`, these follow from the ranks of
`[F E]`, `[F S]`, and `[F E S]`.  This distinguishes a graph over lapse data,
a graph over spatial data, independent sectors, or a more entangled carrier
without guessing a target.  It requires a new prior-art note, protocol commit,
registered targeted verifier, and an adversarial implementation before a
material interpretation is accepted.  The protocol forbids comparison with
119, 120, or any proposed geometric carrier until the target-blind JSON has
been committed.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The last independently reported full baseline was 79/79,
but it is not evidence for later files.  Before stopping, update this file with
the last accepted commit, the current open question, and the exact next test.
