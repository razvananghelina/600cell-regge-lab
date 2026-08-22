# Pre-evaluation correction: explicit branch-B history provenance

Date: 2026-08-22

Original prior-art commit: `eae691a`.

Original protocol commit: `f80e174`.

Repository-prior-art correction commit: `ddb3f6e`.

Status: **CORRECTED BEFORE COMPLETING OR COMMITTING THE NEW VERIFIER AND
BEFORE EXECUTING ANY NEW SECOND-SLAB HESSIAN, PRE-LEGENDRE MATRIX, TANGENT OR
COMPOSITION.**

## Gap found during implementation audit

Section 2 of the frozen protocol requires the independently reconstructed
`q1,h1,r1,m1,pi1,q2,h2,r2,m2,pi2` history to agree with the committed branch-B
history.  The frozen asymptotic-map artifact is necessary later-history
provenance, but its stored `history` begins at `m3,pi3`; it does not contain
`q2,h2,r2`.  Therefore it cannot, by itself, falsify a wrong selection inside
the second-slab branch census.

The bracket `q2 in (31,32)` uniquely fixes the intended root of the frozen
equation, but that is a mathematical selector, not the separately requested
comparison with the committed branch-B record.  Treating the bracket as both
checks would silently drop one preregistered control.

## Frozen missing inputs

The already accepted primary and mechanically different composition
artifacts contain the exact branch-B record:

```text
reproducible/gravity_600cell_finite_height_composition.json
  d4e36141863bd2ae515b96eeeff4f50eb087016cca8cfb6f4b1e3355d6fba447

reproducible/gravity_600cell_finite_height_composition_adversarial.json
  d50e87f736e51585596aa1d7778238febaf7422840d668499878d8bd917f99e9

reproducible/verify_gravity_600cell_finite_height_composition.py
  cb4cf619dc54922d3a64d5e000a6cd0d3c19f71cda2a18b66d69f68271496422

reproducible/verify_gravity_600cell_finite_height_composition_adversarial.py
  8395e921ab1c1f518abb567a114f1eb8bfdf2068be031bff55c8d2f0cff56c2b
```

Require the primary artifact to retain `10/10` and
`FINITE_HEIGHT_TWO_SLAB_NONUNIQUE`.  Require the adversarial artifact to
retain `9/9` and
`FINITE_HEIGHT_TWO_SLAB_NONUNIQUE_ADVERSARIALLY_CORROBORATED`.

## Licensed comparison

After deterministic reconstruction but before any Hessian is assembled:

1. compare `q1,h1,r1` with both the accepted first-tangent background and the
   180-digit adversarial composition record;
2. identify the primary `v=3/2` physical root with `q2 in (31,32)` and the
   adversarial branch explicitly labelled `B`;
3. compare `q2,h2,r2` with both records below the already frozen `1e-65`
   bound;
4. retain the independent equation residual, bracket-width, positivity and
   canonical-junction checks from the original protocol;
5. derive `m1,pi1,m2,pi2` from the frozen recurrence and record them without
   rounding.

No background number from either composition artifact may be used as a root
seed, bracket endpoint, derivative step, tolerance or desired Hessian target.

## What does not change

This correction adds a missing provenance witness.  It changes no equation,
branch bracket, scale-lift identity, derivative step, uncertainty formula,
10/100 classifier, outcome hierarchy or interpretation boundary.  No new
second-slab Hessian or tangent value existed when the correction was made.
