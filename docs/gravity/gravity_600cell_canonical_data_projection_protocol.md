# Protocol: target-blind projection census of compatible canonical data

Date: 2026-08-19

This protocol is frozen before computing any partial augmented rank.  Read it
with `gravity_600cell_canonical_data_projection_prior_art.md` and the frozen
admissibility protocol/artifact.

## Scope and hypotheses

Use exactly the complete variable-face flat-frustum compatibility rows already
constructed by
`reproducible/verify_gravity_600cell_canonical_data_admissibility.py`.
Do not alter the 600-cell, local Jacobian, face transition, right inverse,
representatives, primes, or column definitions.

For each representative `(lambda,tau)=(2,5),(3,11)` and each frozen prime
`p=1000003,1000033`, compute

```text
rF   = rank_p(F)
rFE  = rank_p([F E])
rFS  = rank_p([F S])
rFES = rank_p([F E S])
```

and derive, without constructing or naming a target carrier,

```text
k             = 4440 - rFES
edge_only     = 4320 - rFE
strut_only    = 3720 - rFS
edge_projection  = k - strut_only
strut_projection = k - edge_only
```

Here `edge_only` means compatible data with all strut variations zero, while
`strut_only` means compatible data with all upper-edge variations zero.

## Target blindness

Before the first artifact containing these five dimensions is committed:

- do not compare them with the old static nullity 119;
- do not compare them with 120 vertices or any sum/difference built from it;
- do not propose a scale/lapse, gradient, displacement, cochain, gauge, or
  graviton decomposition;
- do not compute an explicit kernel basis selected to realize such a target;
- do not interpret a modular dimension as a rational theorem.

The artifact commit message must state that no carrier comparison was
performed.  Only a later, separate commit may disclose comparisons or propose
the next exact carrier.

## Construction and convention controls

1. Reproduce the frozen input hashes and the full `rF=3600`, `rFES=4200`
   ranks before accepting a partial rank.
2. Require `rF <= rFE <= rFES` and `rF <= rFS <= rFES`.
3. Require all five derived dimensions to satisfy the corresponding ambient
   bounds and the two rank-nullity identities above exactly.
4. Recompute the census with the alternate exact local right-inverse graph.
5. Recompute it after reversing every face orientation, an odd canonical
   relabelling, and reversing the metric sign.
6. Require equality across both primes, both nonstatic representatives, both
   right-inverse graphs, and all three legitimate convention attacks.
7. Test the quotient accounting on a preregistered synthetic fixture with two
   fixed columns, two edge columns, one strut column, and four equation rows:
   `F=(r0,r1)`, `E=(r0,r2)`, `S=(r2)`.  It must give
   `(k,edge_only,strut_only,edge_projection,strut_projection)=(2,1,0,2,1)`.
   In its negative control replace `S=(r2)` by `S=(r3)`; the cancellation is
   destroyed and the tuple must change to `(1,1,0,1,0)`.  This tests the
   implementation without assuming any unknown 600-cell result.

If the complete convention census is prohibitively expensive, the first
artifact may be labelled intermediate only; it may not claim structural
stability until every preregistered convention has run.

## Outcome hierarchy

- `CANONICAL_DATA_PROJECTION_CONTROL_FAILED`: provenance, frozen ranks,
  identities, bounds, or the zero-data control fails.
- `CANONICAL_DATA_PROJECTION_DISAGREEMENT_OPEN`: valid constructions, primes,
  representatives, or conventions disagree.
- `CANONICAL_DATA_PROJECTION_INTERMEDIATE_MODULAR_OPEN`: the baseline census
  is consistent but at least one preregistered convention is not evaluated.
- `CANONICAL_DATA_PROJECTION_STABLE_MODULAR_OPEN`: every control and
  convention agrees.

Even the stable outcome is **DERIVED (modular)** only.  The exact rational
carrier, its physical interpretation, the action Hessian, symplectic form,
propagation, tick, and speed all remain **OPEN**.

## Required artifact

Record all four ranks and all five derived dimensions for every calculation,
maximum elimination widths, input hashes, source hash, protocol commit,
outcome, and test count.  Register the verifier before its first execution.
Run only this targeted verifier, not the full suite.
