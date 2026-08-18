# Protocol: exact rank of the transported generalized phase intersection

Date: 2026-08-18

Prior-art/framing gate commit: `656b1d7`.

## 1. Frozen inputs and replay

Require exact SHA-256 hashes

```text
phase verifier
  9c4c36b463a8faaa8d40b7db1b6b1852e3c04155c1b6ada4d02fbda747f6fcf3
phase artifact
  45eb9a3e80ead758d9b3c2f8e1eccff44b06e2759251ab00c447aa53e6705743
adversarial corroboration artifact
  c33615ac6d0f3133e53077f46c5ee766b9c633d4d64c32124c24839c9c84c880
```

Replay the phase verifier completely in-process. Require `7/7`, outcome
`GENERALIZED_PHASE_TRANSPORT_REFUTED`, an unchanged phase-artifact hash, 32
high-precision configuration projectors and 16 exact tangent balls. Require
the adversarial artifact to report `8/8` and
`ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED`; do not use its numeric
principal-angle values.

Use `mpmath` at 100 decimal digits and the already reconstructed Flint source
balls at 80 decimal digits. No binary tangent archive may supply the decisive
matrix.

## 2. Complete blind rank census

For each parity, sector `4,5` and all four derivative schedules, form

```text
Q0 = diag(P_old,conjugate(P_old)),
Q1 = diag(P_shifted,conjugate(P_shifted)),
R  = (I-Q1) T_2 Q0.
```

Compute all 60 midpoint singular values at 100 decimal digits. Use the complete
operator error `e` already recorded for the matching full phase residual; it
includes the tangent ball, both projector errors and the frozen arithmetic
floor. Check that the recomputed `||R||` agrees with the committed full
residual norm within `10 e`.

For every singular value `s`, use only

```text
s <= 10 e   SINGULAR_ZERO_CONSISTENT
s > 100 e   SINGULAR_NONZERO_RESOLVED
otherwise   SINGULAR_OPEN
```

Record all labels, values and error units. No singular value or schedule may be
dropped.

## 3. Structural rank control

The exact operator has right factor `Q0` of rank 30, hence exact rank at most
30. Require:

- no cell has more than 30 nonzero-resolved singular values;
- the lower 30 midpoint singular values are all zero-consistent;
- all singular values are finite and sorted nonincreasingly;
- exactly 16 cells and 960 singular-value records exist.

For a cell with 30 nonzero-resolved singular values, Weyl perturbation plus the
structural upper bound certifies

```text
rank(R)=30,
dim(F0 intersection T_2^-1(F1))=0.
```

Fewer than 30 resolved values do not certify a positive intersection and must
remain open.

## 4. Frozen outcome hierarchy

Use the first applicable branch:

1. provenance, replay, finiteness, reconstruction, norm-overlap, ordering,
   census or structural-rank control fails:
   `TRANSPORTED_INTERSECTION_CONTROL_FAILED`;
2. all 16 cells have exactly 30 nonzero-resolved singular values:
   `TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL`;
3. otherwise:
   `TRANSPORTED_INTERSECTION_DIMENSION_OPEN`.

Branch 2 closes the present generalized-fiber phase route but not other Regge
perturbations. Branch 3 authorizes no graph selection; a stronger exact rank
method would need a separate preregistration.

## 5. Deliverable and exclusions

Write a deterministic JSON artifact with all hashes, 16 cell summaries, all
960 singular records, counts, outcome and status ledger. Register the verifier
before its first scientific execution and run twice with a byte-identical
artifact.

Run no full suite and no unrelated verifier. Perform no comparison with a
desired rank and compute no graph, symplectic restriction, propagator, root,
dispersion, mass, inertia or limiting speed in this mission.

## 6. Execution-history addendum: serialized-norm control failure

The first execution of registered verifier commit `6e785d6` terminated `6/7`
with `TRANSPORTED_INTERSECTION_CONTROL_FAILED`.  The SVD norm recomputation and
the committed phase norm differed by `2.56e-30...3.46e-29`, whereas the phase
operator error is approximately `2e-47`.

Cause: the committed phase norm was serialized by `mp.nstr(..., digits=30)`.
It therefore cannot support a replay comparison below its approximately
30-significant-digit text resolution.  The SVD and committed strings agree in
all 30 stored significant digits.  The 480 nonzero-resolved and 480
zero-consistent labels printed by this control-failed execution are diagnostic
only and must not be cited as the scientific result.

Before an admissible rerun, replace the impossible norm-overlap control

```text
|norm_recomputed-norm_committed| <= 10 e
```

by

```text
|norm_recomputed-norm_committed|
    <= 10 e + 1e-29 max(1,|norm_committed|).
```

The additional term bounds only the known 30-digit JSON serialization.  It is
not added to the singular-value error and changes no singular label, resolved
rank, structural upper bound or outcome criterion.  The accepted result still
requires two complete passing runs with byte-identical artifacts.
