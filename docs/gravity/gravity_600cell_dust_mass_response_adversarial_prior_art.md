# Independence gate: invariant-subspace audit of the dust-mass response

Date: 2026-08-18

Primary result artifact commit: `99f855c`

Status: written before evaluating any invariant-subspace leakage.  This is an
adversarial audit of a project result, not a new novelty claim.

## Exact object

The primary calculation reports that the zero-total-mass outgoing response
space `R_d` is separated from both the expanding and contracting invariant
spaces of each frozen tangent block `T_d` in all

```text
2 schedules x 7 sectors x 2 branches = 28
```

comparisons.  Its decisive numerical step used ordered complex Schur spaces
and principal-angle distances.

The independent necessary condition is elementary:

```text
if im(R_d) equals either invariant branch, then
T_d im(R_d) is a subset of im(R_d).
```

Therefore form

```text
L_d = (I-P_R) T_d R_d,
P_R = R_d (R_d* R_d)^-1 R_d*.
```

Any resolved nonzero `L_d` refutes equality with **both** tangent branches
without constructing either branch, ordering eigenvalues or computing a
principal angle.  This is a mechanically different falsification of the
decisive subspace claim.

## KNOWN / CONTROL / OPEN

- **KNOWN:** equality with an invariant subspace implies invariance.  This is
  finite-dimensional linear algebra, not a physics result.
- **CONTROL:** the frozen primary response and tangent matrices, treated as
  exact dyadic inputs; an identity tangent must preserve every response
  space; a synthetic off-block tangent must leak.
- **OPEN:** the leakage rank and margin on the actual matrices; stability
  under precision, column rephasing, phase-block swap and time reversal.
- **OPEN:** source-response construction remains primarily certified by the
  Flint action solve.  This audit changes the decisive comparison, not the
  upstream action differentiation.

## Framing attack

This test can corroborate separation but cannot identify which physical
sector receives the dust response.  Non-invariance also does not imply a
graviton.  It only excludes equality with a complete invariant branch.

No literature search can establish novelty for this project-specific matrix
audit.  The relevant physical prior art remains the inhomogeneous Regge
lattice cosmology and dust references listed in the primary prior-art gate;
the invariant-subspace implication itself needs no empirical assumption.
