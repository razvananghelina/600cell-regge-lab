# Preregistered adversarial protocol: dust-response invariance

Date: 2026-08-18

Independence gate commit: `c0b7091`

Status: frozen before evaluating any actual leakage matrix.

## Inputs and exact-input convention

Require the exact primary JSON/NPZ artifact at commit `99f855c`, its passing
outcome `DUST_MASS_RESPONSE_BOTH_BRANCHES_SEPARATED`, and the frozen tangent
JSON/NPZ.  For each schedule and each of the seven sectors load:

- the operational-primary zero-sum phase response `R`;
- the operational-primary tangent block `T`.

Every binary64 real and imaginary component is converted through
`float.as_integer_ratio()` and represented as the exact dyadic rational in
Flint.  Run the complete audit independently at 100 and 140 decimal digits.
The primary Flint action solve is not rerun.

## Decisive invariant-subspace test

For `R : C^k -> C^n`, form in ball arithmetic

```text
C = (R*R)^-1 R* T R,
L = T R - R C.
```

This avoids constructing an orthogonal complement.  `L=0` if and only if
`im R` is `T`-invariant because the primary rank certificate guarantees full
column rank.  Report

```text
ell = ||L||_F / ||T R||_F.
```

Classification at a given precision:

- `INVARIANCE_REFUTED` if the lower endpoint of the ball for `ell` is
  strictly positive and its midpoint is more than 100 radii from zero;
- `INVARIANT_CONSISTENT` if its upper endpoint is below `10^(-dps+25)`;
- `OPEN` otherwise.

Require the 100- and 140-digit balls to overlap.  A nonzero result in every
one of the 14 schedule/sector cells corroborates both branch separations,
because expanding and contracting spectral branches are invariant spaces.

## Controls

For every distinct `(n,k)` size encountered:

1. use the same `R` with `T=I_n`; it must be `INVARIANT_CONSISTENT`;
2. use `R=[e_1,...,e_k]` and a synthetic tangent with
   `T e_1=e_1+e_(k+1)`; it must be `INVARIANCE_REFUTED`;
3. require the exact Gram determinant ball of every actual `R` to exclude
   zero.

For every actual cell repeat at 100 digits after:

- multiplying source columns by the exact cycle `(1,-1,i,-i,...)`;
- swapping the outgoing `q` and `p` row blocks while applying the same
  similarity to `T`;
- replacing forward evolution by `T^-1`, computed by a Flint solve.

All must retain `INVARIANCE_REFUTED`.  Rephasing and block swap must reproduce
the original `ell` ball; reversal need only reproduce the nonzero label.

As a lower-cost variant control, compute the same leakage ratio with an
independent pivoted-QR complement in SciPy for all four derivative variants.
All 56 ratios must be finite and nonzero; they do not override the exact-input
classification.

## Synthetic control boundary

The actual and separated-control results must have a lower endpoint above
zero at both precisions.  The identity controls must have upper endpoints
below the stated precision-dependent zero gate.  Failure of either control
invalidates the audit rather than changing a threshold.

## Outcome hierarchy

- `ADVERSARIAL_DUST_RESPONSE_SEPARATION_CORROBORATED` if all 14 actual cells
  and all convention variants refute invariance at both precisions, all
  identity/separated controls discriminate, all Gram determinants exclude
  zero, and the precision balls overlap;
- `ADVERSARIAL_DUST_RESPONSE_DISAGREEMENT_OPEN` if controls pass but any
  actual or convention case is open/invariant;
- `ADVERSARIAL_DUST_RESPONSE_CONTROL_FAILED` otherwise.

This audit corroborates only the separation conclusion conditional on the
frozen response matrices.  It does not independently rederive the dust source
from the Regge action and cannot promote the result to a dust clock, scalar
mode or graviton statement.
