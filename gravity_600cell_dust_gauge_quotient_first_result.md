# Gauge-quotient boundary response: first frozen result

Date: 2026-08-13

Prior-art gate: `ff8f404`

Frozen protocol: `25d9ee9`

Implementation commit: `7f0ea2b`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_gauge_quotient.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_gauge_quotient.json`

Targeted run: **13/13 implementation checks passed**.  The full suite was not
run.

## 1. Frozen overall verdict

The preregistered overall verdict is

```text
OPEN NUMERICALLY
```

because the odd schedule parity missed one quotient-spectrum agreement gate.
This outcome is retained before any correction.

The parity labels were:

```text
even: REGULAR_QUOTIENT_29_RESPONSES_PLUS_ONE_BOUNDARY_CONSTRAINT
odd : QUOTIENT_OR_COMPATIBILITY_UNRESOLVED.
```

## 2. Boundary compatibility is resolved in both parities

The independent 90-decimal mixed-action rectangles give exactly the same
reported compatibility row and norm in both parities:

```text
norm(c) = 1.1122241201432517e-3,
epsilon_c = 5.504269839429144e-6.
```

Thus `norm(c) > 100 epsilon_c`, and both parities receive

```text
ONE_LINEAR_BOUNDARY_CONSTRAINT.
```

This is not action noise or an analytic-gradient artefact:

- maximum action/row imaginary contamination is `2.94e-85`;
- all 480 displaced geometries remain Lorentzian and off branch boundaries;
- the old analytic mixed block agrees with the new action-only row to
  `7.67e-15` and `7.68e-15` in the frozen normalization.

Every component lies between

```text
2.030548784e-4 and 2.030676805e-4.
```

## 3. Post-result shape diagnostic

The following was not part of the frozen scientific label and is therefore
**PATTERN** until separately resolved.

The normalized compatibility row has cosine

```text
0.999999999770358
```

with the normalized all-ones boundary vector.  Its component orthogonal to
that vector has relative norm `2.143e-5`; the component spread is `1.280e-8`
around mean `2.030634132e-4`.

This suggests that the single constrained direction is the homogeneous
final-boundary scale, leaving the 29 zero-sum shape perturbations.  The
deviation is smaller than the current mixed-difference error envelope, so the
calculation does not yet decide exact uniformity.  It must not be silently
rounded to an all-ones vector.

## 4. Why the odd quotient remained unresolved

Both projected quotient matrices have 34 singular values above the frozen
absolute `1e-9` threshold.  The even weak quartet agrees with the previously
certified 80-decimal Schur quartet to `0.1395%`, passing the frozen `0.3%`
gate.

The odd weak quartet contains

```text
4.564712224e-8,
4.598541041e-8,
4.598625790e-8,
4.598850024e-8,
```

whereas the independent 80-decimal action reconstruction gives

```text
4.604967055e-8 ... 4.604967056e-8.
```

The resulting `0.8742%` mismatch fails the frozen gate.  This is the same
binary64 weak-cluster resolution problem already seen in the first Jacobian;
it is not evidence for an odd-parity null direction.  Nevertheless, the
protocol required agreement and the odd quotient is formally unresolved
until the high-precision Schur block is inserted by a preregistered block
construction.

## 5. Even-parity response

For the resolved even parity, the deterministic kernel basis of `c` has 29
columns.  The gauge-fixed response has:

```text
rank                    = 29 at 1e-7, 1e-9 and 1e-11,
relative solve residual = 3.746e-11,
operator norm           = 4.750e5,
condition               = 3.612e6.
```

The enormous amplification is physical/numerical information, not a failure
of the linear solve.  It is driven by the four very soft relative
pseudo-constraint modes.  It warns that even tiny allowed boundary
perturbations may leave the linear neighborhood rapidly.

## 6. Status ledger

| Claim | Status |
|---|---|
| One collective internal lapse direction is null | **DERIVED COMPUTATIONAL, upstream exact path** |
| One nonzero linear boundary compatibility condition exists | **DERIVED COMPUTATIONAL, both parities** |
| The condition is exactly homogeneous scale | **PATTERN / unresolved at current mixed step** |
| Even quotient has rank 34 | **DERIVED COMPUTATIONAL** |
| Odd quotient is physically rank deficient | **NOT SUPPORTED** |
| Odd quotient passes the frozen precision gate | **NEGATIVE** |
| All 30 boundary directions propagate | **REFUTED LINEARLY** |
| 29 compatible directions propagate in both parities | **DERIVED only for even; OPEN for odd pending correction** |
| Nonlinear displaced solutions exist | **OPEN** |

## 7. Next correction

Do not loosen the `0.3%` gate.  Reconstruct the quotient through the already
certified block factorization:

```text
A          : recorded regular 30 x 30 block,
S_relative : 80-decimal four-dimensional Schur form,
gauge      : exact analytic collective tangent.
```

This avoids asking the binary64 full Hessian to resolve a `4.6e-8` cluster
under a `2.4e3` leading scale.  In the same correction, add the frozen smaller
mixed step `1.25e-4` and sixth-order extrapolation to decide whether `c` is
exactly uniform or carries a genuine schedule-dependent component.
