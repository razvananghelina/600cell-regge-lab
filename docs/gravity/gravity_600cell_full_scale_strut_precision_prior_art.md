# Prior-art gate: precision audit of the full scale--strut carrier

Date: 2026-08-20

## Exact object and hypotheses

The frozen first artifact contains, for each staircase parity, one real
`1560 x 240` kinematic carrier at the accepted curved background.  Its exact
rank 240 is already established by the pole-identity and unsigned-incidence
block rows.  The open question is narrower:

> Does the `2.85%--2.93%` disagreement between direct binary64 SVD and the
> binary64 eigenvalues of `G^T G` reflect a singular-value defect, or only
> loss of accuracy from forming the normal equations?

Hypotheses are:

```text
the frozen edge order and background coefficients are unchanged;
all matrix entries are finite real numbers;
the direct condition number is approximately 1.10e7;
no rank decision is taken from a floating threshold;
arbitrary-precision arithmetic represents the binary64 matrix entries
exactly before forming its Gram matrix.
```

This mission does not change coordinates, rescale columns, or choose a more
favourable matrix.

## KNOWN

- Forming the normal matrix squares the 2-norm condition number:
  `kappa_2(G^T G)=kappa_2(G)^2` for full-column-rank `G`.
- LAPACK supplies direct SVD drivers and condition estimates rather than
  requiring singular values to be reconstructed from normal equations; see
  the [LAPACK Users' Guide](https://www.netlib.org/lapack/lug/) and the
  documented divide-and-conquer SVD least-squares driver
  [`xGELSD`](https://www.netlib.org/lapack/explore-html/d9/d67/group__gelsd_ga7da8d56f14942ae8cb9e0d681f6c4e20.html).
- At the observed `kappa~=1.10e7`, the elementary scale
  `epsilon_binary64*kappa^2` is about `2.7e-2`, already comparable with the
  observed discrepancy.  This is a prediction from the frozen artifact, not
  a new precision calculation.

## CONTROL

Use the exact-real two-by-two family

```text
A_delta = diag(1,delta) R,
R = 2^(-1/2) [[1,-1],[1,1]],
delta=1e-7.
```

Its singular values are exactly `(1,delta)` in real arithmetic.  In
binary64, direct SVD should remain accurate while the smallest eigenvalue of
`A_delta^T A_delta` is obtained by cancellation at relative scale
`epsilon/delta^2`.  An arbitrary-precision Gram calculation on the exact
binary64 entries must recover the direct result.  This known-answer control
calibrates the diagnostic before it is applied to the carrier.

For each carrier parity compare:

1. LAPACK `GESDD` direct SVD;
2. LAPACK `GESVD` direct SVD;
3. an 80-decimal symmetric eigensolve of the Gram matrix formed from exact
   binary64 entries;
4. the already-frozen inaccurate binary64 Gram route.

## OPEN

- whether the two direct drivers agree on every singular value;
- whether the arbitrary-precision Gram spectrum agrees with direct SVD;
- whether the old discrepancy has the predicted `epsilon*kappa^2` scale;
- whether any nonfinite or nonpositive high-precision eigenvalue appears.

## Scope firewall

This is a numerical-analysis audit only.  Even complete resolution does not
prove the generic geometric formula, does not pass the carrier through an
action, and supplies no gauge, constraint, graviton, clock, tick, `c`, `G`,
Planck or particle-mass information.

