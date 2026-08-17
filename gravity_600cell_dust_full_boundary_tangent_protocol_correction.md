# Protocol correction: binary64 reciprocal-singular floor

Date: 2026-08-17

Original blind protocol: `bc114bf`.

Preserved failed run: `6fbdc08` (`17/19`, outcome
`FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED`).

## What failed

All fourteen operational minimal blocks satisfied the primary calibrated
complex-sector symplectic identity.  Their defect norms were
`6.3e-30`--`1.9e-28`, only about `1e-12` of the permitted calibrated bound.
All 56 Flint balls for the pre-Legendre determinants excluded zero, and all
determinant-modulus controls passed.

The only failed subcontrol was the reciprocal product of binary64 singular
values.  The tangent condition numbers lie between `7.05e10` and `2.66e11`.
Observed reciprocal-product residuals were `2.56e-7`--`2.23e-6`, while the
original code assigned a binary64 floor of only `2.22e-15` because the four
high-precision derivative variants agree much more closely than the
subsequent binary64 SVD can resolve.

## Error in the preregistered proxy

For a backward-stable SVD, the absolute error in each singular value is of
order

```text
delta_sigma <= eps_machine * sigma_max.
```

The error in a reciprocal product is therefore bounded at first order by

```text
delta(sigma_i sigma_j)
 <= delta_sigma * (sigma_i + sigma_j) + delta_sigma^2.
```

For the extreme reciprocal pair this has worst scale

```text
eps_machine * sigma_max^2,
```

not `eps_machine`.  Since a symplectic matrix has
`condition_2(T)=sigma_max^2`, the missing factor is precisely the observed
large condition number.  This is a numerical-analysis error in the original
proxy, not evidence against symplecticity.

## Frozen correction

Do not change any geometry, derivative step, projected Hessian, Flint solve,
tangent entry, singular value or eigenvalue.  Change only the binary64 term
in `epsilon_reciprocal` from

```text
10 eps_machine * max(1,max_abs_reciprocal_product)
```

to

```text
delta_sigma = 10 eps_machine * max(1,sigma_max)
svd_product_floor = delta_sigma * (2 sigma_max) + delta_sigma^2.
```

The factor ten is the same preregistered numerical safety factor already
used elsewhere.  Also compute the operational singular spectrum with both
LAPACK `gesvd` and `gesdd`; add their maximum product-vector difference to
the proxy and report it.  This is an independent binary64 implementation
diagnostic.

The original failure remains committed.  The corrected verifier must be
rerun from scratch.  If the reciprocal products still exceed ten times the
corrected proxy, the canonicality failure stands.

No unit-circle label, schedule comparison, physical interpretation or other
threshold is changed.
