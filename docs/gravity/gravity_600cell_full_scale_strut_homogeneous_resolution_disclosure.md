# Exploratory disclosure: the homogeneous near-null ratio

Date: 2026-08-20  
Status: **PATTERN disclosed before the exact protocol**

## Observed numerical structure

The P100-frozen candidate has five scale coefficients and five strut
coefficients. Within each group their spread is about `1.02e-35`, consistent with
the `h^2=1e-40` finite-difference error after conditioning. The mean scaled ratio
is approximately

```text
x_scale/x_strut = -1.4973543020498256.
```

It is not the simple rational `-3/2`.

After undoing the frozen column scalings, its physical carrier ratio is

```text
sigma/c = -2.3370425140082096e-6.
```

Here `sigma=delta lambda` is the additive upper-vertex scale coordinate and
`c=delta log rho` is the logarithmic pole coordinate.

## Independently predicted ratio

Let `s=log lambda`, `z=log rho`, and let `p_-(s,z)` be the old homogeneous
canonical momentum obtained from the already frozen closed cellular action. At
fixed old momentum,

```text
p_s delta s + p_z delta z = 0.
```

Because `sigma=lambda*delta s` and `c=delta z`, the predicted carrier ratio is

```text
sigma/c = -lambda*p_z/p_s.
```

Direct automatic differentiation of the closed action at the frozen background
gave

```text
p_s = 419653.5961615574293857...,
p_z = 0.9807513514608717338468...,
-lambda*p_z/p_s
    = -2.33704251400820960763253452902063217449...e-6.
```

The endpoint Jacobian stored by the independent homothetic canonical-lapse solve
gives the same value within `5.34e-42`. The mean P100 candidate agrees at the
same finite-difference scale.

## Status and prohibition

This observation makes the next mission target-disclosed. It is not evidence for
an exact null by itself. The exact protocol must derive the line from the action
generating-function identities, close its nullity with certified minors, include
the `sigma=lambda*delta s` convention as a corruption control, and keep the
omitted pole equation outside the intersection verdict.
