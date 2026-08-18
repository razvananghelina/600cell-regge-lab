# Blind enumeration of the full-Hessian primitive spectrum

Date: 2026-08-10

Protocol commit: `3767638`.

Verifier: `reproducible/verify_hopf_full_hessian_spectral_enumeration.py`.
Targeted blind result: `8/8`.

## Provenance boundary

This is STEP 1 only.  No characteristic coefficient was evaluated at a
distinguished `Box_i`; no comparison with `Tr(X^3)`, `C_box`, a signed Hopf
orbit or a desired extremum was performed.  The complete expanded
coefficient data are in
`reproducible/hopf_full_hessian_spectral_enumeration.json`, whose field
`target_comparison_performed` is `false`.

## Exact arena

The full six-label Hessian has the constant vector as an exact universal zero
mode.  Restricting it to the canonical subspace `1^perp` gives a
five-dimensional operator `Hhat_X`.  Coordinates use

```text
X=sum_(a=0)^4 u_a (Box_a-Box_5).
```

The exact norm is

```text
q(X)=Tr(X^2)
    =17280*(sum_a u_a^2 + sum_(a<b) u_a u_b).
```

The restricted operator family is exactly `A5`-equivariant for all 60 group
elements.

## Complete primitive count

With convention

```text
det(lambda I-Hhat_X)
 =lambda^5-e1 lambda^4+e2 lambda^3-e3 lambda^2+e4 lambda-e5,
```

there are exactly

```text
N=5
```

primitive characteristic coefficients.  The five power sums
`Tr(Hhat_X^k)`, `1<=k<=5`, satisfy all Newton identities coefficientwise and
therefore are not additional attempts.

## Blind classification

Exact polynomial comparison gives

```text
e1 = 0,
e2 = -9331200 q.
```

Thus `e1` vanishes and `e2` is constant on the normalized sphere.  The
remaining coefficients have degrees three, four and five:

```text
e3, e4, e5.
```

None is proportional to a norm power times another lower-degree primitive.
Consequently the number of distinct nonconstant primitive spectral
characters on the normalized sphere is exactly

```text
3.
```

Their full integer polynomials, including every monomial coefficient, were
recorded before target comparison in the JSON file rather than abbreviated in
this note.

## STEP 1 ledger

- **DERIVED:** the complete primitive spectral count is `N=5`.
- **DERIVED:** `e1=0` and `e2=-9331200 q`.
- **DERIVED:** precisely three primitive coefficients remain nonconstant on
  the fixed sphere: `e3,e4,e5`.
- **DERIVED:** the characteristic and power-sum descriptions agree by all
  five Newton identities.
- **DERIVED:** the complete operator family is exactly `A5`-equivariant.
- **NOT YET TESTED IN THIS COMMIT:** any relationship to the Hopf target,
  stationarity, extrema, uniqueness, determinant physics or sign.
