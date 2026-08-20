# Symbolic adversarial result: correct generic ideal, unresolved pointwise strata

Date: 2026-08-20

## Provenance

- adversarial protocol: `2b2ee13`;
- registered first implementation: `1fb14f3`;
- factor-classifier first failure: note `1c03caa`;
- classifier-only repair: `94d1bc0`;
- frozen gap artifact: `72cf80f`;
- artifact SHA-256:
  `96a0d91528c7cf69e9286a50bed1f81184ef58adc3a6dd4b8e76b62c6c5c223a`.

The targeted verifier passed `14/14`.  No full suite was run.

## What is exact and positive

The mechanically independent per-vertex construction gives the one-cell
ideal

```text
< A+B-8, C+D-1 >.
```

After the exact lateral face connection is retained and eliminated, the
generic fraction-field Gröbner basis is exactly

```text
B - 2 - 2 tau^2/(lambda-1)^2,
D - lambda/(lambda-1).
```

Thus it yields the disclosed four coefficients.  Time reversal and an odd
local relabelling preserve the basis.  Fixed-frame gluing without the
connection gives the unit ideal `<1>`, so the connection is load-bearing.
The corruption `D -> D+1` leaves sixteen nonzero residuals.  All three finite
global controls are reproduced exactly.

## Why the result is not yet accepted

The preregistered factor collector reported

```text
allowed:
    tau,
    lambda-1,
    (lambda-1)^2-3 tau^2;

additional:
    lambda+tau-1,
    lambda-tau-1,
    (lambda-1)^2+3 tau^2.
```

Consequently the frozen outcome is honestly

```text
FULL_SCALE_STRUT_SYMBOLIC_HYPOTHESIS_GAP.
```

The artifact also isolates what is not responsible:

```text
four vertex-solve determinants:  +/-256 tau,
affine-domain determinant:          16 tau,
lateral-normal norm:
  -((lambda-1)^2-3 tau^2)/(lambda-1)^2.
```

Therefore the additional factors arise among response/connection component
numerators or the chosen component-pivot annihilator, not in the four local
solves or the Lorentz transition's geometric nondegeneracy.

## Framing correction

The first collector treated a zero of any individual rational coefficient as
a rank exception.  That implication is false: a redundant residual component
may vanish while the residual ideal and all relevant ranks remain unchanged.
Likewise, choosing the first nonzero symbolic component of the connection as
an annihilator pivot creates a coordinate chart whose denominator can vanish
even when another connection component remains nonzero.

This diagnosis cannot reclassify the frozen result.  A stronger calculation
is required:

```text
residual belongs to span(connection)
iff c_i r_j-c_j r_i=0 for every i<j.
```

These wedge equations use no component pivot.  Rank-critical factors must be
computed only from determinants/minors, common connection zeros and actual
denominators—not arbitrary response numerators.

## Status

| Claim | Status |
|---|---|
| Fraction-field generic coefficient ideal | **DERIVED EXACT** |
| Three finite global controls | **DERIVED EXACT** |
| Precision of accepted curved matrix | **RESOLVED** |
| Uniform theorem on every frozen real nondegenerate point | **OPEN** |
| Additional factors are genuine geometric strata | **OPEN** |
| Accepted generic kinematic carrier | **OPEN pending gap resolution** |
| Dynamics and physical modes | **NOT EVALUATED** |

