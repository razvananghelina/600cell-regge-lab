# Protocol: coordinate-free resolution of the symbolic hypothesis gap

Date: 2026-08-20

This protocol is committed after the frozen `SYMBOLIC_HYPOTHESIS_GAP` result
and before any wedge residual or special-stratum rank is computed.

## Frozen inputs

| input | SHA-256 |
|---|---|
| gap-result note | `6e81e90427996557e2d9eab62dd1b01d098121d971d66c14f03899d8ed5a017d` |
| first symbolic protocol | `ab1b82702f4acb95998f35eda91b2472870a0b730e46f6c391e73320cc16fea1` |
| first symbolic source | `abf83f50fb5cda2ce6a4820b2d105ce3895da7f91504e886c76b559ad6964f2e` |
| frozen gap artifact | `96a0d91528c7cf69e9286a50bed1f81184ef58adc3a6dd4b8e76b62c6c5c223a` |
| resolved precision artifact | `2a2a79271a92fc2ddde343a9d0651402df6eeb4a90efa2697e26f54cafcdf60f` |

The gap artifact must retain `14/14`, exact generic ideal agreement, and the
three additional reported factors

```text
lambda+tau-1,
lambda-tau-1,
(lambda-1)^2+3 tau^2.
```

## Scope made explicit

The physical Lorentzian variables `lambda,tau` are real.  The uniform theorem
asked for here is on the real domain

```text
lambda != 1,
tau != 0,
(lambda-1)^2-3 tau^2 != 0.
```

The already-proved fraction-field identity remains algebraic over
`Q(lambda,tau)`.  This protocol does not claim a uniform theorem on every
complex specialization.  Complex exceptional strata, if any, are recorded
separately and remain **OPEN** unless directly classified.

## Coordinate-free connection elimination

Rebuild the same two regular tetrahedral frusta from coordinates, without
importing the first symbolic source.  Let `c` be the twelve-component
connection displacement on the shared upper face and let `r` be each
shared-displacement residual column.  Replace the component-pivot
annihilator by all wedge equations

```text
c_i r_j-c_j r_i = 0,   0 <= i < j < 12.
```

These equations are equivalent to `r in span(c)` wherever `c!=0`, use no
component denominator, and remain valid when an individual `c_i` vanishes.

After substituting `A=8-B`, `C=1-D`, the monic Gröbner basis in `(B,D)` must
again equal

```text
B - 2 - 2 tau^2/(lambda-1)^2,
D - lambda/(lambda-1).
```

The fixed-frame unit-ideal and `D+1` corruption controls remain mandatory.

## Correct factor classification

Three classes are recorded separately:

1. **actual denominators:** irreducible factors of denominators in the four
   vertex solves, affine transition, connection, wedge equations and reduced
   Gröbner basis;
2. **rank factors:** irreducible factors of the four solve determinants,
   affine-domain determinant, lateral-normal norm numerator and the common
   nonvanishing certificate for the full connection vector;
3. **coefficient numerators:** recorded diagnostically but never called rank
   exceptions merely because one redundant component vanishes.

For the connection, compute the exact Euclidean component norm

```text
N_c = sum_i c_i^2.
```

Factor its numerator.  On the real domain, any factor of the form
`(lambda-1)^2+3 tau^2` is a positive certificate, not an excluded stratum.
The verifier must show algebraically that it cannot vanish when `tau!=0`.

No actual denominator or rank factor other than

```text
tau,
lambda-1,
(lambda-1)^2-3 tau^2,
(lambda-1)^2+3 tau^2  [positive, not excluded]
```

is allowed.  Any other factor reaches a genuine gap outcome.

## Direct special-stratum attacks

The two linear factors are tested at frozen exact real points:

```text
lambda+tau-1=0:  (lambda,tau)=(-1,2),
lambda-tau-1=0:  (lambda,tau)=( 3,2).
```

Both satisfy every physical nondegeneracy hypothesis.  At each point,
rebuild directly over `Q`:

- all four vertex-solve ranks;
- affine Lorentz transition and connection rank;
- wedge residual ideal in `(B,D)`.

The ideal must equal the specialization of the disclosed formula.  These are
hostile controls specifically capable of showing that the two linear factors
are real exceptional strata.

## Outcome hierarchy

1. `FULL_SCALE_STRUT_GAP_CONTROL_FAILED`: provenance, coordinates,
   pivot-free connection, fixed-frame, corruption or special-point control
   fails.
2. `FULL_SCALE_STRUT_GAP_GENUINE`: an additional actual denominator/rank
   factor or either linear special-point failure remains.
3. `FULL_SCALE_STRUT_GAP_REAL_RESOLVED`: wedge ideal agrees, all real rank
   factors are covered by the frozen hypotheses or positivity, and both
   linear special points pass.
4. `FULL_SCALE_STRUT_GAP_DISAGREEMENT`: the pivot-free generic ideal differs.
5. `FULL_SCALE_STRUT_GAP_OPEN`: any remaining classified case.

Only outcome 3 can combine with the previous finite and precision artifacts
to accept the generic **real kinematic** carrier.  It supplies no action,
canonical stationarity, gauge classification, propagation, tick, `c`, `G`,
Planck scale, graviton or particle mass.

Only the new verifier and static registry guards may be run.

