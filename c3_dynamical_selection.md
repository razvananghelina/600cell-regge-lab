# Dynamical selection in the residual-C3 candidate moduli

Date: 2026-07-24 (evening session)

## Decision

The already-derived physics does **not** provide a new finite selector.

The least-assumption spectral compatibility statement is additivity:

`Tr(D_tot^2)=14880+q_m`,

`(1/2)Tr(D_tot^4)=55920+h_m`,

where `q_m=Tr(D_m^2)` and
`h_m=(1/2)Tr(D_m^4)`.  This is a **DERIVED-CONSTRAINT** for a direct-sum
matter component, but it is bookkeeping rather than a cut.  It leaves the
generic gauge-and-scale moduli at 122 real dimensions.

Requiring instead that the old numerical totals remain exactly
`14880,55920` forces `q_m=0`, hence `D_m=0` by self-adjointness.  The
nonzero projective moduli is then empty.  This is an exact result for the
strict variant, but the strict premise is a **CHOICE**, not an
already-derived constraint: the three moment identities were proved for the
2640-dimensional Kähler--Dirac complex, not for a later enlarged operator.
It therefore cannot honestly close the C3 door.

The founding trace bootstrap has no licensed analogue on `D_m`, and the
derived Galois law relates paired mass-sheet data rather than arbitrary
Dirac coefficients.  Neither may be imposed.  On the known quartic critical
circle, every permitted moment invariant and the canonical metric are
constant.  The two registered CP targets are not distinguished.

Throughout this note, “candidate” means the post-audit
residual-equivariant Krajewski-legal block candidate.  No faithful real
algebra representation, opposite action, matrix-level first-order
double-commutator, or finite spectral-triple Dirac is claimed.

All finite statements below are checked by
`reproducible/verify_c3_dynamical_selection.py`.

## 0. Exact arena

Use the recorded witness

`H+ = A + 2B`, `H- = C + 2B`,

with `A=(2,2)`, `B=(3bar,2)`, `C=(2,1bar)` and C3 weight vectors

`A=(2,1,1)`, `B=(2,2,2)`, `C=(0,1,1)`.

For the upper odd block `T:H+ -> H-`, the complete complex census is

| block | unrestricted | first-order legal |
|---|---:|---:|
| `A -> C` | 2 | 2 |
| `A -> 2B` | 16 | 16 |
| `2B -> C` | 8 | 0 |
| `2B -> 2B` | 48 | 48 |
| total | 74 | 66 |

Self-adjointness sets the reverse block to `T*`.  Therefore the unrestricted
self-adjoint odd arena has real dimension `2(74)=148`, while the legal
linear space has

`d0=2(66)=132`.

These are **STRUCTURAL-CONSTRAINTS** within the explicit candidate model,
with the audit scope above.  The established conditional KO6 sheet relation
determines the second sheet and leaves the same 132 independent real
parameters.

The grading-preserving bimodule commutant has unitary dimension ten.  Its
generic stabilizer is the common scalar `U(1)`, so a generic orbit has
dimension nine.  Removing positive scale gives

`132-9-1=122`.

This is the dimension of the generic principal stratum, not a claim that the
whole quotient is a smooth manifold.

## 1. Consistency cuts

### 1(a). Spectral-moment compatibility

Let `D_g` be the Kähler--Dirac operator on its 2640-dimensional cochain
space.  Its exact finite moments are

`Tr(I_g)=2640`,

`Tr(D_g^2)=14880`,

`(1/2)Tr(D_g^4)=55920`.

For the least-assumption interpretation in the mission, the matter block is
an independent direct-sum component:

`D_tot=D_g direct-sum D_m`.

Thus, exactly,

`Tr(I_tot)=2640+30=2670`,

`Tr(D_tot^2)=14880+Tr(D_m^2)`,

`(1/2)Tr(D_tot^4)=55920+(1/2)Tr(D_m^4)`.

Calling the last two matter terms `q_m,h_m` explicitly parameterizes the
change.  This is the defensible **DERIVED-CONSTRAINT**.  Since `q_m,h_m`
are defined by `D_m` rather than frozen independently, the variety after
this condition is the entire legal space `R^132`; after the generic
gauge-and-scale quotient it remains 122-dimensional.

Three distinct alternatives must not be conflated:

1. **Strict unchanged totals — CHOICE.**  Requiring
   `Tr(D_tot^2)=14880` gives `Tr(D_m^2)=0`.  A self-adjoint `D_m` has real
   eigenvalues, so the trace is a sum of squares and `D_m=0`.  The affine
   survivor is the single point `{0}`; the quartic equation is redundant.
   After the required nonzero scale quotient the survivor is empty.
2. **Chosen nonzero quadratic normalization — CHOICE.**  Fixing
   `q_m=q0>0` gives a positive-definite quadratic sphere `S^131` in the
   132-real-dimensional legal space.  Quotienting its generic nine-dimensional
   gauge orbit gives dimension 122.  It chooses a scale representative but
   does not reduce the projective moduli.
3. **Chosen quadratic and quartic values — CHOICE.**  At a regular,
   algebraically independent nonempty level, the real complete-intersection
   variety
   `{Tr(D_m^2)=q0, (1/2)Tr(D_m^4)=h0}`
   has dimension 130 before gauge and 121 on the generic gauge quotient.
   Exceptional extremal values can be singular.  No derived values `q0,h0`
   exist, so this one-dimensional cut may not be used as a result.

The reduced identity `2(62)^2+1=3(11)(233)` likewise belongs to the
Kähler--Dirac moments after division by 240.  Extending it to the
2670-dimensional direct sum would require a new normalization and is not an
already-derived condition.

### 1(b). Bootstrap self-consistency

No honest analogue exists, so no cut is imposed.

The verified identity is specifically

`Tr((c A_f-A)^3)=3 N_f(c-2)`

on the 120-vertex 600-cell, where the Hopf-fiber adjacency, face count
`N_f=1200`, and vanishing triangle types make the trace polynomial linear.
Setting it equal to `N^2=120^2` uniquely fixes `c=6`.

The 30-dimensional C3 candidate has no derived fiber adjacency, face count,
Coxeter normalization, or one-parameter operator `c A_f-A`.  Replacing
`Box` by `D_m`, choosing a power three because `C3` has order three, or
choosing a right-hand side such as `30^2` would be a new invention, not the
derived bootstrap.  Status: **OPEN**, with no variety cut.

### 1(c). Restriction to the exact critical circle

Let `a,b>0` be the coefficients of the symmetry-breaking action
`b Tr(D^4)-a Tr(D^2)`, and put

`r=sqrt(a/(2b))`.

The exact two-channel family can be gauge-fixed to

`T_AC(theta)=diag(r, r exp(i theta))`, `theta in R/(2 pi Z)`,

with all other legal blocks zero.  Its self-adjoint odd operator is

`D_m(theta) = [[0,T_AC(theta)*],[T_AC(theta),0]]`

on the active four-dimensional subspace and zero on the other 26
single-sheet dimensions.  Its spectrum is

`{+r,+r,-r,-r,0 x 26}`.

Consequently

`Tr(D_m^2)=4r^2`,

`(1/2)Tr(D_m^4)=2r^4`,

independently of `theta`.

The defensible additive compatibility condition therefore leaves the entire
`S1`.  The strict unchanged-total choice has empty intersection for `r>0`.
Any chosen fixed moments either contain the whole circle when they equal the
displayed values or miss it entirely.  No permitted moment cut makes the
intersection finite and nonempty.

### 1(d). Galois covariance

No new equation is licensed.

The derived law `z_b=phi sigma(z_t)` relates two already paired mass
variables.  Its honest analogue would require two separately identified
Dirac blocks and a derived `Q(sqrt(5))` coefficient lattice.  The KO6
sheet-swap covariance already determines the second sheet from the first,
but supplies no arithmetic action on arbitrary complex entries.

In particular, `sigma(D_m)=phi D_m` is not an analogue of the paired law:
applying `sigma` twice gives `D_m=-D_m`, hence zero.  It is an invented fatal
fixed-point equation and is not imposed.  Status: **OPEN**, no cut.

## 2. Physics of the outcomes

Under the defensible derived conditions there is no finite surviving set:
the residual generic moduli dimension is still 122.  Therefore no
gauge-invariant Yukawa/mass list, commutant `Y`, or anomaly-forced charge
tuple is selected.

For completeness, the strict-choice survivor before projectivization is
`D_m=0`.  Its eigenvalue list is thirty zeros, its commutant is the full
unselected

`C_A + M2(C)_B+ + C_C + M2(C)_B-`,

and anomaly forcing has no licensed preferred generator.  Zero-over-zero
mass ratios do not match the frozen exponent list

`(0,5,3,11,16,19,26,17,11)`.

Counting the thirty identical zero eigenvalues as an exponent-zero discovery
would violate the anti-numerology protocol: there are no defined distinct
positive ratios.  Thus the strict result is an excluded zero candidate, not
the discovery of the arc.

## 3. Registered CP phases on the circle

Write `alpha=arctan(sqrt(5))`.  The two registered targets are exact points

`exp(i alpha)=(1+i sqrt(5))/sqrt(6)`,

`exp(3i alpha)=(-7-i sqrt(5))/(3 sqrt(6))`.

Neither is distinguished:

- `Tr(D_m^4)` is constant on the circle.
- The metric induced by the canonical Hilbert--Schmidt inner product is the
  flat constant metric `2r^2 dtheta^2` (up to the fixed convention for the
  real trace pairing), so it has no point extrema.
- The additive moment constraints are constant on the circle.
- A “Galois-pairing overlap” is undefined without the missing arithmetic
  coefficient lattice and paired-block map; inventing one would violate
  Rule Zero.

For look-elsewhere accounting, two registered targets were tested against
the three defined circle diagnostics (quartic trace, canonical metric, and
allowed moment cuts): six comparisons, all null because each diagnostic is
constant.  The proposed Galois overlap was skipped, not counted.  Merely
placing the two target phases on a continuous `S1` is a **PATTERN** with
probability-one availability, not a hit.

## Status ledger

### Strengthened

- **DERIVED-CONSTRAINT:** direct-sum matter moments are additive, and the
  exact contribution is `q_m,h_m`.
- **DERIVED negative:** this least-assumption compatibility cuts none of the
  122 generic projective dimensions.
- **DERIVED conditional on the recorded candidate model:** the critical
  circle has exact spectrum `{+r x2,-r x2,0 x26}` and constant quadratic and
  quartic moments.
- **Exact CHOICE result:** preserving the old total quadratic moment
  unchanged forces `D_m=0` and empties the nonzero projective moduli.

### Downgraded or delimited

- The strict unchanged-total requirement is not derived and cannot close the
  C3 door.
- The founding `Tr(Box^3)=N^2` bootstrap is operator-specific and does not
  transfer to `D_m`.
- Mass-sheet Galois covariance does not define arithmetic covariance of
  arbitrary Dirac coefficients.
- The residual-equivariant Krajewski block remains a candidate, not a
  constructed spectral-triple Dirac.
- The two registered CP values are undistinguished points of a continuous
  critical circle; any association is **PATTERN**.

### Open

- a constructed real finite spectral triple realizing the candidate;
- a derived coupling between the geometric and matter operators rather than
  a direct sum;
- derived nonzero matter moment values or a derived spectral functional;
- a coefficient lattice and paired arithmetic Galois law for matter blocks;
- finite selection of `D_m`, hence masses, `Y`, hypercharge, and anomaly
  forcing.
