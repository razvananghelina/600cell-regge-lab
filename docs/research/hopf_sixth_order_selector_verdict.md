# Sixth-order Hopf selector: exact verdict

Date: 2026-08-10

## Provenance

The candidate moment identities and critical values were already visible in
an exploratory symbolic calculation.  They were disclosed, together with an
exhaustive falsification protocol, in commit `66d89d3`.  This result is
therefore not presented as blind.  Its new content is the exact exhaustion of
the Lagrange ideal and of all three critical fibres.

The registered verifier is
`reproducible/verify_hopf_sixth_order_selector.py`; exact structured output is
stored in `reproducible/hopf_sixth_order_selector.json`.

## Canonical moment hierarchy

Let `P_i`, `i=1,...,6`, be the six equal-status rank-one projectors of the
icosahedral Hopf axes, and define

`S_(2m)(n) = sum_i (n^T P_i n)^m`.

The weights are not fitted: transitivity of the exact six-axis action forces
them to be equal.  Direct arithmetic in `Q(sqrt(5))` gives

```text
S2(n) = 2 |n|^2,
S4(n) = (6/5) |n|^4,
S6(n) is not radial.
```

Thus this particular canonical hierarchy cannot distinguish an axis at
quadratic or quartic order.  Degree six is its first anisotropic member.  This
statement concerns the equal-weight moment hierarchy; the verifier does not
claim a classification of every imaginable non-polynomial invariant.

## Exhaustive critical-point calculation

The twelve exact icosahedron vertices reconstruct the full incidence data
`(V,E,F)=(12,30,20)`.  Their antipodal vertex, edge-centre and face-centre
lines give the three symmetry-axis orbits

```text
orbit                         unoriented lines   signed points
C10 / fivefold / Hopf                 6              12
C4  / twofold                        15              30
C6  / threefold                      10              20
```

On the unit sphere, the exact Lagrange ideal

`grad(S6)-2 lambda n=0,  |n|^2=1`

is zero-dimensional.  Its monic elimination polynomial is

`(5 lambda-12)(15 lambda-34)(25 lambda-78)/1875`.

Euler homogeneity gives `lambda=3 S6`, so the only possible critical values
are

```text
C10/Hopf: S6 = 26/25, lambda = 78/25,
C4:       S6 =  4/5,  lambda = 12/5,
C6:       S6 = 34/45, lambda = 34/15.
```

For the three specialized ideals, the exact quotient-ring dimensions are
respectively `12`, `30` and `20`.  Each equals the number of explicitly
exhibited distinct real signed axis points.  Hence the listed orbits exhaust
the complex critical fibres, with no hidden generic or complex solutions and
no multiplicity deficit.

Compactness of `S^2` then fixes the global extrema:

- the six unoriented Hopf axes are exactly the global maxima, `26/25`;
- the ten unoriented threefold axes are exactly the global minima, `34/45`;
- the fifteen twofold axes form the intermediate critical orbit, `4/5`.

The intermediate orbit is intentionally not labelled a saddle here because
the constrained Hessian was not needed for, and was not included in, the
acceptance boundary.

## What is selected

Mathematically, the invariant potential

`V(n) = -g S6(n),  g>0`

has exactly six unoriented degenerate global minima, the six Hopf axes.  This
is a **DERIVED CONDITIONAL symmetry-breaking selector**: if such a term with
that sign occurs in the action, the continuous choice of a direction collapses
to six symmetry-related vacua.

The opposite sign selects the ten `C6` axes instead.  Therefore the sign is
load-bearing; it cannot be hidden in normalization or dismissed as a
convention.

## What is not selected

- **OPEN:** no existing theory operator has yet been shown to generate this
  sixth-order potential or the sign `-g`.
- **OPEN:** no magnitude for `g` is derived.
- **OPEN:** the order parameter `n` has not yet been identified with a
  physical field, mass condensate, clock or inertial structure.
- **OPEN:** the same six-axis construction occurs in both handed families
  `qH` and `Hq`.  The sixth-order potential selects an axis inside a chosen
  handed sector; it does not by itself choose the chirality.
- **DERIVED NEGATIVE:** the result still does not select a quadratic kinetic
  anisotropy `r`, a Lorentzian sign or a propagation speed.

## Relation to `a1=5`

The exact chain now contains

```text
Tr(P_i P_j)=1/5  (i != j),
S4=(6/5)|n|^4,
six equiangular axes = one chosen axis plus five alternatives.
```

These are **DERIVED icosahedral identities**.  They give the integer five a
clean geometric home: it is the equiangular-frame coherence denominator and
the number of alternative axes.  Equating that role with the physical symbol
`a1=5` remains **STRUCTURAL**, because there is still no derived map from this
order-parameter potential to the repository's coupling or speed formulae.

## Next gate

The next non-arbitrary question is no longer whether a symmetry-breaking
potential can be written: it can, canonically.  It is whether the theory's own
spectral/inner-fluctuation action produces an order parameter in one of the
two handed three-dimensional spaces and a sixth-order term with the required
negative sign.  If it does not, the selector remains mathematically available
but physically external.
