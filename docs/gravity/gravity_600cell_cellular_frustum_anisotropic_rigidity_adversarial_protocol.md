# Adversarial protocol: polynomial fixed-bottom frustum rigidity

Date: 2026-08-19

The primary exact outcome
`CELLULAR_FRUSTUM_SIX_SHAPES_UNDERDETERMINED` is disclosed.  This protocol is
committed before evaluating the independent polynomial Jacobian.

## Frozen inputs

| input | SHA-256 |
|---|---|
| primary verifier | `2f766503296aa43f6192d2cce6ce44faac3b7fb57ba131ba0fbf393a2da80f60` |
| primary artifact | `c55f98313121018ff5ca1fc834260e8f2f075248a21fd7b99a356d89b2d18255` |
| primary protocol | `ac20410bced8408c9cc8ec609653c3036a029b8e1d439a84c3acc3d5960eb1e8` |
| prior-art gate | `92c88042e8233a542b9f21e96a99bc0d09cf13cff89a8e243354f97984baaaab` |

The primary source may be hashed and parsed for provenance but none of its
rigidity, point, isometry or diagonal helper functions may be imported.

## Independent geometry

Fix the bottom tetrahedron pointwise, eliminating global isometries before
the calculation.  Use the irregular equal-radius rational tetrahedron

```text
p0=(5,0,0)   p1=(0,5,0)
p2=(0,0,5)   p3=(3,4,0).
```

It is nondegenerate and every point has norm five.  For

```text
(lambda,tau)=(1,7),(2,7),(3,13)
```

evaluate the top vertices at `qi=(lambda*pi,tau)`.  Equal bottom norms imply
equal timelike struts at every representative.

Create 16 independent symbolic variables for the four top coordinates.  Do
not use a rigidity-row formula.  Instead form the ten squared-length
polynomials directly:

```text
six eta(qi-qj,qi-qj),
four eta(qi-pi,qi-pi),
```

differentiate them symbolically and substitute the rational representative.
All ranks, nullspaces and determinants are exact SymPy results.

## Frozen tests

For each representative:

1. the ten-by-sixteen polynomial Jacobian has exact rank ten and nullity six;
2. for every one of the 24 total colour orders, add the corresponding six
   cross-length polynomials and require the complete sixteen-by-sixteen
   Jacobian determinant to be nonzero;
3. require the six cross differentials to have rank six on the original
   six-dimensional nullspace;
4. record the complete multiset of absolute determinants without interpreting
   their magnitude.

As an unpredicted adversarial census, also enumerate all `2^6=64` independent
choices of one diagonal on each quadrilateral face and report the exact rank
histogram of their complete Jacobians.  No outcome depends on a particular
histogram outside the 24 staircase choices; the census tests whether local
metric completion is even less selective than global staircase conformity.

## Outcome hierarchy

1. `ADVERSARIAL_CELLULAR_FRUSTUM_CONTROL_FAILED` if provenance, equal radius,
   timelike strut, bottom nondegeneracy, polynomial construction or census
   completeness fails.
2. `ADVERSARIAL_CELLULAR_FRUSTUM_SIX_SHAPES_CORROBORATED` if all three base
   Jacobians have rank ten/nullity six and every staircase completion has
   determinant nonzero and rank-six action on the flex kernel.
3. `ADVERSARIAL_CELLULAR_FRUSTUM_DISAGREEMENT_OPEN` otherwise.

## Interpretation firewall

Corroboration confirms only local underdetermination of the current cellular
metric data.  It does not choose additional lengths or establish a physical
ensemble.  The 64-choice census is descriptive and may not be used to select
a completion after inspection.

Only this verifier and static registry guards may be run.

