# Blind direct Regge coefficients on projected 600-cell refinements

Date: 2026-08-17

Prior-art commit: `81db1ec`  
Protocol commit: `23157e2`  
Unit-sphere correction: `05d5685`  
Preserved first failure: `0a57607`  
Quadratic-lapse correction: `e0fcda4` / `8d9983d`  
Preserved second failure: `c67af57`  
Conditioned lapse-step correction: `917717c` / `877a345`

Targeted verifier:

```text
reproducible/verify_gravity_600cell_projected_refinement_acceleration_blind.py
```

Frozen Stage-A artifact:

```text
reproducible/gravity_600cell_projected_refinement_acceleration_blind.json
SHA-256 640bc0dd3d6f1ae727f8113bf29514878874effffd14f539f5a43e3c3b18d069
```

## Blind verdict

**DERIVED NUMERICAL:** the targeted verifier passes `9/9` and returns

```text
PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_DERIVED.
```

No refined coefficient was compared with the continuum target, and no
variant was ranked by such a distance before this artifact was committed.
The known exact level-zero coefficient was used only as the preregistered
implementation calibration.

## Coefficients printed before target comparison

The coefficient `a` is defined by

```text
log(R_1/R_0)=a*(tau/R_0)^2+O((tau/R_0)^4),
```

with unit volume radius, conserved total dust, and the contracting
time-orientation branch.

```text
variant                         level 0          level 1          level 2
common / legacy              -0.5394897453    -0.5095671441    -0.5023722301
first_tie_rank_0             -0.5394897453    -0.5095672089    -0.5023719473
first_tie_rank_1             -0.5394897453    -0.5095671992    -0.5023719533
first_tie_rank_2             -0.5394897453    -0.5095671723    -0.5023720102
```

The complete regulator spread is

```text
level 1: 6.47813e-8,
level 2: 2.82836e-7.
```

Those spreads are reported without interpreting their distance from any
continuum value.

## Geometry and static dust control

For the legacy tower, the unit-sphere chordal volume, scale required for unit
volume radius, spatial Regge curvature, and selected total mass are

```text
level   Vbar          s0             C              M=C/(8*pi)
0     16.69252677   1.057472835   60.41423150      2.403805894
1     18.85946501   1.015313414   59.54077906      2.369052326
2     19.51116056   1.003880949   59.29968342      2.359459435
```

All four variants have the frozen counts

```text
(V,T)=(120,600), (840,4800), (6480,38400),
```

closed face incidence two, Euler characteristic zero, positive tetrahedral
volumes, and Lorentzian frustum inertia `(3,1)` at every audited state.

The exact first-refinement octahedral diagonal tie is quantitatively visible:
the raw floating candidate spread is at most `3.192e-11`, and the four rules
produce four distinct selected-pair digests.  The coefficient agreement is
therefore not caused by accidentally evaluating one identical triangulation.

## Numerical and algebraic controls

Across every carrier, the worst recorded diagnostics are

```text
exact regular closed-action relative error       7.405e-14
static cellular-curvature relative error         7.708e-13
static total-action relative residual            7.708e-13
static lapse relative residual                   1.745e-8
tetrahedron-relabel action relative error         1.153e-13
seam coefficient extrapolation difference        4.337e-7
seam derivative-step difference                  3.001e-7
lapse derivative-step difference                 1.327e-6
independent lapse/seam coefficient difference    2.011e-6
held-out seam affine residual                     2.749e-7
held-out lapse quadratic residual                1.110e-6
static lapse-root residual                        1.874e-7
```

Every number lies inside its unchanged frozen gate.

## What the two preserved failures taught us

The first failure was substantive: the lapse leading equation is quadratic,
not affine.  It has an exact static root and a distinct dynamic root.  The
seam equation remains affine.  That correction is now part of the stated
mathematics rather than a hidden implementation detail.

The second failure was numerical: a real-log derivative step that was too
small amplified the curvature/dust cancellation on one level-two carrier.
The failed artifact and the step audit were committed before replacing only
those steps.  No carrier, eta value, coefficient sentinel, equation,
tolerance or target comparison changed.

## Interpretation boundary before disclosure

**DERIVED NUMERICAL:** the direct non-averaged irregular cellular action has
a stable, uniquely reconstructed homogeneous acceleration coefficient on all
registered finite carriers.

**DERIVED NUMERICAL:** the coefficient is insensitive at the displayed scale
to four disclosed resolutions of the exact first-level diagonal ambiguity.

**OPEN at this Stage-A commit:** whether the coefficients approach the
continuum Friedmann target.  That comparison is deliberately deferred to the
next commit.

**OPEN:** an infinite refinement theorem, exhaustive independence from all
`3^600` diagonal choices, local refined dust, anisotropic stability, tensor
modes, a limiting speed, a selected tick, and Planck or particle scales.
