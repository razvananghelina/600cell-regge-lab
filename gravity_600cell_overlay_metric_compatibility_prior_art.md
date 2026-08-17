# Prior-art gate: metric compatibility of the universal staircase overlay

Date: 2026-08-17

Status: written before evaluating any dynamic metric comparison on the 148
overlay chambers.

## 1. Exact object and complete hypotheses

The parameter carrier is the tetrahedral prism

```text
P = Delta^3 x I,
lambda_i >= 0,  sum_i lambda_i = 1,  0 <= t <= 1.
```

For every order `o=(v_0,v_1,v_2,v_3)`, use the already derived staircase
triangulation with four simplices

```text
S_(o,k) = conv(bottom v_0,...,bottom v_k,
               top v_k,...,top v_3),       k=0,...,3.
```

Let `u_i` be the four unit 600-cell vertex directions of one spatial
tetrahedron.  They satisfy

```text
u_i . u_i = 1,
u_i . u_j = phi/2  for i != j.
```

The inherited homothetic boundary geometry is

```text
X(bottom i) = (R_- u_i, 0),
X(top i)    = (R_+ u_i, T),
R_+-R_- arbitrary,
T^2 = rho+(R_+-R_-)^2,
R_minus = phi L_minus,
R_plus  = phi L_plus,  rho > 0.
```

Thus each same-vertex strut has Lorentzian square `-rho`.  On each staircase,
`X_o` is the unique continuous map that is affine on every `S_(o,k)` and has
these eight boundary values.  Its pullback of the ambient Minkowski form is
the piecewise-flat Regge metric `g_o`.

The already certified 148-chamber arrangement is a common **parametric**
refinement of all 24 staircases.  The present question is whether it is also
a common **metric** refinement:

```text
for every full-dimensional overlay chamber C,
are all restrictions g_o|C exactly equal?
```

Equality of the embeddings `X_o` is a useful stronger diagnostic but is not
the acceptance criterion.  Gravity depends on the intrinsic pullback metric,
so the calculation must compare the complete exact Gram matrices on every
chamber.  No Regge action, dust action, continuum target or fitted parameter
is used in this mission.

## 2. Why a common subdivision is not enough

A common refinement is a statement about domains: each refined cell lies in
a cell of every coarse triangulation.  A piecewise-affine field on a coarse
triangulation restricts to the refinement, but boundary values at the common
coarse vertices do not in general force the piecewise-affine fields belonging
to two different triangulations to coincide.

There is an elementary obstruction here.  A single affine map on the whole
prism would have a top-minus-bottom displacement independent of the spatial
vertex.  The prescribed displacement is

```text
(R_+-R_-) u_i,
```

which depends on `i` unless `R_+=R_-`.  Therefore the static case has a global
affine realization, while a genuinely homothetic dynamic boundary does not.
This observation does not by itself decide whether the 24 induced *metrics*
could nevertheless agree; that is the exact exhaustive test being gated.

## 3. Primary prior art

- Francisco Santos, [*The Cayley trick and triangulations of products of
  simplices*](https://arxiv.org/abs/math/0312069), supplies the established
  product-of-simplices and staircase-triangulation setting.  It does not
  identify the 24 boundary interpolants as one metric.
- Bianca Dittrich and Sebastian Steinhaus,
  [*Path integral measure and triangulation independence in discrete
  gravity*](https://arxiv.org/abs/1110.6866), explicitly study the dependence
  of classical Regge calculus on triangulation in three and four dimensions.
  Triangulation independence is therefore not automatic in 4D Regge gravity.
- Benjamin Bahr and Bianca Dittrich,
  [*Improved and Perfect Actions in Discrete
  Gravity*](https://arxiv.org/abs/0907.4323), show that recovering
  discretization-independent continuum dynamics requires an improved/perfect
  action construction.  A combinatorial common carrier alone is not such a
  construction.
- A. De Felice and E. Fabri,
  [*The Friedmann universe of dust by Regge Calculus: study of its ending
  point*](https://arxiv.org/abs/gr-qc/0009093), is the primary 600-cell
  Regge--dust construction inherited by this repository.  Its five-dimensional
  Lorentzian embedding motivates the homothetic boundary data, but it does not
  supply a schedule-independent PL interpolation on the universal overlay.

No primary source located in this gate gives the exact chamber-by-chamber
comparison of all 24 staircase pullback metrics for this homothetic
`Delta^3 x I` boundary problem.  Absence from this search is not proof of
novelty; external novelty remains **OPEN**.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- The 24 staircase triangulations and their common 148-chamber parametric
  overlay are already certified exactly.
- Each staircase plus the eight boundary images determines one continuous PL
  embedding and one piecewise-flat pullback metric.
- At `R_+=R_-`, all boundary data extend to the same global affine map, so all
  schedules must agree.  This is the static control.
- Classical four-dimensional Regge quantities can depend on triangulation.

### CONTROL

- Reconstruct every chamber/order assignment from the frozen overlay artifact.
- Derive each affine Jacobian independently from its five simplex vertices;
  do not insert a guessed interpolation formula as the primary computation.
- Compare exact symbolic pullback Gram matrices, not floating coordinates or
  sampled lengths.
- Verify the closed-form staircase interpolant against the independent vertex
  solve and verify continuity across all staircase internal facets.
- Require complete agreement for all schedules in the static limit.
- Check the census under the full `S4 x C2` action.

### OPEN

- How many dynamic pullback metrics occur among the 24 schedules on each of
  the 148 chambers.
- Whether any chamber is metric-compatible across all 24 schedules.
- Whether distinct PL embeddings could still give the same pullback metric.
- Whether a new symmetric interpolation, averaging rule, curved-simplex
  action or perfect action can select a unique metric after a negative result.
- Any Regge equation, canonical momentum, tick, clock, `c` or Planck scale on
  the universal overlay.

## 5. Framing attack and verdict scope

Comparing vertex images alone would be too strong and could falsely reject
intrinsically identical metrics.  The load-bearing comparison is therefore
the exact Lorentzian Gram tensor on every full-dimensional overlay chamber.

Conversely, agreement only after permuting spatial labels is too weak.  The
24 schedules are alternative triangulations of the same labelled prism; a
canonical metric must agree at the same parameter point without choosing a
schedule-dependent relabelling.

If the dynamic metrics disagree, the honest verdict is narrow but decisive:
the universal overlay removes the combinatorial schedule choice but does not
inherit a unique Regge metric from the old staircases.  It would not kill
Regge gravity or the homothetic 600-cell model.  It would kill the claim that
common refinement by itself reconciles their dynamics, and any subsequent
metric would need an explicit additional selection principle.
