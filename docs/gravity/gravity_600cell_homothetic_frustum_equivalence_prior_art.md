# Prior-art gate: homothetic 600-cell prism as one Lorentzian 4-frustum

Date: 2026-08-17

Status: written after the identity-pullback result `b56811c`, but before any
evaluation of physical frustum coplanarity, signature, volume or geometric
triangulation.

## 1. Framing correction and exact question

The preceding audit compared the 24 tensors

```text
g_o = X_o^* eta
```

at the same parameter point of a fixed copy of `Delta^3 x I`.  It proved that
they agree under the identity parametrization only when `R_plus=R_minus`.
That comparison is mathematically correct but may be physically too strong.
If every `X_o` is a PL homeomorphism onto the same flat Lorentzian polytope
`Q`, then

```text
h_(o,o') = X_o'^{-1} o X_o
```

is an isometry between the pullback metrics even when they are unequal under
the identity.  In gravity, this distinction between a tensor and its chosen
coordinates is load-bearing.

The new question is therefore:

> Do the 24 homothetic staircase realizations triangulate one and the same
> flat Lorentzian tetrahedral 4-frustum, with the induced isometries fixing
> the two physical time boundaries?

This mission tests geometric equivalence, not equality under a preselected
interior point identification.

## 2. Complete hypotheses

Let `u_0,...,u_3` be the unit directions of one regular tetrahedral cell of
the 600-cell, with

```text
u_i.u_i=1,
u_i.u_j=c=phi/2  for i != j,
phi=(1+sqrt(5))/2.
```

For `R_minus,R_plus,rho>0`, put

```text
p_i^- = (R_minus u_i,0),
p_i^+ = (R_plus  u_i,T),
Delta R = R_plus-R_minus,
T^2 = rho+(Delta R)^2.
```

The target carries signature `(4,1)`, positive on spatial coordinates and
negative on time.  The candidate physical cell is

```text
Q = conv{p_i^-,p_i^+ : i=0,...,3}.
```

For each of the 24 orders and four splits, use exactly the same five-vertex
staircase simplices as in the committed overlay audit.  No averaging,
continuum scale factor, Regge equation, dust action or fitted coordinate is
introduced.

## 3. Analytic structure to be checked

Set

```text
s = sum_i u_i,
D = 1+3c,
n = s/D,
```

so `n.u_i=1` and `n.n=4/D`.  Every one of the eight vertices is predicted to
lie in the affine hyperplane

```text
n.x_space-(Delta R/T) x_time = R_minus.
```

Its Minkowski normal has square

```text
N.N = 4/D-(Delta R)^2/T^2.
```

Since `T^2=rho+(Delta R)^2`, the second term is strictly below one.  For the
600-cell tetrahedron `4/D>1`, so the normal is spacelike and the hyperplane is
Lorentzian.  This would make `Q` a flat four-dimensional spacetime cell, not a
curved image requiring a fitted interior interpolation.

The verifier must rederive all of this exactly and independently check rank,
signature and the geometry of all staircase simplices.  The argument is not
accepted merely because the formula looks plausible.

## 4. Primary prior art

This general route is **KNOWN**, and that substantially changes the novelty
ledger.

- Tsuda and Fujiwara, [*Oscillating 4-Polytopal Universe in Regge
  Calculus*](https://arxiv.org/abs/2011.04120), use the Collins--Williams
  formalism with a tetrahedral cell evolving into a scaled tetrahedral cell;
  the intervening four-dimensional frustum is their fundamental block.  Their
  regular-polytopal cases include the 600-cell.
- De Felice and Fabri, [*The Friedmann universe of dust by Regge Calculus:
  study of its ending point*](https://arxiv.org/abs/gr-qc/0009093), use the
  five-dimensional Lorentzian embedding inherited by this repository for
  600-cell dust evolution.
- Jercher and Steinhaus, [*Cosmology in Lorentzian Regge calculus: causality
  violations, massless scalar field and discrete
  dynamics*](https://arxiv.org/abs/2312.11639), formulate modern Lorentzian
  Regge cosmology directly with flat 4-frusta and emphasize causal character
  of their subcells.  Their explicit model is cubical/spatially flat, not the
  present tetrahedral closed 600-cell.
- Gentle and Miller, [*A fully (3+1)-D Regge calculus model of the Kasner
  cosmology*](https://arxiv.org/abs/gr-qc/9706034), explicitly exhibit
  diffeomorphism freedom in a four-dimensional Regge evolution.  Equality of
  coordinate components is therefore not the correct general equivalence
  criterion.
- Santos, [*The Cayley trick and triangulations of products of
  simplices*](https://arxiv.org/abs/math/0312069), supplies the standard
  staircase triangulations of `Delta^3 x I`.

Thus using a tetrahedral 4-frustum as a cosmological block is not new.  What
remains **OPEN** externally is only the exact reconciliation of this
repository's 24 staircase schedules and the correction of its universal
overlay interpretation.  A search cannot prove even that narrower novelty.

## 5. KNOWN / CONTROL / OPEN

### KNOWN

- Four-dimensional Regge cosmology can be formulated with flat frusta rather
  than treating a simplex schedule as fundamental.
- The homothetic 600-cell embedding fixes the same-vertex proper strut and the
  cross-diagonal square `L_minus L_plus-rho`.
- The repository's homothetic action and canonical roots already agree between
  its independently derived schedule parities to their certified precision.
- The identity-pullback metrics differ for every `R_plus != R_minus`; that
  exact statement remains true.

### CONTROL

- Derive the unique affine hull of the eight vertices twice: from the exact
  Gram data and from an independent rank/nullspace calculation.
- Prove the induced hyperplane signature is `(3,1)` for every `rho>0`, not at a
  numerical sample.
- Exhibit intrinsic affine coordinates in which
  `Q={0<=z<=1, y_i>=0, sum_i y_i=R_minus+z Delta R}` and certify its eight
  vertices and six facets.
- For all 96 labelled four-simplices, compute exact oriented volume
  polynomials and require nondegeneracy for all positive scales.
- For every schedule, require coherent orientation, the correct boundary and
  total simplex volume equal to the exact frustum volume.
- Independently verify the projective/PL equivalence to the standard prism,
  or give an equally strong injectivity and non-overlap proof.
- Show mechanically why `X_o'^{-1} o X_o` is an isometry and why it fixes both
  full spatial time boundaries pointwise.

### OPEN

- Whether all 24 staircase sets are geometric triangulations of the same `Q`
  for every positive pair of scales.
- Whether the isometries glue globally across all 600 tetrahedral frusta for
  the repository's actual global schedules.
- Whether a cellular frustum Regge action exactly equals the schedule actions
  on the homogeneous sector, including boundary terms and dust.
- Whether frustum equivalence survives anisotropic boundary perturbations.
- Any refinement limit, physical clock, causal propagation speed or Planck
  scale.

## 6. Framing attack and possible outcomes

Coplanarity alone is insufficient: a set of nondegenerate simplices could
fold or overlap rather than triangulate the convex hull.  Matching total
absolute volume alone is also insufficient without coherent orientation and
boundary incidence.  All these gates are required.

If every schedule triangulates the same `Q`, then the prior result must be
relabelled as a **coordinate-identity obstruction**, not a dynamic metric
no-go.  The common physical cell is selected before triangulation, and the
naive product/order-complex carrier becomes relevant again as a canonical
barycentric subdivision of that frustum.

If the simplices do not give the same geometric `Q`, the attempted correction
fails and the identity-pullback no-go remains load-bearing.
