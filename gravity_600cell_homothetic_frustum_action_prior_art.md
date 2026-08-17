# Prior-art gate: subdivision invariance of the homogeneous 600-cell frustum action

Date: 2026-08-17

Status: written after the geometric-frustum result `6f9cb78` and before any
new evaluation of a cellular or staircase Regge action.

## 1. Exact object and question

The already certified object is one Lorentzian tetrahedral frustum

```text
Q(L_minus,L_plus,rho)
 = conv{(R_minus u_i,0),(R_plus u_i,T): i=0,...,3},
R_minus=phi L_minus,
R_plus =phi L_plus,
T^2=rho+(R_plus-R_minus)^2,
L_minus,L_plus,rho>0.
```

Its induced four-metric has exact signature `(3,1)`.  Every one of the 24
labelled staircase schedules is a geometric triangulation of the same `Q`,
fixing the lower and upper tetrahedra pointwise.

The new question is deliberately narrower than full triangulation
independence of four-dimensional Regge calculus:

> On the regular homothetic 600-cell slab, does the complete Regge action
> with its boundary term and conserved-dust term depend only on the common
> flat frusta, or does it retain a dependence on the staircase subdivision?

If the reduced actions are identical as functions of the same physical
boundary data, then their reduced pre/post momenta and homogeneous equations
are identical by differentiation.  This does **not** by itself establish
equality of the unrestricted anisotropic actions.

## 2. Complete hypotheses

1. The spatial carrier is the regular 600-cell boundary with
   `(N0,N1,N2,N3)=(120,720,1200,600)`.
2. Consecutive spatial boundaries are homothetic, with common squared spatial
   edge lengths `L_minus^2` and `L_plus^2`.
3. Corresponding vertices are connected by equal timelike struts of squared
   length `-rho`, with `rho>0`.
4. Every tetrahedral world-tube is the flat frustum `Q` above.
5. The gravitational action is the ordinary Lorentzian Regge action with the
   additive Regge boundary term; no cosmological term is included.
6. Dust follows the 120 corresponding-vertex struts and has conserved total
   mass normalization already fixed in the committed dust model.
7. The comparison uses the repository's two independently derived global
   staircase parities and a direct cellular-frustum description.  No action
   coefficient, branch, scale or lapse is fitted.
8. Equality is claimed only on a connected causal domain on which all
   relevant areas and dihedral-angle branches remain defined.

## 3. Primary prior art

### 3.1 Cellular frusta are the Collins--Williams construction

Tsuda and Fujiwara, [*Oscillating 4-Polytopal Universe in Regge
Calculus*](https://doi.org/10.1093/ptep/ptab074), replace a spherical Cauchy
surface by a regular 4-polytope and take four-dimensional frusta as the
fundamental spacetime blocks.  Their Table 2 includes the 600-cell with the
same f-vector used here.  For a tetrahedral cell, the two hinge species are
exactly:

- spatial triangular faces on the Cauchy surfaces;
- lateral isosceles trapezoids between consecutive surfaces.

Their Eq. (13) is the cellular Regge action summed over those hinge types.
Thus the direct frustum action is **KNOWN**, not invented by this project.

### 3.2 Flat-block subdivision

The same paper states immediately after its Regge action that non-simplicial
flat blocks may be fully triangulated by adding hinges with vanishing deficit
angles without changing the Regge action.  The geometric mechanism is:

1. newly introduced hinges strictly inside a flat block have zero deficit;
2. a pre-existing polygonal hinge may be split into coplanar triangles with
   a common deficit/exterior angle, while the signed areas add;
3. block volumes add under a genuine triangulation.

This is a statement about subdivision of the *same flat geometry*.  It is not
a theorem that generic four-dimensional Regge solutions or path integrals
are triangulation independent.

Dittrich and Steinhaus, [*Path integral measure and triangulation
independence in discrete gravity*](https://arxiv.org/abs/1110.6866), show why
the broader claim would be false or at least highly nontrivial in four
dimensions: classical and quantum triangulation dependence must be studied
under Pachner moves, and exact independence is special.  The present flat
frustum is therefore a restricted control case, not evidence for general
discretization independence.

### 3.3 Boundary terms and canonical momenta

Dittrich and Hoehn, [*From covariant to canonical formulations of discrete
gravity*](https://arxiv.org/abs/0912.1817) and [*Canonical simplicial
gravity*](https://arxiv.org/abs/1108.1974), formulate the Regge action with
the boundary term so that it is additive under gluing.  The on-shell action,
or the appropriate one-step discrete action, is the generating function for
canonical evolution; derivatives with respect to old and new boundary edge
lengths give pre- and post-momenta with opposite sign conventions.

Consequently, equality of complete reduced actions on an open domain is
stronger than agreement at a finite set of ticks and automatically implies
equality of the corresponding reduced momenta.  Omitting or mismatching the
boundary term would invalidate that inference.

### 3.4 Lorentzian qualification

Tsuda and Fujiwara derive their displayed angle formulae first in Euclidean
signature and then Wick rotate.  Modern explicitly Lorentzian frustum
cosmology is treated by Jercher and Steinhaus, [*Cosmology in Lorentzian
Regge calculus: causality violations, massless scalar field and discrete
dynamics*](https://arxiv.org/abs/2312.11639).  Lorentzian areas and angles can
change branch when hinge causal type changes.  Therefore the repository must
check the causal domain and angle branches rather than importing a Euclidean
formula formally.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- The physical homothetic cell is one flat Lorentzian tetrahedral frustum.
- All 24 local staircases triangulate that frustum without folds or overlaps.
- Direct regular-polytopal Regge cosmology uses exactly the spatial-triangle
  and lateral-trapezoid hinge description.
- A complete Regge boundary term makes the action additive and its boundary
  derivatives define canonical momenta.
- The two committed global parities agree numerically on the already accepted
  homogeneous ticks; that observation is a control, not an identity proof.

### CONTROL

- Reconstruct both committed global schedule carriers without altering their
  incidence or action conventions.
- Certify the global 600-cell counts and the local-to-global frustum gluing.
- Pair the `1,440` staircase vertical triangles into the `720` lateral
  trapezoids and verify coplanarity, common curvature and signed-area
  additivity.
- Identify every subdivision-only hinge and require its total deficit to be
  exactly zero on the homothetic frustum.
- Check spatial boundary hinges, their exterior-angle convention and the dust
  struts separately.
- Compare the complete reduced action as a function, not only at previously
  selected roots.
- Independently compare derivatives/pre-post momenta at preregistered control
  points to catch a symbolic convention error.

### OPEN

- Exact equality of the cellular and staircase actions in this repository's
  Lorentzian sign and boundary conventions.
- Whether equality holds for all 24 local schedules after conforming global
  gluing, beyond the two derived global parities.
- Equality of unrestricted anisotropic gradients or Hessians.
- Whether generic nonhomothetic boundary data admit a flat cellular frustum
  at all.
- A refinement limit, propagating graviton modes, a selected physical lapse,
  `c`, Planck time or Planck mass.
- External novelty of the narrow 24-schedule reconciliation; literature
  search cannot prove novelty.

## 5. Framing attack

The expected equality is close to a standard subdivision theorem.  A
positive result would remove an important artifact objection, but it would
not be new gravitational dynamics.  The physical information would remain
in the cellular frustum action and its matter coupling, not in choosing one
of its triangulations.

Conversely, a mismatch has three possible meanings which must not be merged:

1. an implementation or Lorentzian branch error;
2. an omitted/mismatched Regge boundary term;
3. a genuine subdivision artifact in the repository's simplicial action.

Only after the first two are excluded can the third be reported.

Finally, equality of `S(L_minus,L_plus,rho)` proves only the homogeneous
collective equations and momenta.  It cannot be cited as evidence that the
30 boundary-edge modes or 35 internal orbit directions have the same
anisotropic dynamics, because the homothetic restriction could erase their
differences.

## 6. Proposed difference from the literature

The proposed contribution is not the existence of frustum cosmology or the
flat-block subdivision principle.  It is the exact audit that the
repository's historically constructed 600-cell staircase action really is a
subdivision of that known cellular action, including its particular
Lorentzian branches, boundary signs, dust normalization and canonical
momenta.  Whether this narrow reconciliation is externally new is **OPEN**.
