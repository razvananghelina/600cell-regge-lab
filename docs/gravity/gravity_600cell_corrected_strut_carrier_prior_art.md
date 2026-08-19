# Prior-art gate: geometry-selected strut carrier on the dynamic 600-cell slab

Date: 2026-08-19

## Exact object and hypotheses

Fix the already accepted non-static Lorentzian Regge--dust slab, separately
for both order-24 staircase parities.  Its lower and upper tetrahedral cells
are homothetic with positive edge-length ratio

```text
lambda = L_plus/L_minus != 1,
```

and each lateral quadrilateral is a planar Lorentzian trapezoid.  Keep all
upper spatial edge data fixed and vary independently the signed squared
strut lengths `s_v` at the 120 vertices.  The complete variable-face gluing
calculation has already proved that these 120 pure-strut directions are a
subspace of the exact 240-dimensional admissible data carrier.

For an oriented staircase diagonal from lower vertex `u` to upper vertex
`v`, with `{u,v}` a boundary edge, the new object is the induced first-order
signed-squared diagonal variation.  It must be derived from the trapezoid,
not selected from a singular vector or Hessian.

Let the bottom edge lie along a unit vector `e`, let the displacement from
the lower to upper copy at `u` be `a`, and suppress the common edge scale.
Then

```text
s_u = a^2,
s_v = (a + (lambda-1)e)^2,
d_uv = (a + lambda e)^2.
```

At fixed boundary scales, eliminating `delta(a.e)` gives, for
`lambda != 1`,

```text
delta d_uv = (-delta s_u + lambda delta s_v)/(lambda-1).       (1)
```

This is the disclosed geometric target.  It is singular on the static
stratum and no limit through `lambda=1` is assumed.

The accepted action uses logarithmic signed-squared edge variables.  If
`rho>0` is the common timelike strut magnitude and

```text
q_diag = lambda L_minus^2-rho > 0,
```

normalize column `c_v` by `delta log(rho_v)=c_v`.  Equation (1) then gives

```text
delta log(q_(u,v+120))
  = rho/((lambda-1) q_diag) * (c_u-lambda c_v),
delta log(rho_v) = c_v,
delta log(q_new) = 0.                                         (2)
```

The resulting `1560 x 120` matrix is the corrected geometric strut carrier
in the pre-Legendre input coordinates `(840 internal,720 new boundary)`.

## Why the old geometric lapse carrier is not the same object

The earlier full-lapse and hyperbolic-alignment audits used a frozen local
normal-displacement ansatz: a column at `v` varied the pole at `v` and only
the staircase diagonals ending at `v`.  The later complete face-gluing audit
refuted that local representative for non-uniform vertex data.  It selected
instead a unique local Poincare correction while leaving the natural upper
edge and strut data unchanged.

Equation (2) contains both endpoints of every oriented diagonal.  Its
uniform column sum reproduces the old collective lapse exactly because

```text
(1-lambda)/(lambda-1) = -1,
```

but its 119 relative directions are different.  Therefore the earlier
negative alignment theorem remains valid for its frozen carrier and does
not decide the corrected carrier.

## Primary-literature boundary

- Dittrich and Hoehn, [*Canonical simplicial gravity*,
  arXiv:1108.1974](https://arxiv.org/abs/1108.1974), derive discrete
  evolution from the action and emphasize that initially free data can be
  fixed later by pre/post constraints.
- Hoehn, [*Canonical linearized Regge Calculus: counting lattice gravitons
  with Pachner moves*, arXiv:1411.5672](https://arxiv.org/abs/1411.5672),
  separates vertex lapse/shift displacements from curvature-carrying lattice
  gravitons on flat backgrounds.
- Bahr and Dittrich, [*(Broken) Gauge Symmetries and Constraints in Regge
  Calculus*, arXiv:0905.1670](https://arxiv.org/abs/0905.1670), explain why
  curved Regge backgrounds generically replace exact gauge generators by
  background-dependent pseudo-constraints.
- Jercher and Steinhaus, [*Cosmology in Lorentzian Regge calculus*,
  arXiv:2312.11639](https://arxiv.org/abs/2312.11639), treat Lorentzian
  frustum heights as dynamical variables and distinguish static
  time-reparametrization freedom from matter-supported evolution.

These works make vertex displacement, lapse-like data, frustum heights and
pseudo-constraints **KNOWN** mechanisms.  None of the located primary
sources supplies equation (2) on this oriented 600-cell staircase or compares
its 120-dimensional image with this repository's Schur/tangent sectors.
Search absence is not a novelty proof; external novelty remains **OPEN**.

## KNOWN / CONTROL / OPEN

- **KNOWN:** planar trapezoid length identities and action-generated
  simplicial pre/post constraints.
- **DERIVED INPUT:** pure struts are admitted by the exact complete
  variable-face carrier, and the old local normal representative fails
  complete gluing for non-uniform data.
- **CONTROL:** equation (1) must be reproduced independently by direct
  differentiation and by both committed rational local response blocks.
- **CONTROL:** the 120 columns must have rank 120, be equivariant under the
  schedule stabilizer, and sum exactly to the previously frozen collective
  lapse column.
- **OPEN:** whether the corrected carrier equals the canonical pole-Schur
  lift or either hyperbolic invariant subspace.
- **OPEN:** whether any such alignment is gauge, a pseudo-constraint, or a
  physical mode; curvature response is still required.

## Anti-fitting separation

The corrected carrier will first be constructed and committed without
loading any Hessian, Schur, tangent, singular-vector or eigenvector artifact.
Only after that commit may a separate protocol compare it with the frozen
dynamic targets.  The first phase cannot report an alignment or physical
interpretation.

