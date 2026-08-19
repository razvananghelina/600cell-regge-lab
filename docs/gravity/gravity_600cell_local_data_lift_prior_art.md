# Prior-art and framing gate: local star-supported data lift

Date: 2026-08-19

This gate precedes comparison of the exact rational lift's support with
600-cell vertex stars or construction of a universal local response block.

## Proposed object

For each tetrahedral cell with ordered vertices `(v0,v1,v2,v3)`, propose one
exact 6-by-8 block mapping

```text
(sigma_v0,...,sigma_v3, strut_v0,...,strut_v3)
```

to its six cell-flex coordinates.  The same block, in local vertex order, must
satisfy every complete face equation.  If it exists and all 48 entries are
nonzero, the support is exactly the four vertex-scale and four strut data in
that cell, while each datum is supported on its 20-cell vertex star.

## What is already known externally

Vertex-associated scale factors for edge lengths are established discrete-
conformal prior art.  Bobenko, Pinkall, and Springborn formulate discrete
conformal equivalence using vertex scale factors in
[Discrete conformal maps and ideal hyperbolic polyhedra](https://arxiv.org/abs/1005.2698).
Glickenstein develops conformal variations on piecewise-flat two- and
three-manifolds and their Regge-curvature variations in
[Discrete conformal variations and scalar curvature on piecewise flat two and
three dimensional manifolds](https://arxiv.org/abs/0906.1560).

Local vertex-displacement and lapse/shift structures in flat linearized Regge
calculus are also known; see Hoehn,
[Canonical linearized Regge Calculus: counting lattice gravitons with Pachner
moves](https://arxiv.org/abs/1411.5672).  A local star-supported response is
therefore structurally plausible and is not, by itself, new physics.

The 600-cell has already been evolved with a generalized set of free
variables using a Sorkin scheme; see De Felice and Fabri,
[Singularities of the closed RW metric in Regge Calculus: a generalized
evolution of the 600-cell](https://arxiv.org/abs/gr-qc/0106077).  The present
question is narrower: the exact local response block for this project's
complete variable-face flat-frustum linearization.

## Framing attacks fixed before calculation

1. Matching support counts is insufficient.  Membership in every vertex star
   must be checked exactly.
2. A local 6-by-8 block is not accepted because it can be fitted to one cell.
   Its 48 coefficients must be determined once and satisfy all 6000 global
   face equations.
3. The new computation must not reuse the 3600-variable global elimination.
   It must solve the much smaller universal-block ansatz directly from the
   face equations and then substitute it globally.
4. The old local block remains a frozen negative.  A positive new block must
   be compared coefficient-by-coefficient and the difference explained as a
   length-preserving local Poincare correction, not silently substituted for
   the failed formula.
5. Locality is kinematic.  It neither identifies physical degrees of freedom
   nor supplies an action, clock, propagation speed, or coupling constant.

## Novelty status

The bounded primary-literature search found the general local and conformal
mechanisms but no exact 6-by-8 block for this complete 600-cell construction.
This absence is not proof of novelty.

**OPEN:** external novelty of the explicit block.

**STRUCTURAL:** a positive result would be an exact internal simplification
and an essential input to the later boundary-Hessian calculation.

