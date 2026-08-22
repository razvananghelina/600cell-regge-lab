# Prior-art and framing gate: finite-height internal carrier rank

Date: 2026-08-22.

Status: **FROZEN BEFORE EVALUATING ANY RANK, SINGULAR VALUE OR NULL VECTOR OF
THE FINITE-HEIGHT INTERNAL-EQUATION MAP.**

## Exact next question

Fix the first positive-height homogeneous 600-cell dust slab generated from
the incoming state `v=3/2`, with the same hypotheses and conventions as the
accepted quadratic parity result.  For staircase parity `p`, split the active
action Hessian into its `840` internal equation rows and all `1560` active
columns, and let

```text
G_p : R^240 -> R^1560
```

be the exact geometry-selected scale-plus-strut tangent carrier.  Define

```text
R_p = H_p[internal,active] G_p : R^240 -> R^840.
```

The background already satisfies all internal equations.  Hence

```text
R_p x = 0
```

is exactly the necessary first-order condition that the kinematic response
`G_p x` remain tangent to the internal stationary locus with the old boundary
held fixed.

The `840` rows must also be split geometrically into

```text
720 oriented cross-diagonal equations,
120 same-vertex pole/strut equations.
```

Report the kernels of the diagonal-only map and the complete map separately.
This distinguishes directions killed by local diagonal stationarity from
directions killed only when the pole/lapse equations are imposed.

No desired nullity, continuum representation, mode spectrum, limiting speed
or fitted source combination may be inspected before this rank classification
is frozen.

## Why this precedes a nonlinear carrier calculation

The phrase "does the infinitesimal carrier integrate?" is not by itself a
falsifiable question.  Let `U` be the open set of logarithmic signed squared
edge data around the nondegenerate Lorentzian background.  For any rank-240
linear map `G`,

```text
Phi(x)=y0+G x
```

maps a sufficiently small neighbourhood of zero smoothly into `U` and has
derivative `G`.  Thus bare local existence of some nonlinear extension is
automatic and carries no geometric evidence.

Nor is an arbitrary second jet invariant.  Under a source-coordinate change

```text
x -> x + (1/2) B(x,x)+O(x^3),
```

the second jet changes by `G B`.  Only its class modulo `im G`, equivalently a
normal second fundamental form after an independently fixed nonlinear
submanifold, can be geometrically meaningful.

A stronger future question can fix exact vertex scaling on the upper
boundary, exact logarithmic strut scaling and exact nonlinear local-frame/
shared-face gluing, then ask whether those data select a unique normal second
jet.  But if `ker R_p=0`, no nonzero first-order kinematic direction reaches
that question while satisfying the action.  The internal rank is therefore
the logically earlier kill gate.

## Relation to the accepted quadratic form

The accepted result compares

```text
Q_p=G_p^T H_p G_p.
```

Equality of `Q_even` and `Q_odd` does not imply a nonzero kernel of `R_p`.
A bilinear form can agree between schedules even when the carrier is
transverse to every internal equation.  Conversely, a nonzero `ker R_p`
would select action-compatible candidates but would not yet prove gauge
reduction, propagation or a physical mode.

The Hessian coordinate caveat remains important.  For an actual nonlinear
embedding with second jet `K`,

```text
d^2(S o Phi)=G^T H G+grad(S) dot K.
```

The accepted parity comparison is valid because the schedule-dependent
internal gradients vanish and the physical-boundary gradients and embedding
are common.  A future physical restricted action, rather than a parity
difference, will still require a selected nonlinear boundary convention or
its equivalent canonical formulation.

## Repository non-duplication

The repository contains three nearby but distinct results.

1. **Older complete-carrier canonical intersection.**  At a different,
   near-static accepted slab, all nonhomogeneous minimal sectors of the
   scale-plus-strut carrier have zero intersection with an action-generated
   strong-equation graph.  A homogeneous weak-pole line survives the reduced
   test and is then removed by the complete pole equation.  That is a
   load-bearing warning, but its background, matrix object and normalization
   differ from `R_p` above.
2. **Older nonlinear boundary covariance no-go.**  On the order-24 quotient
   at the same older background, even and odd nonlinear canonical maps agree
   to first order and split quadratically on all 32 frozen rays.  It covers
   neither the present finite-height state nor the full rank-240 carrier.
3. **Current finite-height quadratic parity result.**  The two schedule
   parities induce the same one-sided pulled-back quadratic form within a
   180-digit calibrated bound.  That verifier did not compute or classify
   `R_p`.

An exhaustive repository search found no rank census of the present
`840 x 240` finite-height maps.  This establishes internal non-duplication,
not external novelty.

## Primary literature

### KNOWN

- Dittrich and Höhn derive action-generated discrete canonical evolution and
  show that higher-order Regge dynamics can turn linear constraints into
  background-dependent pseudo-constraints in
  [arXiv:0912.1817](https://arxiv.org/abs/0912.1817).
- Bahr and Dittrich show that curved Regge backgrounds generically break exact
  discrete gauge symmetry in
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670).
- Dittrich and Höhn formulate the general constrained canonical evolution of
  simplicial systems in
  [arXiv:1108.1974](https://arxiv.org/abs/1108.1974).
- Höhn separates gauge and propagating lattice-graviton variables only after
  the action Hessian and constraint structure have been analysed in
  [arXiv:1411.5672](https://arxiv.org/abs/1411.5672).
- Glickenstein treats vertex-based conformal variations of piecewise-flat
  two- and three-manifolds and the Regge scalar-curvature functional in
  [arXiv:0906.1560](https://arxiv.org/abs/0906.1560).  Bobenko, Pinkall and
  Springborn give the exact exponential vertex-scaling notion for
  triangulated surfaces in
  [arXiv:1005.2698](https://arxiv.org/abs/1005.2698).

These sources make the general mechanisms known.  They do not supply the
present finite-height `600`-cell internal rank.

### CONTROL

- The primary finite-height verifier must reproduce the internal-gradient
  stationarity and exact carrier rank/support before forming `R_p`.
- The two parities must be reconstructed separately.
- A zero matrix and an embedded identity must reproduce their known
  nullities under the same classifier.
- A frozen carrier coefficient corruption and a frozen internal-Hessian row
  corruption must change the classified matrix above the calibrated error.
- Rank must be certified by precision/derivative stability and a
  mechanically different method; one binary64 SVD threshold is insufficient.

### OPEN

- the diagonal-only and complete ranks/nullities for both parities;
- whether their kernel projectors agree in the common 240 data coordinates;
- external novelty of the exact finite computation;
- the physical classification of any survivor.

## Falsifiable outcomes

### `FINITE_HEIGHT_INTERNAL_CARRIER_FULL_COLUMN_RANK`

Use only if both parity maps have certified complete rank `240`.  Then

```text
ker R_even = ker R_odd = 0.
```

This is a **DERIVED NEGATIVE** for nonzero action-compatible perturbations on
the present geometry-selected carrier at this slab.  It closes this carrier
route at first order; it does not refute Regge calculus or more general
boundary data.

### `FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_SELECTED`

Use only if both parities have a resolved positive nullity and their kernel
projectors agree within the preregistered error.  This selects candidates for
the next canonical/gauge calculation; it does not yet identify gravitons.

### `FINITE_HEIGHT_INTERNAL_CARRIER_SCHEDULE_DEPENDENT`

Use if the resolved nullities differ, or equal-dimensional kernel projectors
are resolved inequivalent.  Then equality of the scalar quadratic forms was
not enough to remove schedule ambiguity from the internal equations.

### `CONTROL_FAILED` or `NUMERICALLY_OPEN`

Use for any failed prerequisite or unresolved singular-value/minor band.  No
physical conclusion follows.

No outcome derives a tick, `c`, `G`, a Planck scale, particle mass or
Standard-Model sector.
