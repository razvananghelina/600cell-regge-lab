# Prior-art gate: product-lapse null line and boundary compatibility

Date: 2026-08-21

Status: written after the primary internal-singularity result and before any
new null-vector or boundary-coupling evaluation.

## Exact object and complete hypotheses

Use the same internally on-shell static product over
`K0=P(sd K_600)`, all 24 staircase schedules, `tau0=0.0102`, exact projected
rank geometry and the independently selected curvature masses.  Boundary and
internal total-orbit log squared-edge coordinates retain the frozen `12+10`
ordering.

The primary audit found one compatible zero eigenvalue of every internal
Hessian `H_ii`.  This gate asks whether that line is:

1. the analytically induced product-lapse tangent;
2. a null direction of the complete Hessian, or instead couples to boundary
   perturbations through `c=H_bi n`;
3. schedule-independent after the fixed time-reversal identification.

No pseudoinverse or effective Hessian is in scope.

## Linear-algebra distinction

For a symmetric block Hessian

```text
H=[[H_bb,H_bi],
   [H_ib,H_ii]],
```

and `H_ii n=0`, the linearized internal equation is

```text
H_ib delta_b + H_ii delta_i=0.
```

Multiplying by `n^T` gives the necessary boundary compatibility condition

```text
c^T delta_b=0,  where c=H_bi n.                  (1)
```

- If `c=0`, the null line decouples and may be a gauge direction; elimination
  can be made independent of a gauge fixing.
- If `c!=0`, the internal equations are soluble only on the codimension-one
  boundary subspace `ker(c^T)`.  An unconstrained Moore--Penrose Schur
  complement would be mathematically and physically misleading.

This is finite-dimensional linear algebra, not a physical interpretation.

## Geometry-derived tangent and disclosed post-result prediction

Let `u=log(tau^2)`.  Along the exact static product family, at fixed spatial
boundary edges,

```text
rho_r=tau^2,
q_cross,rs=q_spatial,rs-tau^2.
```

The unnormalised internal tangent is therefore fixed without an eigensolver:

```text
n_cross,rs=-tau^2/q_cross,rs,
n_rho,r=1.                                        (2)
```

The local curvature-mass identity holds for every positive `tau` on this
product branch, so differentiating the ten vanishing internal equations
predicts `H_ii n=0`.

After seeing the primary singularity, the already accepted boundary-hinge
formula supplies a sharper analytic prediction.  Every boundary action
gradient is proportional to `tau`, while (2) differentiates with respect to
`u=log(tau^2)`.  Hence

```text
c=H_bi n=(1/2)*g_boundary,                         (3)
```

where in the frozen canonical convention

```text
g_boundary=(-P_pre, P_post)=(-P_pre,-P_pre).
```

Equation (3) is a **POST-PRIMARY DERIVED PREDICTION**, not a blind prediction.
It is nevertheless falsifiable component by component and contains no fitted
coefficient: the factor `1/2` is `d log(tau)/d log(tau^2)`.

Because the accepted boundary covector is nonzero, (3) predicts that the
internal null line is not a null line of the complete Hessian.  It predicts
one common rank-one family of compatibility covectors across all schedules.

## KNOWN from primary literature

- Exact Hessian null vectors on flat Regge backgrounds encode vertex
  displacement gauge symmetries and canonical constraints: Dittrich and
  Hoehn, [arXiv:0912.1817](https://arxiv.org/abs/0912.1817), DOI
  `10.1088/0264-9381/27/15/155001`.
- Linearized canonical Regge calculus distinguishes lapse/shift gauge
  variables, constraints and propagating lattice-graviton variables: Hoehn,
  [arXiv:1411.5672](https://arxiv.org/abs/1411.5672), DOI
  `10.1103/PhysRevD.91.124034`.
- In an effective Regge Hessian, a gauge-fixed internal inverse is independent
  of the gauge parameters only when the internal null vectors also decouple
  from the boundary block.  An explicit area-Regge example is given by
  Dittrich and Kogios,
  [arXiv:2203.02409](https://arxiv.org/abs/2203.02409), DOI
  `10.1088/1361-6382/acc5d9`.

These sources establish the distinction, not equations (2)--(3) for this
projected 600-cell.  External novelty remains **OPEN**.

## Framing attack

- Calling every lapse null direction "gauge" is circular.  The complete
  coupling test (1) decides whether gauge-independent elimination is allowed.
- A nonzero `c` is not yet a propagating Hamiltonian constraint of continuum
  GR.  It is a linear compatibility covector in the invariant finite system.
- Rank one is expected once one common product parameter is identified; its
  dimension alone carries no evidential weight.  The signal is the exact
  componentwise relation (3), independently reconstructed.
- This `H4` sector still cannot see nonhomogeneous gravitons or determine a
  causal speed.

## Next admissible calculation

Preregister a directional high-precision verifier that differentiates the
complete 22-component gradient only along (2), rather than rebuilding the
full Hessian.  Test `H_ii n=0`, equation (3), schedule/time-reversal covariance
and rank of the 24 compatibility rows.  Then perform a mechanically distinct
analytic boundary-hinge reconstruction before accepting the result.

