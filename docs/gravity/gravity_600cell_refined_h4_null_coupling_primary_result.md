# Primary result: the product-lapse null line imposes one boundary compatibility condition

Date: 2026-08-21

Status: primary result; not accepted before the frozen adversarial
boundary-hinge reconstruction.

## Provenance and reproducibility

| stage | commit |
|---|---|
| prior-art gate | `32854ba` |
| frozen protocol | `752891b` |
| verifier registered before execution | `fab5f4d` |

The targeted directional verifier passed `16/16` twice and produced the
byte-identical artifact

```text
reproducible/gravity_600cell_refined_h4_null_coupling.json
SHA-256 6b6fbd95b07f365b3fcac332fa3546021e8d756a510af0184bc974e52d5efa79.
```

No full Hessian, pseudoinverse, Schur complement, root search, nested census,
spectrum or full suite was run.

## Finite product-family control

For every one of the 24 schedules and each

```text
tau/tau0 in {1/2,1,2},
```

the exact product cross diagonals and lapse squares remain internally
stationary under the fixed curvature-selected masses.  Across all 72 points,

```text
maximum internal residual      1.19e-76
maximum angle-identity error   6.32e-138
maximum imaginary curvature    1.65e-137.
```

This confirms a finite one-parameter stationary family, not merely one small
eigenvalue at `tau0`.

## The exact internal null line

With `u=log(tau^2)`, the geometry-derived tangent is

```text
n_cross,rs=-tau0^2/q_cross,rs,
n_rho,r=1.
```

Directional differentiation of the complete 22-component gradient gives

```text
H_ii n=0
```

for all 24 schedules.  The maximum displayed internal image is `4.11e-44`,
inside a maximum propagated envelope `9.74e-41`.  Combined with the frozen
`(9,1,0)` internal inertia, this identifies the full one-dimensional internal
kernel without choosing an eigensolver sign.

## The line is not a full-Hessian gauge direction

Before loading the frozen boundary vector, the computation produced

```text
c=H_bi n=
(0.0182753971831513704...,
 0.00121804711711642587...,
 0.0000327518273255248438...,
 0.0366773216171774596...,
 0.000494989505267934317...,
 0.0191118520096157871...,
 0.0182753971831513704...,
 0.00121804711711642587...,
 0.0000327518273255248438...,
 0.0366773216171774596...,
 0.000494989505267934317...,
 0.0191118520096157871...).
```

Componentwise,

```text
c=(1/2)*g_boundary=(1/2)*(-P_pre,P_post),
```

with maximum direct and frozen-vector disagreement `6.50e-44`.  The minimum
ratio of a resolved coupling component scale to its numerical envelope is
`3.76e38`; the nonzero decision is not marginal.

All 24 unnormalised rows agree within their envelopes, so their row rank is
one.  The maximum schedule difference is `9.18e-44`, and the maximum fixed
old/new time-reversal difference is `2.19e-125`.

The frozen outcome is

```text
REFINED_H4_NULL_COUPLING_COMPATIBILITY_CONFIRMED.
```

## Exact meaning

The linearized internal equation

```text
H_ib delta_b+H_ii delta_i=0
```

has the necessary compatibility condition

```text
c^T delta_b=0.                                    (1)
```

- **PRIMARY DERIVED COMPUTATIONAL:** the product duration is a true flat
  internal direction on this static curvature-matched branch.
- **PRIMARY DERIVED COMPUTATIONAL:** it is not a null direction of the full
  Hessian.  It changes the boundary momenta, exactly at half their log-lapse
  scaling rate.
- **PRIMARY DERIVED STRUCTURAL:** equation (1) is one schedule-independent
  linear boundary constraint in the invariant 12-coordinate system.
- **DERIVED NEGATIVE:** an ordinary unconstrained Schur complement and an
  arbitrary Moore--Penrose inverse remain forbidden.
- **OPEN:** the quadratic form reduced simultaneously by the internal null
  line and the boundary compatibility hyperplane.
- **NOT ESTABLISHED:** that (1) is the continuum Hamiltonian constraint, a
  physical tick, propagation, `c`, `G`, Planck units or particle physics.

The factor `1/2` was derived after the singularity was seen, from
`d log(tau)/d log(tau^2)`.  Its provenance is therefore
**POST-PRIMARY DERIVED PREDICTION**, not a blind hit.

## Adversarial gate

Do not differentiate the action again.  Rebuild all actual spatial edge
curvatures and use the independent local hinge identity

```text
g_old,rs=g_new,rs=tau0*C_rs/4,
c_rs=tau0*C_rs/8.
```

Require all twelve components, the internal product-family identity, the
rank-one schedule statement and corruption controls before acceptance.

