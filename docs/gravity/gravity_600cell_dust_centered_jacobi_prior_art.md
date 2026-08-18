# Prior-art gate: centered mass--drift--stiffness decomposition

Date: 2026-08-18

## Exact object, carrier and hypotheses

Start only from the blindly committed three-slice Regge Jacobi balls

```text
K_- delta q_0 + K_0 delta q_1 + K_+ delta q_2 = 0
```

on the complete `720`-position boundary carrier, in the same seven minimal
binary-tetrahedral sectors, four derivative variants and two staircase
schedules.

Define the centered coefficients uniquely by

```text
M = (K_- + K_+)/2,
N = (K_+ - K_-)/2,
V =  K_- + K_0 + K_+.
```

Then the exact algebraic identity is

```text
M (delta q_2 - 2 delta q_1 + delta q_0)
+ N (delta q_2 - delta q_0)
+ V delta q_1 = 0.
```

If `M` is regular, define only as finite matrices

```text
Gamma = M^-1 N,
Omega = M^-1 V.
```

No proper-time division is made.  Thus `Gamma` and `Omega` are respectively
per-tick and per-tick-squared coefficients only in the declared dimensionless
coordinate convention; they are not yet physical damping or frequency
operators.

## Primary prior art

- [Marsden--West](https://doi.org/10.1017/S096249290100006X) give the general
  discrete variational framework in which a regular discrete Lagrangian
  yields a three-point linearized recurrence.
- [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974) formulate the analogous
  action-generated canonical evolution for simplicial gravity.
- [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817) show that the linearized
  coefficients and pseudo-constraints depend on the Regge background.
- [Hoehn](https://arxiv.org/abs/1411.5672) identifies physical lattice
  gravitons only after curvature and gauge reduction, not from a raw
  recurrence spectrum.
- [Rostworowski](https://arxiv.org/abs/1902.05090) supplies a continuum FLRW
  control: after constraint reduction, metric perturbations obey wave-type
  master equations on a changing background, while a matter mode obeys a
  transport equation.
- [Christiansen](https://arxiv.org/abs/1106.4266) relates the quadratic
  three-dimensional Regge action to a continuum `curl^T curl` operator and
  proves convergence in its stated Euclidean setting.  It does not identify
  the present four-dimensional Lorentzian dust stencil.

Finite-difference decomposition into second difference, first difference and
zero-frequency stiffness is therefore **KNOWN ALGEBRA**.  Neither the formula
nor the words “mass/drift/stiffness” are a discovery.

No located primary source prints the complete centered coefficients for this
fixed dust 600-cell carrier.  External novelty is **OPEN**, not inferred from
search failure.

## KNOWN / CONTROL / OPEN

- **KNOWN:** the centered decomposition is an exact change of variables in a
  three-point stencil.
- **CONTROL:** the source Jacobi artifact is target-free, regular,
  variational, schedule robust and byte-reproducible.
- **CONTROL:** the raw future/past asymmetry is already resolved at about
  `0.22--0.24`; a first-difference term must not be dropped by assumption.
- **OPEN:** whether every centered `M` is invertible.
- **OPEN:** inertia of the Hermitian part of `M` and whether it is definite.
- **OPEN:** non-Hermiticity of `M`, `N`, `V` under the literal boundary
  identification.
- **OPEN:** reality and schedule robustness of the finite generalized
  stiffness spectrum `Omega`.
- **OPEN:** whether any part of `V` is an intrinsic spatial tensor operator.
- **OPEN:** proper-time normalization, continuum refinement, dispersion and
  limiting speed.

## Framing attack

On a time-dependent background, the three slice tangent spaces are distinct.
Identifying their edge labels is geometrically natural on this fixed carrier,
but it is not a derived parallel transport in superspace.  Consequently:

- `M` need not be Hermitian;
- the positivity of its Hermitian part is only a **STRUCTURAL NECESSARY
  DIAGNOSTIC**, not a coordinate-free no-ghost theorem;
- eigenvalues of `Omega=M^-1 V` are not continuum squared frequencies;
- the decomposition can be algebraically exact even when no physical wave
  interpretation survives.

The globally meaningful self-adjoint object remains the full two-slab action
Hessian on `(q_0,q_1,q_2)`.  Its off-diagonal blocks connect different time
fibres and need not become one self-adjoint spatial operator under the literal
identification.

The repository's Whitney/Hodge Laplacian acts on oriented differential
cochains, whereas `delta q` is an unoriented Regge edge-metric variation.
Comparing their eigenvalue lists without a geometry-selected intertwiner
would be a category error.  The same warning applies to importing a desired
continuum Lichnerowicz spectrum.

This mission therefore stops at a blind finite-operator census.  A later
spatial comparison is licensed only if its carrier map and symmetry
intertwiner are independently derived and committed first.

