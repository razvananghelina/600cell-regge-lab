# The continuum Kähler--Dirac operator sees every Hopf spin-two direction

Date: 2026-08-12  
Protocol commit: `31ecea7`  
Registered verifier: `reproducible/verify_hopf_kahler_induced_gravity.py`  
Machine-readable result: `reproducible/hopf_kahler_induced_gravity.json`

## Headline

**STRUCTURAL INDUCED-GRAVITY ADVANCE.**  On the unit round `S3`, the ordinary
heat trace of the continuum de Rham Kähler--Dirac operator has a
nondegenerate curvature response on the complete five-dimensional Hopf
spin-two carrier.

For a fixed-volume homogeneous tracefree metric perturbation

```text
H in Sym^2_0(R^3),
```

the exact quadratic response of the first curvature heat coefficient is

```text
delta^2 A2(H,H) = (8/3) Vol(S3) Tr(H^2),
```

apart from the conventional common `(4 pi)^(-3/2)` factor.  Its inertia on
the five-dimensional carrier is therefore

```text
(positive,zero,negative) = (5,0,0).
```

This is more than the existence of a correctly typed tensor: the theory's
continuum geometric operator actually sees that tensor through curvature.
It is still not a theory-selected gravitational action or a propagating
graviton.

## 1. Complete scope

The calculation assumes the unit round `S3=SU(2)` and the continuum operator

```text
D_g=d+d*_g
```

on the full exterior algebra.  The tested perturbations have constant
components in a left-invariant frame.  They are the five homogeneous TT
fields already identified by the Hopf projectors.  Inversion is an isometry
of the bi-invariant round metric and exchanges left and right, so the same
calculation holds for the opposite-handed homogeneous space.  These are two
spaces of global fields, not ten independent components in one local tensor
fibre; the local fibre dimension remains five.

The metric is varied only for this continuum response calculation.  The
finite incidence operator and fixed Regge metric in the repository have not
thereby been promoted to a variable-metric theory.

No Einstein equation, Regge action, Lorentzian time, source, measured target,
Newton constant or Planck scale enters.

## 2. Curvature derived directly from `SU(2)`

Let `G` be a generic positive symmetric Gram matrix in a fixed
left-invariant Lie-algebra frame with

```text
[e_a,e_b]=2 epsilon_abc e_c.
```

The verifier constructs the Levi--Civita connection from the Koszul formula,
checks metric compatibility and zero torsion, and contracts its Riemann
tensor.  It obtains exactly

```text
R(G)=2 ((Tr G)^2-2 Tr(G^2))/det(G).
```

At the unit round point,

```text
Ric=2g,   R=6.
```

To remove the pure volume direction, use the normalized homogeneous
Einstein functional per fixed coordinate-volume factor,

```text
Y(G)=R(G) det(G)^(1/3).
```

This equals the scalar-curvature integral up to a constant on determinant-one
paths.  The round point is stationary, the scale direction is exactly null,
and on an orthonormal Frobenius basis of `Sym^2_0(R^3)` its Hessian is

```text
Hess(Y)|_STF = -4 I_5.
```

Equivalently,

```text
delta^2Y(H,H)=-4 Tr(H^2).
```

All six centered Hopf projectors have `Tr(T_i^2)=2/3`, span rank five, and
therefore have the same nonzero response `-8/3`.

## 3. The actual smooth heat coefficient

This step must not be confused with the finite moments
`(2640,14880,55920)`.  For the smooth continuum operator,

```text
Tr exp(-t D_g^2)
 ~ (4 pi t)^(-3/2) [A0+t A2+...].
```

For the Hodge Laplacian on `p`-forms, the rank and trace of the Weitzenböck
curvature term in dimension three give:

| degree `p` | rank | `Tr(q_p)/R` | curvature multiplier in `A2` |
|---:|---:|---:|---:|
| 0 | 1 | 0 | `1/6` |
| 1 | 3 | 1 | `-1/2` |
| 2 | 3 | 1 | `-1/2` |
| 3 | 1 | 0 | `1/6` |

Thus the **ordinary** trace on the full exterior algebra has

```text
A0 = 8 Vol,
A2 = -(2/3) integral R dVol.
```

The graded supertrace coefficient cancels exactly.  This distinction is
load-bearing: an index supertrace has no such metric response, whereas the
ordinary trace used by a spectral action does.

On determinant-one Hopf metric paths the volume term is constant.  Combining
the exact multiplier `-2/3` with the curvature Hessian `-4 I_5` gives

```text
Hess(A2)|_Hopf = (8/3) I_5.
```

No relative tensor coefficient or phenomenological sign has been fitted.

## 4. Attack on the positive interpretation

The coefficient `-2/3` is universal for the full de Rham complex in three
dimensions.  It is not special arithmetic selected by `a1=5`.  Any suitable
three-dimensional Riemannian manifold with this operator has the analogous
curvature heat term.

The theory-specific content is narrower:

1. the 600-cell/Hopf construction selected six tensors forming an exact tight
   frame of the required five-dimensional local spin-two type;
2. those tensors lift canonically to homogeneous TT fields on its round
   quaternionic continuum carrier;
3. the continuum limit of the already used geometric operator has a
   nonzero, full-rank curvature response on that exact carrier.

This proves compatibility of the pieces and removes a possible decoupling
obstruction.  It does not prove that the finite theory selects metric
variation or the spectral functional.

There is also a dimension boundary.  This is a three-dimensional Riemannian
curvature functional.  Pure three-dimensional Einstein gravity has no local
propagating graviton, and the five fields tested here are homogeneous global
anisotropies.  A later `3+1` interpretation requires a selected time
direction and kinetic/gauge constraints; neither is supplied by this
calculation.

## 5. What the repository actually selects

The source audit confirms three separate facts.

- The certified 2640-state incidence Kähler--Dirac operator has a fixed
  Euclidean cochain inner product and no continuous metric variable.
- The Whitney operator genuinely derives its adjoint from exact metric mass
  matrices, but those matrices are evaluated on the fixed piecewise-flat
  geometry.
- The finite spectral-action verifier certifies ordinary moments of one fixed
  matrix.  It does not define a variable-metric configuration space or select
  a cutoff function and scale.

Therefore the positive continuum coefficient is not already an action on
the finite configuration space.  Promoting the metric to a field and choosing
the spectral functional remain physical axioms unless a refinement-natural
selection theorem is found.

## 6. Status ledger

| Claim | Status |
|---|---|
| Direct `SU(2)` scalar-curvature formula | **DERIVED** |
| Round normalized curvature Hessian is `-4 I5` on STF | **DERIVED** |
| Full de Rham ordinary heat coefficient is `-(2/3) integral R` | **DERIVED CONTINUUM** |
| Heat-coefficient Hessian is `+(8/3) I5` on Hopf STF | **DERIVED CONTINUUM** |
| Graded/index trace supplies the same response | **DERIVED NEGATIVE** |
| The Hopf carrier is not decoupled from continuum curvature | **STRUCTURAL POSITIVE** |
| The coefficient is special to `a1=5` | **REFUTED** |
| The finite/refined theory already varies the metric | **DERIVED NEGATIVE as currently defined** |
| A spectral functional and normalization are selected | **OPEN** |
| These homogeneous modes are propagating gravitons | **OPEN / not established** |
| Lorentzian gauge constraints and universal source coupling | **OPEN** |
| Newton/Planck scale | **OPEN** |

## 7. Next strict gate

### Subsequent metric-selection audit

The proposed transfer gate has now been attacked.  The exact Whitney target
metric `g_R` and the unit round metric `g_0` are distinct but both pass full
600-cell symmetry, radial naturality, face compatibility, uniform equivalence
and exact Whitney refinement.  Indeed the whole family
`g_u=(1-u)g_R+u g_0` passes.  Therefore those kinematic conditions do not
select a baseline, and norm equivalence does not transfer the round heat
Hessian to the fixed-Regge theory.  See
`hopf_whitney_metric_selection_result.md`.

The remaining strict gate is consequently dynamical rather than merely a
mass-matrix construction:

An already licensed spectral functional would have to select one metric from
the complete admissible family before inspecting the tensor spectrum and
must answer all of the following:

1. Does the functional select the metric baseline and Hopf perturbation
   without an inserted interpolation or endpoint choice?
2. Do the resulting mass matrices and ordinary spectral response converge
   to the continuum response above in a jointly controlled `h,t` regime?
3. Does the construction itself select the cutoff/normalization, rather than
   merely permit them?

Passing the first two would transfer the structural coupling to the finite
tower.  Failure of the third would still leave induced gravity compatible
but unselected.  Lorentzian propagation and sources remain later gates.

## Reproduction

Run only the targeted verifier:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_hopf_kahler_induced_gravity.py
```

Result: `19/19` checks passed.  The full suite was not run under the user's
current instruction.
