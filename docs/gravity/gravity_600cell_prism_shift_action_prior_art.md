# Prior-art gate: Regge action on the global prism-shift family

Date: 2026-08-19

Status: **completed before protocol and before evaluating any shifted
cellular action or Hessian**.

## 1. Exact object and hypotheses

Use the regular 600-cell boundary with all spatial edge lengths equal to
`L>0`.  Its bottom and top copies are identical.  Every corresponding vertex
pair has timelike squared length `-rho`, `rho>0`.

The preceding exact theorem supplies the complete shape-matched equal-scale
prism family.  A potential

```text
phi in C^0(K;R)/constants
```

gives, on every tetrahedron `sigma=(v0,v1,v2,v3)`,

```text
a_sigma,i=phi(vi)-phi(v0),
H_sigma=[G_sigma a_sigma; a_sigma^T -rho].
```

Every `H_sigma` is Lorentzian, and all shared lateral triangular-prism faces
match.  The family has 119 parameters and identical bottom, top and vertical
natural edge lengths.

The proposed action is the ordinary zero-cosmological-constant Lorentzian
Regge action for this flat polytopal complex, including the additive spatial
boundary term.  Dust follows the 120 vertical edges with their already fixed
conserved masses, so its action depends only on `sqrt(rho)` and is constant
on the present `phi` family.

This mission is local at the static point `phi=constant`.  It does not vary
`L`, `rho`, the dust masses or either spatial boundary.

## 2. Cellular Regge action is known

Tsuda and Fujiwara formulate the Collins--Williams four-polytopal Regge
action directly in terms of the two codimension-two hinge species: spatial
polygons and lateral trapezoids.  Their Eq. (13) is

```text
sum_hinge area * deficit/exterior angle
```

plus an optional cosmological-volume term:

- R. Tsuda and T. Fujiwara, *Oscillating 4-Polytopal Universe in Regge
  Calculus*, PTEP 2021, 083E01, arXiv:`2011.04120`, DOI
  `10.1093/ptep/ptab079`, Sections 2--3.

The repository has already reconciled this cellular action, its Lorentzian
branch and its boundary term with both homogeneous staircase actions.  The
present extension changes the shapes of the flat prisms rather than adding a
new action coefficient.

General flat-polytope Regge data require full shape matching, which the
previous mission now supplies:

- P. Donà, M. Fanizza, G. Sarno and S. Speziale, *SU(2) graph invariants,
  Regge actions and polytopes*, arXiv:`1708.01727`, Section 5.

## 3. Schläfli reduction and boundary term

For a flat four-cell, the Schläfli identity removes the angle variations from
the first variation of the Regge action.  With the complete boundary term,
the remaining boundary angle variations cancel as well.  The discrete
one-step action then generates the canonical pre/post data:

- B. Dittrich and P. A. Hoehn, *From covariant to canonical formulations of
  discrete gravity*, arXiv:`0912.1817`, especially the Regge action and
  boundary variation in Section 2;
- B. Dittrich and P. A. Hoehn, *Canonical simplicial gravity*,
  arXiv:`1108.1974`.

The repository independently checked the same complex Lorentzian Schläfli
identity simplex by simplex before using it in earlier Hessians.  That is a
control on conventions, not permission to omit a new polytopal branch test.

For an equal-scale translated prism, top and bottom facets are parallel and
isometric.  Their oriented exterior angles are opposite.  Since their
spatial triangle areas are equal, the total top-plus-bottom boundary action
cancels throughout the shift family.  This claim must be tested directly in
the registered implementation.

## 4. Prediction frozen before evaluation

For an oriented spatial edge `e=(u,v)`, write

```text
x_e=phi(v)-phi(u).
```

Its lateral quadrilateral is a Lorentzian parallelogram spanned by the
spatial edge and the cell translation.  Its real area magnitude is

```text
A_e(phi)=sqrt(rho L^2+x_e^2).
```

At `phi=constant`, every cell angle around `e` is `acos(1/3)`, so

```text
epsilon=2*pi-5*acos(1/3)>0.
```

The complete Schläfli-reduced first variation is predicted to be

```text
dS=sum_e epsilon_e dA_e.
```

Since `dA_e=0` at `x_e=0`, the static shift is stationary.  Its second
variation is predicted exactly as

```text
delta^2 S
 = epsilon/(L*sqrt(rho))
   * sum_(u,v in edges) (delta phi(v)-delta phi(u))^2.
```

Thus, in the standard Euclidean coordinate inner product on vertex
potentials,

```text
H_phi = kappa * Delta_0,
kappa = epsilon/(L*sqrt(rho)),
Delta_0 = 12 I-A_600.
```

This formula is a preregistered analytic prediction.  It has not been
inferred from a shifted action evaluation.

The exact unscaled Laplacian spectrum is already certified independently:

| eigenvalue | multiplicity |
|---:|---:|
| `0` | 1 |
| `12-6*varphi` | 4 |
| `12-4*varphi` | 9 |
| `9` | 16 |
| `12` | 25 |
| `14` | 36 |
| `8+4*varphi` | 9 |
| `15` | 16 |
| `6+6*varphi` | 4 |

where `varphi=(1+sqrt(5))/2`.  If the prediction is correct, the constant
potential is the sole null and all 119 quotient modes are positive because
`epsilon>0`.

## 5. Framing attack

Several stronger interpretations would be false even if the matrix identity
passes.

1. A graph Laplacian is expected whenever an action is a sum of equal
   edge-difference squares.  The matrix form is structurally natural; it is
   not on its own exotic new physics.
2. The variables `phi` have dimensions of length squared because they are
   mixed Gram entries, not scalar matter fields.
3. These 119 longitudinal modes are not three arbitrary shift functions per
   vertex and are not tensor gravitons.
4. A positive Hessian means that the static action selects zero shift locally
   within this family.  It does not establish exact diffeomorphism gauge
   symmetry; indeed it would refute exact null-gauge behavior here.
5. A fixed coarse graph spectrum is not a dispersion relation.  Refinement
   and a time recurrence would still be required.

Flat-background vertex-displacement gauge modes and their possible breaking
on curved Regge backgrounds are known:

- P. A. Hoehn, arXiv:`1411.5672`, DOI
  `10.1103/PhysRevD.91.124034`;
- B. Bahr and B. Dittrich, arXiv:`0905.1670`, DOI
  `10.1088/0264-9381/26/22/225011`.

The 600-cell spatial slice is intrinsically curved, so a nonzero stiffness
would be compatible with the known pseudo-constraint mechanism.

## 6. KNOWN / CONTROL / OPEN

### KNOWN

- cellular Regge action on flat four-polytopal frusta;
- Schläfli reduction of its first variation with a complete boundary term;
- the exact 600-cell vertex Laplacian and spectrum;
- flat-background Regge lapse/shift gauge generators and their possible
  curvature breaking.

### CONTROL

- the static `phi=constant` action must reproduce
  `720 L sqrt(rho) epsilon` before dust;
- all 600 cell metrics must have signature `(3,1)` near the origin;
- every lateral hinge is incident on exactly five cells;
- top and bottom boundary actions must cancel directly;
- finite differences of the complete area-angle action must agree with the
  Schläfli-reduced gradient before a Hessian is interpreted;
- the constant potential must be an exact null by construction.

### OPEN

- the displayed exact Hessian identity;
- whether any nonconstant mode is null after the full boundary action;
- whether the 119 modes are pseudo-constraints, multipliers or physical
  internal variables in the complete canonical theory;
- the corresponding result on the projected rank-edgewise refinement;
- any time propagation, tensor sector, limiting speed or continuum limit;
- external novelty of this exact 600-cell application.

## 7. Decision boundary

If the direct action, boundary or Schläfli controls fail, the Hessian verdict
is **OPEN** and no graph-Laplacian claim may be made.

If the Hessian is zero on the 119-dimensional quotient, these directions
survive as candidates for exact discrete gauge.

If it is exactly `kappa*Delta_0`, the action locally selects constant `phi`
and the modes are stiff rather than null.  This is a derived
gravity-to-spectral-operator bridge, but only in the longitudinal
shift-potential sector.

If a different nonzero operator appears, that operator and the failed
prediction are the result; no coefficients may be fitted afterward.

