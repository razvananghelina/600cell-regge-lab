# Result: relative prism shifts are canonically fixed, not freely propagated

Date: 2026-08-19

## Headline

**DERIVED COMPUTATIONAL NEGATIVE FOR THE FREE-SHIFT ROUTE.**  At the first
accepted unequal-scale dust--Regge slab, with old boundary geometry and
incoming momentum fixed, the complete canonical equations eliminate all
`119` relative prism-shift/strut directions.  An independent geometric audit
confirms that the exact polytopal extension of these directions changes only
the `720` internal cross diagonals and the `120` internal poles; it has zero
boundary support in both staircase schedules.

Thus the previously derived static operator

```text
H_phi = [2*pi-5*acos(1/3)]/(L*sqrt(rho)) * Delta_600
```

is real but is not the spatial half of a freely propagating scalar field.
The most defensible present reading of `phi` is an auxiliary or
pseudo-constraint-like longitudinal variable.  That reading remains
**STRUCTURAL**, because the calculation is on a curved finite background
where exact discrete gauge symmetry need not survive.

This closes one shortcut; it does not block the gravity programme.  The
action-selected `600`-position shape subsystem remains the carrier on which
physical curvature/tensor dynamics must be sought.

## Provenance ledger

| stage | commit |
|---|---|
| prior-art and framing gate | `d90a44b` |
| primary composition protocol | `e411f93` |
| registered primary verifier | `7df66ba` |
| frozen primary `13/13` artifact | `78fa42c` |
| adversarial geometry protocol | `d01217b` |
| registered adversarial verifier | `2e28678` |
| frozen adversarial `14/14` artifact | `33d909b` |

Artifacts:

```text
reproducible/gravity_600cell_prism_shift_canonical_elimination.json
SHA-256 b9e31d56670c397232937ae5f2e7e002632cc715807c221b2b98b47e20dde332

reproducible/gravity_600cell_prism_shift_canonical_elimination_adversarial.json
SHA-256 1821e02966c6b4f97ee499c00c6366f5a63a5b688852bdd1056c70961f2e708c
```

Only the two mission-specific verifiers and static guards are run in the
final check.  The full suite is excluded by user instruction.

## Complete hypotheses

The result requires all of the following:

1. the fixed regular 600-cell and its two frozen staircase triangulations;
2. the accepted nonstatic, fixed-total-mass dust--Regge background;
3. the complete `2280` edge carrier, split into `720` old-boundary, `840`
   internal and `720` new-boundary edges;
4. the frozen complete pre-Legendre Jacobian and its `1440` strong plus `120`
   weak canonical partition;
5. invertibility, with certified determinant and singular-value errors, of
   every strong block `A` and every pole Schur block
   `S=D-C A^-1 B` in all seven binary-tetrahedral sectors and both schedules;
6. the exact relative-pole map `R e_i=e_i-e_119`;
7. unequal scale ratio `q != 1` for the polytopal reconstruction;
8. homogeneous linearized equations at fixed old geometry and incoming
   momentum.

It is not a theorem for arbitrary Regge backgrounds, nonlinear finite
perturbations, a continuum limit or sourced boundary data.

## Primary canonical implication

The frozen complete canonical linearization has the block equations

```text
A x + B z = 0,
C x + D z = 0,
```

where `z` contains the `120` pole variables and `x` contains the `1440`
strong internal/new-boundary variables.  Since every `A` and every Schur
complement

```text
S = D-C A^-1 B
```

is invertible,

```text
x=-A^-1 B z,
S z=0,
z=0.
```

The exact relative map `R` has rank `119` over both `F_101` and
`F_1000003`; its all-ones complement supplies the collective pole/lapse.
Consequently `S R` is injective and no relative homogeneous direction
survives.

The smallest independently recomputed Schur singular value is

```text
4.24456181727093006e-9,
```

with minimum margin about `4.08e12` over the stored `100*epsilon` boundary.
The rank decision is therefore not marginal inside the declared error model.

## Why the polytopal diagonal graph does not change the Schur result

A geometric realization generally changes strong cross diagonals together
with the poles.  Write that graph as `x=Gz+y`.  The transformed blocks are

```text
B' = A G+B,
D' = C G+D.
```

Exactly,

```text
D'-C A^-1 B'
 = (C G+D)-C A^-1(A G+B)
 = D-C A^-1 B.
```

The primary verifier checks this cancellation on a nonsymmetric rational
matrix control.  Therefore the effective weak equation is invariant under a
legitimate strong-coordinate graph.  The adversarial audit was needed to
establish that the actual 600-cell polytopal map really is such an internal
graph rather than a hidden boundary variation.

## Independent geometric audit

The audit rebuilt the binary-icosahedral multiplication, the five
binary-tetrahedral right cosets, all `600` spatial tetrahedra and both
staircase slabs without importing the Regge slab constructor.  Each schedule
gave exactly

```text
2400 4-simplexes,
720 old + 720 new + 720 cross + 120 pole edges = 2280.
```

For every spatial edge, the schedule selects exactly the diagonal from the
bottom endpoint of larger color to the top copy of the endpoint of smaller
color.  The resulting `720` cross diagonals and `120` poles are all internal;
their intersection with both boundary edge sets is empty.

For a centered regular tetrahedron, direct exact coordinate differentiation
in two nonsymmetric rational controls rederived

```text
delta D_ij = [q delta r_j-delta r_i]/(q-1),

delta log D_ij
  = rho*(z_i-q z_j)/[(q-1)D_ij].
```

All direct residuals are exactly zero.  Reversing the endpoint convention
produces nonzero residuals in all `24/24` directed controls.

For `q=3/2` and `q=5/3`, in both schedules, the `720 x 120` diagonal graph
has rank `120` over both `F_101` and `F_1000003`; its pullback to the relative
hyperplane has rank `119`.  Deleting every diagonal incident to one vertex
drops the rank to `119`, so the rank test has demonstrated falsification
power.  On the collective vector the formula reduces exactly to the frozen
geometric lapse coefficient `-rho/D`.

## Status ledger

- **DERIVED:** the static equal-scale `phi` Hessian is a positive multiple of
  the 600-cell graph Laplacian.
- **DERIVED:** unequal homogeneous scale and common strut length forbid a
  nonzero independent local shift.
- **DERIVED COMPUTATIONAL:** nonuniform pole data extend to a complete
  internal polytopal graph of rank `120`, with relative rank `119`.
- **DERIVED COMPUTATIONAL:** the complete fixed-data canonical equations
  eliminate all `119` relative directions and the collective pole direction.
- **STRUCTURAL:** `phi` is auxiliary/pseudo-constraint-like rather than a
  propagating scalar.
- **OPEN:** the nonzero response induced on new-boundary data when old data
  are perturbed.
- **OPEN:** the physical scalar/vector/tensor and constraint decomposition of
  the `600`-position action-selected shape carrier.
- **OPEN:** multi-tick evolution, continuum convergence, dispersion, a
  limiting speed, Planck units and particle masses.

## Relation to prior art

The interpretation is consistent with, but not proved by, the canonical
Regge literature:

- [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974) show that simplicial
  evolution may introduce initially free data which later constraints fix;
- [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817) derive
  background-dependent pseudo-constraints once nonlinear effects break the
  flat-background symmetry;
- [Bahr--Dittrich](https://arxiv.org/abs/0905.1670) explain that curved Regge
  solutions generically lack exact discrete gauge symmetries;
- [Hoehn](https://arxiv.org/abs/1411.5672) separates vertex lapse/shift
  variables from curvature-carrying lattice gravitons in canonical
  linearized Regge calculus around flat backgrounds.

No located source states the present exact 600-cell graph formula or this
complete-sector Schur composition.  Search absence is not novelty proof;
external novelty is **OPEN**.

## Next load-bearing gate

Do not attempt to extract `c`, inertia or a scalar wave equation from `phi`.
The live route is already selected by the action:

```text
720 spatial edge perturbations
= 120 conformal directions + 600 shape directions.
```

The next physical calculation should operate on the `600`-position shape
recurrence and derive a curvature/tensor intertwiner before comparing any
spectrum with continuum harmonics.  In particular, the existing
`30`-position autonomous negative-stiffness subsystem must first be
classified as scalar, vector or tensor (and tested under refinement); its
literal polynomial root count is `15/0/15`, but the project-wide `100x`
safety margin remains open.  This route is harder than the discarded
longitudinal shortcut, but it is the one not eliminated by the canonical
equations.

