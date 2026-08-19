# Adversarial protocol: prism-shift Regge Hessian

Date: 2026-08-19

Primary result commit: `1cfb400`.

Primary artifact SHA-256:

```text
63c9fe41ea4b4de2457f1308a91689786e3871d09ffe8be9008912300e6a4260
```

The primary verifier reports

```text
H_phi = [2*pi-5*acos(1/3)]/(L*sqrt(rho)) * Delta_0.
```

This audit attacks that identity without reusing the primary full-action
evaluator, its finite differences or the Schläfli-reduced gradient.

## 1. Exact local differentiation, no Schläfli assumption

Work on one regular edge-one tetrahedral prism with symbolic independent
local potentials

```text
(p0,p1,p2,p3)=(0,x1,x2,x3),
rho=1.
```

Build its inverse metric directly from

```text
G_ii=1, G_ij=1/2,
H=[G,a;a^T,-1].
```

For all six spatial edges compute the lateral dihedral cosine from the two
opposite facet conormals.  Do not call the primary angle function.  Differentiate
each cosine exactly through second order at the origin.

At the origin all six cosines must equal `1/3`.  Their first derivatives
must vanish.  Convert each cosine Hessian to the corresponding angle Hessian
using the exact derivative of `acos` at `1/3`, then require

```text
sum_over_6_edges Hessian(theta_edge)=0
```

entrywise in `Q(sqrt(2))`.  This is a direct symbolic check of the local
angle cancellation; the Schläfli identity is the conclusion, not an input.

## 2. Direct local wedge Hessian

For every local edge `(i,j)` use

```text
A_ij=sqrt(1+(p_j-p_i)^2).
```

The local wedge functional is

```text
W=sum_(i,j) A_ij theta_ij.
```

Using the independently computed angle Hessians, assemble `Hess(W)` at the
origin.  Require exact equality

```text
Hess(W)=acos(1/3)*K4_reduced,
```

where `K4_reduced` is the `3 x 3` quadratic form obtained from
`sum_(i,j)(delta p_j-delta p_i)^2` after `p0=0`.

No finite-difference tolerance is involved in this decisive check.

## 3. Independent global assembly

Reconstruct the 600-cell edge and tetrahedron incidences.  Assemble the
global Hessian directly as

```text
2*pi * sum_global_edge_area_Hessians
-sum_600_cells local_wedge_Hessian.
```

Do not insert the number five by hand.  The cell assembly must discover it
from the six local edges of every tetrahedron.  Compare the resulting
integer/symbolic matrix entrywise with

```text
[2*pi-5*acos(1/3)]*(12I-A600).
```

## 4. Boundary audit

For every lateral facet, bottom and top conormals differ by sign.  Verify
algebraically at the cosine level that

```text
cos(theta_top)=-cos(theta_bottom).
```

On the branch connected to the static anchor this implies

```text
theta_top=pi-theta_bottom.
```

With identical bottom/top triangle areas and exactly two incident cells per
triangle, the full boundary Hessian must vanish.  Derive this cancellation
symbolically; do not set the boundary contribution to zero by declaration.

## 5. Negative controls

1. Replace the discovered five-cell edge incidence by four in a shadow
   assembly.  It must produce

```text
[2*pi-4*acos(1/3)]*Delta_0,
```

and fail the primary coefficient.
2. Delete the angle-Hessian contribution before local summation.  At least
   one local symbolic residual must become nonzero before the exact six-angle
   cancellation is applied.  This guards against an angle routine returning
   constant `1/3` for every potential.
3. The constant global potential must remain an exact null, while a one-
   vertex potential must have strictly positive quadratic form.

## 6. Spectral checksum

Instead of numerical diagonalization, evaluate the assembled Hessian on the
already certified exact adjacency spectral polynomial

```text
prod_lambda (Delta_0-lambda I)=0
```

using the nine frozen eigenvalues.  At minimum, require exact trace and
trace-square agreement with the multiplicity table and rank `119` from the
known connected graph.

## 7. Verdict

Return

```text
SHIFT_LAPLACIAN_HESSIAN_CORROBORATED
```

only if the local symbolic derivatives, boundary audit, direct global
assembly, negative controls and spectral checks all pass.

Return `PRIMARY_HESSIAN_REFUTED` on any exact disagreement.  There is no
precision-open category for the decisive symbolic matrix identities.

A passing audit remains restricted to the longitudinal shift-potential
sector.  It does not establish gravitons, propagation or an observable
speed.

Only this targeted audit and static guards may run; the full suite is
excluded.

