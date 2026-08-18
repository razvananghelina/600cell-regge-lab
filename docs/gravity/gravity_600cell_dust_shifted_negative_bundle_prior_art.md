# Prior-art gate: shifted negative shape bundle

Date: 2026-08-18

## Exact proposed object

The fixed carrier is the `720`-dimensional logarithmic signed-squared edge
space of one labelled regular 600-cell slice.  Consecutive homothetic slices
use the literal edge-label identification already certified by the two- and
three-slab tangent verifiers.

At each of the first two nonstationary middle slices, and separately in both
staircase schedules, all seven frozen binary-tetrahedral sectors and all four
derivative variants, the action fixes

```text
M_j = (K^-_j + K^+_j)/2,
V_j = K^-_j + K^0_j + K^+_j,
K_j = im C,
S_j = ker(C* H(M_j)),
A_j = - W_j* H(V_j) W_j,
```

where `C` is the literal unsigned vertex--edge incidence, `H(X)` is the
Hermitian part and the columns of `W_j` span `S_j`.  The proposed object is
the ambient spectral projector

```text
P^-_j = W_j 1_{(-infinity,0)}(A_j) W_j*.
```

The new calculation first asks blindly for the complete shifted inertia.
Only afterwards does it compare the ranks, sector support and projector
distance `||P^-_2-P^-_1||_2`.  A finite two-step recurrence restricted to the
negative carrier is permitted only if one common carrier is certified within
the propagated numerical enclosure.  No Procrustes, polar-decomposition or
overlap-based rotation is allowed to manufacture a temporal connection.

## Primary literature checked

1. Dittrich and Höhn, [*Canonical simplicial gravity*](https://arxiv.org/abs/1108.1974),
   derive discrete canonical evolution from Hamilton's principal function and
   allow evolving phase spaces and constraints.  This supports using the
   action's mixed Hessian blocks as the canonical temporal data; it does not
   identify a negative spectral bundle for this 600-cell dust background.
2. Dittrich and Höhn, [*From covariant to canonical formulations of discrete
   gravity*](https://arxiv.org/abs/0912.1817), derive linearized Regge dynamics
   and explain how curvature and nonlinear order turn exact constraints into
   background-dependent pseudo-constraints.  This makes persistence on a
   curved nonstationary background a calculation rather than a symmetry
   theorem.
3. Bahr and Dittrich, [*(Broken) Gauge Symmetries and Constraints in Regge
   Calculus*](https://arxiv.org/abs/0905.1670), show that curved Regge solutions
   generically lack exact discrete gauge symmetry.  Consequently a change of
   rank or sector support cannot automatically be dismissed as gauge.
4. Dittrich, Freidel and Speziale, [*Linearized dynamics from the 4-simplex
   Regge action*](https://arxiv.org/abs/0707.4513), compute an action Hessian
   and relate its null structure to linearized gravity in a different carrier
   and background.  It is a control for interpreting Hessians, not a result
   about the present negative modes.
5. Davis and Kahan, [*The Rotation of Eigenvectors by a Perturbation. III*](https://doi.org/10.1137/0707001),
   give basis-independent eigenspace perturbation bounds in terms of spectral
   gaps.  This supports comparing projectors rather than individual
   eigenvectors.  It does not supply a physical connection between unequal
   fibers.

## Gate ledger

- **KNOWN:** an action Hessian/principal function canonically determines the
  linearized discrete evolution; invariant spectral subspaces are compared
  basis-independently by their orthogonal projectors and eigengaps.
- **CONTROL:** the repository already certifies the first centered Jacobi
  operator, its action-relative conformal/shape split, and `30` resolved
  negative shape directions supported in two frozen sectors.  It separately
  certifies the shifted centered Jacobi operator, but no shifted shape
  stiffness or negative projector has yet been evaluated.
- **OPEN:** whether the shifted action has the same negative rank and sector
  support; whether the two ambient projectors agree within error or rotate;
  whether the action itself selects a temporal connection if they rotate.
- **OPEN external novelty:** the search found no primary source computing this
  complete third-slab, full-600-cell dust-Regge negative bundle.  A literature
  search is not proof of novelty.

## Framing attack

The statement "the same 30 modes persist" has three inequivalent meanings:

1. equal rank and sector labels only;
2. a resolved but rotating rank-30 spectral bundle;
3. one common ambient subspace under the literal edge identification.

Only (3) permits a reduced two-step product without adding transport data.
Case (2) remains mathematically structured, but an arbitrary unitary matching
the fibers would be fitting.  Case (1) alone is only a multiplicity pattern.

No result below may call these modes gravitons, physical instabilities or
propagating polarizations before the constraint quotient, physical clock and
continuum/refinement controls exist.
