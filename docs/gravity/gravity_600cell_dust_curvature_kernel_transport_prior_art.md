# Prior-art gate: transport of the curvature-kernel line between two slabs

Date: 2026-08-17

## Exact object, operator, carrier and hypotheses

Use the two already accepted, canonically glued, fixed-mass homogeneous dust
slabs on the same 600-cell staircase carrier:

```text
slab 1: (a_0,r_0) -> (a_1,r_1),
slab 2: (a_1,r_1) -> (a_2,r_2).
```

Here `a_n=log(L_n/L_0)` and `r_n=log(rho_n/rho_0)`.  The shared boundary
momenta satisfy the complete committed canonical seam in all thirty quotient
components.  The total dust mass, Lorentzian angle branch, two schedule
parities and logarithmic canonical coordinates are unchanged.

For each slab `n`, reconstruct independently:

```text
F_n : C^60 -> C^160,
K_n = ker F_n,
T_n : C^60 -> C^60.
```

`F_n` is the derivative of all internal triangle-deficit orbit values after
solving that slab's linearized internal/final pre-Legendre equations.  `T_n`
is the action-generated canonical boundary tangent.  The only comparison in
this mission is made at the common intermediate phase space:

```text
T_1 K_1  ?=  K_2.
```

The final boundary of slab 1 and old boundary of slab 2 are identified by the
literal layer-shift orbit map.  No fitted phase-space transformation or
eigenbasis is allowed.  Lines are compared projectively, so their arbitrary
overall normalization is irrelevant.

## Primary prior art

- Dittrich and Hoehn formulate canonical simplicial evolution using the
  action as Hamilton's principal function.  Pre/post constraints and free
  data can change with the evolution move and may be fixed by later moves:
  <https://arxiv.org/abs/1108.1974>.
- Their covariant-to-canonical analysis shows exact propagation of
  vertex-displacement constraints for flat-background linearized Regge
  calculus, while higher-order curvature produces background-dependent
  pseudo-constraints: <https://arxiv.org/abs/0912.1817>.
- Hoehn explicitly tracks vertex-displacement generators and lattice
  curvature observables through Pachner moves on flat backgrounds:
  <https://arxiv.org/abs/1411.5672>.
- Bahr and Dittrich show that curved Regge solutions generically lack exact
  discrete gauge symmetry and replace constraints by pseudo-constraints:
  <https://arxiv.org/abs/0905.1670>.
- Gambini and Pullin obtain a constraint-free canonical transformation in a
  different consistent-discretization framework:
  <https://arxiv.org/abs/gr-qc/0511096>.

Thus transport/preservation of canonical constraints is **KNOWN** in the
appropriate discrete frameworks, and its failure on curved Regge backgrounds
is also expected in general.  No located primary source defines the present
160-orbit curvature-response kernels or tests their transport on this
fixed-mass 600-cell dust trajectory.  External novelty is **OPEN**.

## KNOWN / CONTROL / OPEN

- **KNOWN:** the two committed slabs share one canonically matched boundary.
- **KNOWN:** an action-generated tangent maps the first slab's initial phase
  tangent space into its final phase tangent space.
- **CONTROL:** `K_1` is a unique homogeneous line and the complete homogeneous
  plane is invariant under `T_1`.
- **CONTROL:** `K_1` is not a fixed eigenline: `T_1 K_1 != K_1` is already a
  resolved result.
- **CONTROL:** both slab backgrounds solve all internal equations and remain
  on the same certified Lorentzian branch.
- **OPEN:** whether `F_2` again has one-dimensional kernel.
- **OPEN:** whether that kernel lies in the same geometry-fixed homogeneous
  phase plane.
- **OPEN:** whether `T_1 K_1=K_2` after the literal seam identification.
- **OPEN:** schedule robustness, later-tick recurrence and refinement.

## Framing attack

`K_n` is the kernel of an observable response: it is a line of *tangent
vectors*.  A Hamiltonian constraint is a function/covector whose zero set
restricts phase space.  Therefore even exact transport of `K_n` would not by
itself prove a first-class constraint, gauge symmetry or lapse.  It would
prove only a background-dependent invariant line distribution selected by
internal curvature blindness.

Conversely, failed transport would reject this specific candidate
distribution, not every possible Regge pseudo-constraint.  The curvature
response is a physically motivated diagnostic but is not the complete
linearized constraint operator.

There is a second hazard: comparing `T_1 K_1` with `K_1` silently treats an
evolving background as stationary.  The accepted trajectory changes both
scale and lapse, so `K_2` is the correct target at the shared boundary.  This
mission repairs that framing without weakening the falsification criterion.

Only the new targeted verifier will be run.  The full suite remains excluded
by instruction.
