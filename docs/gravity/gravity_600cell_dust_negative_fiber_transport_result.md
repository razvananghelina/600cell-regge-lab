# Result: common configuration fiber, naive phase closure refuted

Date: 2026-08-18

## Headline

The old and shifted rank-`15` negative-stiffness projectors are geometrically
the same within certified error in all `16/16` cells.  However, their natural
unweighted phase lift `E + E*` is not invariant under the action-derived
second-slab tangent.

The preregistered outcome is

```text
NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED.
```

Thus the result is stronger than persistent rank but weaker than an autonomous
dynamical sector: there is a common rank-`30` configuration subsystem, while
the naive momentum companion is wrong.

## Provenance ledger

| stage | commit |
|---|---|
| prior-art and canonicity gate | `42dc0e2` |
| frozen transport protocol | `07de221` |
| registered verifier before first execution | `95c86bf` |
| pre-target integration fix (`I`) | `a493255` |
| deterministic result artifact | `da5325b` |

The first execution stopped before any projector comparison because an
imported sector helper lacked its global imaginary unit.  The one-line fix was
committed before rerunning and changed no mathematical definition.

The final verifier source has SHA-256

```text
f462e507500d7f02ecf799f0d4b320e05795216a36a0d10eb908d6dc67b48181
```

and the byte-identical artifact from two executions has SHA-256

```text
d630bf07066f88c35eee5a62a80ec1f43399a95ea882a43528289220c67f4599.
```

Both executions reported `8/8` checks passed.  Only the targeted verifier was
run; the full suite was not run.

## Result A: the configuration fiber is common

All `32` old/shifted midpoint projectors have a separated `15/10` spectral
gap with conservative shape-carrier and Davis--Kahan-type error propagation.
After the unique boundary map,

```text
16/16  COMMON_FIBER_RESOLVED.
```

The observed projector distance is only about `4.35e-4` of its already
conservative error bound.  This is not merely equal dimension: the
configuration subspaces themselves agree within the certified enclosure.

Because the comparison uses the unique seam edge ordering and no optimized
unitary, this result is **GEOMETRIC STRUCTURAL**.

## Result B: the full phase lift is not closed

For the canonical tangent `T_2=[A B; C D]`, the complete `64`-block leakage
census is

| tangent block | zero-consistent | nonzero-resolved | error-unit range |
|---|---:|---:|---:|
| `A` | 16 | 0 | about `3.58e-4` |
| `B` | 0 | 16 | about `914.10` |
| `C` | 16 | 0 | about `3.66e-4` |
| `D` | 0 | 16 | about `655.08` |

The pattern is exact across both parities, both sectors and all four derivative
schedules:

```text
32 LEAKAGE_ZERO_CONSISTENT
32 LEAKAGE_NONZERO_RESOLVED.
```

Position perturbations in the negative fiber are carried back into the common
configuration/cotangent fibers by `A,C`.  Arbitrary momentum perturbations in
the conjugate fiber leak through `B,D`.  Hence the direct sum `E + E*` is not
an invariant phase subsystem.

## Framing correction

The stiffness inertia of `A=-V_S` is a genuine, persistent quadratic-form
property, but its Euclidean negative projector is not by itself the normal-mode
projector of the dynamics.  The kinetic form `B=-M_S` participates in the
equations of motion.  A dynamically meaningful phase fiber may therefore
require one of the following action-selected objects:

1. the generalized Hermitian pencil `(A,B)`;
2. a time-dependent Lagrangian graph `p=R q` propagated by the canonical
   tangent/Riccati relation;
3. a separately derived constraint-reduced phase lift.

Choosing an `R` or an alignment to minimize leakage would be fitting and is
forbidden.  The next test must let the action select it before comparison.

## Status ledger

- **DERIVED COMPUTATIONAL:** sectors `4,5` retain rank `15` negative stiffness
  at the shifted tick under the frozen derivative family.
- **GEOMETRIC STRUCTURAL:** all old/shifted negative configuration projectors
  are common within conservative certified error.
- **DERIVED COMPUTATIONAL NEGATIVE:** the naive complexified cotangent lift
  `E+E*` is not invariant under the canonical second-slab tangent.
- **OPEN:** whether the generalized kinetic--stiffness pencil selects a closed
  fiber.
- **OPEN:** an action-derived Riccati/Lagrangian graph, constraint reduction,
  longer-time transport and refinement.
- **OPEN:** particle inertia, mass, graviton/wave interpretation, dispersion,
  limiting speed and continuum physics.

The negative result does not erase the common configuration fiber.  It says
precisely that “same shape direction” is not yet “same phase-space mode.”
