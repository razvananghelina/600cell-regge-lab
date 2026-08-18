# Prior-art and framing gate: generalized-fiber canonical phase transport

Date: 2026-08-18

## Exact object, carrier and hypotheses

The fixed background is the committed three-slab homogeneous Lorentzian
600-cell Regge--dust trajectory.  At the two middle boundary slices, the
residual-certified calculation supplies distinct rank-15 generalized
configuration fibers

```text
E_old, E_shifted subset C^30
```

in each of two symmetry sectors, both schedule parities and all four disclosed
finite-difference schedules.  Their Hermitian projectors are `P_old` and
`P_shifted`.

The fixed canonical second-slab tangent, derived from the discrete Regge
action and the unique seam ordering, is

```text
T_2 = [A B; C D] : (q_old,p_old) -> (q_new,p_new)
```

on the 60-dimensional complexified phase carrier of each one-dimensional
symmetry sector.

This mission asks whether the natural full cotangent phase lift

```text
F_t = E_t direct-sum E_t^*,
Q_t = diag(P_t, conjugate(P_t))
```

is transported by the action-generated tangent:

```text
(I-Q_shifted) T_2 Q_old = 0.
```

No fitted alignment, graph, Riccati solution or selected momentum relation is
allowed in this first test.

## Primary literature

**KNOWN.** Marsden and West derive discrete Legendre transforms and symplectic
evolution directly from a discrete variational principle; this supports using
the action-generated tangent rather than an independently chosen alignment:
<https://authors.library.caltech.edu/records/1h96d-ymc40>.

**KNOWN.** Dittrich and Höhn use Hamilton's principal function to generate
canonical discrete evolution in simplicial gravity, including constraints and
possibly changing phase-space dimensions:
<https://arxiv.org/abs/1108.1974>.

**KNOWN.** Their covariant-to-canonical construction derives the linearized
Regge dynamics from the action and warns that curved discrete backgrounds
generically replace exact gauge constraints by background-dependent
pseudo-constraints:
<https://arxiv.org/abs/0912.1817>.

**KNOWN.** Constrained variational integrators require a separately derived
discrete null-space/constraint reduction; an arbitrary restriction is not a
substitute:
<https://authors.library.caltech.edu/records/cjm35-0e711/latest>.

These sources fix the canonical object and the interpretation of a failed
phase lift.  They do not imply closure of the present generalized fibers.
External novelty remains **OPEN**.

## Canonicity and the substantive limitation

For a configuration subspace `E` in the already fixed unitary sector
coordinates, its full cotangent dual has projector `conjugate(P)`.  Therefore
`Q=diag(P,conjugate(P))` is a basis-independent phase subspace selected by the
configuration fiber and canonical `q,p` pairing.  The unique seam ordering
and `T_2` are already derived; no basis matching remains.

However, closure of the whole `F_t` is a sufficient, not necessary, condition
for mode dynamics.  A discrete solution may occupy a lower-dimensional
Lagrangian graph `p=R q`, or a constraint-reduced phase subspace, inside
`F_t`.  Such a graph is not supplied by the generalized configuration pencil
alone.  If the full lift fails, the honest result is:

```text
the unrestricted cotangent lift is not transported,
```

not “the modes have no dynamics.”

A subsequent test may form the canonical intersection

```text
F_old intersect T_2^{-1}(F_shifted)
```

and audit its rank, graph property and symplectic character.  It must be
preregistered separately; choosing a graph after minimizing leakage would be
fitting.

## Relation to the previous negative

The earlier tangent audit applied the same phase construction to the
Euclidean negative-stiffness fiber and resolved leakage in blocks `B,D`.  The
generalized Hermitian-definite fiber is a different, kinetic-dependent
configuration object.  Its success is not implied by recurrence closure, and
its failure is not implied by the Euclidean result.  Reusing the old numbers
would therefore be invalid.

## KNOWN, CONTROL, OPEN

- **DERIVED UPSTREAM:** all 32 generalized configuration projectors have a
  `15+10` split and all old/new schedule pairs resolve a time-dependent
  rotation.
- **DERIVED UPSTREAM:** the exact second-slab tangent is regular, symplectic
  and uses the unique identity seam ordering.
- **CONTROL:** reconstruct both projectors and `T_2` at high precision from
  their action source; no binary archive may be the sole scientific input.
- **CONTROL:** require the four block residuals `A,B,C,D` and the combined
  phase residual to agree.
- **OPEN:** transport of the full generalized cotangent phase lift.
- **OPEN:** the dimension and graph/Lagrangian nature of the canonical
  transported intersection if full closure fails.
- **OPEN:** constraint reduction, longer-time transport and refinement.
- **FORBIDDEN HERE:** fitted phase graphs, polar/Procrustes transport, reduced
  roots, dispersion, mass, `c`, graviton or particle-inertia claims.
