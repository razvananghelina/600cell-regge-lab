# Protocol: weak first-class completion of the Whitney Hamiltonian

Date: 2026-08-12

This protocol is committed before the definitive computation.  A preliminary
reconnaissance on the boundary of a 4-simplex found that strict endpoint-star
support fails in the five top-form rows; that observation is disclosed here
and is not a blind result.  No particle, mass, coupling, clock, `a1=5`, Planck
or other phenomenological target is used.

## Question and fixed hypotheses

The previous minimal first-class audit imposed the strong condition `H K=0`.
The present audit asks whether that condition was unnecessarily restrictive.

Use the already derived occurrence carrier, exact block-local Whitney metric
`M>0`, local weak Kaehler--Dirac form `A=A*`, canonical complete copy
constraint `C`, and reducible auxiliary Poisson tensor

\[
 G=C M^{-1}C^*.
\]

On `z=(u,eta)` set

\[
 L=(C,I),\qquad \Phi=Lz,\qquad
 K=\binom{-M^{-1}C^*}{G}.
\]

The top-left block of the Hermitian quadratic Hamiltonian is fixed to the
geometrically derived `A`:

\[
 H=\begin{pmatrix}A&B\\B^*&D\end{pmatrix},\qquad D=D^*.
\]

No coefficient in `A`, `C`, `M` or `G` may be changed or fitted.

## Correct weak first-class condition

First-class invariance only requires the Hamiltonian variation to vanish on
the constraint surface:

\[
 HK\in\operatorname{im}L^*.
\]

Since `ker L` is parametrised by

\[
 S=\binom{I}{-C},
\]

this is equivalent to the exact linear identity

\[
 S^*HK=
 -A Q+B G+C^*B^*Q-C^*DG=0,
 \qquad Q=M^{-1}C^*.
\]

The verifier must establish this equivalence and may not substitute the
stronger condition `HK=0`.

## Frozen top-degree obstruction

The complex is a closed triangulated 3-manifold.  In degree two every global
triangle has exactly two tetrahedral occurrences, so the canonical
neighbour constraint `C_2` has one disjoint difference row per triangle and
has full row rank.  Hence

\[
 G_2=C_2M_2^{-1}C_2^*>0
\]

is invertible.  There are no copy constraints in degree three, so `C_3=0`.
Taking the `(u_3,gauge_2)` block of the weak identity removes both terms
left-multiplied by `C^*` and freezes

\[
 B_{32}G_2=A_{32}M_2^{-1}C_2^*,
\]

therefore

\[
 B_{32}=A_{32}M_2^{-1}C_2^*G_2^{-1}.
\]

The audit must prove this uniqueness exactly.  It follows before any support
calculation and means that `D` and every other block of `B` are irrelevant to
the top-degree locality question.

## Frozen locality criterion and controls

For a row belonging to a tetrahedron and a constraint column belonging to a
shared triangle, endpoint-star locality means that the tetrahedron is one of
the two endpoints of that triangle.  A nonzero outside those endpoints is
called remote.  This is a necessary, deliberately generous element-locality
criterion: failure closes strict element-local completion, while passing
would not by itself establish a causal physical tick.

Perform the following target-free checks in order:

1. exact rational boundary-of-a-4-simplex control;
2. independently assembled numerical base 600-cell control, reporting a
   residual tolerance and the maximum dual-graph distance reached by entries
   above a frozen relative threshold of `1e-10`;
3. a known local negative control obtained by replacing each consistent
   degree-two Whitney block by its diagonal row-sum lumping (and rebuilding
   `Q_2,G_2` from that metric); this makes the distinct triangle constraints
   orthogonal.  It checks that the support census is capable of returning
   zero remote entries when appropriate.  The local coboundary and top-form
   mass are unchanged.

The exact control decides nonzero versus zero.  The base computation is a
support-depth audit, not an exact proof of numerical nonzeros.

## Attack on the framing

The conclusion is scoped to:

- the fixed canonical occurrence carrier and auxiliary bracket `+iG`;
- linear Hermitian quadratic Hamiltonians;
- fixed top-left Whitney block `A`;
- strict element/endpoint-star locality.

It is not a theorem against nonlinear Hamiltonians, a different auxiliary
carrier, non-Hermitian evolution, or a separately derived Lorentzian notion
of locality.  Conversely, merely finding some weak completion is not a
physical result unless the fixed top-degree block is local and the physical
Hamiltonian is selected without fitted coefficients.

## Decision boundary

- **DERIVED SCOPED NEGATIVE:** if the unique exact `B_32` has any remote
  endpoint-star entry, weak first-classness does not rescue a strictly local
  Whitney Hamiltonian in the stated class.
- **STRUCTURAL/PATTERN only:** a finite numerical support depth without the
  exact control.
- **OPEN:** if the top block is local, solve the remaining Hermitian linear
  equations for `B,D` without enlarging support or choosing fitted
  coefficients.

No result here derives physical time, causality, inertia, mass, the speed of
light or Planck units.
