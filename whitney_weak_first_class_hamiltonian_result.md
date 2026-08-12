# Weak first-classness does not restore a local Whitney Hamiltonian

Date: 2026-08-12

Protocol commits: `50a5f79`, corrected negative control `bb02023`.

## Outcome

> **DERIVED SCOPED NEGATIVE.**  Replacing strong invariance `HK=0` by the
> correct weak first-class condition does not rescue a strictly
> endpoint-local quadratic Hamiltonian on the fixed reducible Whitney--BRST
> carrier.

The obstruction is smaller and sharper than the earlier full inverse
calculation.  It is already forced by the block mapping degree-two auxiliary
constraints into degree-three occurrence forms.  That block is unique, is
independent of every remaining Hamiltonian block, and contains an exact
nonzero entry at the largest possible tetrahedron-to-triangle distance on the
base 600-cell.

This does **not** remove the already derived Whitney/Galerkin evolution on the
assembled quotient.  It closes one proposed *strictly local first-class
realisation* of that evolution.  Time, a causal cone, mass, `c` and Planck
units remain unconstructed.

## 1. Complete hypotheses

The statement assumes all of the following:

- the canonical duplicated tetrahedral Whitney carrier;
- the exact block-local consistent Whitney metric `M`;
- the already derived local weak Kaehler--Dirac form `A=A*`;
- every canonical neighbour copy constraint `C`, including its
  reducibilities;
- the fixed auxiliary Poisson tensor `+iG`, where
  `G=C M^-1 C*`;
- a linear Hermitian quadratic Hamiltonian with top-left block fixed to `A`;
- endpoint-star locality: a tetrahedron row may couple to a triangle
  constraint only when it is one of that triangle's two parent tetrahedra.

No phenomenological target or fitted coefficient enters the calculation.

## 2. The framing correction was real

For

\[
 L=(C,I),\qquad
 K=\binom{-M^{-1}C^*}{G},\qquad
 H=\begin{pmatrix}A&B\\B^*&D\end{pmatrix},
\]

strong invariance `HK=0` is more than first-class dynamics requires.  The
correct condition is

\[
 HK\in\operatorname{im}L^*.
\]

The constraint surface is parametrised by

\[
 S=\binom{I}{-C},\qquad LS=0.
\]

Because `S` has full column rank and `L` has full row rank,
`im S=ker L`.  Thus the weak condition is exactly

\[
 S^*HK=0,
\]

or, with `Q=M^-1 C*`,

\[
 -AQ+BG+C^*B^*Q-C^*DG=0.                 \tag{1}
\]

The exact control uses the full canonical redundant constraint matrix:
`L` is `70 x 145`, `S` is `145 x 75`, and the identities and ranks hold
exactly.  The audit therefore does not quietly fall back to strong
invariance or to an independent constraint basis.

## 3. Why the top block cannot be adjusted

On a closed triangulated 3-manifold each triangle has exactly two
tetrahedral occurrences.  Consequently the degree-two copy constraints
`C2` have one disjoint difference row per triangle and full row rank.  Hence

\[
 G_2=C_2M_2^{-1}C_2^*>0
\]

is invertible.  Degree-three forms have one occurrence per tetrahedron, so
there are no degree-three copy constraints: `C3=0`.

Take the `(u3,gauge2)` block of (1).  Both terms left-multiplied by `C*`
vanish there.  The result is

\[
 B_{32}G_2=A_{32}M_2^{-1}C_2^*,
\]

and therefore the unique solution is

\[
 B_{32}=A_{32}M_2^{-1}C_2^*G_2^{-1}.     \tag{2}
\]

Neither `D` nor any other block of `B` can change (2).  Thus weak
first-classness really was tested, but its additional freedom is absent in
the load-bearing top-degree block.

## 4. Exact reduction

For the regular tetrahedral Whitney metric the verifier derives, rather than
assumes,

\[
 M_2^{-1}=6I+s s^*,\qquad s=(1,-1,1,-1)^T,
\]

and

\[
 A_{32}M_2^{-1}=\frac{15}{4}s^*.
\]

Let `N=C2 diag(s,...,s)`.  Since the rows of `C2` have disjoint occurrence
support,

\[
 G_2=12I+NN^*.
\]

Equation (2) reduces from the triangle space to the tetrahedron space:

\[
 \boxed{B_{32}=\frac{15}{4}(12I+N^*N)^{-1}N^*.}       \tag{3}
\]

This formula has no Schur coefficient to choose.  Changes of orientation
only multiply rows or columns by signs and therefore cannot change its zero
pattern.

## 5. Results

| control | exact/numerical | result |
|---|---:|---|
| boundary of a 4-simplex | exact rational | 20 incident and 30 remote nonzeros |
| diagonal row-sum mass lumping | exact rational | zero remote nonzeros |
| base 600-cell, relative threshold `1e-10` | numerical, residual `6.79e-16` | 436800 remote entries, reaching dual distance 8 |
| base 600-cell, lexicographic maximum-distance pair | exact modulo `1000003` | nonzero at endpoint distance 14 |

The mass-lumped negative control matters: it returns

\[
 B_{32}^{\rm lump}=\frac{3}{16}N^*,

which is exactly endpoint-local.  The locality census therefore does not
declare every construction nonlocal by definition.  The remote support is a
consequence of the consistent Whitney metric.

For the base 600-cell, the exact certificate selects the lexicographically
first tetrahedron/triangle pair at maximum endpoint distance using topology
alone.  The pair is row `0`, triangle `1196`, whose endpoints are tetrahedra
`579` and `599`.  The distance is `14`.  The reduced integer matrix has full
rank `600` modulo `1000003`, and the corresponding value of `B32` is
`979391` modulo that prime.  Hence the rational entry cannot be zero.  This
is stronger than interpreting a very small floating-point tail.

## 6. Attack on the conclusion

What is closed:

- a linear Hermitian quadratic completion;
- on this fixed canonical auxiliary carrier;
- with the exact Whitney `A` kept as its physical top block;
- under strict element/endpoint-star locality.

What is not closed:

- nonlinear Hamiltonians;
- a differently derived auxiliary carrier or Poisson structure;
- a Lorentzian construction with another independently justified locality
  notion;
- the already valid, but spatially assembled, Whitney/Galerkin evolution.

Allowing an arbitrary wider stencil would evade the binary endpoint test,
but would not solve the problem: equation (3) selects support at the maximum
base dual distance.  Calling such a block “local” would change the physical
criterion rather than satisfy it.

## Status ledger

- **DERIVED:** weak first-classness is `HK in im L*`, not necessarily
  `HK=0`.
- **DERIVED:** equation (1) is the exact weak Hamiltonian condition.
- **DERIVED:** `C3=0` and invertibility of `G2` uniquely freeze `B32`.
- **DERIVED:** formula (3) follows from the exact regular Whitney metric.
- **DERIVED CONTROL:** the exact 4-simplex boundary has 30 remote entries.
- **DERIVED NEGATIVE CONTROL:** row-sum lumping has zero remote entries.
- **DERIVED BASE CERTIFICATE:** an exact nonzero reaches maximum endpoint
  distance 14 on the 600-cell.
- **DERIVED SCOPED NEGATIVE:** weak first-classness does not yield a strict
  local Hamiltonian in the stated class.
- **STILL DERIVED:** the assembled Whitney equation gives nontrivial unitary
  mathematical evolution.
- **OPEN:** a selected local physical Hamiltonian from a genuinely different
  construction.
- **NOT CLAIMED:** physical time, causality, inertia, mass, `c`, `hbar`,
  Newton's `G` or Planck units.

## Reproduction

Only the verifier relevant to the current mission was run, following the
user's instruction not to run the full suite:

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_weak_first_class_hamiltonian.py
```

Result: `15/15` checks passed.
