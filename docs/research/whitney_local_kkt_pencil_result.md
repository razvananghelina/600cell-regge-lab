# Exact local Whitney spectral pencil, but no unitary tick

Date: 2026-08-11

Preregistration commit: `ba6035e`

Targeted verifier:
`reproducible/verify_whitney_local_kkt_pencil.py`

Targeted result: **13/13 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Main result

Let (M_{\rm loc}) and (A_{\rm loc}) be the exact block-local metric and
symmetric weak Kähler--Dirac matrices on the 9,000 tetrahedron-local Whitney
cochains.  Let (J:\mathbb C^{2640}\to\mathbb C^{9000}) be conforming
assembly and

\[
Q=I-J(J^TJ)^{-1}J^T.
\]

On the multiplier range ​(operatorname{ran}Q), consider

\[
\mathcal K(z)=
\begin{pmatrix}
A_{\rm loc}-zM_{\rm loc} & Q\\
Q&0
\end{pmatrix}.
\]

The verifier certifies coefficientwise, as exact integer matrix identities,

\[
\exists\mu\in\operatorname{ran}Q:\
\mathcal K(z)\binom{Jx}{\mu}=0
\quad\Longleftrightarrow\quad
(A-zM)x=0,
\]

where

\[
A=J^TA_{\rm loc}J,
\qquad
M=J^TM_{\rm loc}J.
\]

> **DERIVED LOCAL SPECTRAL PENCIL:** the exact Whitney generalized spectrum
> can be represented without forming (M^{-1}), using element-local weak
> blocks plus equality constraints.

This is a genuine way around the exploding inverse-polynomial degree for
spectral calculations.  It is not an alternative formula for the same dense
inverse: the inverse never appears in the pencil coefficients.

## Exact constraint census

The occurrence multiplicities are

\[
(20,5,2,1).
\]

After scaling (Q) by 20, the verifier checks exact integer identities

\[
(20Q)^2=20(20Q),
\qquad
(20Q)J=0,
\qquad
J^T(20Q)=0.
\]

Its exact rank is 6,360 and its kernel dimension is 2,640, exactly the
conforming space.  Every nonzero off-diagonal constraint entry couples copies
of the same global simplex only.  The largest block has 20 copies.

The complete constrained carrier has dimension

\[
9000+6360=15360
\]

when the multiplier is restricted to its nonredundant range.  Using a full
9,000-component multiplier would introduce an additional 2,640-dimensional
gauge kernel and must not be counted as physical states.

## Independent spectral calibration

On the boundary of a 4-simplex, a small triangulated 3-sphere with five
tetrahedra, an independent full-row-rank difference basis was used rather
than the projector formula.  Its descriptor pencil has exactly 30 finite
roots, equal to the 30-dimensional assembled Whitney generalized spectrum.

The maximum discrepancy is

\[
6.22\times10^{-15}.
\]

This numerical control follows the exact full-complex proof; it is not the
evidence on which equivalence rests.

## Why this is not yet time evolution

The natural descriptor metric is

\[
\mathcal M=
\begin{pmatrix}
M_{\rm loc}&0\\
0&0
\end{pmatrix}.
\]

It has rank 9,000 and exactly 6,360 null directions on the canonical
multiplier range.  Consequently the system is a constrained
differential-algebraic equation, not an ordinary Schrödinger equation on a
positive-definite Hilbert space.

Eliminating the multiplier and solving for an autonomous derivative returns
the metric projection and hence the same assembled inverse.  Keeping the
multiplier avoids the inverse only because the constraint remains algebraic.

There is also a general obstruction on the literal physical summand.  For a
self-adjoint block Hamiltonian

\[
H=\begin{pmatrix}H_{PP}&B\\B^*&H_{QQ}\end{pmatrix},
\]

exact invariance of the first summand forces (B^*=0), hence (B=0).  A
nontrivial auxiliary coupling can work only with a different embedded
physical subspace, nonautonomous/reduced dynamics, or an approximation.  Its
embedding must be constructed and tested for locality.

> **DERIVED NEGATIVE:** the canonical KKT pencil is not a positive-metric
> unitary Whitney tick.

## Refinement caveat found by hostile audit

The current (Q) is canonical but only **star-local**: it averages all copies
of one simplex.  On the base complex its largest block is 20.  At the first
barycentric refinement an original vertex belongs to 120 fine tetrahedra, so
complete-star averaging is not a bounded-degree microscopic rule.

This does not invalidate spectral equivalence.  It prevents promotion to a
uniformly local refinement mechanism.

There is a sharper candidate left to test: impose pairwise equality only
between tetrahedron copies sharing a codimension-one face that contains the
simplex.  In a 3-manifold the resulting star-dual constraint graphs have
degrees at most 3, 2 and 1 for vertex, edge and triangle copies.  If they are
connected, their kernel is the same conforming space while locality remains
bounded under subdivision.  This is the next gate; it is not assumed here.

## What this advances

The result separates two roles that had been conflated:

- the spectral action needs the generalized eigenvalues, for which a local
  constrained pencil is enough;
- causal unitary evolution needs a positive-metric generator, which remains
  unsolved.

Thus the refinement explosion of (M^{-1}) does not kill local **spectral
geometry**, but it does still obstruct a local **physical tick**.

## Status ledger

- **DERIVED:** exact element/star-local constrained Whitney pencil.
- **DERIVED:** complete equality constraint rank `6360`, kernel `2640`.
- **DERIVED:** exact generalized spectral equivalence for the full 600-cell.
- **DERIVED:** independent 30-root control on ​(partial\Delta^4).
- **DERIVED NEGATIVE:** singular descriptor metric, hence no unitary tick.
- **STRUCTURAL:** interpreting the constrained pencil as the local object
  underlying the spectral action.
- **OPEN:** bounded-degree neighbour constraints under refinement.
- **OPEN:** positive-metric local dilation with a selected embedded physical
  subspace.
- **NOT CLAIMED:** Lorentzian time, inertia, mass or (c).

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_local_kkt_pencil.py
```

Expected result: `13/13`.
