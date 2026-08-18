# Preregistration: a local constrained pencil versus a unitary Whitney tick

Date: 2026-08-11

## Motivation and hostile framing

The assembled Whitney generator is (D=M^{-1}A), where (A=MD) is the
symmetric weak operator.  Exact inverse-polynomial depth grows sharply at the
first barycentric refinement.  A standard mixed formulation may avoid
forming (M^{-1}), but it can also merely hide the same global solve in a
Schur complement.

This protocol distinguishes two claims before constructing the block system:

1. **spectral locality:** the exact generalized spectrum is represented by a
   sparse local matrix pencil with algebraic constraints;
2. **unitary dynamical locality:** an ordinary self-adjoint Hamiltonian with a
   positive-definite Hilbert metric generates the same autonomous physical
   evolution.

Passing claim 1 must not be reported as passing claim 2.

## Fixed element carrier

Use the exact 9,000-dimensional direct sum of all 600 regular tetrahedron
Whitney cochain spaces from protocol `162ce61`.  Let

\[
M_{\rm loc}=\bigoplus_t M_t,
\]

and let (A_{\rm loc}) be the corresponding symmetric weak Kähler--Dirac
matrix.  Both are block-local inside one tetrahedron.  Let

\[
J:\mathbb C^{2640}\longrightarrow\mathbb C^{9000}

\]

duplicate each global oriented cochain into all incident tetrahedra.  Define

\[
M=J^TM_{\rm loc}J,
\qquad
A=J^TA_{\rm loc}J.
\]

No fitted coefficient enters these definitions.

## Canonical constraint

The earlier anchor-difference matrix has the correct kernel but depends on a
chosen anchor occurrence.  It is unsuitable as the definition of a canonical
operator.

Instead use the Euclidean equality projector

\[
P_{\rm eq}=J(J^TJ)^{-1}J^T,
\qquad
Q=I-P_{\rm eq}.
\]

Because every degree has uniform occurrence multiplicity `(20,5,2,1)`, (Q)
is a direct sum of complete-star mean-subtraction blocks.  It is selected by
copy equality, is symmetric and satisfies exactly

\[
Q^2=Q,qquad QJ=0,qquad\ker Q=\operatorname{im}J.
\]

It couples only copies belonging to the same global simplex star.  Record the
largest block size and do not call this nearest-tetrahedron locality; it is
**star-local**.

## Frozen KKT pencil

For spectral parameter (z), on
(mathcal H_{\rm loc}\oplus\operatorname{ran}Q) define

\[
\mathcal K(z)=
\begin{pmatrix}
A_{\rm loc}-zM_{\rm loc} & Q\\
Q & 0
\end{pmatrix}.
\]

The multiplier is restricted to ​(operatorname{ran}Q); using the full
9,000-dimensional multiplier space would add the gauge kernel
(ker Q) and must be reported as redundant rather than physical.

The exact equivalence to be certified is:

\[
\exists\mu\in\operatorname{ran}Q:\
\mathcal K(z)\binom{Jx}{\mu}=0
\quad\Longleftrightarrow\quad
(A-zM)x=0.
\]

Forward implication follows by multiplying the first row by (J^T).  For the
reverse implication, the residual
((A_{\rm loc}-zM_{\rm loc})Jx) is Euclidean-orthogonal to
(operatorname{im}J), hence lies in ​(operatorname{ran}Q), and fixes
(mu) uniquely in that range.

## Exact and calibrated checks

1. Rebuild all local Whitney masses and the full 600-cell assembly exactly.
2. Clear uniform multiplicity denominators and certify (Q^2=Q), (QJ=0),
   `rank(Q)=6360`, and `ker(Q)=im(J)` degree by degree.
3. Certify exact symmetry and block locality of (M_{\rm loc}) and
   (A_{\rm loc}), and exact assembly of (M,A).
4. Verify the two directions of the pencil equivalence as complete symbolic
   matrix identities in the coefficients of (z), not at selected
   eigenvalues.
5. On the boundary of a 4-simplex (five tetrahedra, a small triangulated
   3-sphere), construct an independent full-rank constraint basis and compare
   the finite generalized KKT eigenvalues with the assembled Whitney
   generalized spectrum.  This is a numerical calibration after the exact
   identities, not the proof.

## Decision boundary

- **DERIVED LOCAL SPECTRAL PENCIL:** every exact identity passes and the small
  calibration agrees.  Then the exact Whitney spectrum has a coefficient-free
  element/star-local constrained representation without explicitly forming
  the assembled inverse.
- **REFUTED:** any exact identity or calibrated finite spectrum fails.

## Mandatory dynamical negative

The generalized time pencil has multiplier metric

\[
\mathcal M=\begin{pmatrix}M_{\rm loc}&0\\0&0\end{pmatrix},
\]

which is singular by exactly ​(operatorname{rank}Q=6360) directions.
Therefore the KKT system is a descriptor/differential-algebraic system, not an
ordinary Hilbert-space Schrödinger generator.

Also record the general self-adjoint dilation fact.  If

\[
H=\begin{pmatrix}H_{PP}&B\\B^*&H_{QQ}\end{pmatrix}

\]

is self-adjoint and the literal physical subspace (P\mathcal H) is invariant,
then (B^*=0), hence (B=0): any nontrivial auxiliary coupling destroys
exact autonomous invariance.  A different embedded physical subspace remains
possible, but its local isometry must be constructed rather than assumed.

Thus even a positive spectral result must carry:

- **DERIVED NEGATIVE:** this canonical KKT pencil is not a unitary Whitney
  tick;
- **OPEN:** a positive-metric local dilation with a selected nontrivial
  physical embedding;
- **OPEN:** whether its depth remains bounded under refinement.

Only the new targeted verifier will be run.  The full suite remains excluded
by user instruction.
