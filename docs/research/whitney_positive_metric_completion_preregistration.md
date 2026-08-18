# Preregistration: can a positive multiplier metric preserve Whitney dynamics?

Date: 2026-08-11

## Fixed family and scope

The canonical constrained Whitney pencil has independent constraint rank
(r) and descriptor matrices

\[
K=\begin{pmatrix}A_{\rm loc}&C^*\\C&0\end{pmatrix},
\qquad
B_0=\begin{pmatrix}M_{\rm loc}&0\\0&0\end{pmatrix}.
\]

It has the correct (n-r) finite physical eigenvalues but a singular metric.
Test the minimal ordinary-Hilbert completion which changes no operator block
and only gives the independent multipliers a positive kinetic metric:

\[
B_\varepsilon=
\begin{pmatrix}M_{\rm loc}&0\\0&\varepsilon N\end{pmatrix},
\qquad \varepsilon>0,quad N>0.
\]

The theorem must be stated for every positive (N).  The numerical control
uses (N=I) in a recorded independent constraint basis only as a calibration;
that basis is not promoted to canonical geometry.

This protocol does not cover new Hamiltonian blocks, extra Stueckelberg
variables, an indefinite/Krein metric or a different embedded physical
subspace.

## Frozen dimension obstruction

For ​(arepsilon>0), (B_\varepsilon) is positive definite on the complete
((n+r))-dimensional complex carrier.  Therefore the generalized problem has
(n+r) finite eigenvalues counting algebraic multiplicity, whereas the exact
descriptor has only (n-r) physical finite eigenvalues.  The completion adds

\[
(n+r)-(n-r)=2r

\]

extra finite complex spectral slots.  This equals the real second-class
constraint count.

Record the numbers:

- base: (n=9000), (r=6360), physical `2640`, positive-metric dimension
  `15360`, extra `12720`;
- first barycentric: (n=216000), (r=153120), physical `62880`,
  positive-metric dimension `369120`, extra `306240`.

No claim that the extra modes are experimentally excluded is permitted.  The
point is exact spectral non-equivalence unless an additional selection rule
removes them.

## Frozen eigenvalue-shift theorem

For a simple finite descriptor eigenpair

\[
Ky=zB_0y,qquad y=(u,\lambda),qquad u^*M_{\rm loc}u\ne0,

\]

differentiate (Ky(\varepsilon)=z(\varepsilon)B_\varepsilon
y(\varepsilon)) at zero.  Symmetry gives

\[
z'(0)=-z\,
\frac{\lambda^*N\lambda}{u^*M_{\rm loc}u}.
\]

Thus every nonzero simple mode with ​(lambda\ne0) moves at first order for
every (N>0).  In a degenerate nonzero eigenspace the first-order splitting
matrix is ​(-zL^*NL); it vanishes on the entire eigenspace only if every
multiplier vector vanishes.

The prior full-rank leakage theorem shows that local and assembled evolution
do not intertwine on positive-degree conforming inputs.  Do not jump from
that statement to an all-eigenmode claim.  The small control must explicitly
measure each finite eigenspace.

## Boundary-of-4-simplex calibration

Use the same five-tetrahedron 3-sphere and full-row-rank 45-constraint basis
as the previous exact controls.  It has

\[
n=75,qquad r=45,qquad n-r=30,qquad n+r=120.

\]

1. Compute the 30 finite descriptor eigenvalues and left/right eigenvectors.
2. Group exact/numerical degeneracies before interpreting derivatives.
3. For every nonzero eigenspace, compute the eigenvalues of the frozen
   first-order splitting matrix with (N=I).  Record how many spaces have a
   nonzero shift.
4. For ​(arepsilon=10^{-6}), compute all 120 positive-metric eigenvalues.
   Match only the 30 branches continuously nearest the descriptor finite
   values and compare finite differences with the perturbative splitting.
   This is a calibration, not a fitted ​(arepsilon).
5. Record the 90 additional finite eigenvalues that descend from the
   descriptor's infinite sector.

## Decision boundary

- **DERIVED POSITIVE-METRIC NO-GO FOR THE MINIMAL COMPLETION:** all nonzero
  physical eigenspaces shift and the finite spectral count jumps by (2r).
- **PARTIAL NEGATIVE:** only some nonzero eigenspaces shift.  Report the hit
  fraction and do not generalize beyond it.
- **REFUTED:** the perturbative formula or spectral count fails.

A negative result closes only the block-fixed multiplier-metric completion.
It does not close a new first-class extension or an enlarged Hamiltonian whose
physical subspace is selected by additional geometry.

Only the new targeted verifier will be run.
