# Preregistration: is the Whitney metric inverse a short local polynomial?

Date: 2026-08-11

## Framing attack

For every finite invertible matrix (M), Cayley--Hamilton already implies
that (M^{-1}) is a polynomial in (M) of degree at most (N-1).  Merely
exhibiting a finite polynomial would therefore be content-free and would not
remove the locality problem.

The nontrivial questions frozen here are:

1. what is the **exact minimal-polynomial degree** of each assembled Whitney
   mass block;
2. what is the exact least-degree inverse polynomial selected by that minimal
   polynomial;
3. how far its support propagates in the one-step sparsity graph of the mass;
4. whether the degree is small because of 600-cell symmetry, rather than only
   finite because the carrier is finite.

No phenomenological number or desired answer enters this audit.

## Fixed matrices

Rebuild the full 600-cell simplicial complex and derive the four regular
tetrahedron Whitney masses from the defining form integrals.  Assemble

\[
M_p=J_p^TM_p^{\rm loc}J_p,
\qquad p=0,1,2,3,
\]

then clear the unique rational denominator and divide the integer numerator
by the gcd of all its entries.  Call the resulting primitive symmetric
integer matrix (B_p).  Scaling changes the inverse coefficients but not its
minimal-polynomial degree or support.

Record dimension, nonzero count, maximum absolute row sum and connected
components of the off-diagonal support graph before computing any inverse.

## Exact minimal-polynomial certificate

Numerical eigenvalues may be printed as a diagnostic only.  They cannot
certify the number of distinct eigenvalues.

For every (B_p):

1. Generate a deterministic integer probe vector from a recorded seed.
2. Over several recorded primes, compute the scalar Krylov sequence
   (v^TB_p^kv) and its Berlekamp--Massey recurrence.  A recurrence of degree
   (s) is an exact lower bound on the rational minimal-polynomial degree.
3. Reconstruct the common monic integer polynomial (m_p(x)) by Chinese
   remaindering.  Continue until the modulus exceeds twice the a priori
   coefficient bound obtained from the maximum absolute row sum.
4. Certify (m_p(B_p)=0) on the **whole carrier**, not merely on the probe,
   by modular Horner evaluation on the identity matrix.  Use enough primes
   that their product exceeds twice the corresponding integer residual
   bound.  This is an exact zero certificate, not a probabilistic test.

If the single probe misses a factor, its candidate will fail the full-matrix
annihilation.  Add deterministic probe vectors and take polynomial lcms; do
not accept a polynomial until the whole-matrix certificate passes.  Once a
degree-(s) annihilator passes and a scalar Krylov sequence has exact linear
complexity (s) for at least one prime, minimality is proved.

## Exact inverse and locality

If

\[
m_p(x)=c_0+c_1x+\cdots+c_{s-1}x^{s-1}+x^s
\]

has (c_0\ne0), then the uniquely selected inverse polynomial of degree below
(s) is

\[
B_p^{-1}=-\frac{c_1I+c_2B_p+\cdots+B_p^{s-1}}{c_0}.
\]

Evaluate its integer numerator by modular Horner arithmetic.  For every
entry, a nonzero residue proves nonzero support; entries zero for all primes
are certified as exact zeros only when the modulus product exceeds twice the
entry bound.

Compute shortest-path distance in the off-diagonal support graph of (B_p),
and report:

- minimal-polynomial degree (s_p) and inverse-polynomial degree;
- graph diameter of every connected component;
- maximum distance at which the inverse is nonzero;
- exact nonzero fraction of the inverse, including cross-component zeros;
- the trivial Cayley--Hamilton ratio
  `(inverse degree)/(matrix dimension - 1)`.

## Decision boundaries

- **DERIVED SHORT FINITE ALGEBRA:** the exact degrees are far below (N-1)
  because only a small number of spectral classes occurs.  Report the numbers
  without converting “small” into a continuum claim.
- **DERIVED NONLOCAL ON THE FIXED COMPLEX:** an inverse has nonzero entries at
  the diameter of its connected support graph.  Then it cannot be a one-step
  operator in that graph, even though a finite-depth polynomial exists.
- **DERIVED LOCAL BLOCK:** the inverse has no support beyond one-step
  neighbours in a block.
- **INCONCLUSIVE:** coefficient reconstruction or the exact whole-matrix
  annihilation certificate fails.  Numerical clustering alone is not a
  result.

Even a short exact polynomial on this one complex does not show bounded-depth
locality under refinement.  That requires a preregistered refinement family
and degree scaling.  Nor does a polynomial automatically give a unitary
microscopic tick; a reversible implementation and any ancillas must still be
constructed.

Only the new targeted verifier will be run.  The full suite remains excluded
by user instruction.
