# Iterated barycentric subdivision is provably not shape regular

Date: 2026-08-11

Preregistration commit: `dd1f5c1`

Targeted verifier:
`reproducible/verify_whitney_barycentric_shape.py`

Targeted result: **10/10 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Headline

The refinement carrier itself, not only the Whitney inverse, has a rigorous
asymptotic obstruction:

> **DERIVED NEGATIVE:** iterated barycentric subdivision of the tetrahedral
> carrier is not uniformly shape regular.

This means that the recently observed stiffness factors cannot be interpreted
as a clean physical renormalization flow.  They include compensation for
increasingly anisotropic numerical elements.

The result does not depend on a fitted exponent, finite-depth trend, or
eigensolver.  One exactly constructed nested flag proves the failure.

## Exact repeated-flag theorem

Write the parent edge matrix as

\[
 A=[p_1-p_0\;p_2-p_0\;p_3-p_0].
\]

At every level choose the same barycentric flag `(0,1,2,3)`.  Direct exact
subtraction of its barycentres gives

\[
 A'=AT,
 \qquad
 T=
 \begin{pmatrix}
 1/2&1/3&1/4\\
 0&1/3&1/4\\
 0&0&1/4
 \end{pmatrix}.
\]

The verifier checks this identity exactly.  Since (T) is triangular,

\[
 \sigma_{\rm eig}(T)=\{1/2,1/3,1/4\},
 \qquad
 \det T=1/24.
\]

After (n) repetitions,

\[
 A_n=A_0T^n.
\]

For every matrix, its spectral norm is at least its spectral radius.  Hence

\[
 \lVert T^n\rVert_2\geq(1/2)^n,
 \qquad
 \lVert T^{-n}\rVert_2\geq4^n,
\]

and therefore

\[
 \kappa_2(T^n)\geq2^n.
\]

The regular reference edge matrix has exact Gram eigenvalues

\[
 (16,4,4),
\]

so (kappa_2(A_0)=2).  Submultiplicativity applied to
(T^n=A_0^{-1}A_n) gives

\[
 \kappa_2(A_n)\geq2^{n-1}.
\]

Finally, for (G_n=A_n^*A_n),

\[
 \kappa_2(G_n)=\kappa_2(A_n)^2
 \geq4^{n-1}.
\]

Thus no level-independent shape-regularity constant exists.  The conclusion
is asymptotic and exact.

## Complete exact type census

All 24 child orderings were propagated exactly over the rationals through
four levels.  Identical Gram matrices were merged exactly, not by tolerance.

| level | exact Gram types | total tetrahedra represented | minimum mean-ratio quality | maximum affine condition | maximum (h\lVert D\rVert) |
|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 1.000000 | 2.000 | 10.954 |
| 1 | 1 | 24 | 0.596793 | 7.101 | 20.675 |
| 2 | 24 | 576 | 0.195389 | 23.158 | 71.551 |
| 3 | 576 | 13,824 | 0.053073 | 73.053 | 232.215 |
| 4 | 13,824 | 331,776 | 0.014013 | 243.455 | 736.701 |

The mean-ratio quality is permutation invariant and normalized to one for a
regular tetrahedron.  Its rapid decay confirms that the obstruction is not an
artefact of the ordered affine-coordinate condition number.

The dimensionless quantity (h\lVert D\rVert) removes uniform element size.
Its growth therefore measures shape damage rather than the expected
(1/h) Dirac scaling.

## Repeated-chain numerical control

The exact recurrence was evaluated through eight nested steps.  The Gram
condition numbers are approximately

```text
4,
50.4,
512,
4.14e3,
2.77e4,
1.60e5,
8.33e5,
4.02e6,
1.84e7.
```

Over the same chain:

- mean-ratio quality falls from `1` to `0.000959`;
- (h\lVert D\rVert) rises from `10.95` to approximately `15,480`.

These numbers calibrate the theorem; they are not used to infer an exponent.

The numerical Gram-only Whitney integration independently reproduces the
exact symbolic local Dirac norms with relative residual
(8.88\times10^{-16}).

## Consequence for the stiffness results

The exact trace-stiffness construction remains valid at every fixed finite
level.  What fails is the interpretation of its barycentric refinement
factors as pure scale flow.

The observed growth mixes at least two effects:

1. smaller physical elements, which should produce the ordinary (1/h)
   Dirac scaling;
2. worsening shapes, which make the dimensionless factor
   (h\lVert D\rVert) diverge along a nested subroute.

Therefore the apparent clustering of degreewise stiffness factors is not
evidence for a fundamental constant until it survives on a uniformly
shape-regular tower.

## What this closes

It is no longer legitimate to use indefinite iterated barycentric
subdivision as a neutral continuum limit and attribute all operator growth to
physics.

This closes:

- a level-independent local condition bound on that tower;
- a straightforward physical reading of its fitted or observed refinement
  factors;
- the assumption that more barycentric levels automatically improve the
  geometric approximation uniformly.

## What it does not close

The theorem does not reject:

- the base 600-cell or either fixed finite refinement;
- exact Whitney and trace-stiffness identities on fixed carriers;
- a shape-regular tetrahedral refinement;
- adaptive remeshing selected by additional geometry;
- a non-simplicial local microscopic carrier.

The next route must confront a genuine selection tradeoff: barycentric
subdivision is canonical and symmetric but degenerates; common shape-regular
tetrahedral refinements usually require choices of edge order or octahedral
diagonal.  Any replacement must be selected by the theory rather than chosen
for a favorable spectrum.

## Status ledger

- **DERIVED:** exact repeated-child affine transform.
- **DERIVED:** exact eigenvalues `(1/2,1/3,1/4)` and determinant `1/24`.
- **DERIVED NEGATIVE:** unbounded affine and Gram conditioning.
- **DERIVED NEGATIVE:** no uniform shape regularity for the barycentric tower.
- **DERIVED CONTROL:** complete exact type census through level four.
- **DERIVED CONTROL:** mean-ratio and normalized-Dirac degradation.
- **STRUCTURAL CORRECTION:** barycentric stiffness flow is contaminated by
  mesh shape.
- **OPEN:** a geometry-selected shape-regular refinement.
- **OPEN:** trace stiffness and spectral convergence on that carrier.
- **OPEN:** absolute stiffness, chirality, Lorentzian time, and causal speed.
- **NOT CLAIMED:** mass, inertia, (c), (hbar), Newton's (G), or a Planck
  scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_barycentric_shape.py
```

Expected result: `10/10`.
