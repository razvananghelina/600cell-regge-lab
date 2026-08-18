# Protocol: fixed-Regge FEEC transfer

Date: 2026-08-12

This protocol is committed before the transfer verifier and before changing
the continuum verdict in `whitney_edgewise_continuum_dynamics_result.md`.
It contains no spectral value, particle, mass, coupling, clock, Planck or
Standard-Model target.

## Question

The exact Whitney tower uses the fixed piecewise-Euclidean metric on the
boundary of the 600-cell.  The existing result left Hodge--Laplacian spectral
convergence **OPEN** because the inspected FEEC theorems treated either a
Euclidean domain or a smooth compact Riemannian manifold.

Can that scope gap be closed without replacing the Regge metric or its exact
Whitney mass matrices?

## Frozen carrier identification

Let `P` be the convex regular 600-cell, centred at the origin, with vertices
on the unit 3-sphere, and let

\[
 R:\partial P\longrightarrow S^3,\qquad R(x)=x/\lVert x\rVert .
\]

Convexity makes `R` a bijection.  On each tetrahedral facet it is a smooth
embedding, so it transports every rank-edgewise triangulation to an intrinsic
smooth triangulation of the *smooth carrier* `S^3`.

Crucially, transport the piecewise-flat metric as well:

\[
 g_R=(R^{-1})^*g_{\rm flat}.
\]

Thus `R` is only a coordinate identification.  The physical/discrete metric
remains the exact Regge metric and the Whitney mass matrices do not change.
The round metric `g_0` on `S^3` may be used only as an auxiliary norm for an
analytic theorem.

For a unit-circumradius regular 600-cell facet, adjacent vertices have scalar
product `phi/2`.  Its supporting hyperplane has squared distance

\[
 a^2=\frac{2+3\phi}{8}.
\]

The verifier must derive the pullback of the round metric,

\[
 (R^*g_0)_x(v,w)
 =\frac{v\cdot w}{r^2}
  -\frac{(x\cdot v)(x\cdot w)}{r^4},
 \qquad r=\lVert x\rVert,
\]

and its tangent eigenvalues `1/r^2, 1/r^2, a^2/r^4`.  Since
`a <= r <= 1`, this gives a level-independent equivalence between `g_R` and
`g_0`.

## Frozen analytic transfer lemma

Use the following claim only after proving it in the result note.

> If two Hilbert inner products on every degree of the same closed complex
> are equivalent, then closed range and the compactness property are
> invariant under the change of inner product.

The proof must not assert that the adjoint domains are identical; in general
they are not.  It must instead use the metric-independent quotient
`V^k/ker(d_k)` and the compact minimal right inverses of `d_k`.  Orthogonal
minimal representatives for either inner product differ by bounded
projections onto `ker(d_k)`, so compactness transfers.

## Frozen FEEC hypothesis chain

The acceptance route must establish all of the following.

1. The transported flat metric is a Regge metric: it is smooth on each
   simplex and its tangential--tangential trace agrees across every face.
2. The exact cochain spaces are the intrinsic lowest-order trimmed spaces
   `P^-_1 Lambda^k` for the transported smooth triangulation; the metric enters
   only their `L2` products.
3. The previously certified finite shape set and `h -> 0`, together with the
   fixed radial norm equivalence, give uniform regularity in the auxiliary
   round metric.
4. Licht's smoothed commuting projections for the smooth triangulations are
   uniformly `L2(g_0)`-bounded.  Fixed norm equivalence transfers the same
   projections to uniformly `L2(g_R)`-bounded commuting projections onto the
   *same* Whitney spaces.
5. Approximation in the smooth norm transfers to the Regge norm.  No round
   mass, radial finite-element operator or geometric-error limit may replace
   the exact metric.
6. Compactness of the smooth de Rham Hilbert complex transfers to the fixed
   Regge inner product by the lemma above.
7. Arnold--Falk--Winther's abstract eigenvalue theorem then applies degree by
   degree.  Since the discrete `D_h^2` is the direct sum of those discrete
   Hodge Laplacians, its fixed-index eigenvalues and eigenspaces converge
   without spurious modes.

Primary sources to audit against their actual hypotheses:

- Gawlik and McKee, *Intrinsic Finite Element Error Analysis on Manifolds
  with Regge Metrics*, arXiv:`2410.15579`, especially the Regge `L2/H(d)`
  construction and metric quasi-isometry result;
- Licht, *Smoothed Projections over Manifolds in Finite Element Exterior
  Calculus*, arXiv:`2310.14276`, especially the smooth-triangulation finite
  element spaces and uniformly bounded commuting projections;
- Arnold, Falk and Winther, *Finite Element Exterior Calculus: from Hodge
  Theory to Numerical Stability*, arXiv:`0906.4325`, especially the
  compactness hypothesis and Hodge--Laplacian eigenvalue theorem.

## Decision boundaries

### Acceptance

All seven hypotheses are established without changing `g_R`.  Label the
result **DERIVED/STRUCTURAL: fixed-Regge Hodge--Laplacian spectral
convergence**.  The earlier analytic gap is closed.

This permits the already derived continuum Hodge--Dirac symbol and standard
finite-propagation theorem to be applied to the fixed Regge continuum, but it
still does not derive Lorentzian time, the numerical value of `c`, masses,
`hbar`, Newton's `G`, Planck units, or a fourth dimension.

### Kill / remain open

Retain the **OPEN ANALYTIC GAP** if any of the following occurs:

- radial transport changes the exact mass metric;
- the intrinsic FE spaces differ from the exact Whitney spaces;
- the auxiliary smooth-metric projections are not uniform on the edgewise
  tower or do not transfer under norm equivalence;
- compactness is not invariant under the equivalent inner products;
- the abstract spectral theorem needs an unverified regularity hypothesis.

No finite numerical trend or fitted continuum spectrum can override a failed
hypothesis.
