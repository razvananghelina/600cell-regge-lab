# The fixed Regge continuum gap is closed

Date: 2026-08-12

Preregistered protocol commit: `1682a46`.

Targeted verifier:
`reproducible/verify_whitney_regge_continuum_transfer.py`.

## Headline

The continuum carrier need not be changed to use the smooth-manifold FEEC
theorem.  Radial projection identifies the boundary of the regular 600-cell
with the smooth 3-sphere, while the exact piecewise-flat metric is transported
with it.  A smooth round metric is used only to construct uniformly bounded
commuting projections.  Fixed equivalence of the two Hilbert norms transfers
those projections and the compactness property back to the exact Regge
metric.

> **DERIVED/STRUCTURAL POSITIVE:** the rank-edgewise Whitney
> Hodge--Laplacians converge spectrally, degree by degree, to the Hodge
> Laplacians of the fixed piecewise-flat 600-cell boundary.  Eigenvalue
> multiplicities and eigenspaces converge and no spurious fixed-index modes
> occur.

Consequently the earlier `OPEN ANALYTIC GAP` is closed.  This is a theorem
about the exact three-dimensional Regge continuum.  It does not supply a
fourth dimension or any physical unit.

## 1. Why radial projection does not change the metric

Let `P` be the convex regular 600-cell, centred at the origin and normalized
to unit circumradius.  Every ray from the origin meets its boundary once, so

\[
 R:\partial P\longrightarrow S^3,
 \qquad R(x)=\frac{x}{\lVert x\rVert}
\]

is a bijection.  On every tetrahedral facet it is a smooth embedding.  It
therefore sends every fine simplex of

\[
 K_n=\operatorname{Esd}_{2^n}(\operatorname{sd}K)
\]

to a smooth intrinsic simplex in `S^3`.

The essential point is to push the flat metric forward:

\[
 g_R=(R^{-1})^*g_{\rm flat}.
\]

Then `R` is an isometry by definition.  A flat Whitney form and its
transported form have exactly the same `L2` product.  Thus every consistent
mass matrix and every discrete Hodge adjoint in the existing construction is
unchanged.  No round-metric finite-element operator has been substituted.

The round metric `g_0` on `S^3` is only an auxiliary analytic norm.

## 2. The norm equivalence is exact and uniform in level

For adjacent unit 600-cell vertices,

\[
 u_i\cdot u_j=\frac{\phi}{2},
 \qquad
 \lVert u_i-u_j\rVert^2=\frac1{\phi^2}.
\]

The centroid of four vertices of a facet is normal to that facet.  The exact
calculation gives its squared distance from the origin:

\[
 a^2=\frac{2+3\phi}{8}
     =\frac{7+3\sqrt5}{16}.
\]

The verifier derives, rather than assumes, the differential of `R`:

\[
 (R^*g_0)_x(v,w)
 =\frac{v\cdot w}{r^2}
  -\frac{(x\cdot v)(x\cdot w)}{r^4},
 \qquad r=\lVert x\rVert.
\]

On the three-dimensional tangent space of a facet, the eigenvalues relative
to its flat metric are

\[
 \frac{a^2}{r^4},\quad \frac1{r^2},\quad \frac1{r^2}.
\]

Since `a <= r <= 1`, all lie in the fixed interval

\[
 \frac{7+3\sqrt5}{16}
 \leq\lambda\leq
 28-12\sqrt5.
\]

The same fixed comparison induces degree-dependent but level-independent
equivalence constants for every `L2` space of forms.  It also implies
equivalence of the `H(d)` graph norms because `d` is metric-independent.

Across a shared triangular face, both tetrahedra induce the same equilateral
face Gram matrix.  Radial transport uses the same face map from either side.
Hence the tangential--tangential components of `g_R` agree exactly: `g_R` is
a Regge metric in the sense used by [Gawlik and
McKee](https://arxiv.org/abs/2410.15579).  Their Regge `L2/H(d)` construction
and metric quasi-isometry result support precisely this norm comparison.

That paper does **not** by itself prove our conclusion: its error-analysis
section starts from a smooth target metric approximated by Regge metrics.
Here the Regge metric is the fixed target.  The missing step is supplied
below.

## 3. The finite-element spaces are exactly the same spaces

On a flat fine simplex, a Whitney `k`-form is the pullback of a lowest-order
trimmed polynomial form from the reference simplex.  After applying `R`, it
is still the pullback of that same reference form, now through the smooth
simplex embedding `R` composed with the affine map.

Therefore the transported spaces are exactly

\[
 \mathcal P^-_1\Lambda^k(K_n)
\]

on the intrinsic smooth triangulations of `S^3`.  They are not merely
isomorphic spaces chosen after seeing a spectrum.  The finite-element
degrees of freedom and the exterior derivative are metric-independent; only
the `L2` product uses `g_R`.

The earlier exact certificates already supply:

- a conforming subcomplex in every degree;
- nesting and `h_n=h_0 2^{-n} -> 0`;
- three normalized tetrahedral shape classes at every level;
- exact consistent `L2(g_R)` products.

The radial map and its inverse have fixed derivative bounds on the finite
parent-facet cover.  Its higher derivatives are bounded there as well because
`r >= a > 0`.  Thus the transported smooth triangulations retain uniform
regularity in `g_0`.  Smooth forms are approximated by the canonical Whitney
interpolants, and density extends this approximability to the `H(d)` complex.

## 4. Uniform commuting projections transfer to the Regge norm

[Licht](https://arxiv.org/abs/2310.14276) constructs projections

\[
 \pi_h^k:L^2\Lambda^k(S^3,g_0)
 \longrightarrow \mathcal P^-_1\Lambda^k(K_h)
\]

which commute with `d`, act as the identity on the finite-element space, and
are uniformly `L2(g_0)`-bounded when the mesh regularity parameters are
uniform.  His construction is for intrinsic smooth triangulations, exactly
the role played by the radially transported `K_n`.

Let the fixed norm equivalence be

\[
 m_k\lVert u\rVert_{0}
 \leq\lVert u\rVert_R
 \leq M_k\lVert u\rVert_{0}.
\]

Then the *same linear projections* satisfy

\[
 \lVert\pi_h^ku\rVert_R
 \leq \frac{M_k}{m_k}
       \sup_h\lVert\pi_h^k\rVert_0
       \lVert u\rVert_R.
\]

The new bound has no `h`.  Commutation, projection onto the exact Whitney
space and approximation are algebraic/topological properties and are not
altered by the equivalent norm.

This is why the auxiliary round metric is legitimate: it constructs an
analytic projector but never replaces the Regge inner product used by the
discrete or continuum operator.

## 5. Compactness really is invariant here

This is the delicate step.  Equivalent `L2` inner products do **not** give the
same adjoint `d*`, and the verifier includes an exact finite-dimensional
control where the two minimal right inverses are genuinely different.

Let `W_j^k` be the same vector space with either equivalent inner product
`j=0,R`; let `V^k=dom(d_k)`, `Z^k=ker(d_k)` and
`B^{k+1}=ran(d_k)`.  Equivalence of graph norms makes `V`, `Z`, and `B` the
same topological spaces.  Closed range is therefore invariant.

For either metric, define the minimal right inverse

\[
 S_k^{(j)}:B^{k+1}longrightarrow
 V^k\cap(Z^k)^{\perp_j},
 \qquad d_kS_k^{(j)}=I.
\]

If the smooth complex has the compactness property, `S_k^(0)` is compact:
its image of a bounded set is bounded in the joint `d/d*` graph norm.  Let
`P_Z^(R)` be Regge-orthogonal projection onto `Z^k`.  It is bounded in either
equivalent norm, and

\[
 S_k^{(R)}=(I-P_Z^{(R)})S_k^{(0)}.
\]

Hence `S_k^(R)` is compact.  Apply the same argument at degree `k-1`; its
Regge Hilbert adjoint is compact too.

Every element of
`dom(d_k) intersect dom(d_(k-1)^(*R))` has the Regge Hodge decomposition

\[
 u=h
   +S_k^{(R)}d_ku
   +(S_{k-1}^{(R)})^*d_{k-1}^{*R}u.
\]

The harmonic part is finite-dimensional because its dimension is the
metric-independent cohomology dimension.  The other two terms are images of
bounded sets under compact operators.  The inclusion into `W_R^k` is
therefore compact.  Reversing the two metrics gives the converse.

> **DERIVED LEMMA:** for equivalent Hilbert inner products on the same closed
> complex, closed range and the Hilbert-complex compactness property are
> invariant, even though the adjoints and their domains may change.

The exact verifier control realizes the identity

\[
 S^{(R)}=(I-P_Z^{(R)})S^{(0)}
\]

with rational matrices and explicitly rejects equality of the two minimal
inverses.

## 6. Spectral conclusion

We now have every hypothesis of Theorem 3.19 in [Arnold, Falk and
Winther](https://arxiv.org/abs/0906.4325):

1. a closed Hilbert complex with the compactness property;
2. approximating finite-dimensional subcomplexes;
3. uniformly `L2(g_R)`-bounded commuting cochain projections.

It follows that the discrete Hodge--Laplacian eigenvalues and eigenspaces
converge degree by degree to those of the fixed Regge Hilbert complex, with
multiplicity and without spectral pollution at fixed index.

The exact finite Kaehler--Dirac operator satisfies

\[
 D_h^2=\bigoplus_{k=0}^3\Delta_{h,k}.
\]

Therefore the same fixed-index spectral statement holds for `D_h^2`.  This
does not claim uniform control of modes whose index runs to the lattice
cutoff.

## 7. What follows physically, and what does not

On the Regge continuum, `D=d+d*_R` is the self-adjoint Hodge--Dirac operator
of the closed Hilbert complex.  Away from the measure-zero Regge skeleton its
principal symbol is the already certified Clifford symbol

\[
 \sigma_D(\xi)^2=\lVert\xi\rVert_{g_R}^2 I.
\]

The commutator with a Lipschitz scalar is multiplication by this bounded
symbol.  The first-order finite-propagation theorem of [McIntosh and
Morris](https://arxiv.org/abs/1201.1818), which includes metric-measure
systems, applies to the self-adjoint evolution.  Thus, after supplying a time
parameter and coefficient `c`,

\[
 i\partial_t\psi=(cD+\mu\gamma)\psi
\]

has propagation speed bounded by `|c|` in the fixed Regge length metric; the
zeroth-order mass term does not change the cone.

This is now **DERIVED CONDITIONAL CONTINUUM CAUSALITY**, conditional only in
the physical sense that the theory has not selected time or the numerical
constant `c`.  It is no longer conditional on a missing continuum-convergence
theorem.

The ultraviolet negative is untouched: at every finite mesh the inverse
consistent mass is nonlocal, and the calibrated circle has cutoff velocity
`sqrt(2)c`.  Spectral convergence at fixed index does not manufacture a
strict finite-lattice cone.

## Attack on the framing

The earlier statement that radial projection would necessarily change the
metric was too coarse.  Radial projection changes the metric only if one
silently replaces the pushed-forward flat metric by the round metric.  Used
as a coordinate identification with `g_R` transported along it, it changes
nothing in the exact discrete construction.

The 2026 Regge-FEEC source was also not a turnkey solution: it assumes
commuting projections and analyzes a smooth target metric with Regge
approximants.  Treating its title as the desired theorem would still have
been wrong.  The decisive bridge is the equivalent-norm transfer of Licht's
smooth projections plus the compactness lemma above.

This closes a mathematical existence/convergence question.  It does not show
that nature selects this carrier or this dynamics.

## Status ledger

- **DERIVED:** radial identification of the convex 600-cell boundary with
  `S^3` while preserving the exact metric by pushforward.
- **DERIVED:** exact, level-independent Regge/round norm equivalence.
- **DERIVED:** the transported spaces are the same intrinsic Whitney spaces.
- **STRUCTURAL:** uniformly bounded smooth-metric commuting projections from
  Licht's theorem.
- **DERIVED:** uniform transfer of those projections to the Regge norm.
- **DERIVED:** compactness invariance under equivalent Hilbert metrics.
- **STRUCTURAL:** AFW spectral convergence theorem applies to each degree.
- **DERIVED/STRUCTURAL POSITIVE:** fixed-index spectrum/eigenspaces of
  `D_h^2` converge to the fixed Regge continuum without spurious modes.
- **DERIVED CONDITIONAL:** continuum propagation speed `|c|` once physical
  time and `c` are supplied.
- **DERIVED NEGATIVE:** no strict finite-mesh `c`-cone; cutoff witness remains
  `sqrt(2)c`.
- **NOT DERIVED:** Lorentzian time, numerical `c`, any mass value, a fourth
  dimension, `hbar`, Newton's `G`, Planck length/time/mass, or the Standard
  Model.

## Reproduction

Only the targeted verifier was run:

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_regge_continuum_transfer.py
```

Result: `15/15` checks passed.  The full suite was not run by user request.

The registry was audited without executing the suite: `178` distinct scripts
are listed exactly once, `2` additional files are deliberately skipped with
recorded reasons, and there are no unregistered or missing verifier files.
