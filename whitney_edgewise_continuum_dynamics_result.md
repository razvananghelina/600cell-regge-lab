# Exact edgewise dynamics exists; the Regge continuum theorem is still missing

Date: 2026-08-12

Preregistered protocol commit: `075cb38`.

## Headline

The canonical shape-regular tower

\[
 K_n=\operatorname{Esd}_{2^n}(\operatorname{sd}K)
\]

supports the complete consistent Whitney--Kaehler--Dirac construction in all
cochain degrees.

> **DERIVED EDGEWISE GALERKIN-INDUCTIVE DYNAMICS.**  The inclusions commute
> with `d`, preserve the exact `L2` metric, compress the weak Dirac form, and
> preserve form parity exactly.  The weak generator is simplex-star local and
> the finite dynamics is nontrivial.

However, a hostile reading of the primary FEEC sources found a remaining
scope gap:

> **OPEN ANALYTIC GAP.**  The cited convergence theorems explicitly treat a
> Lipschitz domain in Euclidean space or a smooth compact manifold.  The exact
> theory uses a *closed piecewise-flat Regge manifold*, singular along
> codimension-two cells.  None of the inspected statements explicitly proves
> spectral convergence on that carrier.

Thus the old mesh-degeneration blocker is removed, but continuum convergence
is not yet promoted to a theorem by citation alone.

## 1. Exact finite and all-level construction

The verifier independently reconstructs the Edelsbrunner--Grayson color
schemes on one rank-ordered barycentric orthoscheme.  The exact local carriers
are

| carrier | f-vector |
|---|---:|
| `Esd_1` | `(4,6,4,1)` |
| `Esd_2` | `(10,25,24,8)` |

For every degree `p=0,1,2,3`, it derives the consistent Whitney mass by exact
affine integration and constructs the inclusion by barycentric determinants.
It then proves

\[
 d_fP_p=P_{p+1}d_c,
 \qquad
 P_p^*M_{f,p}P_p=M_{c,p}.
\]

With `A=M D` and the form-parity grading `gamma`, it also proves

\[
 P^*A_fP=A_c,
 \qquad
 \gamma_fP=P\gamma_c.
\]

The calculation is local and affine-natural.  Combined with the already
proved face conformity and composition of edgewise subdivision, it applies
on every level and on the full 600-cell assembly.  The number of top cells is

\[
 14400\,8^n,
\]

the mesh scale is `h_n=h_0/2^n`, and only three normalized tetrahedral shape
classes occur.  The construction is therefore an actual infinite inductive
system, not a finite sequence guessed from a trend.  The edgewise composition
and finite-shape statements are the content of the primary construction of
[Edelsbrunner and Grayson](https://doi.org/10.1007/s004540010063).

## 2. The dynamics is nontrivial, but induction is Galerkin

On the exact `Esd_2` control:

- the smallest audited mass eigenvalue is `1.117113e-3`;
- the weak Dirac stencil has `308` directed nonzeros and zero simplex-star
  locality violations;
- the Dirac kernel has dimension one, the Betti sum of a tetrahedron;
- the spectral radius is `22.409846...`;
- the grading-forced `+/-` pairing residual is `3.55e-14`.

Strong operator intertwining still fails.  The exact adjoint-leakage ranks are

\[
 (5,4,1)
\]

in degrees `(0,1,2)`.  Every leakage is invisible after Galerkin compression.
This is the same essential distinction as before: exact quadratic-form
induction is true; equality of the strong operators on inherited states is
false.

## 3. Why ordinary mass lumping is not the repair

Row-sum lumping breaks exact refinement isometry with ranks

\[
 (3,6,4,0)
\]

in degrees `(0,1,2,3)`.

A post-protocol framing attack found a more basic problem.  Reversing the
orientation of a cochain basis vector changes a mass matrix by `M -> S M S`.
A geometric metric rule must transform covariantly.  Ordinary row-sum lumping
does not do so for one- and two-forms; it passes only in top degree on this
control.

> **DERIVED POST-PROTOCOL NEGATIVE:** the row-sum operation used as a local
> negative control is not itself a basis-independent geometric metric for the
> middle cochain degrees.

This does not invalidate its use as a diagnostic that the support census can
return a local answer.  It does invalidate promoting that particular lumping
rule to fundamental dynamics without a new orientation-covariant derivation.

## 4. What the FEEC literature does and does not prove here

The finite construction now supplies the familiar FEEC inputs:

- a conforming de Rham subcomplex;
- exact consistent `L2` products;
- nested meshes with `h -> 0`;
- uniform shape regularity;
- a compact fixed cellular carrier.

[Arnold, Falk and Winther](https://arxiv.org/abs/0906.4325) show that a
subcomplex with uniformly bounded commuting projections yields stable
Hodge--Laplacian and eigenvalue approximation.  [Christiansen's finite-element
systems framework](https://arxiv.org/abs/1006.4779) works over cellular
complexes, includes polyhedral grids, constructs stable commuting projections,
and obtains Hodge--Laplacian eigenpair approximation.  His separate
[compact-manifold result](https://arxiv.org/abs/1007.1120) proves discrete
Poincare and Rellich compactness.

But the detailed scope matters:

- the analytic convergence section of `1006.4779` specializes to a domain
  `S` in Euclidean space;
- `1007.1120` specializes to smooth compact manifolds;
- the exact 600-cell metric is a closed piecewise-flat metric with Regge cone
  singularities.

It is plausible that the cellular proof extends by a finite atlas/gluing
argument, but the repository does not yet contain that proof and the inspected
sources do not state it in this exact form.  Therefore spectral convergence on
the theory's own continuum carrier remains **OPEN**, not **STRUCTURAL**.

No smooth radial projection may be substituted silently: it changes the
metric and hence changes the exact mass matrices being studied.

## 5. The causal symbol is exact, conditionally

On the eight-dimensional exterior algebra of a three-dimensional cotangent
space, the verifier constructs exterior multiplication and contraction and
proves

\[
 \sigma_D(\xi)
 =i\bigl(\epsilon(\xi)-\iota(\xi)\bigr),
 \qquad
 \sigma_D(\xi)^2=|\xi|^2I.
\]

With form parity,

\[
 \bigl(c\,\sigma_D(\xi)+\mu\gamma\bigr)^2
   =\bigl(c^2|\xi|^2+\mu^2\bigr)I.
\]

Therefore, **if** the edgewise family is proved to converge to the intended
self-adjoint Hodge--Dirac continuum, the equation

\[
 i\partial_t\psi=(cD+\mu\gamma)\psi
\]

has characteristic speed `|c|`; the mass is zeroth order and does not alter
the causal cone.  Finite propagation for self-adjoint first-order systems on
Riemannian/metric-measure carriers is established by
[McIntosh and Morris](https://arxiv.org/abs/1201.1818).

This is an exact principal-symbol result, not a derivation of physical time or
of the numerical value of `c`.

## 6. The ultraviolet problem is not hidden

For the calibrated consistent-Whitney circle dispersion, the velocity at
dimensionless momentum `q=kh` has expansion

\[
 \frac{v}{c}=1+\frac{q^2}{8}+\frac{q^4}{384}+O(q^5).
\]

Thus fixed physical modes approach speed `c` as `h -> 0`.  But the finite
cutoff witness remains

\[
 \frac{v}{c}=\sqrt2
 \quad\text{at}\quad q=\frac{2\pi}{3}.
\]

The correct interpretation is therefore:

- a continuum limit, if established, has the Dirac characteristic speed;
- no finite consistent-mass lattice in this family has a strict speed-`c`
  cone.

## Status ledger

- **DERIVED:** canonical nested shape-regular edgewise carrier at all levels.
- **DERIVED:** exact all-degree Whitney inclusions on `Esd_1 -> Esd_2`.
- **DERIVED:** exact `d`, metric, weak-Dirac and grading compatibility.
- **DERIVED:** nontrivial finite metric-unitary Galerkin dynamics.
- **DERIVED:** simplex-star locality of the weak generator.
- **DERIVED NEGATIVE:** strong leakage ranks `(5,4,1)`.
- **DERIVED NEGATIVE:** lumped isometry residual ranks `(3,6,4,0)`.
- **DERIVED POST-PROTOCOL NEGATIVE:** row-sum lumping is not
  orientation-covariant in middle degrees.
- **DERIVED:** three-dimensional continuum Hodge--Dirac principal symbol and
  conditional speed `|c|`.
- **DERIVED NEGATIVE:** the finite cutoff still permits `sqrt(2)c`.
- **OPEN ANALYTIC GAP:** spectral/strong-resolvent convergence on the closed
  Regge carrier.
- **CONDITIONAL:** continuum finite propagation, once that gap is closed.
- **NOT DERIVED:** time, numerical `c`, mass values, a fourth dimension,
  `hbar`, Newton's `G` or Planck units.

## Reproduction

Only the targeted verifier was run:

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_edgewise_continuum_dynamics.py
```

Result: `22/22` checks passed.  The full suite was not run by user request.
