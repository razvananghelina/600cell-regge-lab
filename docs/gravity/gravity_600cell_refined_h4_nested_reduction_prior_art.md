# Prior-art gate: nested 6+3+1 reduction of the refined H4 equations

Date: 2026-08-20

Status: completed before solving any new reduced equation.

## Exact object and complete hypotheses

Use the same equal-boundary `K0=P(sd K_600)` slab, exact rank geometry,
corrected complex Lorentzian Regge-plus-dust action, selected total mass,
conditional `P1` weights and 12 schedule/time-reversal classes as the accepted
bounded ten-equation search.

Write the internal variables as

```text
x in R^6 : cross-diagonal log coordinates,
z in R^4 : rank-lapse log coordinates.
```

First solve the six cross equations `F(x,z)=0` for `x=x(z)`. Then decompose

```text
z=t*1_4 + Q*u,
```

where the three columns of `Q` are the normalized standard Helmert contrasts.
At fixed `t`, solve the three contrast equations

```text
Q^T R(t,u)=0,
R(t,u)=G_rho(x(t,u),t*1+Q*u),
```

and inspect the remaining common scalar

```text
g(t)=1_4^T R(t,u(t))/2.
```

A zero of all nested equations is exactly a zero of the original ten
equations on the connected solved branch. No common-lapse restriction is
imposed because `u` is solved, not set to zero.

## KNOWN mathematics and primary Regge context

The implicit-function theorem and finite-dimensional
Lyapunov--Schmidt/Schur elimination are standard: an invertible partial
Jacobian locally licenses elimination of those variables, while the reduced
equations retain the zeros of the full system on that branch. This is a
method, not a new physical principle.

For Regge evolution, internal equations must be eliminated before a boundary
Hamilton principal function or canonical map is assigned, and curved
discretizations can turn continuum lapse constraints into
background-dependent pseudo-constraints:

- Dittrich and Hoehn,
  [Canonical simplicial gravity](https://arxiv.org/abs/1108.1974);
- Dittrich and Hoehn,
  [From covariant to canonical formulations of discrete gravity](https://arxiv.org/abs/0912.1817);
- Bahr and Dittrich,
  [(Broken) Gauge Symmetries and Constraints in Regge Calculus](https://arxiv.org/abs/0905.1670).

A post-result search using `Lyapunov--Schmidt`, `Schur complement`,
`internal Regge equations`, `lapse pseudo-constraint` and `stationary
simplicial slab` located no primary source computing this projected
600-cell system. Search absence is not a novelty proof; external novelty is
**OPEN**.

## DERIVED controls already available

From the committed 12 high-precision Hessian classes:

- every `6x6` cross block has rank six;
- its smallest singular value is `5.93908093...` and condition number is
  `2048.651...` in every class;
- after eliminating it, every `4x4` lapse Schur complement has rank four,
  inertia `(2,0,2)` and condition number `1.0913...`;
- its eigenvalues are approximately
  `(-0.018369,-0.018212,+0.017603,+0.019210)`;
- all 12 lapse Schur complements coincide at the inherited fill even though
  the full ten-dimensional Hessians form 12 classes.

Thus the nested reduction is locally legitimate and much better conditioned
than the original ten-dimensional merit problem. These facts do not prove
that the branch continues globally.

## Framing attack

Scanning only `z=t*1_4` would not test the full system. The four rank orbits
are not permuted by the retained symmetry, and the original search allowed
four independent lapse coordinates. Such a one-dimensional restriction could
not falsify asymmetric roots and is forbidden as exclusion evidence.

The Helmert contrasts repair this: they are solved at each `t`. Even then, a
single continuation follows only one connected transverse branch. It cannot
exclude disconnected cross/contrast branches without additional seeds or an
interval/root-count certificate.

A sign-definite sampled `g(t)` is not a continuous sign theorem. Brackets may
license localized roots, while exclusion requires interval enclosures or
another rigorous continuous certificate.

## OPEN

- existence and uniqueness of the six-equation cross branch over the lapse
  domain;
- existence and uniqueness of the three contrast solutions;
- zeros and signs of the final scalar;
- equality of the nonlinear reduced branches across schedule classes;
- disconnected branches and global root nonexistence.

## Next admissible calculation

Preregister a target-free nested continuation with frozen lapse grid, seeds,
branch gates, solver tolerances, scalar bracketing and independent
high-precision substitution into all ten original equations. Report every
unresolved branch point; never interpolate across a gap.
