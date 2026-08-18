# The connected three-bond walk is strongly anisotropic after one period

Date: 2026-08-11  
Preregistration commit: `6811eee`

## Result

The three-bond schedule repairs support connectivity, but the fixed published
coin does not produce an isotropic local step on the first barycentric
600-cell geometry.

> **DERIVED NUMERICAL ONE-PERIOD LOCAL ANISOTROPY:** none of the six temporal
> orders passes the preregistered zero-drift and tangent-covariance gate.

\[
\text{isotropy hits}=0/6.
\]

For the designated rank-forward order ((01)\to(12)\to(23)),

\[
\|\mu\|=0.030868454,
\qquad
\frac{\lambda_{\max}}{\lambda_{\min}}=11.692106950,
\qquad
A=0.542433387.
\]

This is a large failure, not numerical noise around one.

The targeted verifier passes `9/9` in about one second.  No full suite was
run.

## Calibrated observable

The test starts from one chamber with the maximally mixed active spin state
(I_4/4), applies the fixed published coin after each of the three
translations, and averages the four evolved basis-state probabilities.

Final chamber positions are mapped to the tangent space of the initial point
by the exact spherical logarithm on (S^3).  The verifier then measures the
mean displacement and the three eigenvalues of the centred tangent
covariance.

The known-answer control uses the four equally weighted directions from the
centre of a regular tetrahedron.  It returns

\[
\mu=0,
\qquad R=1,
\qquad A=1.92\times10^{-16}.
\]

Thus the estimator recognizes isotropy when it is present.

## Geometry already breaks equality of the four microsteps

The four colour-neighbour geodesic lengths are constant across chambers, as
required by (H_4) transitivity, but unequal between colours:

\[
(\ell_0,\ell_1,\ell_2,\ell_3)
=
(0.15706954,\ 0.08251489,\ 0.05599560,\ 0.06858481).
\]

One tick per chamber facet therefore does not mean one equal physical
distance per tick.  This is the concrete metric reason that a strict graph
light cone is not yet a physical light cone.

## All six preregistered orders

| schedule | support chambers | drift norm | covariance ratio (R) | residual (A) |
|---|---:|---:|---:|---:|
| 01→12→23 | 7 | 0.030868 | 11.6921 | 0.5424 |
| 01→23→12 | 10 | 0.042329 | 7.2486 | 0.4978 |
| 12→01→23 | 10 | 0.030355 | 8.2556 | 0.4476 |
| 12→23→01 | 10 | 0.039469 | 3.9941 | 0.2747 |
| 23→01→12 | 10 | 0.029205 | 8.6579 | 0.4152 |
| 23→12→01 | 7 | 0.030868 | 11.6921 | 0.5424 |

The fourth row is numerically the least anisotropic of the six, but choosing
it for that reason would be target-driven fitting.  No schedule is promoted.

## Physical interpretation

The current candidate now has:

- exact unitarity;
- a strict finite graph cone;
- global support connectivity;
- **no local metric isotropy**.

Therefore it still does not derive a universal light speed or the massless
Dirac equation on this geometry.  Graph locality alone was insufficient.

This does not yet close all continuum behaviour.  It is a one-period result
at one finite refinement level.  However, the repository already contains an
exact repeated-flag witness showing that unmodified iterated barycentric
refinement is not shape-regular and develops unbounded anisotropy.  It is
therefore not defensible simply to assume that the present mismatch washes
out.

The next honest route is metric-aware: ask whether the four derived tangent
steps admit canonical positive weights or a local tetrad coin that cancels
drift and yields a tight frame.  Such weights must be derived from the
geometry/Whitney metric and preregistered, not selected by minimizing the
table above.

## Status ledger

- **DERIVED NUMERICAL:** all four colour step lengths are individually
  constant but mutually unequal.
- **DERIVED NUMERICAL:** designated one-period drift and covariance are
  strongly anisotropic.
- **DERIVED NUMERICAL:** zero isotropy hits out of six complete schedules.
- **DERIVED CONTROL:** the estimator returns isotropy on a regular
  tetrahedral tight frame.
- **STRUCTURAL:** connected combinatorial propagation survives.
- **OPEN:** existence and uniqueness of a geometry-derived weighted tight
  frame / tetrad coin.
- **OPEN:** amplitude-level multistep and controlled refinement behaviour.
- **NOT CLAIMED:** a physical light cone, Dirac limit, mass or Planck scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_h4_three_bond_local_isotropy.py
```

Expected result: `9/9`.
