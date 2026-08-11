# Scalar weights cannot isotropize the local H4 chamber

Date: 2026-08-11  
Preregistration commit: `9beb49c`

## Result

The four metric neighbour steps admit a unique positive choice of weights
that cancels the local drift.  Those same weights do **not** make the second
moment isotropic.  In fact, the complete unconstrained linear system is
inconsistent before positivity is imposed.

> **DERIVED NUMERICAL SCALAR NO-GO AT FIRST SCALE:** neither the literal
> geodesic steps nor their unit directions admit positive scalar weights
> satisfying zero drift and tangent-space isotropy simultaneously.

The targeted verifier passes `14/14` in about 1.4 seconds.  It reconstructs
the geometry independently at 80 decimal digits from the exact 600-cell
coordinate alphabet.  No full suite was run.

## Complete preregistered census

The preregistered attempt count was exactly two.

| variant | unique zero-drift weights | covariance eigenvalue ratio | normalized traceless residual | simultaneous-system residual |
|---|---|---:|---:|---:|
| literal geodesic steps | `(0.055054, 0.204799, 0.444934, 0.295213)` | 3.40584 | 0.251585 | 0.00139230 |
| unit directions | `(0.122296, 0.238997, 0.352357, 0.286349)` | 6.07313 | 0.367250 | 0.275381 |

For both variants,

\[
\operatorname{rank}(A)=5,
\qquad
\operatorname{rank}([A\mid b])=6.
\]

Thus the equations are inconsistent even over unrestricted real weights.
The positivity restriction is not what causes the failure.

The literal-system least-squares residual looks numerically smaller because
the moment equations contain squared short geodesic lengths whereas the
normalization equation has unit scale.  It is not a near hit: the
scale-independent covariance ratio is `3.41`, its normalized traceless part
is `0.252`, and the rank mismatch is discrete.

## Independent geometric obstruction

There is a useful four-vector lemma.  Let `V` contain four spanning vectors
in a three-dimensional tangent space and suppose positive normalized weights
`p` satisfy

\[
Vp=0,
\qquad
V\operatorname{diag}(p)V^T=cI_3.
\]

Set

\[
W=V\operatorname{diag}(\sqrt p).
\]

Then `WW^T=cI_3`, while `W sqrt(p)=0`.  Because `W` has a one-dimensional
kernel,

\[
W^TW=c\left(I_4-\sqrt p\sqrt p^{,T}\right).
\]

Consequently every off-diagonal inner product must have the same value:

\[
v_i\mathbin{\cdot}v_j=-c
\qquad(i\ne j).
\]

The regular-tetrahedron control has all six values exactly `-1/3`.  The H4
data do not:

- literal-step off-diagonal spread: `0.00653995`;
- unit-direction off-diagonal spread: `0.809355`.

This independently explains the linear-system failure: the four local H4
directions do not form the required simplex tight frame.  The obstruction is
stable between double precision and the 80-decimal reconstruction.

## Calibration and numerical robustness

The regular-tetrahedron known-answer control returns unique weights

\[
p_0=p_1=p_2=p_3=\frac14,
\]

zero drift, covariance `I_3/3`, and an 80-decimal simultaneous residual of
`2.85e-84`.

For H4, the double- and 80-decimal simultaneous residuals agree to the shown
digits:

\[
1.392298144\times10^{-3}
\]

for literal steps and

\[
2.753806046\times10^{-1}
\]

for unit directions.  The relevant condition numbers are at most `30.4`, so
the negative is not an ill-conditioning artefact.  The exact-alphabet and
ordinary geometry step coordinates agree within `5.5e-12`.

The labelled four-step Gram matrix was also checked over all 14,400 chambers;
its maximum spread is `3.75e-12`.  Chamber 0 was therefore not an exceptional
or selected location.

## Physical meaning

This closes the cheapest repair of the anisotropic connected walk:

- scalar probabilities cannot simultaneously cancel drift and yield a local
  isotropic diffusion tensor;
- normalizing away the unequal step lengths makes the result worse, not
  better;
- adding a stay-put probability only rescales the same tensor and cannot
  repair its shape.

The result does **not** close quantum dynamics on the 600-cell.  A scalar
probability vector is not automatically realizable as the outgoing
distribution of a unitary four-state coin acting on a maximally mixed state;
that state remains maximally mixed under every such coin.  Conversely,
matrix-valued, position-dependent coins can carry directional Clifford/tetrad
data that this scalar feasibility problem does not contain.

The next admissible route is therefore a genuine matrix-valued discrete
tetrad/coin.  Its coefficients must be selected by the chamber metric and H4
equivariance, not fitted to isotropy or to a desired speed.  Before building
one, its complete equivariant parameter space and uniqueness must be counted;
otherwise the extra matrices merely hide the same fitting freedom.

## Status ledger

- **DERIVED CONTROL:** both double and 80-decimal solvers recover the regular
  tetrahedral tight frame.
- **DERIVED NUMERICAL:** the zero-drift weights are unique and positive for
  both preregistered H4 variants.
- **DERIVED NUMERICAL:** neither set of weights is isotropic.
- **DERIVED NUMERICAL:** both complete unconstrained systems have augmented
  rank larger than coefficient rank.
- **DERIVED STRUCTURAL LEMMA:** a positive four-vector tight frame with zero
  drift requires all six off-diagonal Gram entries to agree.
- **DERIVED NUMERICAL:** both H4 variants violate that necessary equality.
- **CLOSED:** first-scale scalar reweighting of the four fixed directions.
- **OPEN:** existence and selection of an H4-equivariant matrix tetrad/coin.
- **NOT CLAIMED:** a Dirac continuum limit, physical light speed, mass or
  Planck scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_h4_local_tight_frame.py
```

Expected result: `14/14`.
