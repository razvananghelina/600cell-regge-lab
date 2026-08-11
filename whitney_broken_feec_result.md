# Broken FEEC repairs topology, but does not select the finite operator

Date: 2026-08-11

Preregistration commit: `4452a9d`

Targeted verifier:
`reproducible/verify_whitney_broken_feec.py`

Targeted result: **12/12 PASS**.  The verifier is registered exactly once.
The full suite was not run, by explicit user request.

## Headline

Both preregistered local conforming projections give a valid broken-FEEC
(CONGA) Hilbert complex on the three-dimensional duplicated Whitney carrier.
For every positive stabilization they recover precisely the harmonic kernel
of the conforming triangulated 3-sphere:

\[
 (b_0,b_1,b_2,b_3)=(1,0,0,1).
\]

> **DERIVED APPLICATION OF A KNOWN THEOREM:** topology is robust under the two
> frozen projection choices.

> **STRUCTURAL AMBIGUITY:** the projections, broken differentials and positive
> finite spectra differ after nontrivial refinement.  The present theory still
> does not select a unique finite operator.

This is a sharper result than either “the penalty fails” or “the flux solves
everything.”  A projection-based differential plus symmetric stabilization
solves the topology problem; it does not solve the physical selection problem.

## Prior-art correction

The mechanism is established numerical mathematics, not a newly invented
physical principle.  Campos Pinto and Güçlü define the broken differential

\[
 d_h=d_{\rm pw}P,
\]

its adjoint in the broken metric, and the stabilized Hodge Laplacian.  Their
kernel theorem gives the conforming harmonic fields for every positive
stabilization.  The primary reference is:

- Martin Campos Pinto and Yaman Güçlü, *Broken-FEEC discretizations and Hodge
  Laplace problems*, Mathematics of Computation,
  [DOI 10.1090/mcom/4085](https://doi.org/10.1090/mcom/4085),
  [arXiv:2109.02553v3](https://arxiv.org/abs/2109.02553).

The earlier circle calculation independently reconstructed the lowest-order
one-dimensional instance: equal conforming projection, projected derivative,
local broken adjoint and mismatch stabilization.  Its exact circle numbers
remain useful calibration, but the method itself must not be presented as a
new discovery.

## Operators tested

For both the counting and diagonal-Whitney recoveries, let

\[
 P_p=J_pL_p,
 \qquad
 d_{h,p}=D_{{\rm pw},p}P_p,
\]

and

\[
 d_{h,p}^*=M_p^{-1}d_{h,p}^TM_{p+1}.
\]

At the frozen numerical control `alpha=1`, the degree-`p` weak pencil is

\[
 \begin{split}
 K_p={}&d_{h,p}^TM_{p+1}d_{h,p}
 +(d_{h,p-1}^*)^TM_{p-1}d_{h,p-1}^*\\
 &+(I-P_p)^TM_p(I-P_p).
 \end{split}
\]

Every term is a positive Gram factor.  The value one is not a physical
coupling: every positive coefficient has the same exact kernel.

## Exact complex structure

For both projections, on

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}\partial\Delta^4),
 \qquad k=1,2,
\]

the verifier establishes:

\[
 P_p^2=P_p,
 \qquad
 P_pJ_p=J_p,
 \qquad
 L_pJ_p=I,
\]

and the oriented intertwiner

\[
 D_{{\rm pw},p}J_p=J_{p+1}d_p.
\]

Consequently,

\[
 P_{p+1}d_{h,p}=d_{h,p},
 \qquad
 d_{h,p+1}d_{h,p}=0.
\]

The largest observed output/nilpotency residual is `4.17e-17`.  The largest
broken-metric adjoint residual is `5.69e-14`, well inside the preregistered
`1e-11` gate.  All sixteen weak pencils are symmetric and positive
semidefinite.

Unlike the orthogonal conforming recovery audited previously, the adjoint
here uses the inverse **broken** mass.  That mass is block diagonal by
tetrahedron, so its inverse is local.  This is the precise mechanism missed by
the earlier global-adjoint framing.

## Why the kernel is exact

The stabilization term is zero exactly when

\[
 (I-P_p)u=0,
\]

so `u` lies in `im(P_p)=im(J_p)`, the conforming subspace.  On that subspace,
the verified intertwiner reduces the other two positive factors to the
ordinary conforming Whitney differential and coderivative.  Their common
kernel is the conforming harmonic space.

Both controls are subdivisions of the boundary of a 4-simplex and hence are
triangulations of `S^3`.  Their Betti numbers are exactly

\[
 (1,0,0,1).
\]

This maps every hypothesis of the published kernel theorem.  The numerical
low-spectrum calculation is a control, not the logical source of the kernel
claim.

| edgewise level | projection | observed degreewise nullities |
|---:|---:|---:|
| 1 | counting | `(1,0,0,1)` |
| 1 | diagonal Whitney | `(1,0,0,1)` |
| 2 | counting | `(1,0,0,1)` |
| 2 | diagonal Whitney | `(1,0,0,1)` |

## The ambiguity survives away from zero

At `k=1`, the two projections coincide because all copies relevant to one
global simplex have equal diagonal Whitney weights.  At `k=2`, they differ:

| differential degree | nonzero coefficient differences | maximum absolute difference |
|---:|---:|---:|
| 0 -> 1 | 0 | 0 |
| 1 -> 2 | 31,680 | 0.0666667 = `1/15` |
| 2 -> 3 | 3,840 | 0.15625 = `5/32` |

The first positive generalized Hodge eigenvalues at `k=2` are

| form degree | counting | diagonal Whitney | relative difference |
|---:|---:|---:|---:|
| 0 | 0.889571245449 | 0.889571245449 | below `1e-14` |
| 1 | 0.939262420044 | 0.935331908666 | 0.4185% |
| 2 | 0.875727775400 | 0.876595261866 | 0.0990% |
| 3 | 4.374103306795 | 4.369517911725 | 0.1048% |

The differences are small on this control, but one refinement cannot support
a universality or convergence claim.  They are not zero, and no candidate may
be selected because its values look preferable.

## What has and has not been repaired

Repaired:

- the duplicated local carrier can preserve the correct Betti modes at every
  positive stabilization;
- the differential, its Hilbert adjoint and the Hodge Laplacian remain local
  because only block-local mass inverses occur;
- a singular `kappa h -> infinity` limit is unnecessary for the exact kernel.

Not repaired:

- counting and diagonal-Whitney projections remain different natural inputs;
- the stabilized Hodge Laplacian is not shown to be the square of a uniquely
  derived odd first-order microscopic evolution;
- the positive finite spectrum is projection-dependent;
- no physical value or scale for the stabilization has been selected.

## Attack on the positive framing

The published paper calls the dual commuting sequence canonical once a
conforming projection `P` has been supplied.  It also describes equal
averaging as “a simple conforming projection.”  This does not prove that equal
averaging is the unique projection allowed by the geometry.  The exact
counterexample in `whitney_3d_flux_canonicity_result.md` remains decisive on
that narrower question.

Thus literature rescues locality and topology, not uniqueness.  Reporting
CONGA as a derivation of the universe's dynamics would still be an
overstatement.

## Status ledger

- **DERIVED:** both frozen `P` maps are conforming projections.
- **DERIVED:** both projected derivatives form complexes and have local
  broken-metric adjoints.
- **DERIVED APPLICATION OF KNOWN THEOREM:** both stabilized kernels have Betti
  multiplicities `(1,0,0,1)`.
- **DERIVED NUMERICAL:** the complete low-spectrum controls reproduce those
  nullities at `k=1,2`.
- **DERIVED NEGATIVE:** the two positive finite pencils are not identical.
- **STRUCTURAL POSITIVE:** topology is robust to the demonstrated projection
  ambiguity.
- **STRUCTURAL NEGATIVE:** finite-cutoff dynamics remains unselected.
- **OPEN:** whether the candidate differences vanish in a common continuum
  universality class.
- **OPEN:** a canonical odd first-order realization of the stabilization.
- **NOT CLAIMED:** Lorentzian time, causal speed, inertia, mass or Planck
  units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_broken_feec.py
```

Expected result: `12/12`.

