# Connected canonical continuation of the 600-cell dust slab: result

Date: 2026-08-16

Prior-art gate: `52a6d50`

Frozen protocol and pre-evaluation clarifications: `393e528`, `c2e942d`,
`cf00b38`

Implementation registered before evaluation: `1a6f796`

Reproduction-only calibration failure: `1ea960b`

Upstream-uncertainty amendment and implementation: `8cab22c`, `c80e434`

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_canonical_continuation.py`

Artifact:
`reproducible/gravity_600cell_dust_canonical_continuation.json`

Artifact SHA-256:
`808b5080648428eee221fcbc6399a8f182c062fcb5d8db07708eebf23489e7a2`

Result/artifact commit: `fcfe0c9`

## Verdict

**DERIVED COMPUTATIONAL:** the frozen verifier passes `8/8`.  Both derived
schedule parities return

```text
CANONICAL_CONTINUATION_NUMERICALLY_OPEN
```

The connected branch was followed through all coarse points up to
`lambda=20/64=0.3125`.  The next coarse target `lambda=21/64=0.328125`
reached the preregistered 30-iteration Newton limit.  Exactly 20 mandated
bisections then placed the last accepted root at

```text
lambda = 0.32812498509883880615234375,
upper unresolved target = 0.328125,
gap = 1.490116119384765625e-8.
```

No `lambda=1` endpoint was reached.  Therefore no canonically selected next
frame, expansion or contraction has been established.

## Accepted branch

**DERIVED COMPUTATIONAL:** all 41 accepted records (reproduction, coarse
continuation and bisections) have the same spatial boundary, within the
registered arbitrary-precision controls.  Across the entire accepted path,

```text
max |Delta log spatial length| = 3.134e-53.
```

At the last accepted point,

```text
rho/rho0                    = 0.1181640829890975297181739071
Delta log spatial length   = 1.120e-60
reduced residual infinity  = 5.374e-56
full residual infinity     = 9.543e-40
minimum leading minor      = 1.2293791194e-5
minimum angle argument     = 0.9984530285
reduced Jacobian s_min     = 1.7240946191e-10
combined Jacobian error    = 6.1743452304e-21
s_min / error              = 2.792e10
```

Thus the unresolved endpoint is not accompanied by a certified Jacobian
rank loss, Lorentzian branch loss, vanishing principal minor or angle cusp.
The frozen outcome is numerical openness, not physical branch termination.

**DERIVED COMPUTATIONAL:** the even and odd states, momenta, accepted
`lambda` values and `rho/rho0` values agree exactly at the stored precision.
Their largest reported scalar diagnostic difference is below `5.6e-96`.

## Simple lapse law

**PATTERN:** on every accepted point the solution obeys

```text
q_new(lambda) = q_new(0),
rho(lambda)/rho(0) = (1 - 2 lambda)^2
```

to the inherited reproduction uncertainty.  The largest absolute deviation
from the quadratic law is `6.543e-21`, occurring at the reproduction point;
the last-point deviation is `7.732e-22`.  The explicit lapse-relation
residual and spatial-scale displacement are at most about `3.14e-53` after
reproduction.

This is strong evidence that the homotopy changes the common lapse while
leaving the spatial slice fixed.  It is not yet an analytic identity.  In
particular, `lambda` is a momentum homotopy parameter, not physical time.
Extrapolating the observed law to zero lapse at `lambda=1/2` is **OPEN** and
is not licensed by this run, which stopped near `lambda=0.328125`.

## Why the failed coarse target is not a no-go

At the failed coarse attempt the residual fell from `3.465e-5` to
`7.115e-11`, but the fixed Armijo search settled at damping `2^-9` and did
not meet `1e-50` within 30 accepted iterations.  The subsequent frozen
bisections converged in two Newton iterations arbitrarily close to that
target.  This behavior, together with the resolved Jacobian and healthy
branch margins, gives no evidence of a geometric endpoint at `21/64`.

The protocol correctly prevents a post-result retry of the upper target
from a better seed from being called the preregistered result.  A later
method may diagnose or cross the numerical boundary, but it must be frozen
separately.

## Framing correction

**DERIVED COMPUTATIONAL:** the canonical inverse response refutes the claim
that the homogeneous seam itself asks for one Newton step in spatial scale.
On the connected solution found here the new spatial boundary remains fixed;
the varying coordinate is the collective lapse.  The seam is a momentum
covector, and applying the inverse pre-Legendre Jacobian is essential before
interpreting its direction in configuration space.

Forcing the lapse to remain fixed would add a clock/gauge condition and
define a different problem.  Adding a scale impulse would change the initial
canonical datum.  Neither operation can be described as the evolution of the
same preregistered state without a new physical selection rule.

## Post-result primary-source audit

The result remains inside known canonical-Regge structure:

- Dittrich and Hoehn formulate action-generated canonical simplicial
  evolution and a priori/a posteriori data in arXiv:1108.1974.
- Dittrich and Hoehn show that nonlinear curvature replaces exact
  constraints by background-dependent pseudo-constraints in
  arXiv:0912.1817.
- Bahr and Dittrich relate broken discrete gauge symmetry to
  pseudo-constraints in arXiv:0905.1670.
- Gambini and Pullin describe consistent discretizations in which lapse-like
  variables can be determined by the discrete equations in
  arXiv:gr-qc/0511096.

No located primary source supplies the exact quadratic lapse law or the
connected branch for this order-24 carrier, action, schedule pair and
published datum.  External novelty remains **OPEN**; the general phenomenon
that discretization can determine would-be gauge data is not new.

## Next falsifiable step

The next clean mission is analytic, not a looser continuation run:

1. derive the regular-lapse family directly from the exact symmetric Regge
   action;
2. compute its pre-momentum as a function of `rho`;
3. prove or refute `p_pre(rho) proportional sqrt(rho)` with all branch and
   sign hypotheses stated;
4. only if that identity is established, determine whether the target
   homotopy reaches zero lapse at `lambda=1/2` and whether an opposite
   temporal-orientation convention supplies a nonzero forward datum.

Expansion, a dust proper-time clock, a Friedmann continuum limit, shape-mode
propagation, `c`, Planck units and particle masses remain **OPEN / NOT
TESTED**.
