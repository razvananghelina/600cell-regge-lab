# Prior-art gate: regular-lapse momentum identity for the 600-cell dust slab

Date: 2026-08-16

Status: **completed before evaluating the candidate identity away from the
already committed continuation points**.

Upstream continuation result: `fcfe0c9` (provenance update `cc52549`).

This gate maps primary literature; it does not establish external novelty.

## 1. Exact object and complete hypotheses

Retain the derived 600-cell boundary, the two order-24 schedule parities, the
staircase triangulation of the slab, the corrected Lorentzian angle branch,
the Schlaefli gradient and the De Felice--Fabri dust term used by the committed
canonical-rank and continuation verifiers.

Let `L=ARB_L0`, keep every old and new spatial edge square equal to `L^2`, and
define the one-parameter regular family for `rho>0` by

```text
old spatial square       = L^2,
new spatial square       = L^2,
staircase diagonal square= L^2-rho,
pole edge square         = -rho.
```

The implementation stores the positive pole magnitude `rho`; its edge
Jacobian inserts the Lorentzian minus sign.  Write `tau=sqrt(rho)>0`.

The candidate identity is the conjunction

```text
all 35 internal Regge+dust equations vanish on the family,
p_pre(rho)  = sqrt(rho/rho0) p_pre(rho0),
P p_post(rho) = -p_pre(rho),
```

for both derived schedule parities on a precisely certified connected
Lorentzian domain.  `P` is the already derived old-to-final orbit map.  The
claim concerns canonical momentum versus proper lapse.  It does not make the
continuation parameter `lambda` physical time and does not assert expansion.

## 2. KNOWN

Sorkin evolution of the 600-cell and dust-filled regular cosmologies predates
this repository.  De Felice and Fabri analyze the dust 600-cell evolution and
its causal endpoint:

- <https://arxiv.org/abs/gr-qc/0009093>;
- <https://arxiv.org/abs/gr-qc/0106077>.

Collins--Williams-type regular-polytopal reductions and their Hamiltonian and
evolution equations are known.  Liu and Williams analyze closed Regge FLRW
models and report abrupt endpoints when timelike struts become null:

- <https://arxiv.org/abs/1501.07614>.

Tsuda and Fujiwara derive a continuum-time regular-polytopal Regge action with
time-reparameterization invariance; strut variation gives a Hamiltonian
constraint:

- <https://arxiv.org/abs/2109.01075>.

The action-generated pre/post momentum framework and the replacement of exact
constraints by pseudo-constraints in curved discrete gravity are known:

- <https://arxiv.org/abs/1108.1974>;
- <https://arxiv.org/abs/0912.1817>;
- <https://arxiv.org/abs/0905.1670>.

Therefore a regular 600-cell cosmology, a lapse/strut equation, a null-strut
endpoint and action-generated canonical momenta are not new.

## 3. CONTROL

**DERIVED UPSTREAM:** the published regular datum satisfies the complete
internal equations and has a calibrated rank-65 pre-Legendre Jacobian in both
schedule parities.

**DERIVED UPSTREAM:** the frozen nonlinear continuation accepted 41 points in
each parity, with unchanged spatial boundary and no detected Lorentzian branch
loss.

**PATTERN UPSTREAM:** those points obey

```text
rho(lambda)/rho0 = (1-2 lambda)^2
```

to the inherited target uncertainty.  This finite numerical pattern is the
motivation for the present mission, not evidence for its analytic conclusion.

## 4. OPEN difference

No located source gives the exact finite-slab identity above for this
staircase carrier, corrected angle convention, dust normalization, two
schedule parities and canonical boundary derivative.  Whether the entire
regular family is an exact stationary lapse family and whether its boundary
momentum is exactly linear in `tau` are **OPEN**.

External novelty remains **OPEN** even if the identity is proved.  A proof
could be a special instance of a more general static-prism factorization not
located by this targeted search.

## 5. Framing attack

Even a successful identity would not be a cosmological tick.  It would show
that the observed continuation moves along a lapse family of one fixed
spatial geometry.  A zero-lapse endpoint would then be a degenerating
canonical parametrization, not a Big Bang, maximum expansion or emergent
clock.

Conversely, failure of the identity at generic `rho` would not erase the
committed local continuation.  It would demote the quadratic law to a local
numerical pattern and require identifying the first nonzero correction.

