# Prior-art reconciliation: two inequivalent zero-strut mechanisms

Date: 2026-08-20

Status: **source-level audit; no new numerical evaluation**.

This note attacks a potentially misleading comparison between the exact
zero-lapse boundary derived in this repository and the null-strut endpoints
reported in older 600-cell Regge cosmologies.  Similar terminology is not an
identity of mechanisms.

## 1. Repository object and complete hypotheses

The repository result
`gravity_600cell_dust_regular_lapse_identity_result.md` concerns the exact
regular **static product** family

```text
L_old = L_new = L,
tau = sqrt(rho) > 0,
spatial squares = L^2,
same-vertex pole square = -tau^2,
different-vertex cross square = L^2-tau^2,
M = (90/pi) epsilon3 L.
```

On this family all unrestricted internal equations vanish and

```text
p_pre = -epsilon3 L tau/4,
p_post = +epsilon3 L tau/4.
```

The momentum homotopy therefore forces

```text
tau(lambda) = tau0 (1-2 lambda),   0 <= lambda < 1/2.
```

At `lambda=1/2` the temporal separation itself vanishes while the spatial
scale does not change.  The homotopy parameter is not physical time.

## 2. De Felice--Fabri dust endpoint

De Felice and Fabri evolve a genuinely changing spatial 600-cell with a
Sorkin algorithm.  Their stopping diagnostic is a spatial contraction rate:

```text
Delta l / tau approximately 1,
```

and, in their generalized calculation, the ratio grows beyond one.  In their
five-dimensional embedding the endpoint is characterized by
`|dR/dv| -> 1`; extending farther would change the induced metric from
Lorentzian to positive definite.  This is a dynamical causality-breaking
boundary reached at nonzero spatial radius.

Primary sources:

- A. De Felice and E. Fabri, *The Friedmann universe of dust by Regge
  Calculus: study of its ending point*,
  <https://arxiv.org/abs/gr-qc/0009093>.
- A. De Felice and E. Fabri, *Singularities of the closed RW metric in Regge
  Calculus: a generalized evolution of the 600-cell*,
  <https://arxiv.org/abs/gr-qc/0106077>, especially Sections 2, 3 and 5.

## 3. Liu--Williams vacuum-Lambda endpoint

For a Collins--Williams block with lower/upper equilateral edge lengths
`l_i,l_(i+1)` and coordinate separation `delta t_i`, Liu and Williams derive

```text
m_i^2 = ((3/8) dot(l)_i^2 - 1) delta(t)_i^2,
dot(l)_i = (l_(i+1)-l_i)/delta(t)_i.
```

Their strut becomes null at

```text
dot(l)^2 = 8/3,
```

while the model is dynamically expanding.  The endpoint is delayed as the
spatial carrier is refined.  Their matter content is a positive cosmological
constant, not the fixed point-dust action used in the repository.

Primary source: R. G. Liu and R. M. Williams, *Regge calculus models of the
closed vacuum Lambda-FLRW universe*,
<https://arxiv.org/abs/1501.07614v2>, Section III A (the equations labelled
`strut` and `diag`) and the parent/child evolution discussion.

## 4. Decisive invariant comparison

On every nondegenerate member of the repository's exact family,

```text
Delta L / tau = 0,
pole square = -tau^2 < 0.
```

The older dynamical endpoints instead require a nonzero scale velocity that
drives the strut norm to zero.  Consequently:

- **DERIVED:** the repository's `lambda=1/2` boundary and the published
  dynamical null-strut endpoints are not the same mechanism under the stated
  hypotheses.
- **STRUCTURAL:** both are degenerations of a temporal edge, which explains
  the verbal resemblance but supplies no physical equivalence.
- **OPEN:** a future nonstatic repository trajectory could encounter the
  published velocity-driven null condition.  It must be tested using its
  actual `Delta L/tau`; it cannot be inferred from the static momentum
  homotopy.

## 5. Consequence for the programme

The old papers do not provide a no-go theorem for Regge evolution.  They show
that fixed coarse carriers can develop causal artifacts and that refinement
changes or delays them.  They therefore strengthen, rather than replace, the
current requirement: compare the unrestricted canonical/Jacobi dynamics on
canonically related spatial refinements before assigning physical meaning to
any fixed-600-cell soft mode, endpoint or speed.

No claim of external novelty follows from this audit.
