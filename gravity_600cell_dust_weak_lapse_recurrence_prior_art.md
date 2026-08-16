# Prior-art gate: weak-lapse recurrence of the canonical dust map

Date: 2026-08-16

## Exact object and hypotheses

The carrier, action, fixed dust mass, branch convention and order-24 schedule
parities are exactly those of the three accepted homothetic ticks.  The
published static identity holds for every admitted positive proper duration
`tau`, so introduce a control parameter

```text
tau_static(lambda) = lambda*tau0,
rho_static(lambda) = lambda^2*rho0,
lambda in {1/2,1/4,1/8}.
```

Keep `L0` and

```text
M=(90/pi)*epsilon3*L0
```

fixed.  At each `lambda`, start from the exact static output momentum
`+lambda*k`, where `k=epsilon3*L0*tau0/4`, and apply the same reduced canonical
map for three successive slabs:

```text
(a_n,p_n) -> (a_(n+1),p_(n+1)),
```

with the pole equation determining the slab lapse.  No tick-four target or
experimental datum is used.

The object measured is the **small-lapse asymptotic law of the fixed-carrier
canonical map**, not a refinement of the spatial triangulation.

## Why lapse, not momentum, is the expansion parameter

On the exact static family, boundary momentum is proportional to `tau`, while
the first three scale and relative-lapse logarithms are of order `tau^2`.
Therefore a Taylor expansion at fixed published lapse and varying momentum is
not the relevant continuum-like jet.  The natural dimensionless expansion is
in `lambda`, with geometric corrections expected in even powers beginning at
`lambda^2`.

## Primary prior art

- Gambini and Pullin formulate consistent discretized Regge dynamics as a
  canonical transformation.  This supplies the broad setting for studying a
  small-step canonical map:
  <https://arxiv.org/abs/gr-qc/0511096>.
- Dittrich and Höhn derive action-generated pre/post evolution and explain how
  discrete consistency conditions replace continuum constraints:
  <https://arxiv.org/abs/1108.1974>.
- De Felice and Fabri evolve the dust-filled 600-cell and compare it with a
  continuum Friedmann description:
  <https://arxiv.org/abs/gr-qc/0009093>.
- Dittrich, Gielen and Schander show that dust-filled Lorentzian Regge shells
  can approximate closed continuum cosmology, with discretization dependence:
  <https://arxiv.org/abs/2109.00875>.

## KNOWN / CONTROL / OPEN

- **KNOWN:** Regge actions can generate discrete canonical transformations,
  and 600-cell/dust cosmologies have multi-step evolutions.
- **CONTROL:** the exact static all-lapse identity supplies the initial state
  for every `lambda` without fitting.
- **CONTROL:** the accepted `lambda=1` three-tick trajectory fixes a connected
  branch and provides zero-order continuation seeds, not acceptance targets.
- **CONTROL:** all scaled solves must pass complete 35-equation, 30-momentum,
  branch, Jacobian and parity gates.
- **OPEN:** whether the normalized scale increments tend to `(1,2,3)`, lapse
  increments to `(1,3,5)`, and cumulative values to triangular numbers and
  squares as `lambda -> 0`.
- **OPEN:** the convergence order and leading coefficients.
- **OPEN:** any spatial-refinement or physical-continuum statement.

## Framing attack

Even a clean integer limit would not be a new fundamental law by itself.  A
constant-force/constant-acceleration symplectic step generically produces odd
momentum increments and triangular positions.  The result may explain the
observed sequence as ordinary local Friedmann kinematics rather than expose
new microscopic physics.

Shrinking `tau` on one fixed triangulation is not spatial refinement and does
not prove convergence to Einstein gravity.  It only identifies the temporal
jet of this reduced discrete map.  Likewise, `lambda` is an externally chosen
control and does not make time emergent.

A source search cannot prove novelty.  Only the new targeted verifier will be
run; the full suite will not be run.
