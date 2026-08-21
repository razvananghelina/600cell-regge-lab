# Prior-art gate: generic-velocity relational composition

Date: 2026-08-21

Status: written after the adversarially corroborated turning-point half-step
no-go, before expanding the action at nonzero velocity and before evaluating
any mass/velocity branch.

## 1. Exact object and hypotheses

Retain the certified homogeneous cellular 600-cell action, positive
Lorentzian branch, zero cosmological constant and conserved dust mass.  Use
global scale covariance to set the incoming spatial scale to one, but do not
freeze the dimensionless mass ratio to its static value.  Write

```text
mu = M/L_minus > 0,
tau = s*e,
log(L_plus/L_minus)=s*v*e+O(e^2),
rho=tau^2[1+O(e)].
```

Here `v` is a fixed nonzero dimensionless logarithmic velocity and `s=1` or
`1/2` labels a coarse or nominal half interval.  The same physical state
means the same incoming scale, canonical momentum, mass ratio and velocity;
only the interval factor `s` changes.

The first question is exact and leading-order only: does the full cellular
action reduce to a reparametrization-covariant principal function

```text
S = s*e*L0(v,mu)+O(e^2),
```

whose lapse constraint and incoming momentum are independent of `s` when
the endpoint displacement is scaled by the same `s`?  This is a necessary,
not sufficient, condition for a generic-velocity same-state half-step.

## 2. Why this is not the preceding test

At the static mass ratio the physical continuum velocity is zero.  The
accepted coarse boundary momentum is then `O(e)`, and the preceding exact
test showed that its half-step lapse and momentum equations disagree at
`O(e)`.  A nonzero-velocity state instead has endpoint displacement `O(e)`
and canonical momentum `O(1)`.  It is a different asymptotic stratum and is
not covered by the `log L=O(e^2)` no-go.

Choosing one convenient numerical mass ratio would introduce a hidden
look-elsewhere parameter.  The new calculation must retain symbolic `mu` and
`v`, eliminate one only through the exact lapse constraint, and state the
complete real-domain hypotheses before evaluating controls.

## 3. Primary prior art

- Reparametrized discrete systems and improved/perfect actions are developed
  by Bahr and Dittrich, arXiv:`0907.4323`; their continuum/perfect limit
  restores the time-reparametrization symmetry broken by a generic finite
  discretization: <https://arxiv.org/abs/0907.4323>.
- Exact discrete Lagrangians generate fixed-time Hamiltonian flow and are the
  reference for variational error analysis in Marrero, Martin de Diego and
  Martinez, arXiv:`1608.01586`:
  <https://arxiv.org/abs/1608.01586>.
- Canonical simplicial pre/post data and action composition are established
  by Dittrich and Hoehn, arXiv:`1108.1974`:
  <https://arxiv.org/abs/1108.1974>.
- Brown--Kuchar dust supplies an explicit relational proper-time coordinate
  and conjugate momentum in the continuum theory, arXiv:`gr-qc/9409001`:
  <https://arxiv.org/abs/gr-qc/9409001>.
- Regular 600-cell dust cosmology and its continuum comparison are known from
  the Collins--Williams and De Felice--Fabri line, including
  <https://arxiv.org/abs/gr-qc/0009093>.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- The exact cellular action is homogeneous of length degree two and common
  to all homogeneous staircase schedules.
- Its static mass-balanced line has arbitrary lapse and zero physical
  velocity.
- The turning-point same-state analytic half-step is absent.
- Continuum generally covariant dynamics is locally reparametrization
  invariant; generic finite Regge discretizations can break that symmetry.

### CONTROL

- Derive the leading action directly for symbolic `(s,v,mu)`.
- Verify its `s` dependence exactly, not at sampled values.
- Differentiate the leading action using the correct fixed-endpoint rule for
  the lapse equation; differentiating at fixed `v` would be wrong because
  `v=Delta/tau` changes under a lapse variation.
- Derive the canonical pre-momentum independently from the complete action.
- Recover the static normalization only as the `v->0` control after the
  generic formulas are committed.

### OPEN

- Whether the leading constraint and momentum are interval-factor
  independent.
- The exact mass--velocity constraint and its real branches.
- Next-order same-state composition, where finite-discretization error first
  appears.
- Generic anisotropic evolution, spatial refinement and external novelty.

## 5. Framing attack

A positive leading theorem is expected continuum kinematics, not a new law
of gravity and not evidence for a fundamental tick.  It only shows that the
turning-point obstruction does not occur one order earlier on the generic-
velocity stratum.  The decisive physics remains the next-order
one-versus-two comparison.

A negative leading theorem would be more serious: the cellular action would
fail even the local reparametrization-covariant principal-function form away
from the turning point.  It would direct the programme immediately to an
improved/perfect or explicitly deparameterized dust-clock action.

External novelty is **OPEN** regardless of outcome.

