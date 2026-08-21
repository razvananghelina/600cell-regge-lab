# First direct-action adversarial result: two auxiliary controls fail

Date: 2026-08-21

Status: frozen formal control failure; no correction has been applied.

Protocol commits: `8c16996`, `b4916f8`, `a9ee74b`.  Verifier registration:
`7de8f9e`.

The first complete direct-action run produced

```text
reproducible/gravity_600cell_refined_h4_constrained_response_adversarial.json
SHA-256 a23ef4cc23d08ad8768f1df66789aa900cdb95a7f3529486df80697a53b1fe81
tests 15/17
outcome ADVERSARIAL_REFINED_H4_CONSTRAINED_RESPONSE_CONTROL_FAILED
```

No full suite or deferred nonlinear census was run.

## What passed

The mechanically independent route reconstructed a complete `20 x 20`
restricted second variation for every schedule from 210 scalar-action
directions and 55,610 scalar-action evaluations in total.  It did not use a
primary Hessian block or primary internal lift.

- **TENTATIVE DERIVED COMPUTATIONAL:** all 24 direct reductions lie in one
  class at direct response envelopes between `1.3984066e-54` and
  `1.3984072e-54`.
- **TENTATIVE DERIVED COMPUTATIONAL:** direct time reversal is covariant; its
  maximum used fraction is `1.1563e-16`.
- **TENTATIVE DERIVED COMPUTATIONAL:** all 24 direct matrices match the frozen
  analytic-gradient primary matrices; the maximum used cross-method fraction
  is `6.6667e-5`.
- The two complete 220-digit repeats agree with the 180-digit route; maximum
  used fraction `6.8530e-91`.
- All direct internal blocks are positive with global minimum eigenvalue
  `1.3780099e-5`; the maximum normalized solve residual is `1.7954e-181`.
- The known polynomial and wrong-polarization controls pass.  A deliberate
  response corruption is resolved by `1.9085e24` comparison envelopes.

These are diagnostic facts only.  Under the frozen outcome ordering they do
not license acceptance while either auxiliary control fails.

## Failure 1: direct stationarity threshold

The preregistered two-step fourth-order scalar-action derivative required all
240 internal components below `1e-60`.  The maximum was

```text
5.16616600924884e-54.
```

This is seven orders above the frozen gate.  A plausible explanation is the
`O(h^4)` truncation of the Richardson derivative at `h=(1e-15,5e-16)`, with a
large fifth derivative.  That explanation has not yet been tested and is
therefore **OPEN**.  It must be checked using a preregistered multi-level
scalar-action ladder; the threshold may not be relaxed retroactively.

## Failure 2: overstrong off-shell curvature reality gate

Across all perturbed action evaluations the recorded extrema were

```text
angle identity residual          9.0444e-178
relative imaginary total action  1.1259e-176
minimum angle argument            0.9122203
imaginary individual curvature   2.0750e-8
```

Only the last quantity violates the frozen `<1e-80` branch gate.  Individual
Lorentzian deficit-angle expressions can become complex under arbitrary
off-shell log-edge perturbations even when the complete action stays real and
the analytic angle identity stays on a safe branch.  This is a plausible
framing error in the control, not yet a result: it is **OPEN** until a
step-halving and conjugation/cancellation diagnostic shows that the imaginary
individual terms vanish continuously at the real background while the total
action remains real.  If they do not, the direct audit remains killed by a
branch ambiguity.

## Frozen verdict and next action

**FORMAL DERIVED NEGATIVE FOR THIS EXECUTION:** the first adversarial protocol
did not pass and the primary single-class result is not yet accepted.

Preregister a small diagnostic before inspecting finer stationarity or
off-shell-curvature ladders.  It must test the expected derivative convergence
and the scaling/cancellation of individual imaginary curvatures.  Only if
that diagnostic identifies both failures as control-framing/truncation errors
may a separately preregistered corrected adversarial verifier be run.

