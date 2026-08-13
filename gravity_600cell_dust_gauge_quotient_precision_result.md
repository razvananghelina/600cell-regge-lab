# Precision-corrected gauge quotient: frozen result

Date: 2026-08-13

Prior-art gate: `44bd4cf`

Frozen protocol: `da34272`

Implementation commit: `979e518`

Mechanical typo correction: `511fce8`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_gauge_quotient_precision.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_gauge_quotient_precision.json`

Targeted run: **14/14 implementation checks passed**.  Per the current work
rule, the full suite was not run.

## 1. Result

Both five-stage schedule parities receive the frozen outcome

```text
REGULAR_QUOTIENT_29_ZERO_SUM_RESPONSES.
```

Thus the previous odd-parity `OPEN NUMERICALLY` label is resolved by the
separately preregistered precision correction; it is not retroactively
changed.  On the restricted order-24 carrier, at this stationary dust
sandwich and to linear order:

- one exact collective lapse direction is gauge;
- the quotient internal Hessian has rank 34;
- one final-boundary compatibility equation removes the homogeneous scale;
- the remaining 29 zero-sum boundary shape directions have a rank-29 linear
  response in both schedules.

This is **DERIVED COMPUTATIONAL LINEAR**.  The identification of the row with
the homogeneous all-ones scale direction has **PATTERN / cross-resolution
confirmed provenance**, because that target was learned from the first run.
Finite numerical equality is not an analytic proof of component equality.

## 2. Independent boundary-row refinement

The new calculation evaluated 360 complete-action points per schedule at 100
decimal digits, adding the preregistered step `1.25e-4` to the old `5e-4` and
`2.5e-4` levels and forming the frozen sixth-order extrapolation.

The two parities give the same reported values:

```text
norm(c6)                 = 1.112196291297975e-3
sixth-order error        = 2.289244841840344e-9
norm(c6)/error           = 4.858354471179548e5
mean component           = 2.030583323724958e-4
component spread         = 1.103354449363651e-11
norm(nonuniform part)    = 2.0543233456e-11
nonuniform/error         = 8.9738035e-3
cosine with all-ones     = 1 within binary64 rounding
```

The nonuniform component is only `0.009` of the independent extrapolation
error, so the frozen classification is

```text
ONE_BOUNDARY_CONSTRAINT,
UNIFORM_WITHIN_FROZEN_ERROR.
```

This is a much stronger cross-resolution result than the first fourth-order
row, whose apparent nonuniformity was `2.143e-5` relative.  It remains honest
to say “uniform within the frozen numerical error,” not “analytically equal.”

All 720 displaced configurations remained Lorentzian and away from branch
boundaries.  Maximum imaginary contamination was `7.45e-95`.

## 3. The precision correction is perturbatively small

The correction used no fitted eigenvector or Schur coefficient.  It combined:

1. the robust recorded `30 x 30` regular block;
2. the independently reconstructed 80-decimal relative five-pole Schur form;
3. the exact collective tangent from the constant-action lapse family;
4. unique minimum-Frobenius updates enforcing the two exact compatibility
   identities.

For even/odd schedules respectively:

```text
relative coupling correction  = 4.424e-8 / 4.681e-8
relative full-H correction     = 7.536e-12 / 7.974e-12
relative mixed-B correction    = 7.844e-12 / 7.844e-12
normalized exact-null residual = 2.818e-21 / 1.789e-21
```

Both corrected full Hessians have exactly 34 singular values above the
absolute `1e-9` gate and one collective null.  Their four weak quotient
singular values are all about `4.60497e-8`, agreeing with the independent
Schur reconstruction.  This refutes the concern that resolving the odd
schedule required a large or target-driven matrix alteration.

## 4. Linear response and its warning

For each schedule, the exact zero-sum boundary basis has 29 columns.  The
response has rank 29 at all frozen relative thresholds `1e-7`, `1e-9` and
`1e-11`.

```text
                         even              odd
response norm            4.742934e5         4.742931e5
condition                 3.606609e6         3.619108e6
quotient residual         3.363e-11          3.848e-11
full corrected residual   8.997e-11          9.153e-11
```

Almost the entire response norm lies in the four relative-pole sector.  This
large amplification is a physical/numerical warning: the linear derivative
exists, but the four pseudo-constraint curvatures are so soft that a small
boundary perturbation can require a very large internal displacement.  The
linear result therefore does not by itself supply a useful finite evolution
step.

## 5. What is and is not new

The broad mechanism is **STRUCTURAL / KNOWN**.  Homogeneous variation and a
Hamiltonian or initial-value constraint are standard in closed Regge--FLRW
models.  The present result identifies that structure explicitly in the
dust-filled, locally varied 600-cell staircase.

No located primary source prints this exact mixed row, precision-corrected
quotient or 29-direction response on this carrier.  External novelty of that
explicit realization is **OPEN**, not claimed.

## 6. Status ledger

| Claim | Status |
|---|---|
| Collective lapse is an exact null family | **DERIVED COMPUTATIONAL, upstream** |
| Four relative pole modes are nonzero | **DERIVED COMPUTATIONAL** |
| Both quotient Hessians have rank 34 | **DERIVED COMPUTATIONAL LINEAR** |
| Exactly one nonzero boundary compatibility equation occurs | **DERIVED COMPUTATIONAL LINEAR** |
| Its direction is homogeneous scale within the frozen error | **PATTERN / cross-resolution confirmed** |
| Its thirty components are analytically identical | **OPEN** |
| The compatible tangent space is the 29-dimensional zero-sum space | **DERIVED COMPUTATIONAL under the numerical-uniformity gate** |
| The response has rank 29 in both schedules | **DERIVED COMPUTATIONAL LINEAR** |
| Small finite boundary perturbations admit nearby nonlinear roots | **OPEN** |
| The response describes physical gravitons | **NOT ESTABLISHED** |
| The result covers the full 840-edge carrier | **NOT TESTED** |
| A clock, Planck scale or multi-tick dynamics is selected | **NOT ESTABLISHED** |

## 7. Next falsification test

The next meaningful step is not another Hessian spectrum.  It is a
preregistered nonlinear continuation along deterministic zero-sum boundary
directions.  For each amplitude, solve the 34 quotient equations with one
collective lapse gauge condition and then evaluate the omitted compatibility
equation independently.

The route advances only if complete 35-equation residuals converge and the
solution scales linearly as the amplitude tends to zero.  Because the linear
response norm is about `4.74e5`, amplitudes must be fixed from the certified
linear-neighborhood scale before any solve; choosing them after convergence
would be fitting.  Failure would show that the present tangent response is a
singular linear artifact rather than usable local evolution.

## 8. Post-result coordinate audit

The first post-result base audit compared quantities in inconsistent
coordinates: raw derivatives from the published residual table were projected
onto the logarithmic Hessian without multiplying by their squared-length
coordinates.  Its claimed odd correction `8.675e-3` is **RETRACTED**.

The coordinate-consistent calculation is:

```text
                         even            odd
norm(log residual)       2.576e-13       2.729e-12
norm(transverse part)    2.576e-13       2.729e-12
linear base correction   2.521e-10       9.019e-7
```

Both corrections lie below the preregistered `1e-5` weak-scale tolerance.
Consequently the printed point remains an adequate base for the frozen linear
matrix theorem.  A smaller-step action-only audit is still required before
claiming a stationary family to precision below the four soft curvatures, but
there is no evidence here for an order-`1e-2` displacement of the odd base.
