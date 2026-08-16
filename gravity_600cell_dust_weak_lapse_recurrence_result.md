# Weak-lapse recurrence of the canonical 600-cell dust map

Date: 2026-08-16

## Provenance

- prior-art gate: `90757a6`;
- frozen protocol: `1c26fdb`;
- registered implementation, committed before its first evaluation: `d7e5d6c`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_weak_lapse_recurrence.py`;
- artifact:
  `reproducible/gravity_600cell_dust_weak_lapse_recurrence.json`;
- artifact SHA-256:
  `500be1c4e2d7ec4104b9773bc1cfc71065c9d930607119eb616367d18fa5d8f9`.

Only this targeted verifier was run.  It returned **5/5**.  The full suite was
not run.

The artifact records `tick4_target_parsed=false`.  No fourth-tick target or
datum entered the calculation.

## Verdict

```text
WEAK_LAPSE_QUADRATIC_INTEGER_LAW
```

**DERIVED COMPUTATIONAL:** on the frozen fixed-mass contracting branch, all
18 non-static solves pass the complete action, momentum, branch, Jacobian and
parity gates.  For three consecutive steps and

```text
lambda in {1/2, 1/4, 1/8},
tau(lambda) = lambda*tau0,
```

all eleven preregistered normalized observables converge quadratically to
their frozen integer limits.  Both observed orders for every observable lie
within `[1.8,2.2]`; in fact they lie between `1.999996` and `2.000002`.
Both independently computed Richardson intercepts are consistent with the
same preregistered integers.

## Leading recurrence

Let `a_n` be the absolute log scale, `r_n` the absolute relative-lapse log,
and

```text
u_n = a_n-a_(n-1),
v_n = r_n-r_(n-1).
```

Let `k_lambda=lambda*epsilon3*L0*tau0/4` be the exact static boundary
momentum.  The computation establishes the empirical small-lapse asymptotic
law

```text
u_n/u_1                 = n     + O(lambda^2),
v_n/v_1                 = 2n-1  + O(lambda^2),
p_post,n/k_lambda       = 2n+1  + O(lambda^2),

a_n/u_1                 = n(n+1)/2 + O(lambda^2),
r_n/v_1                 = n^2       + O(lambda^2),
```

for `n=1,2,3` on the tested branch.  Thus the previously observed triangular
scale logs and square lapse logs are the leading weak-lapse jet of the
canonical map, rather than accidental near-integers at the published lapse.

This is an asymptotic computational certificate on three steps, not an exact
all-order recurrence theorem.

## Numerical evidence

The three columns below correspond to `lambda=(1/2,1/4,1/8)`.

| observable | target | computed values | observed orders |
|---|---:|---|---|
| `u2/u1` | 2 | `2.000000668216`, `2.000000167054`, `2.000000041764` | `2.0000006`, `2.0000002` |
| `u3/u1` | 3 | `3.000002672867`, `3.000000668216`, `3.000000167054` | `2.0000012`, `2.0000003` |
| `a2/u1` | 3 | `3.000000668216`, `3.000000167054`, `3.000000041764` | `2.0000006`, `2.0000002` |
| `a3/u1` | 6 | `6.000003341083`, `6.000000835270`, `6.000000208818` | `2.0000010`, `2.0000003` |
| `v2/v1` | 3 | `3.000001054386`, `3.000000263596`, `3.000000065899` | `2.0000005`, `2.0000003` |
| `v3/v1` | 5 | `5.000005271936`, `5.000001317983`, `5.000000329496` | `2.0000015`, `2.0000002` |
| `r2/v1` | 4 | `4.000001054386`, `4.000000263596`, `4.000000065899` | `2.0000005`, `2.0000003` |
| `r3/v1` | 9 | `9.000006326322`, `9.000001581579`, `9.000000395395` | `2.0000013`, `2.0000002` |
| `p_post,1/k` | 3 | `2.999995994129`, `2.999998998532`, `2.999999749633` | `1.9999992`, `1.9999998` |
| `p_post,2/k` | 5 | `4.999979970668`, `4.999994992660`, `4.999998748165` | `1.9999980`, `1.9999995` |
| `p_post,3/k` | 7 | `6.999943917963`, `6.999985979454`, `6.999996494861` | `1.9999962`, `1.9999990` |

The target-independent leading coefficients also stabilize quadratically:

```text
u1/lambda^2 -> approximately -3.11605992e-6,
v1/lambda^2 -> approximately -3.55925586e-6.
```

Even and odd order-24 schedules agree far inside every frozen tolerance.  The
largest recorded parity discrepancies are numerical noise, not resolved
schedule dependence.  Every Newton solve converges in at most three accepted
full steps, and all evaluated action, derivative and trial states remain on
the Lorentzian branch.

## What the result means

- **DERIVED COMPUTATIONAL:** the canonical fixed-carrier map has a resolved
  temporal small-lapse jet with constant leading momentum change and linearly
  changing log-scale increments over the first three steps.
- **STRUCTURAL:** odd momenta and triangular positions are the generic form of
  constant-force/constant-acceleration stepping.  The integer sequence is
  therefore evidence that the three accepted ticks share one coherent local
  dynamics, not by itself evidence for new microscopic physics.
- **STRUCTURAL:** shrinking the lapse on one 600-cell carrier is temporal
  asymptotics, not spatial refinement.
- **OPEN:** an analytic derivation of the Taylor coefficients directly from
  the reduced action.
- **OPEN:** whether the recurrence holds beyond the three tested steps or at
  finite lapse to all orders.
- **OPEN:** convergence under spatial refinement, anisotropic stability and
  recovery of continuum Einstein/Friedmann dynamics.
- **OPEN:** emergent time and absolute units.  `lambda` and `tau0` remain
  external inputs, so this result does not derive `c`, Planck time or Planck
  mass.

## Post-result primary-source audit

The broad setting is established: consistent Regge discretizations can define
canonical transformations (<https://arxiv.org/abs/gr-qc/0511096>), and
action-generated discrete evolution has pre/post data and consistency
conditions (<https://arxiv.org/abs/1108.1974>).  Multi-step dust evolution of
the 600-cell is also known (<https://arxiv.org/abs/gr-qc/0009093> and
<https://arxiv.org/abs/gr-qc/0106077>).

Variational-integrator theory relates the order of a discrete canonical map to
the accuracy of its generating function
(<https://arxiv.org/abs/1609.02309>).  This supports the conservative reading
of the observed quadratic corrections as an ordinary local integration law.

No located primary source states the present eleven-ratio
triangular/square/odd-momentum certificate for this fixed staircase action.
External novelty remains **OPEN**; a search cannot establish novelty.

## Next decisive test

Do not calculate tick four and then reinterpret it.  First derive a numerical
tick-four prediction solely from this committed weak-lapse jet, including an
uncertainty rule fixed from the Richardson remainders.  Commit that prediction
before evaluating any fourth slab.  Tick four is then a genuine out-of-sample
test of whether the local recurrence continues beyond the three steps used to
identify it.
