# Preserved second blind failure: lapse differencing was too small on level 2

Date: 2026-08-17

Quadratic-lapse protocol correction: `e0fcda4`  
Corrected implementation: `8d9983d`

Second artifact:

```text
reproducible/gravity_600cell_projected_refinement_acceleration_blind.json
SHA-256 32a534a57224b722682d41f5aa1bbef2711ac55e64146ee8107d591dbd9dbc92
```

## Outcome

The second targeted run again exited `7/9` with

```text
PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_OPEN.
```

The quadratic lapse correction worked on the coarse carrier and on seven of
the eight refined carriers.  Exactly one gate failed:

```text
carrier: first_tie_rank_0 / level2
held-out lapse quadratic residual: 3.03059e-6
frozen threshold:                 2.00000e-6
```

Its dynamic lapse root still differed from its independently computed seam
root by only `9.30e-6`, inside the unchanged `5e-5` agreement gate.  No
continuum distance was computed by the blind verifier.

## Internal conditioning audit

The failure occurs after subtracting `O(eta)` curvature and dust terms and
then dividing the remainder by `eta^3`.  On the single failed carrier, the
same held-out quadratic identity and dynamic root were recomputed at four
real-log derivative base steps:

```text
h          dynamic root          held-out quadratic residual
8e-3      -0.5023727903                 4.21e-7
4e-3      -0.5023722989                 2.11e-8
2e-3      -0.5023812437                 3.03e-6   <- registered failure
1e-3      -0.5023624287                 5.42e-6
```

The loss worsens as the subtraction step becomes smaller, while the two
larger steps agree in the root to `4.92e-7`.  This is the signature of
floating-point cancellation, not derivative truncation or a nonquadratic
physical equation.  The comparison uses the held-out algebraic identity,
not closeness to the continuum target.

## Frozen correction required

Before rerunning, replace only the lapse real-log base steps by

```text
primary h_F=4e-3,
repeat  h_F=8e-3.
```

Keep the seam estimator, eta values, sentinels, carriers, coefficients,
quadratic threshold and lapse/seam threshold unchanged.  Because the step
choice was audited after inspecting one refined internal residual, it is a
transparent numerical-conditioning correction, not a prospective blind
choice.  The second failed artifact remains in git history.
