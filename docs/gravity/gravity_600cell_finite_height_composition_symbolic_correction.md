# Frozen correction: local junction identity certificate

Date: 2026-08-21.

Protocol commit: `386466e`.

Registered first implementation: `c7bffc3`.

Status: frozen after interrupting the first execution and before modifying or
rerunning the verifier.

## Preserved failure

The first execution passed provenance and exact degree-two homogeneity, then
spent approximately 180 seconds inside a global `sympy.simplify` call on

```text
(L1/2)*d(S1+S2)/dL1-(p_post,1-p_pre,2).
```

It was manually interrupted before any outgoing momentum, root or physical
outcome was evaluated.  The traceback and status are recorded in

```text
reproducible/gravity_600cell_finite_height_composition_first_timeout.json.
```

## Allowed correction

The registered physical and algebraic condition is unchanged.  Construct

```text
d1=dS1/dL1,
d2=dS2/dL1,

shared=(L1/2)*(d1+d2),
p_post,1=(L1/2)*d1,
p_pre,2=-(L1/2)*d2.
```

Then certify the local polynomial identity

```text
shared-(p_post,1-p_pre,2)=0
```

before substituting the large Regge derivatives.  This tests the same sign
convention and follows from distributivity; it does not ask a general-purpose
simplifier to rediscover linearity after expanding two nested inverse-
trigonometric actions.

No root bracket, state, tolerance, action, momentum convention or outcome may
change.  If the local certificate or later direct shared-slice residuals fail,
the composition result remains `OPEN`.
