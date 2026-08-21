# Correction protocol: adversarial composition convergence

Date: 2026-08-21

Registered adversarial implementation commit: `bc7e1d0`.

Preserved first adversarial artifact:

```text
reproducible/gravity_600cell_generic_velocity_next_order_adversarial_first_failure.json
SHA-256 2fde9147a9632248662d26195ca78c24f3375d79e0a57735a7495fcd3707ceba
```

The first adversarial run returned `9/11` and outcome
`GENERIC_NEXT_ORDER_ADVERSARIAL_DISAGREEMENT`.

Every exact gate passed:

```text
direct derivative-first C1 and P1 match the primary exactly,
K and B factorizations match,
the monotonic exceptional-set theorem passes,
x_star agrees beyond 80 decimals,
all three hostile controls pass,
all 12 new one-slab numerical controls pass.
```

Seventeen of the eighteen stationary composition convergence controls passed.
The sole failure was

```text
v=11/5, action defect,
h=(1/500,1/1000),
residuals=(0.6670230668,0.4500156570),
order=0.5677614563.
```

This velocity lies below but relatively near the exceptional value and has
the large common coefficient `a=-29.2727906269...`; the frozen pair may be
pre-asymptotic.  This explanation is not accepted without a new registered
diagnostic.

Retain the original failed control and evaluate the same unexpanded action
defect at 100 decimals for

```text
v=11/5,
h in {1/2000,1/4000,1/8000,1/16000}.
```

Adjudicate it as pre-asymptotic only if:

```text
- the original failure set contains exactly (11/5,action);
- all four new absolute residuals decrease strictly;
- all three new halving orders lie in [0.8,1.2].
```

Otherwise retain adversarial disagreement.  Do not change the exact route,
root, endpoints, original points, scaling powers, thresholds, hostile
controls or outcome hierarchy.

