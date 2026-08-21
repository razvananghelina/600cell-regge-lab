# First execution failure: refined boundary cotangent verifier

Date: 2026-08-21

The first registered execution returned

```text
REFINED_BOUNDARY_COTANGENT_CONTROL_FAILED
15/16 PASS
```

The preserved artifact has SHA-256

```text
8a0f7d8f1a6c141e6b8e739abf9a45e3dd9a7e807f5c5dcb5d6a2e195262bf50.
```

No result is accepted from this run.

## Failure diagnosis

The synthetic kernel perturbation was built inside a 100-decimal context as

```text
delta=(1e-6,-1e-6,0,0,0,0).
```

The artifact proves that its pullback change was exactly displayed as zero
and its infinity distance was exactly displayed as
`0.0000010000...`.  Nevertheless the boolean control compared that stored
high-precision `mpf` outside the context with a newly constructed
default-precision `mpf("1e-6")`.  Those two binary approximations are not
bit-identical, so the equality test returned false.

Every scientific gate passed in the failed run, including:

- all 24 complete covectors were schedule independent to `1.91e-96`;
- the action-derived seed was internally stationary;
- Richardson derivatives agreed to `1.73e-42` relative;
- the refined pullback identity agreed to `4.32e-76`;
- the raw fixed-radius ratio was `0.98419037738852...`;
- the mass-normalized ratio was `1.00000000000000...`.

These values are retained only as diagnostic output until the control passes.

## Allowed correction

Evaluate both the pullback and distance assertions inside a 100-decimal
context, constructing `1e-6` and the comparison envelope there.  Do not alter
the perturbation, action, schedules, tolerances, physical criteria or outcome
hierarchy.
