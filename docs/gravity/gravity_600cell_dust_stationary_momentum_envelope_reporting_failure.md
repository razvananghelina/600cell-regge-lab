# Stationary momentum envelope: base-bracket reporting failure

Date: 2026-08-16

## Provenance

- prior-art gate: `dedcbc6`;
- frozen protocol: `ed1cd6a`;
- implementation before evaluation: `0bd3fc4`;
- first artifact SHA-256:
  `92d79110507cb639a1feacaee55670af3fcbed29df24660d09d17cb4eb89b846`.

Only the targeted verifier was run.  It exited nonzero with **4/5** and

```text
MOMENTUM_ENVELOPE_BASE_BRACKET_FAILED.
```

## What is established

The target firewall and accepted-geometry hash passed.  The frozen symmetric
base-bracket algorithm did not accept a sign bracket in its thirteen allowed
widths.  No curve or target comparison was evaluated.

## What is not established

The first implementation failed to serialize the thirteen bracket attempts.
Consequently the artifact does not reveal whether:

- the endpoint values of `G` retained the same sign;
- one or both endpoints left the Lorentzian branch;
- a sign change existed but failed the branch requirement.

Therefore this run does **not** establish absence of an internally stationary
base point.  Its scientific status is **OPEN / reporting failure**.

## Permitted correction

A reporting-only patch may serialize, for every already frozen expansion:

```text
expansion index, half-width, left/right b,
left/right G, branch flags and branch margins.
```

It may also make the test harness treat the preregistered mechanical negative
as a classified passing execution rather than an implementation failure.  It
must not change the center, widths, number of expansions, branch gate, signs,
equations, target firewall or outcome hierarchy.

The identical targeted verifier must then be rerun.  Any new search strategy
requires a separate protocol after the corrected record is committed.
