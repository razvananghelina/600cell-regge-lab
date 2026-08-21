# Correction protocol: exact reduction in the adversarial half-step audit

Date: 2026-08-21

Adversarial verifier commit: `18f69cd`.

Preserved first-failure artifact:

```text
reproducible/gravity_600cell_relational_step_composition_adversarial_first_failure.json
SHA-256 02d86a8c1413c7a03b192dd3702380cdbf072645c14df06a085a2b9bbef388fc
```

## Failure

The first adversarial run returned `7/9`.  Its exact-action-first limit
printed

```text
P = 30[-12*pi*A+10*sqrt(2)*A+30*A*acos(1/3)
       -45*acos(1/3)+18*pi].
```

The independently expected expression was

```text
180[(5*sqrt(2)/3-epsilon)A+3*epsilon/2],
epsilon=2*pi-5*acos(1/3).
```

Expanding the second expression gives the first term by term.  The check used
generic `simplify`, which retained the zero as an unreduced linear
combination.  This is an implementation-control failure, not a mathematical
disagreement.  The same run independently passed:

- the nonzero exact resultant and both root obstructions;
- the changed-state hostile control;
- all three direct 100-decimal vanishing residuals with order two;
- nonzero obstruction signal-to-drift ratios `563740` and `155073`.

## Frozen correction

Change only the exact equality predicate for the two limit differences to

```text
cancel(expand(left-right), extension=sqrt(2)) == 0.
```

Do not alter the action, derivatives, substitutions, roots, precision,
sample values, thresholds, outcomes or stored first-failure artifact.  Rerun
only the targeted adversarial verifier.

