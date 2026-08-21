# Correction protocol: adversarial generic-velocity verifier

Date: 2026-08-21

Registered adversarial implementation commit: `74653d6`.

Preserved first-run artifact:

```text
reproducible/gravity_600cell_generic_velocity_composition_adversarial_first_failure.json
SHA-256 cdd84af7b9554227ed925259aed339e781d9d8e810fbd40e87c78b013233c605
```

The first run returned `7/9` and outcome
`GENERIC_VELOCITY_ADVERSARIAL_DISAGREEMENT`.  The independently derived
action, constraint, momentum and mass formulas all matched the frozen primary
artifact exactly.  The symmetry, hostile and numerical gates also passed.

## Isolated failure

The branch gate tested the coefficient of `mu` by Python structural equality:

```text
mu_coefficient == -8*pi
```

The direct constraint retains the equivalent unfactored coefficient

```text
-8*pi*(v^2+4)/(v^2+4),
```

so the structural comparison is false, while the independently evaluated
residual is exactly

```text
simplify(mu_coefficient + 8*pi) = 0.
```

All other branch predicates are individually exact:

```text
simplify(C(v,mu(v))) = 0,
mass_match = True,
both exact cosine-bound residuals = 0.
```

Replace only the structural coefficient comparison by

```text
simplify(mu_coefficient + 8*pi) == 0.
```

This does not change the action, the constraint, the branch, any registered
point, any tolerance or the outcome hierarchy.  Rerun only the targeted
adversarial verifier.

