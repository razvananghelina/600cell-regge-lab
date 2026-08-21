# Flag-incidence result: orbit reduction is exact; discriminator state differed

Date: 2026-08-21

## Provenance

- prior-art gate: `c14b5ac`;
- frozen protocol: `c90833d`;
- registered implementation before evaluation: `d4210c8`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_orbit_action_incidence.py`;
- artifact:
  `reproducible/gravity_600cell_orbit_action_incidence.json`;
- artifact SHA-256:
  `3c7073e937309ef6670b68d81b463cfe4393bfd1d73e390f315fd7f3d5e306e6`.

Only the targeted verifier was run.  It returned **9/12** and exited nonzero.
The failed mechanical outcome is retained.

## Mechanical outcome

```text
ORBIT_ACTION_DIRECT_CONSTRUCTION_SUSPECT
```

That label is not accepted as the scientific diagnosis.  The exact incidence
stage falsified the proposed multiplicity bug, and subsequent source inspection
found that the compared action states were not identical.

## Exact incidence result

For each parity the order-24 action gives

```text
triangle orbits = 260,
simplex orbits  = 100,
flag orbits     = 1000,
labelled flags  = 24000.
```

For every flag orbit the exact coefficient `|F_a|/|T_i|` equals the coefficient
implicitly used by the representative-simplex shortcut.  The number of
mismatching flag orbits is exactly **zero** in both parities.  All double-count
and per-triangle row-sum identities pass.

**DERIVED EXACT NEGATIVE:** the direct/orbit disagreement is not caused by a
triangle--simplex flag-incidence multiplicity error.

## Action controls

The exact flag reduction and the shortcut agree on the regular published
control to about `1e-95`; both reproduce the stored regular action within its
50--60 digit serialization.  The exact flag action also reproduces the stored
80-decimal direct off-shell action to `1e-72`--`1e-73`.

It does not reproduce the primary orbit artifact at that supposed same state,
which caused the preregistered discriminator gates to fail.

## Source-level cause found after the frozen result

The original protocol and primary verifier define three separate perturbations:

```text
old boundary : modulus 7, center 3,
internal     : modulus 5, center 2,
new boundary : modulus 11, center 5.
```

The adversarial verifier instead concatenated all 35 internal and 30 final
variables and applied the internal `modulus 5, center 2` rule to all 65.  The
precision adjudicator and this incidence audit inherited that noncompliant
state.  Thus the direct and primary actions were evaluated on different final
boundaries.  Both states are individually stabilizer-invariant, so the
invariance gate could not detect this error.

This diagnosis is exact from the frozen source text, but its numerical repair
must be preregistered and executed separately.

## Status

- **DERIVED EXACT:** carrier maps agree and the orbit flag reduction is exact.
- **DERIVED IMPLEMENTATION FAILURE:** the adversarial state did not implement
  the frozen `q_new` perturbation.
- **RETRACTED AS DIAGNOSIS:** the earlier implication that the orbit action or
  direct action was wrong.
- **OPEN:** direct/orbit agreement on the genuinely identical frozen state.
- **UNCHANGED:** the exact dimensional scale-covariance argument; both wrongly
  and correctly perturbed states are homogeneous under global rescaling.

## Required correction

Keep all failed artifacts.  A new registered verifier must reconstruct
`old`, `internal`, and `new` separately from literal formulas, compare every
labelled edge value between the primary and direct carriers, then rerun the
direct binary64 raw-gradient audit and arbitrary-precision direct action on the
correct state.  Only that result may close the adversarial gate.

