# Reproducible gravity verifiers

This directory contains only the current 600-cell Regge programme:

- explicit registered `verify_gravity_*.py` verifiers;
- frozen `gravity_*.json` and `gravity_*.npz` artifacts;
- the coverage-guarded registry in `run_all.py`;
- one documentation-layout guard.

Legacy particle-mass, coupling, E8 and fitted-parameter scripts are absent from
the current tree.

## Environment

Use the project interpreter:

```bash
/home/razvan/science/.venv/bin/python
```

The retained verifiers use packages already installed in that environment,
including NumPy, SciPy, SymPy, mpmath and, for selected certified calculations,
python-flint.

## Latest accepted pair

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_asymptotic_map.py

/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_asymptotic_map_adversarial.py
```

Expected summaries:

```text
primary:     11/11 PASS
adversarial: 12/12 PASS
```

Accepted artifact hashes:

```text
gravity_600cell_finite_height_asymptotic_map.json
a93837d2bbec340ddbac528c0be4da52aefe45c8f0d4310496eb1aef6a7b19b6

gravity_600cell_finite_height_asymptotic_map_adversarial.json
5215b2f07140be44f9e864b2688afa5e8e522b310a33ee5f7efa6cfccebc7405
```

## Registry policy

Every `verify_*.py` file must occur exactly once in the explicit `scripts`
list in `run_all.py`, unless it occurs in `DELIBERATELY_SKIPPED` with a
nonempty reason. The guard exits with status 2 for unregistered files,
duplicates, missing files, invalid exclusions or registration/exclusion
overlap.

The full suite is intentionally not run during focused missions. Run only the
targeted verifier and its independent replication unless a full-suite audit is
explicitly requested.

## Scientific scope

A passing verifier proves only the claim and hypotheses printed by that
verifier and its corresponding result note. It does not by itself establish
physical uniqueness, continuum general relativity, an absolute tick, a
limiting speed, Planck units or particle physics.
