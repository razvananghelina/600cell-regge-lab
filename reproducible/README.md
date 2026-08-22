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
  reproducible/verify_gravity_600cell_finite_height_fourth_slab.py

/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_fourth_slab_adversarial.py
```

Expected summaries:

```text
primary:     7/7 PASS
adversarial: 9/9 PASS
```

Accepted artifact hashes:

```text
gravity_600cell_finite_height_fourth_slab.json
cf322cf0d60668d8f3f58e251425c9ad6bf43b112f22f9f3aebbc28f86212468

gravity_600cell_finite_height_fourth_slab_adversarial.json
ac1ed19fd72549cf7cd054107d921e2819580704391fbc294a55e106a8f7a1bd
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
