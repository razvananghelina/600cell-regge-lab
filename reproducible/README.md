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
  reproducible/verify_gravity_600cell_finite_height_invariant_region.py

/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_invariant_region_adversarial_resolution.py
```

Expected summaries:

```text
primary:     14/14 PASS
adversarial: 12/12 PASS
```

Accepted artifact hashes:

```text
gravity_600cell_finite_height_invariant_region.json
9b6a473c462e7d23af50878cdd4d849bb66c69068c3178b82235e5d0e39926b9

gravity_600cell_finite_height_invariant_region_adversarial_resolution.json
813e05bd66b47cc3ae1cd35d0a2eddb9c645a850d84abeaad37d15b14a6a380f
```

The first direct-quotient adversarial route is deliberately retained as an
`OPEN` result; its interval dependency failure is documented in the result
note and is not counted as corroboration.

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
