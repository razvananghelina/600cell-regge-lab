# Reproducible gravity verifiers

This directory contains only the current 600-cell Regge programme:

- explicit registered `verify_gravity_*.py` verifiers;
- frozen `gravity_*.json` and `gravity_*.npz` artifacts;
- the coverage-guarded registry in `run_all.py`;
- two governance guards for the public layout and authoritative theory map.

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

## Theory-map guard

The authoritative route and no-go registry is checked without running a
scientific calculation:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_theory_map.py
```

The expected summary is `14/14 PASS`. It requires unique route IDs, exactly
one active gate, existing evidence paths, an acyclic dependency graph,
explicit reopening conditions for bounded no-go results, and exactly one
registration of the verifier itself.

## Latest bounded no-go

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_constraint_quotient_no_go.py
```

Expected outcome:

```text
FINITE_HEIGHT_TWO_STEP_EXACT_CONSTRAINT_QUOTIENT_BOUNDED_NO_GO
14/14 PASS
```

Accepted artifact SHA-256:

```text
bb3735fd475c08068f8240ed1f85a0b205a529228cb46f2d6d8a788114b9e49f
```

The verifier checks the accepted first/second Legendre rank margins, two
independent no-constraint arguments, a truly singular positive constraint
control, an exact `10^-100` anti-threshold control and a later-singular-move
scope control. It computes no mode spectrum and runs no full suite.

## Latest refinement feasibility result

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_nested_vertex_displacement.py
```

Expected outcome:

```text
NESTED_TANGENTIAL_VERTEX_CARRIER_DERIVED
14/14 PASS
```

Accepted artifact SHA-256:

```text
ef7f565fd54487885e9de459d97fbe61af8dcc6e9add30768ae2fe94ca7d4250
```

The verifier reconstructs the two carriers by an edgewise-facet route and an
independent old-plus-edge route. It certifies parent independence, an exact
old-vertex left inverse of rank 7920, complete fine-edge differentiation and
`O(4)` covariance. The preserved first run is `13/14`: its binary64 centered
edge derivative failed the frozen convergence control before the
preregistered extended-precision correction. The accepted result is a
tangential configuration prolongation, not a momentum lift or restored
constraint.

## Earlier internal-carrier pair

The latest targeted internal-carrier rank pair is

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_internal_carrier_rank.py

/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_internal_carrier_rank_adversarial.py
```

Expected summaries are `25/25 PASS` and `19/19 PASS`.  They select one common
homogeneous internal-constraint tangent and no nonhomogeneous survivor.  They
do not by themselves impose the fixed incoming canonical momentum equation.

The exact reconciliation is

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_internal_kernel_canonical_reconciliation.py
```

Its expected summary is `14/14 PASS`.  It proves that the surviving line is
the homogeneous lapse-constraint tangent and that fixed incoming momentum
removes it.  The next object is a forced response to varying incoming data,
not another kernel search.

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

## Latest incoming-domain discovery

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_incoming_basin_discovery.py
```

The frozen targeted execution returned `8/8` with outcome
`INCOMING_BASIN_CANDIDATE_SKELETON_FROZEN`. Its lossless compressed artifact
is

```text
gravity_600cell_finite_height_incoming_basin_discovery.json.gz
f492f50cfcaa8e171fb6faa21524d824b4d11b3701b7d635ce483500aaffeb8d
```

The archive expands to JSON SHA-256
`146a6a1426044e2065c66a2a8974bd94ff24deba4ef0527acf90c5ba459dee58`.
This is a finite candidate skeleton, not a continuum basin theorem.

## Latest local theorem

The primary and exact final resolver are

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_local_signature.py

/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_local_signature_identity_resolution.py
```

Expected summaries are `10/10 PASS` and `7/7 PASS`. Accepted artifact hashes
are

```text
gravity_600cell_finite_height_local_signature.json
9f524cc22df8cfb5083f372481b3efd19868252b85551d56378327eea7a6d613

gravity_600cell_finite_height_local_signature_identity_resolution.json
ccedea4620f7cd485381f8002a8fa29b39a7842a94867e430a630024e6e7eb60
```

The mechanically different direct-bracket chain deliberately preserves a
`4/11 OPEN` raw-interval attempt and an `8/10 OPEN` monotone-factor attempt.
The latter certifies every root and terminal; the `7/7` exact resolver closes
its sole auxiliary identity-width failure algebraically. See the consolidated
result note before interpreting individual exit codes.

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
