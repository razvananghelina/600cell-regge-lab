# Incoming-state basin discovery: finite candidate skeleton

Date: 2026-08-22.

Status: **DERIVED finite census / OPEN continuum basin**.

## Provenance

- prior-art gate: `edd5db0`;
- frozen discovery protocol: `5da21c0`;
- verifier registered before first execution: `7be5a04`;
- first implementation failure preserved: `e6cdb71`;
- initial diagonal-root correction: `d2c791d`;
- second implementation failure preserved: `e40d0ae`;
- stationary-free tail correction used by the successful run: `5dae367`.

Neither correction changed the frozen incoming intervals, node rule, node
count, branch depth, precision, physical gates, terminal labels or success
criteria. The failures and their permitted corrections are recorded in
`gravity_600cell_finite_height_incoming_basin_first_run_failure.md`.

Only the targeted verifier was run. No full-suite claim is made.

## Fixed hypotheses and scope

The census concerns only the fixed homogeneous tetrahedral-frustum 600-cell
model, zero cosmological constant, conserved global dust, the committed
canonical-momentum convention, and the special incoming curve

```text
(m,pi)=(mu(v),p(v)).
```

It does not cover the full two-dimensional canonical state space, any
nonhomogeneous perturbation, a local evolution law, an absolute tick, `c`,
`G`, a Planck scale or particle physics.

## Targeted result

The successful execution reports

```text
RESULT: 8/8 checks passed
OUTCOME: INCOMING_BASIN_CANDIDATE_SKELETON_FROZEN
CANDIDATE_CELLS: 50
DISTINCT_SIGNATURES: 36
```

All 3072 frozen Gauss--Chebyshev inputs had a complete all-real root census
and a complete physical tree through slab four. There were no unresolved
nodes. The largest tree used 8 nodes, below the frozen budget of 256.

The known controls were evaluated only after the blind skeleton existed:

```text
v=3/2: two second-slab branches, terminal labels DEAD and ENTERED_D;
v=3:   two second-slab branches;
v=20:  no physical second-slab successor.
```

The altered post-momentum scale, reversed sign and dust-mass reset all changed
the branch signature. One representative of each of the 36 distinct
signatures satisfied the independently redifferentiated action residual test.

## What the finite census says

The node counts below are deterministic diagnostics, not a probability or a
measure on incoming states:

| Incoming component | Frozen nodes | Diagnostic terminal behaviour |
|---|---:|---|
| `(v_A,v_star)` | 1024 | 227 `ENTERED_D` only; 97 `DEAD+ENTERED_D`; 700 contain at least one branch still live outside `D` at depth 4 |
| `(v_star,v_M)` | 1024 | 644 `DEAD` only; 380 contain at least one branch still live outside `D` at depth 4 |
| `(v_M,v_C)` | 1024 | all 1024 are `DEAD` only |

Across all components, 1080 inputs have at least one
`LIVE_OUTSIDE_D_AT_DEPTH_4` terminal. Their trees contain 4632 such terminal
branches in total. Therefore the frozen rule in the protocol applies:

> the complete incoming basin remains **OPEN**, and the discovery depth must
> not be increased in this mission.

The `DEAD+ENTERED_D` behaviour at `v=3/2` is not an isolated sampled point:
97 adjacent frozen nodes have the same terminal-label multiset. However, this
does not establish a continuum interval. It is surrounded by other branch
types, and the full domain contains 36 distinct signatures in 43 contiguous
signature runs.

The discovery located 50 adjacent candidate cells: 40 contain a combinatorial
signature change and 10 contain an `m*q-125` invariant-entry sign change.
Together with the four inherited boundaries `v_A`, `v_star`, `v_M`, `v_C`,
these form the frozen input for a later interval proof. A signature-change
cell retains both possible intrinsic mechanisms, branch birth/merger and
zero endpoint, until an exact interval calculation distinguishes them.

## Evidential ledger

- **DERIVED (finite census):** all 3072 frozen branch trees, all-real roots,
  physical gates, terminal labels, 36 signatures, 50 candidate cells and the
  three delayed controls.
- **PATTERN:** the sampled `v=3/2` behaviour persists over adjacent nodes, and
  the whole sampled component `(v_M,v_C)` dies after its first slab.
- **OPEN:** every continuum branch count, every exact transition value and the
  complete basin on the incoming curve. The 1080 live-at-depth-four inputs
  prevent a global basin verdict from this mission.
- **STRUCTURAL:** selection by complete forward extendibility remains a
  global-in-time consistency criterion. It has not been derived as a local
  dynamical law.
- **NOT TESTED:** the full `(m,pi)` state plane and nonhomogeneous/local
  degrees of freedom.

This result removes the claim that the whole story was established by one
convenient value `v=3/2`, but it also refutes the stronger hope that its simple
two-branch pattern represents the entire physical incoming domain.

## Artifact integrity

The complete artifact is stored losslessly as
`reproducible/gravity_600cell_finite_height_incoming_basin_discovery.json.gz`.

```text
uncompressed JSON SHA-256:
146a6a1426044e2065c66a2a8974bd94ff24deba4ef0527acf90c5ba459dee58

deterministic gzip SHA-256:
f492f50cfcaa8e171fb6faa21524d824b4d11b3701b7d635ce483500aaffeb8d
```

Decompressing the committed archive reproduces the uncompressed hash exactly.

## Next admissible gate

Do not run more homogeneous slabs merely to see what happens. The next
mathematical step is to preregister an interval-cover proof over the 43 runs
and 50 candidate cells, using the complete branch paths already serialized.
Its admissible outcomes are:

1. exact continuum subintervals whose paths die or enter `D`;
2. a rigorously surviving branch outside `D`, which refutes the proposed
   complete selector on the incoming curve;
3. an unresolved candidate cell, leaving the basin open.

Only if a nontrivial open basin is proved should this homogeneous background
be promoted to the next physical gate: a preregistered nonhomogeneous
perturbation and its linearized evolution spectrum.

