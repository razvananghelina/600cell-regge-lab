# First result: the homogeneous refined momentum lift has five free directions

Date: 2026-08-21

Status: primary exact result; adversarial replication still required before
consolidation.

## Provenance

| stage | commit |
|---|---|
| prior-art and framing gate | `3e188e7` |
| target-free protocol | `3aefd29` |
| registered verifier before execution | `26e6ccc` |

Verifier:

```text
reproducible/verify_gravity_600cell_refined_homogeneous_cotangent_lift.py
```

Artifact:

```text
reproducible/gravity_600cell_refined_homogeneous_cotangent_lift.json
SHA-256 93dd857bff3b406e86d41a8a4b05d6441cb0e3e1c11e4f53d098555b1218924b
```

Two targeted executions were byte-identical and each passed `12/12`.  The
full suite and the nested refined root census were not run.

## Primary result

For refined orbit-total momenta, preservation of the canonical one-form on
the homothetic configuration line gives

```text
p_s = 2(P_01+P_02+P_03+P_12+P_13+P_23).
```

The exact pullback row has rank one and nullity five.  Hence, for every fixed
coarse homogeneous momentum `p_s`, the compatible refined momenta form an
affine five-parameter family.

For common per-edge momenta the row is

```text
2(1440,3600,2400,3600,3600,2400).
```

It also has exact rank one and nullity five.  The positive population
diagonal has determinant

```text
386983526400000000000
```

and exactly intertwines the two conventions.  The result is therefore not a
choice between per-edge and orbit-total normalization.

The one-orbit synthetic control has nullity zero and unique unit-momentum
lift `1/2`.  Log-length versus log-squared-length coordinates and reversal of
the orbit order preserve the actual rank and nullity.  Zeroing one synthetic
population correctly destroys convention invertibility.

Frozen primary outcome:

```text
REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNDERDETERMINED
```

## Provisional interpretation

**STRUCTURAL, pending adversarial replication:** projected barycentric
geometry, full spatial `H4` invariance and the canonical pairing do not select
a unique homogeneous coarse-to-refined momentum lift.  They define the
pullback from refined covectors to the coarse covector, whose inverse fiber
has dimension five.

This does not refute a lift selected by the refined on-shell action, an
independently derived supermetric or perfect-action coarse graining.  It does
refute using a Euclidean pseudoinverse, population weighting or any other
preferred point without declaring and deriving that extra structure.

No refined tick, tensor mode, dispersion, `c`, `G` or Planck scale follows.

## Required adversarial gate

Rebuild the result without a matrix-rank or nullspace routine:

1. construct five explicit independent differences of compatible lifts;
2. use exact rational satisfiability to exhibit distinct lifts for the same
   coarse momentum;
3. prove the one-orbit control admits no distinct pair;
4. show that imposing an additional full `S6` permutation symmetry makes the
   lift unique, while the actual rank-coloured `H4` action does not supply
   that symmetry.

Until this passes, the primary verdict is not consolidated under Rule 4.
