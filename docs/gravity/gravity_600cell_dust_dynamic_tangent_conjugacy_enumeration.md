# Preregistered geometric conjugacy enumeration

Date: 2026-08-17

This file records Stage A before either tangent matrix is parsed or any
intertwining residual is evaluated.

## Provenance

- prior-art gate: `0129053`
- frozen protocol: `e006b8c`
- registered enumeration verifier before first execution: `0efa93a`

Artifact:
`reproducible/gravity_600cell_dust_dynamic_tangent_conjugacy_enumeration.json`

Artifact SHA-256:
`51b52457eba84ca1e41926b6e4fb1c51032f788b70bde916a3fb755d0323cb3e`

The targeted enumeration passed `7/7`.  It records
`tangent_matrices_parsed=false` and `spectral_target_parsed=false`.  The full
suite was not run.

## Frozen counts

Across all 14,400 `H4` actions:

```text
complete direct even-to-odd slab actions       0
complete time-reversed even-to-odd actions     0
distinct complete direct candidates            0
distinct complete reversed candidates          0

boundary-partition-preserving H4 actions     1440
distinct boundary permutations                 60
```

The physical-edge-set identification is one of the 60 boundary permutations
and is supported by 24 `H4` actions.  Because the independently sorted even
and odd orbit coordinates differ, it is not the numerical identity
permutation.

## Status before comparison

- **DERIVED NEGATIVE:** no `H4` action is an isomorphism between the two
  complete ordered slabs, either preserving or reversing their layers.
- **DERIVED ENUMERATION:** the admissible boundary-only family has exactly 60
  distinct permutations.
- **STRUCTURAL:** Stage B therefore has exactly 120 weaker attempts, one direct
  and one reversed phase-space lift per boundary permutation.
- **OPEN:** whether any of those 120 boundary identities intertwines the
  dynamic tangent maps.
- **OPEN:** the origin of the calibrated dynamic isospectrality.

No complete-slab covariance hit is possible in Stage B because the
geometrically canonical candidate family is empty.  A boundary hit, if any,
cannot be promoted to a four-dimensional slab isomorphism.
