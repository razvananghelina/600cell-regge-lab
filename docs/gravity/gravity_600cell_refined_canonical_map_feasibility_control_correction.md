# Preregistered control correction: refined map feasibility census

Date: 2026-08-20

First-failure record commit: `fc92015`.

This correction is frozen after the first infrastructure-only failure and
before the first carrier reconstruction or schedule enumeration.

## Allowed changes

Exactly these three changes are allowed in
`verify_gravity_600cell_refined_canonical_map_feasibility.py`:

1. replace the malformed 62-character expected digest for
   `verify_gravity_600cell_projected_rank_edgewise_carrier.py` by its directly
   measured 64-character SHA-256
   `50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23`;
2. accept the carrier artifact only when its actual scalar-schema values obey
   `passed == tests == 16` and its outcome is exactly
   `PROJECTED_RANK_EDGEWISE_CARRIER_DERIVED`;
3. keep the balanced-slab control in its actual nested schema and require
   `tests.passed == tests.total == 15`,
   `selection.existence_passes == true`, and
   `selection.ordered_slab_alternatives == 24`.

## Frozen exclusions

No carrier construction, colour convention, staircase, schedule, count,
formula, threshold, outcome ordering, target exclusion or output field may be
changed.  In particular, this correction supplies no permission to select or
average a schedule and no permission to compute an action, Hessian, spectrum
or continuum target.

The corrected verifier must still be run twice and its JSON artifacts must
have byte-identical SHA-256 digests.
