# First failure: refined canonical-map feasibility census

Date: 2026-08-20

Protocol commit: `b0d42b8`  
Implementation/registration commit: `12319fe`

## Complete scope of the failed run

The first targeted execution of
`reproducible/verify_gravity_600cell_refined_canonical_map_feasibility.py`
stopped in the frozen-input control, before reconstructing either refined
carrier and before enumerating any temporal schedule.  It therefore produced
no scientific count, no temporal-carrier verdict and no JSON artifact.

The visible output began with

```text
[FAIL] the frozen carrier and slab inputs have exact provenance
```

and then stopped with

```text
TypeError: 'int' object is not subscriptable
```

at the attempted access `carrier_artifact["tests"]["passed"]`.

## DERIVED diagnosis

Three control-only defects are present in the registered implementation.

1. The expected SHA-256 for
   `verify_gravity_600cell_projected_rank_edgewise_carrier.py` has length 62.
   It omitted the two hexadecimal characters `f2` after position 20.  The
   file's actual 64-character SHA-256 is
   `50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23`.
   The other four frozen hashes match exactly.
2. The carrier artifact uses the older scalar schema
   `"passed": 16, "tests": 16`, not a nested `tests` object.
3. The carrier artifact's certified outcome is
   `PROJECTED_RANK_EDGEWISE_CARRIER_DERIVED`, not the string
   `CANONICAL_PROJECTED_CARRIER_EXISTS` assumed by the new verifier.  The
   balanced-slab artifact independently uses the nested `15/15` schema and
   certifies `selection.existence_passes = true` with 24 ordered alternatives.

## Verdict

**DERIVED CONTROL FAILURE:** this is an implementation/provenance-parser
failure and cannot support either a positive or a negative statement about
the 24 schedule carriers.  The frozen enumeration, count formulae, outcome
hierarchy and forbidden-target rules were never reached and remain unchanged.

Any correction must be preregistered and restricted to the three controls
above.  Changing a geometric construction, threshold, schedule convention or
outcome after this failure is forbidden.
