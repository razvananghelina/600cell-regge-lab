# Result: complete scale+strut carrier/full-equation intersection is zero

Date: 2026-08-20  
Status: **DERIVED LOGICAL/COMPUTATIONAL; route kill boundary**

## Mechanical verdict

After preserving a first `5/6` coverage-parser failure, the corrected registered
verifier reports

```text
FULL_SCALE_STRUT_FULL_EQUATION_INTERSECTION_ZERO
6/6 tests passed
```

Corrected artifact SHA-256:

```text
964e993fd9078387eab7064537b5f496d46abfcfd77182671bbc0903ec6e29a4
```

Verifier source SHA-256:

```text
5a1c2dd44abc4b84f55f2e03b14d4207b0b4b910a5916880bc68bec3da91951a
```

Only the targeted logical verifier and static registry audit ran.  The full
suite did not run.

## Complete certificate composition

The fixed carrier has seven binary-symmetry sectors in each of two parities.
The complete ledger is:

```text
12 nonhomogeneous cells:
  weak carrier/canonical intersection = 0,
  supported by 24 aggregate D/K records containing 48/48
  cross-precision direct-minor certificates;

2 homogeneous cells:
  weak carrier/canonical intersection dimension = 1,
  unique generator transverse to the complete pole equation,
  full-equation intersection dimension = 0.
```

Therefore all fourteen full-equation intersection dimensions are zero.  The
pole-null negative control leaves exactly two homogeneous lines, so the
verifier does not return zero independently of the pole result.

## Preserved first failure

The first consolidation incorrectly required 48 JSON records.  The upstream
artifact actually stores 24 records, each containing two cross-precision
certificates.  Its `5/6 FULL_EQUATION_CARRIER_COVERAGE_OPEN` artifact is
preserved in commit `a1aed02`.  The frozen correction checks all nested
certificates and writes a distinct artifact; it changes no scientific matrix or
dimension.

## Kill boundary and scope

**DERIVED NEGATIVE:** on the accepted nonstatic background and under all frozen
carrier, action, symmetry, coordinate, branch and fixed-input hypotheses, the
complete geometric scale+strut carrier contains no nonzero tangent satisfying
the full canonical equations.

This closes the route that attempted to select physical perturbations as the
full-equation intersection of this exact 240-column carrier with the canonical
graph.

It does **not** establish any of the following:

- that the accepted four-step homogeneous root sequence is false;
- that the unrestricted 1,440-dimensional action-generated tangent does not
  exist;
- that Regge gravity has no perturbations;
- that every possible independently derived carrier must fail.

The result says the present carrier is not a closed physical phase-tangent
space.  Enlarging it after this no-go without an independent geometric
selection rule would be fitting.

## Consequence for the next mission

No further spectrum, inertia, `c` or mass interpretation should be extracted
from this carrier.  The defensible continuation is a spatial-refinement test of
the unrestricted canonical/Jacobi map, tracking which soft
pseudo-constraint-like and curvature-carrying sectors persist or converge.
Constructing a different carrier is allowed only if geometry fixes it before
comparison with desired physical modes.

