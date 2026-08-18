# Protocol implementation correction: OPEN is a valid passing result

Date: 2026-08-17

Original protocol: `e50a0ea`.

Preserved first run: `2eb8b40` (`12/14`).

The first run reconstructed every object, labelled all 56 principal-angle
comparisons and correctly assigned

```text
HYPERBOLIC_EXTREME_SUBSPACE_OPEN.
```

The homogeneous sector's fifth expanding/contracting pair has gap `1.006`,
below the frozen `>2` gate.  That scientific result is unchanged.

The implementation error is procedural.  Section 7 of the protocol lists
`HYPERBOLIC_EXTREME_SUBSPACE_OPEN` as a valid mechanical outcome, and states
that the verifier passes when it reconstructs and classifies the object.
Nevertheless, the code used the truth of the physical gap predicate itself
as two `check(...)` conditions, causing a valid OPEN result to exit nonzero.

Correct only the test semantics:

- require every fixed-count Schur selection to return the requested `5d`
  dimension and a finite positive gap;
- record separately whether every gap exceeds two;
- use the latter only in the already frozen outcome hierarchy.

Do not change the gap threshold, any matrix, response, subspace, principal
angle, uncertainty, label or outcome.  Also remove the misleading printed
count from the check label: there are 112 variant/branch selections, while
56 is the number of final candidate comparisons.
