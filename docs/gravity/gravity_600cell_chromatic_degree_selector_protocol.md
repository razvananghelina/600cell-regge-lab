# Preregistration: exact chromatic-degree selector

Date: 2026-08-17

Prior-art gate: `38cce14`.

Status: frozen before the first degree computation.

## 1. Inputs and exclusions

Rebuild the 600-cell and its fixed left-coset binary-tetrahedral cover using
only:

```text
commons/cell600.py
SHA-256 ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

The previous orientation-census source may be compared as an implementation
control only:

```text
reproducible/verify_gravity_600cell_staircase_orientation_selector.py
SHA-256 4885fd9c69ecc82c2d0aa31b5cde72b123999ef01792f1aebc8435ba063dc90e
```

Do not parse its result artifact, the nonlinear result, a Regge action,
schedule quality, continuum value or desired sign.

Require the frozen carrier counts `120/720/1200/600`, five disjoint colour
classes of size 24, and four distinct colours on every tetrahedron.

## 2. Source and target fundamental chains

For each sorted spatial tetrahedron `t=(v0,v1,v2,v3)`, define

```text
c_K(t) = sign det[v0;v1;v2;v3]
```

from the four stored vertex-coordinate rows.  Require absolute determinant
greater than `1e-10` for all 600 tetrahedra and exact integral cancellation of
the complete simplicial boundary of

```text
C_K = sum_t c_K(t)[t].
```

Use the target boundary chain

```text
C_target = sum_i (-1)^i [0,...,omit i,...,4].
```

## 3. Exact pushforward degree

For each source tetrahedron, list its four colours in source vertex order and
sort them into the canonical target facet missing colour `i`.  Multiply
`c_K(t)` by the sign of this sorting permutation.  Sum these integer
coefficients separately for all five missing colours, producing `P_i`.

Define five degree candidates

```text
d_i = (-1)^i P_i.
```

They must be identical.  Independently recompute one `P_i` by direct signed
preimage enumeration and require the same integer.

No normalization by the number of tetrahedra or colour classes is allowed.

## 4. Complete 120-order census

For every total order `sigma` of the original colours, relabel an old colour
`c` by its rank in `sigma` and recompute all five `d_i`.  Require:

```text
d_sigma = sign(sigma) * d_identity
```

for all 120 orders.  Report the exact degree multiset separately on the two
already defined `A5` parity classes, but do not load their Regge outputs.

Rebuild the 60 induced even cover permutations from the exact `H4` setwise
stabilizer and require that they preserve `d_identity`.

## 5. Mechanical outcome

- If every control passes and `d_identity = 0`, return
  `CHROMATIC_DEGREE_ZERO`.
- If every control passes, `d_identity != 0`, and the two 60-element order
  classes have opposite nonzero degrees, return
  `CHROMATIC_ORIENTATION_LINE_DERIVED`.
- Otherwise return `OPEN_CONTROL_FAILURE`.

## 6. Interpretation boundary

A nonzero degree derives an orientation of the abstract colour simplex from
the oriented 600-cell colouring.  It permits a target-blind convention such as
"positive chromatic degree", but no current Regge, causal or matter axiom has
been shown to require that convention.

Therefore a nonzero result is **DERIVED MATHEMATICAL / STRUCTURAL**, while its
use as a physical time-schedule selector remains **OPEN**.  It must not be used
post hoc to choose whichever nonlinear schedule gives a preferred physical
answer.

Register the verifier before its first execution.  Run only that verifier.
