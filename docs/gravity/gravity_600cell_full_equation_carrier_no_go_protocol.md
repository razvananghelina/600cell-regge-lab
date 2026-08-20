# Protocol: complete full-equation scale+strut carrier no-go

Date: 2026-08-20  
Status: **preregistered logical/computational consolidation**

## Inputs and no-recomputation rule

Pin by SHA-256:

- the prior-art/framing gate and this protocol;
- the adversarial nonhomogeneous direct-minor artifact;
- the primary and repaired-adversarial homogeneous-line artifacts;
- the pole-transversality artifact.

Do not rebuild any Hessian, carrier, symmetry basis, root or candidate.  This
verifier tests the logical coverage and outcome composition only.

## Required ledger

Require exactly:

```text
2 parities x 6 nonhomogeneous sectors = 12 zero weak intersections,
48/48 direct nonhomogeneous rank certificates,
2 parities x 1 homogeneous sector = 2 one-dimensional weak intersections,
2/2 homogeneous generators transverse to the pole equation.
```

For every nonhomogeneous cell set the full-equation dimension to zero because
adding equations cannot enlarge an intersection.  For every homogeneous cell
set it to zero because its unique weak generator has nonzero pole derivative.
Require a total of fourteen zero full-equation cells.

Negative control: replacing the two pole-transverse flags by pole-null flags
must leave two one-dimensional homogeneous full intersections.  This prevents a
verifier that returns zero independently of the pole result.

## Outcomes

1. `FULL_EQUATION_CARRIER_NO_GO_CONTROL_FAILED`;
2. `FULL_EQUATION_CARRIER_COVERAGE_OPEN`;
3. `FULL_EQUATION_CARRIER_INTERSECTION_OPEN`;
4. `FULL_SCALE_STRUT_FULL_EQUATION_INTERSECTION_ZERO`.

Outcome 4 is **DERIVED LOGICAL/COMPUTATIONAL** under the complete hypothesis
list.  It is a kill boundary for this carrier-intersection selection route, not
for the canonical map or Regge gravity generally.

Only the new targeted verifier and static registry check may run.  No full
suite.

