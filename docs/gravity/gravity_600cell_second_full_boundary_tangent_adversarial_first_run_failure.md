# First dense second-tangent adversarial run: control failure

Date: 2026-08-22

Status: **PRESERVED BEFORE ANY IMPLEMENTATION CORRECTION.**

Protocol commit: `eceda30`  
Registry commit: `4ce8ee9`  
Implementation commit: `977f3ba`  
Verifier SHA-256:
`c6575c5f74bee291fd0dfff796f171ef6367ca7f13080deedd1b6ad09766d570`  
Failed artifact SHA-256:
`2338fd58ab50ff309b50b30bc7a9beeac5596bc362cd6d07732e23c3a125acc9`

Only the targeted verifier was run. It ended with `27/28 PASS` and

```text
SECOND_FULL_TANGENT_DENSE_CONTROL_FAILED
```

## What passed

- both complete carriers, Lorentzian branch controls and raw reciprocity
  gates;
- all fourteen raw/Richardson Hessian scale comparisons, with
  `SCALE_LIFT_CONFIRMED` in both parities;
- all eighteen full pre-Legendre rank classifications and dense tangent
  constructions;
- exact reproduction of all six previously accepted dense first-slab tangent
  hashes;
- the direct `c=r1^2` tangent conjugacy and rejection of the identity, `r1`
  and omitted-`K_NO` hostile tangent constructions;
- direct second-slab schedule agreement;
- canonical classification of all four dense two-step products and agreement
  of all six schedule pairs.

## What failed

The common `1e-3` product corruption was classified `OPEN`, with normalized
distance `5.554945944616681e-09` and uncertainty
`2.9756217108579626e-10`. The identity and `r1` product lifts were still
refuted. The product uncertainty used the product of the two pre-Legendre
condition estimates; whether this is a justified propagation or an overly
loose implementation is **OPEN** until derived, not adjusted to obtain a
desired label.

The delayed entrywise comparison produced 24 `AGREES` and 32 `REFUTED`
labels. The split was exact by minimal-sector dimension:

| sector dimensions | first | second | products |
|---|---:|---:|---:|
| `3,2,2,2` | 8 refuted | 8 refuted | 16 refuted |
| `1,1,1` | 6 agree | 6 agree | 12 agree |

The order-one discrepancies in every non-scalar sector, with agreement in
every scalar sector, are evidence of a possible left/right or conjugate-basis
error in the delayed closure map. They are not evidence for that diagnosis:
the primary result remains **OPEN** until a mechanically derived basis
convention either removes or confirms the disagreement.

No tangent spectrum, continuum target or full suite was run. No scientific
positive is accepted from this execution.
