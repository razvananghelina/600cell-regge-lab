# Result: refined canonical-map feasibility census

Date: 2026-08-20

Prior-art commit: `883b4e7`  
Protocol commit: `b0d42b8`  
First failed execution recorded at commit: `fc92015`  
Control correction protocol commit: `70b13cb`  
Corrected implementation commit: `969bafd`

## Complete hypotheses

The census uses only the two already certified spatial carriers

```text
K0 = P(sd K_600),
K1 = P(Esd_2(sd K_600)),
```

their face-dimension-derived proper four-colourings, fixed time orientation,
and the standard four-pentachoron staircase for each of all `4! = 24` linear
orders of the colour classes.  Layer vertices remain labelled.  An internal
edge is an edge not contained entirely in either spatial boundary.

No Regge derivative, Hessian value, eigenvalue, continuum target, schedule
average, schedule selection or coarse/fine transport enters this census.

## Reproducibility controls

The corrected targeted verifier passed `8/8` twice.  Both executions wrote a
30,887-byte JSON artifact with SHA-256

```text
ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e
```

The first execution before the disclosed correction reached no enumeration
and is retained in
`gravity_600cell_refined_canonical_map_feasibility_first_failure.md`.

## DERIVED exact census

| quantity | `K0` | `K1` |
|---|---:|---:|
| spatial `(V,E,F,T)` | `(2640,17040,28800,14400)` | `(19680,134880,230400,115200)` |
| schedules enumerated | 24 | 24 |
| slab pentachora `4T` | 57,600 | 460,800 |
| boundary edges `2E` | 34,080 | 269,760 |
| internal edges `V+E` | 19,680 | 154,560 |
| total slab edges `V+3E` | 53,760 | 424,320 |
| pre-Legendre dimension `V+2E` | 36,720 | 289,440 |
| boundary phase-map dimension `2E` | 34,080 | 269,760 |
| distinct internal-edge sets | **24** | **24** |
| cross-diagonal intersection | **0** | **0** |
| cross-diagonal union `2E` | 34,080 | 269,760 |
| schedule pairs checked | 276 | 276 |
| minimum positive weighted distance | 1,440 | 18,240 |
| maximum weighted distance | 17,040 | 134,880 |
| local upper Hessian-incidence count `55(4T)` | 3,168,000 | 25,344,000 |

Here the weighted distance is the exact number of spatial edges whose selected
cross-layer diagonal changes.  Direct edge-set comparison agrees for every
one of the 276 schedule pairs with the independent formula obtained by
summing the colour-pair edge populations over inverted colour pairs.

The colour-pair populations are

```text
K0: 01=1440, 02=3600, 03=2400, 12=3600, 13=3600, 23=2400
K1: 01=18240, 02=28800, 03=20400,
    12=18240, 13=28800, 23=20400.
```

The complete weighted-distance histograms are frozen in the JSON artifact.
Their extrema already prove that neither refinement hides a repeated schedule:
the smallest change affects an entire colour-pair orbit, while reversing an
order changes the diagonal on every spatial edge.

## STRUCTURAL size ledger

Dense `float64` materialization would require

| array | `K0` | `K1` |
|---|---:|---:|
| slab Hessian | 21.533 GiB | 1,341.458 GiB |
| pre-Legendre Jacobian | 10.046 GiB | 624.176 GiB |
| boundary phase map | 8.653 GiB | 542.182 GiB |

These are storage facts, not a physical no-go.  The local upper-incidence
counts show that sparse or matrix-free operations remain admissible.  They do
close direct dense construction on `K1` as a practical implementation route.

## Verdict

**DERIVED:** the frozen geometry has 24 distinct legitimate internal temporal
edge sets at each refinement.  There is no schedule-independent cross-diagonal
core and no geometry-selected single simplicial slab.

Therefore the preregistered outcome is

```text
REFINED_MAP_SCHEDULE_ELIMINATION_REQUIRED
```

This is not yet evidence that the physical boundary dynamics depends on the
schedule.  Distinct bulk triangulations can still yield equal effective
boundary quadratic forms after their internal variables are eliminated.  It
would be circular to choose the cheapest or most attractive schedule before
testing that statement.

## Exact next falsification boundary

On `K0`, construct the action Hessian for all 24 schedules at one frozen
homogeneous Lorentzian background, eliminate the schedule-specific internal
edge variables, identify their common labelled boundary edge space, and test
equality/covariance of the resulting effective boundary quadratic operators.
Use sparse or matrix-free operations; do not materialize a dense phase map.

- If all 24 effective boundary operators agree under the already derived
  colour relabellings, temporal-schedule ambiguity is removed and a refined
  dynamics calculation is licensed.
- If inequivalent boundary operators remain and no independent geometric
  principle selects or sums them, the simplicial 600-cell dynamics is not
  canonical in the stated construction and this route reaches its kill
  boundary.

Even a positive result would leave the independent coarse/fine phase-space
transport problem open.  This census derives no graviton, physical tick,
effective `c`, `G` or Planck scale.

