# Primary result: local curvature selects the refined rank masses

Date: 2026-08-21

Status: primary post-hoc confirmation; not accepted before the adversarial
replication gate.

## Provenance

| stage | commit |
|---|---|
| prior-art gate and exploratory disclosure | `f12f56c` |
| frozen protocol | `bb80f28` |
| verifier registered before execution | `e76f756` |
| preserved `13/15` first failure | `ead1a1c` |
| frozen control correction | `2db5457` |
| corrected implementation | `1fa48a2` |

The corrected verifier passed `15/15` twice and wrote byte-identical artifact

```text
reproducible/gravity_600cell_refined_local_curvature_mass.json
SHA-256 180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091.
```

The failed artifact and both disclosed implementation defects remain in the
history.  No full suite or nested root census was run.

## Complete hypotheses

Use `K0=P(sd K_600)`, its exact projected rank metric, the induced static
Minkowski-product fill at supplied `tau0=0.0102`, the corrected Lorentzian
Regge action with boundary terms and each of the 24 colour-ordered staircase
triangulations.  Restrict by spatial `H4`, giving four vertical rank orbits.

For a spatial edge `e`, put `C_e=l_e epsilon_e` and split it equally between
its endpoints:

```text
K_v=(1/2) sum_(e incident on v) C_e,
K_r=sum_(rank(v)=r) K_v.
```

The matter model in this calculation consists of four conserved rank mass
totals `mu_r` with action `-8*pi sum_r mu_r sqrt(rho_r)`.

## Computed identity

For all `24*4=96` rank-lapse equations, the independently reconstructed
spatial curvature and the committed action gradients satisfy

```text
dS_grav/dlog(rho_r) = tau0*K_r/2.
```

The maximum absolute disagreement is

```text
1.4016460e-76,
```

and the spread among all 24 schedules is zero at the stored precision.  The
four mass derivatives have Jacobian `-4*pi*tau0*I_4`, hence exact rank four.
Stationarity therefore selects uniquely

```text
mu_r=K_r/(8*pi).
```

The selected mass fractions are

```text
(0.1287831657723389984083291990898...,
 0.3657000761313201399462782272255...,
 0.3759856918014127686434122686089...,
 0.1295310662949280930019803050758...).
```

All are positive and sum to the previously selected global mass.  The total
rank masses are

```text
(0.3046755523027747147650360218964...,
 0.8651742019562024879531654941241...,
 0.8895079385611822780763795593009...,
 0.3064449373262038343639094971899...).
```

## What this does and does not repair

**DERIVED COMPUTATIONAL / STRUCTURAL, POST-HOC PRIMARY:** the four local
lapse equations do not demand an arbitrary fit.  On the fixed static product
geometry they demand the endpoint-half Regge curvature mass at each rank.
This is the local refinement of the older global balance
`M=K/(8*pi)`.

**DERIVED NEGATIVE:** the conditional affine-exact `P1` density puts mass
fraction `1/4` in every rank and is not locally stationary on this finite
carrier.  Keeping matter only on original vertices also fails.  The visually
near `(1,3,3,1)/8` pattern is not exact; its maximum fraction error is
`0.00929992386868...`.

The curvature-selected density is not homogeneous relative to the `P1` dual
volumes.  Its four density ratios are approximately

```text
(0.5151, 1.4628, 1.5039, 0.5181),
```

with maximum/minimum ratio `2.919525...`.  Thus the result supplies consistent
time-symmetric discrete initial matter, not a derivation of homogeneous dust.
Whether this rank contrast tends to one under refinement is **OPEN** and is
the physically relevant convergence question.

The result does not yet supply an unequal-boundary tick, a physical lapse,
tensor propagation, `c`, `G`, Planck units or particle masses.

## Acceptance gate

The exploratory equality was seen before the formal protocol, so byte-level
reproduction of the primary verifier is insufficient.  Before acceptance,
the adversarial route must:

1. reconstruct the actual 2,640 vertices and 17,040 spatial edges and sum
   curvature locally rather than using only six orbit formulas;
2. evaluate a dust-free Lorentzian action derivative directly, rather than
   recovering it by subtracting the old `P1` term;
3. include known-pass and deliberate-failure controls.

Until those agree, the headline remains provisional.
