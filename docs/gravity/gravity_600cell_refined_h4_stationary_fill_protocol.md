# Protocol: all-schedule H4 stationary-fill census

Date: 2026-08-20

Prior-art gate commit: `5518fa7`.

This protocol is frozen before computing any of the 240 internal residuals.

## 1. Frozen inputs

Hash and use only:

```text
commons/cell600.py
reproducible/verify_gravity_global_regge_orbits.py
reproducible/gravity_600cell_refined_canonical_map_feasibility.json
reproducible/verify_gravity_600cell_refined_canonical_map_feasibility.py
reproducible/gravity_600cell_projected_rank_edgewise_acceleration_blind.json
reproducible/verify_gravity_600cell_projected_rank_edgewise_acceleration_blind.py
reproducible/gravity_600cell_projected_rank_edgewise_local_dust.json
reproducible/verify_gravity_600cell_projected_rank_edgewise_local_dust.py
reproducible/gravity_600cell_projected_rank_edgewise_balanced_slab.json
```

Their expected SHA-256 values, in that order, are

```text
ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf
ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e
36fba835048e6e0f0676b749192a9d882406932770a00ba1396929bbc4d04a32
2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2
496ee770ad06cf4a7f0bca79153042cca7e2821c179bdb5b27f1fdb9f393ba2b
53463e5271301ae41eb26564875d26991ddea8024a9e09ae3c302d428ad39779
e1064380d8580e458ebfcb990285181b2ca5b092f7738a6729b0feec528df2e0
0a9e9e796cd671c82f2e428bfa21ba63ccb07fe76867e4553979c3c54b22a0d5
```

Only `K0=P(sd K_600)` is tested.  No result from `K1` is inferred.

## 2. Exact rank geometry

Do not use rounded projected coordinates to determine the six chamber edge
lengths.  Put

```text
phi=(1+sqrt(5))/2,  c=phi/2,
N_k=k(1+(k-1)c).
```

For nested faces of sizes `a<b` in one regular 600-cell tetrahedron, use

```text
u_a . u_b = a(1+(b-1)c)/sqrt(N_a N_b),
l_ab^2 = 2-2(u_a . u_b).
```

Ranks are `r=a-1`.  Reconstruct the chamber tetrahedral volume, total chordal
volume, spatial Regge curvature and

```text
s0=(2*pi^2/V_chord)^(1/3),
M=R_Regge/(8*pi)
```

at high precision.  Compare `V_chord` and `M` with the frozen blind artifact.
Because every barycentric chamber contains one vertex of each rank and all
chambers are congruent, the `P1` dust mass in each rank class must be `M/4`;
also reconstruct and check this against the frozen local weights.

## 3. Frozen variables and fill

For every one of the 24 colour orders use these 22 positive coordinates:

```text
q^-_rs, q^+_rs  for 0<=r<s<=3,      12 boundary coordinates;
d_rs             for 0<=r<s<=3,       6 cross diagonals;
rho_r            for r=0,1,2,3,        4 vertical coordinates.
```

The Lorentzian squared edge assigned to a vertical coordinate is `-rho_r`.
At the inherited induced flat static fill,

```text
q^-_rs=q^+_rs=s0^2 l_rs^2,
rho_r=tau0^2,              tau0=0.0102,
d_rs=s0^2 l_rs^2-tau0^2.
```

No coordinate is fitted or solved in this census.

## 4. Action and exact combinatorial reduction

Use the corrected complex Lorentzian angle branch and triangle-area
derivatives frozen in `verify_gravity_global_regge_orbits.py`.  Boundary
triangles start with `pi`; internal triangles start with `2*pi`.  Include

```text
S_dust=-8*pi*sum_v m_v sqrt(rho_rank(v)).
```

For every schedule explicitly enumerate its 57,600 pentachora and all triangle
incidences.  A state type is `(rank,layer)`.  The reduction is admissible only
if every triangle with the same state type has exactly the same multiset of
abstract incident simplex/hinge types.  Record all type counts and reject any
mixed type.

The reduced action and analytic log-coordinate gradient follow from those
integer multiplicities and the Schlaefli identity.  Check the action and all
22 analytic derivatives independently by centred finite differences on the
lexicographically first order and its reverse.  Use 100 decimal digits,
steps `1e-15` and `5e-16`, and Richardson extrapolation.  The relative
analytic/Richardson disagreement must be below `1e-24`.  Repeat the analytic
residual census at 140 decimal digits.

## 5. Residual certificates

For each internal orbit report both the total and per-edge log residual.  Let
`g100` and `g140` be the 100- and 140-digit per-edge values and let

```text
scale=max(1, all absolute per-edge residuals),
epsilon=100*max_abs(g100-g140)+1e-60*scale.
```

Classify an entry as certified nonzero only when

```text
abs(g140) > epsilon.
```

Otherwise label it unresolved/zero-compatible; do not silently round it to
zero.  Report separately the six cross-diagonal and four vertical results,
the full multiset of 24 ten-vectors, and the number of distinct vectors under
the same envelope.

Mandatory collective controls:

1. all 24 total actions agree at the induced fill within the precision
   envelope;
2. reversing an order and exchanging the two equal boundaries gives the same
   internal residual vector;
3. the common induced-lapse directional derivative

   ```text
   sum_r dS/dlog(rho_r)
   - sum_rs (rho/d_rs) dS/dlog(d_rs)
   ```

   is zero-compatible for every schedule and reproduces the frozen static
   lapse equation.  This control prevents mistaking a sign or dust
   normalization error for extra equations.

## 6. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_H4_STATIONARY_FILL_CONTROL_FAILED` if provenance, topology,
   exact-geometry, branch, combinatorial reduction, finite-difference,
   time-reversal or collective-lapse controls fail.
2. `REFINED_H4_INDUCED_FILL_OFF_SHELL` if at least one of the 240 internal
   entries is certified nonzero.
3. `REFINED_H4_INDUCED_FILL_STATIONARY_CANDIDATE` only if every internal
   entry is zero-compatible at both precisions.

Outcome 2 licenses a separately preregistered ten-variable stationary-fill
solve for all 24 schedules; it does not kill the 600-cell route.  It does
invalidate an effective Hessian evaluated at the inherited fill.

Outcome 3 licenses an all-schedule `H4` effective-Hessian comparison, but it
does not establish full non-invariant stationarity or schedule independence.

## 7. Exclusions and deliverables

Forbidden: a root solve, Hessian, eigenvalue, continuum comparison, schedule
selection/average, physical-mode count, `c`, `G` or Planck interpretation.

Deliver a registered verifier, deterministic JSON artifact, result note,
targeted rerun and static registry audit only.  Do not run the full suite.
