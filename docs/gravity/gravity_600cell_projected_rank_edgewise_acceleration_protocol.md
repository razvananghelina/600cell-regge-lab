# Protocol: blind homogeneous acceleration on the canonical projected carrier

Date: 2026-08-19

Prior-art commit: `1921519`

## Separation from the comparison target

Stage A computes and commits the two finite-carrier coefficients without
loading or printing:

- the continuum closed-dust coefficient;
- any projected-red coefficient;
- any error, ratio, ranking or convergence verdict relative to either.

Only after the Stage-A JSON is committed may Stage B load the continuum
target and compare the two frozen values.  This ordering prevents a numerical
choice made while looking at the desired answer.

## Frozen inputs

```text
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f

reproducible/verify_gravity_600cell_projected_rank_edgewise_carrier.py
  50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23

reproducible/gravity_600cell_projected_rank_edgewise_carrier.json
  b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84

reproducible/verify_gravity_600cell_projected_refinement_acceleration_blind.py
  e88111adaeb333abf80b68e06e23d7840ef14399238ada9d0f3cd722d7934e50
```

The last file freezes the already audited action, branch and coefficient
extraction formula.  The new verifier may copy those formulas but must build
the rank-edgewise carriers independently rather than import and execute the
old projected-red verifier.

## Frozen carriers

Reconstruct exactly

```text
K_0 = P(sd K_600),
K_1 = P(Esd_2(sd K_600)).
```

The builder must recover the frozen f-vectors and unit-sphere chordal-volume
scalars:

```text
K_0: (2640,17040,28800,14400),  V_bar=19.147932918312847
K_1: (19680,134880,230400,115200), V_bar=19.583480465413963.
```

The comparison tolerance for each volume is `5e-11` absolute.  All
tetrahedral volumes must be positive, every triangle must have incidence two,
and each Euler characteristic must be zero.

## Frozen action and normalization

For each unit-sphere carrier:

1. set the scale for unit volume radius to

   ```text
   s0=(2*pi^2/V_bar)^(1/3);
   ```

2. compute the direct three-dimensional Regge curvature

   ```text
   C=s0*(2*pi*sum_edges(l_e)-sum_tet,local_edge(l*theta));
   ```

3. select the total global dust mass once as

   ```text
   M=C/(8*pi);
   ```

4. evaluate the complete cellular Lorentzian frustum action, including both
   spatial-boundary contributions, plus `-8*pi*M*sqrt(rho)`;
5. use the contracting time orientation and define

   ```text
   s_plus=s0*exp(a*eta^2), rho=eta^2.
   ```

No local dust distribution or local lapse is introduced.

## Frozen extraction

Use exactly:

```text
seam eta values       (0.04,0.02,0.01,0.005)
lapse eta values      (0.04,0.02,0.01)
coefficient sentinels (0,-1,-2,-3)
seam log steps        2e-5 and 1e-5
lapse real-log steps  4e-3 and 8e-3.
```

For the seam residual, apply two levels of fourth-order Richardson
extrapolation and solve the frozen affine equation in `a`.  For the lapse
residual, apply one fourth-order Richardson extrapolation, remove the exact
static root and solve the remaining quadratic dynamic root.  The verifier
must record both routes and both derivative-step variants.

## Preregistered controls and thresholds

For each carrier require:

- static gravitational-action relative error `<5e-9`;
- static total-action relative residual `<5e-8`;
- static lapse relative residual `<5e-8`;
- bottom and top `pi/2` residuals `<5e-9`;
- tetrahedron-relabel action relative error `<5e-9`;
- seam extrapolation difference `<2e-6`;
- seam derivative-step difference `<2e-6`;
- lapse derivative-step difference `<5e-5`;
- seam/lapse coefficient difference `<5e-5`;
- seam affine residual `<2e-6`;
- lapse quadratic residual `<2e-6`;
- lapse static-root residual `<2e-6`;
- maximum imaginary residual `<5e-8` on real states;
- Lorentzian inertia `(3,1)` for every sampled frustum.

These are the unchanged successful gates of the earlier direct-action
calculation.  They are fixed before either new coefficient is evaluated.

## Stage-A output and outcome hierarchy

Write

```text
reproducible/gravity_600cell_projected_rank_edgewise_acceleration_blind.json
```

with provenance, carrier diagnostics, selected masses, the complete seam and
lapse audit, and the two coefficient values.

Exactly one Stage-A outcome is allowed:

1. `CANONICAL_CARRIER_ACCELERATION_INTERNAL_FAILURE` if any frozen control
   fails;
2. `CANONICAL_CARRIER_ACCELERATION_COEFFICIENTS_DERIVED` otherwise.

Stage A makes no claim about improvement, convergence or physics.

## Stage-B acceptance and kill boundary

After the Stage-A artifact is committed, compare both coefficients with the
continuum value `-1/2`.

- **CALIBRATION PASS:** the absolute fine-level error is strictly smaller
  than the base-level error and both Stage-A carrier audits passed.
- **CALIBRATION FAILURE:** otherwise.  Do not proceed to the large
  inhomogeneous Hessian until the failure is understood.

Even a pass is only **DERIVED NUMERICAL CONTROL**.  It does not establish an
infinite limit, gravitational waves, a limiting speed, a physical tick,
inertia, Planck scales or particle masses.
