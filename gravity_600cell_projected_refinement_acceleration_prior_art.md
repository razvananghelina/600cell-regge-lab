# Prior-art gate: direct Regge acceleration on projected 600-cell refinements

Date: 2026-08-17

## Exact question and hypotheses

Let `K_0` be the regular 600-cell boundary and let `K_1,K_2` be the two
already frozen projected red refinements from
`smooth_hopf_red_refinement_preregistration.md`.  Every mesh edge has one
shared midpoint, the midpoint is projected radially to the unit `S3`, and
each tetrahedron is split into four corner tetrahedra and four tetrahedra
along the geometrically shortest diagonal of its central octahedron.  The
three carriers have respectively

```text
level        vertices        tetrahedra
0                 120                600
1                 840              4,800
2               6,480             38,400
```

The question is:

> For the direct cellular Lorentzian Regge action of a homothetic one-slab
> evolution of each complete, irregular tetrahedral carrier, does the
> weak-lapse acceleration at a time-symmetric closed-dust turning point move
> toward the continuum closed-FLRW value under spatial refinement?

The hypotheses are fixed as follows.

1. Every spatial vertex `p_v` lies on the unit sphere in `R4`.  A homogeneous
   spatial slice has vertices `s p_v`; only the common scale `s` evolves.
2. Each tetrahedron sweeps out its flat Lorentzian tetrahedral frustum.  The
   gravitational action is the direct cellular Regge action: timelike
   trapezoids over spatial edges plus the lower and upper spatial-boundary
   triangle terms.  No averaged Schlaefli symbol or averaged angle is used.
3. Every vertex strut has the same squared proper length `rho>0`.  This is a
   global homogeneous lapse.  Local lapse and anisotropic edge equations are
   outside this mission.
4. Dust has one conserved total mass and contributes
   `-8*pi*M*sqrt(rho)`.  This tests the same global dust minisuperspace used by
   the exact fixed-600-cell theorem; it does not specify a local dust
   distribution on the refined vertices.
5. The comparison radius is the volume radius, not the embedding
   circumradius and not a selected edge.  If `V_j` is the total Euclidean
   tetrahedral volume for the unit-sphere positions, define

   ```text
   R = s * (V_j/(2*pi^2))^(1/3).
   ```

   Thus `R=1` means that the piecewise-flat slice has the volume `2*pi^2` of
   the unit round `S3`.  On the regular carrier this reduces exactly to the
   repository's already frozen map
   `L/R=(pi^2*sqrt(2)/50)^(1/3)`.
6. At each refinement level separately, the conserved mass is fixed before
   the dynamic coefficient by the exact static lapse constraint at `R=1`:

   ```text
   M_j = (1/(8*pi)) * sum_edges l_e delta_e,
   ```

   where `delta_e` is the three-dimensional Regge deficit of the spatial
   slice.  Holding the coarse mass fixed would compare different turning
   radii and is therefore not the stated test.
7. The same time-symmetric half-step convention as in
   `gravity_600cell_cellular_weak_lapse_protocol.md` is used.  With
   `eta=tau/R`, the continuum target is

   ```text
   log(R_1/R_0) = -(1/2)*eta^2 + O(eta^4).
   ```

No coefficient from levels 1 or 2 has been evaluated while writing this
gate.

## Primary prior art

### Collins--Williams and direct finite-resolution models

Collins and Williams introduced regular 5-, 16- and 600-cell Cauchy surfaces
and flat world-tube/frustum blocks for closed Friedmann dynamics:
[Dynamics of the Friedmann Universe Using Regge
Calculus](https://doi.org/10.1103/PhysRevD.7.965).

De Felice and Fabri explicitly evolved a dust-filled 600-cell with a Sorkin
scheme and found the familiar finite-carrier stopping/causality issue:
[The Friedmann universe of dust by Regge
Calculus](https://arxiv.org/abs/gr-qc/0009093).

Liu and Williams compared regular Collins--Williams models with Brewin-type
subdivided carriers.  They distinguish global from local variation and find
resolution-dependent deviations from continuum closed FLRW:
[Regge calculus models of the closed vacuum FLRW
universe](https://doi.org/10.1103/PhysRevD.93.024032).

### Projected 600-cell geodesic domes

Tsuda and Fujiwara subdivide cells of a regular 600-cell, project the new
vertices onto its circumsphere, and call the resulting carrier a geodesic
4-dome.  For frequency two, their untriangulated cell consists of four corner
tetrahedra and one central octahedron.  This is closely related to one level
of the present red refinement before the central-octahedron diagonal is
chosen.  Their published dynamics replaces the increasingly complicated
direct irregular action by a *pseudo-regular* model with averaged angular
data.  That averaged model tends to continuum FLRW at infinite frequency:
[Oscillating 4-Polytopal Universe in Regge
Calculus](https://doi.org/10.1093/ptep/ptab079).

The paper explicitly warns that a naive averaged Schlaefli symbol does not
even have the right smooth limit and repairs it by a different angular
averaging prescription.  Consequently, its convergence theorem is evidence
for the broad geodesic-dome programme, but it is not a calculation of the
direct non-averaged irregular action proposed here.

### Canonical and refinement framing

Dittrich and Hoehn show that a complete one-step Regge action is Hamilton's
principal function for the corresponding discrete canonical update, while
refinement moves can add data and constraints:
[Canonical simplicial gravity](https://arxiv.org/abs/1108.1974).  This is why
the complete boundary term and not merely a spatial curvature sum is required
for the acceleration calculation.

## KNOWN / CONTROL / OPEN

### KNOWN

- Regular-polytopal and dust 600-cell Friedmann dynamics in Regge calculus is
  established prior art.
- Projecting subdivided 600-cell vertices to the circumsphere is established
  geodesic-dome prior art.
- Pseudo-regular angular averages for increasing frequency approach FLRW in
  the cited Tsuda--Fujiwara model.
- Complete one-step actions generate discrete canonical momenta.
- The repository already proves that the selected projected red meshes are
  closed, shape-controlled through two levels, and converge on specified
  round-`S3` spectral calibration modes.

### CONTROL

- Level zero must reproduce the existing exact regular cellular action and
  its coefficient ratio `1.078979468041351...`.
- At every level the static action must reduce independently to
  `tau*sum_e l_e delta_e`, the spatial-boundary angles must be `pi/2`, and the
  selected mass must cancel the static lapse equation.
- The action must be real on the chosen Lorentzian branch, invariant under
  tetrahedron vertex relabelling, and insensitive to the central-octahedron
  diagonal only to the extent implied by the already fixed spatial mesh.
- The volume-radius normalization must reproduce the exact regular
  `L/R` conversion at level zero.

### OPEN

- The direct, non-averaged weak-lapse coefficient on `K_1` and `K_2`.
- Whether those two coefficients move monotonically toward `-1/2`.
- Whether the present repeated red tower has a uniform infinite refinement
  limit; two refined levels cannot prove it.
- Local refined-lapse equations, a dust discretization per vertex,
  anisotropic stability, tensor modes and a full GR continuum limit.
- External novelty of the narrow direct irregular-dust audit.  A primary
  search found no identical coefficient calculation, but absence from a
  search is not proof of novelty.

## Framing attack

A positive trend is not a new theory of gravity.  It would be a finite
convergence control for a known Regge approximation.  Moreover, the
homogeneous global variation can only test the Friedmann scale mode; symmetry
removes the local gravitational degrees of freedom.  Even exact convergence
of this coefficient would not establish gravitational waves or general
Einstein dynamics.

Conversely, failure to converge would refute this selected projected-red,
global-dust implementation.  It would not refute Regge calculus or every
possible refinement.  The distinction is part of the acceptance boundary,
not a qualification to be added after seeing the numbers.
