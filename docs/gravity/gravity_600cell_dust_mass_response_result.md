# Result: conserved dust-mass response is full rank but is not a tangent branch

Date: 2026-08-18

## Headline

On each frozen five-stage 600-cell schedule, the complete zero-total-mass
source carrier has dimension `119` and its forced response has:

```text
outgoing boundary-phase rank = 119,
internal Regge-curvature rank = 119.
```

Nevertheless, that response space is neither the full expanding nor the full
contracting invariant branch of the canonical tangent map.  The primary
target-disclosed comparison classified all

```text
2 schedules x 7 minimal sectors x 2 branches = 28/28
```

comparisons as `SEPARATED`.  A mechanically different audit then refuted
invariance of the response image in all `14/14` schedule/sector cells, at
both `100` and `140` decimal digits.  Since either tangent branch is
invariant, this independently excludes equality with both branches.

The consolidated classification is **DERIVED COMPUTATIONAL,
ADVERSARIALLY CORROBORATED NEGATIVE**.  It closes the simplest proposed
identification of the known `119` strong tangent pairs with the response to
`120-1` conserved comoving mass contrasts.  It does not close matter
coupling, full Regge dynamics or a future extension with dynamical dust
fields.

## Complete hypotheses

The statement above is conditional on all of the following:

- the fixed 2,400-simplex Lorentzian Regge slab and its accepted angle branch;
- the two frozen five-stage schedules;
- 720 old-boundary, 840 internal and 720 new-boundary logarithmic
  squared-edge variables;
- the certified regular pre-Legendre Jacobian at the homogeneous dust
  solution;
- 120 fixed comoving pole world-lines with conserved weights
  `m_v=(M/120)(1+eta_v)` and `sum_v eta_v=0`;
- the unchanged point-particle action
  `S_dust=-8 pi sum_v m_v sqrt(rho_v)`;
- fixed incoming canonical data `(q_old,p_old)`;
- the already certified 1,440-dimensional outgoing tangent map;
- comparison with the preregistered largest- and smallest-modulus invariant
  tangent branches, with no post-result choice of dimension or sector.

The `eta_v` are parameters of the action, not canonical dust variables.
There are no perturbed dust positions, velocities, proper-time fields or
conjugate momenta in this calculation.

## Provenance ledger

| stage | commit |
|---|---|
| primary prior-art and framing gate | `6b587ea` |
| target-disclosed primary protocol | `a0b253d` |
| registered primary verifier | `c1c635d` |
| primary reproducible artifacts | `99f855c` |
| adversarial independence gate | `c0b7091` |
| adversarial protocol | `df924c7` |
| registered adversarial verifier | `5534de2` |
| disclosed archive-array harness correction | `61beabc` |
| adversarial reproducible artifact | `e05f2c0` |

The primary protocol explicitly disclosed the numerical coincidence
`120-1=119` and fixed all `28` comparisons before constructing the response.
The adversarial protocol was committed before evaluating its leakage
matrices.  The only post-launch code correction changed the expected NPZ
array count from `112` to its correct value `224`; the initial run stopped at
that provenance gate before any scientific verdict, and no operator,
threshold or outcome rule changed.

## Result A: action-derived source response

For `z_v=log rho_v`, the only nonzero mixed source derivative is selected by
pole incidence and was derived analytically as

```text
partial^2 S_dust / (partial z_v partial eta_v)
    = -4 pi (M/120) sqrt(rho_v).
```

The complete source has exactly 120 nonzero incidence-selected entries.  The
uniform source is retained as a control, while the physical comparison uses
its fixed 119-dimensional zero-sum complement.  In each binary-tetrahedral
sector the forced response is the solution

```text
Y_m = -J^-1 B_m,
R_m = (delta q_new, delta p_new).
```

The solved-variable rank is only an algebraic control: regular `J` maps an
injective source to an injective solution.  The contentful results are that
projection to outgoing phase does not lose a direction and that applying the
independently reconstructed Jacobian of all 3,840 internal causal deficits
also does not lose a direction.

Restored over the regular representation and both schedules, the censuses
are:

```text
phase:     238 NONZERO_RESOLVED, 0 ZERO_CONSISTENT, 0 OPEN,
curvature: 238 NONZERO_RESOLVED, 0 ZERO_CONSISTENT, 0 OPEN.
```

Thus each schedule has phase rank `119` and curvature rank `119`.

The phase singular spectra are schedule-robust in all seven minimal sectors.
The raw Euclidean singular spectrum of the internal-deficit response is
schedule-robust in one sector and schedule-separated in six.  This does not
change its rank.  The two schedules triangulate the slab interior differently,
so their unweighted deficit row spaces and Euclidean row norms are not a
common physical norm.  Equality of those raw singular values is therefore
not promoted to an invariant requirement.

The primary principal-angle distances range from approximately
`2.3075e-6` to `1.0451e-3`.  All 28 pass the frozen `distance > 100 epsilon`
separation rule.  The weakest margin is only approximately `105.59 epsilon`,
which is why the independent test below is load-bearing.

Two complete primary executions reported `20/20` and wrote byte-identical
artifacts:

```text
JSON SHA-256 48de8f4a9edabb84145d3ce960aab808fa45c13431ff30e8756a4314f9e1ef60
NPZ  SHA-256 ae550de064c7853cba8f5b1375276a6809d8a4631bf3abfaa656ea9a05555af6
```

## Result B: independent invariant-subspace falsification

The audit did not reconstruct, order or compare with either spectral branch.
For each frozen response matrix `R` and tangent block `T`, it instead formed

```text
C = (R* R)^-1 R* T R,
L = T R - R C.
```

If `im R` equalled either invariant branch, necessarily `L=0`.  Every stored
binary64 component of `R` and `T` was converted to its exact dyadic rational,
and the complete calculation was repeated in Flint ball arithmetic at `100`
and `140` decimal digits.

The result is:

```text
14/14 actual cells: INVARIANCE_REFUTED at 100 digits,
14/14 actual cells: INVARIANCE_REFUTED at 140 digits,
14/14 precision balls overlap,
14/14 Gram determinants exclude zero.
```

Exact source-column rephasing, a simultaneous `q/p` block swap and time
reversal through `T^-1` preserve the non-invariance verdict in every cell.
Four synthetic size classes discriminate an exactly invariant identity
control from an explicitly leaking control.  A separate SciPy pivoted-QR
implementation gives finite nonzero leakage in all `56/56` derivative
variants, with ratios from approximately `2.3063e-6` to `4.2621e-6`.

Two complete executions reported `7/7` and wrote the same artifact:

```text
SHA-256 4c1fa9a660b09696bf559ce45aa8da8ba1b1a51dfa925924fa9cd2bb56e9684c
```

This audit is independent at the decisive subspace test.  It is conditional
on the frozen response matrices and does not independently rederive the
dust-source column from the action.  That scope limit is explicit rather than
silently promoted to a second end-to-end derivation.

## Physical verdict

- **DERIVED COMPUTATIONAL:** conserved zero-sum mass contrasts generate a
  119-dimensional outgoing geometric response on each schedule.
- **DERIVED COMPUTATIONAL:** every such source direction changes internal
  Regge curvature at first order.
- **DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED NEGATIVE:** the response
  image is not either complete strong tangent branch.
- **STRUCTURAL:** the equality `119=120-1` is an equality of counts only; it is
  not a subspace identification or a particle/gravity correspondence.
- **CLOSED:** the proposal that the 119 strong pairs are simply the response
  to frozen comoving density weights.
- **OPEN:** which constraint-reduced modes, if any, are continuum scalar,
  vector or tensor perturbations.
- **OPEN:** full dynamical dust, refinement, a continuum limit and external
  novelty of the project-specific response calculation.
- **NOT ESTABLISHED:** gravitons, inertia, a limiting speed, Planck units or
  particle masses.

The negative does not say that matter fails to gravitate: the mass response
is injective and curvature-changing.  It says that a source sensitivity at
fixed incoming data is not itself a closed propagating tangent branch.

## Post-result primary-source check

The refined search used the now-available terms forced canonical response,
mass perturbation, tangent evolution and invariant subspace.  It again found
the single-perturbed-mass lattice-universe calculation of Liu and Williams,
especially Section IV, where local Regge equations are derived for one
unequal mass:

- R. G. Liu and R. M. Williams, *Regge calculus models of closed lattice
  universes*, arXiv:`1502.03000`, DOI `10.1103/PhysRevD.93.023502`.

That is direct prior art for inhomogeneous Regge point masses, but it does not
construct the present complete `119 -> 1440` fixed-input response or compare
its image with invariant branches of a full canonical tangent map.

The canonical interpretation continues to be governed by Sections 2 and 6
of:

- B. Dittrich and P. A. Hoehn, *Canonical simplicial gravity*,
  arXiv:`1108.1974`, DOI `10.1088/0264-9381/29/11/115009`.

There a discrete action generates pre/post momenta and momentum matching,
while physical constraint/gauge directions must be established from the
action.  That framework supports the present caution: a forced response is
not automatically an invariant propagating degree of freedom.

Finally, Brown--Kuchar dust contains proper-time and comoving-label fields
with conjugate momenta:

- J. D. Brown and K. V. Kuchar, *Dust as a Standard of Space and Time in
  Canonical Quantum Gravity*, arXiv:`gr-qc/9409001`, DOI
  `10.1103/PhysRevD.51.5600`.

Adding such fields would be a new model input, not an interpretation of the
119 frozen mass weights already computed.  The refined search located no
primary source with the exact finite response/invariance calculation here,
but a search is not a novelty proof.  External novelty remains **OPEN**.

Only the two targeted mass-response verifiers and the documentation/registry
guards are required for this result.  The full suite is deliberately not
run.
