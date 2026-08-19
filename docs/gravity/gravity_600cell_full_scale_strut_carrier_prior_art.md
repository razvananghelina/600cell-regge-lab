# Prior-art gate: complete scale--strut boundary carrier

Date: 2026-08-19

## Exact mission

On the fixed nonstatic Lorentzian 600-cell slab, construct the linear map
from 120 vertex-scale coordinates and 120 logarithmic strut-magnitude
coordinates into the 840 internal plus 720 upper-boundary log squared edge
coordinates.  Cross-diagonal response must be derived from local
squared-length geometry and shared-face compatibility, not chosen from a
rank or representation target.

The candidate formula is already disclosed in
`gravity_600cell_full_scale_strut_carrier_exploratory_disclosure.md` because
it arose during the coordinate audit before this gate.  The next protocol is
therefore target-disclosed and must attempt to falsify it.

Complete generic hypotheses are:

```text
regular spacelike tetrahedral lower cells;
Lorentz metric diag(1,1,1,-1);
q_i=lambda p_i+tau n with lambda!=1 and tau!=0;
(lambda-1)^2-3 tau^2 != 0;
raw signed squared-length differentiation;
physical shared vertices and exact compatible affine face transitions;
no independent face-connection variable beyond the existing length data.
```

For the curved-slab instantiation additionally require

```text
L0^2>0, rho>0, q_diag=lambda L0^2-rho>0,
the two frozen staircase schedules describe the same accepted background,
and the carrier is expressed in the action's frozen logarithmic coordinates.
```

No `lambda=1` limiting statement belongs to this mission.

## Primary literature

### KNOWN

- Dittrich and Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974), derive discrete canonical
  evolution from Hamilton's principal function and simplex gluing/Pachner
  moves.  Their abstract explicitly notes that data free at one move may be
  fixed later by constraints.  This supports the general mechanism, not our
  carrier coefficients.
- Dittrich and Ryan, [*Simplicity in simplicial phase
  space*](https://arxiv.org/abs/1006.4295), explain how gluing conditions can
  arise from secondary simplicity constraints in a different phase-space
  formulation.  It is relevant to the distinction between metric boundary
  data and extra connection data, but it does not supply the present
  length-response map.
- Jercher and Steinhaus, [*Cosmology in Lorentzian Regge
  calculus*](https://arxiv.org/abs/2312.11639), study homogeneous Lorentzian
  4-frusta with dynamical height and matter and recover expanding/contracting
  branches under stated causality conditions.  Their carrier is
  symmetry-reduced and cuboidal, not the nonhomogeneous tetrahedral
  600-cell carrier here.
- The repository's accepted two-frustum audit proves that, with no new
  connection variable, compatible local strut-preserving Poincare motions
  are diagonal.  The accepted universal local-lift audit proves exact local
  vertex-star support at `(lambda,tau)=(2,5),(3,11)`.

### CONTROL

- one-frustum compatibility must leave exactly the two relations
  `A+B=8` and `C+D=1`, rather than falsely claiming uniqueness;
- two-cell shared-face compatibility must retain the one-dimensional full
  Poincare face-stabilizer control during elimination;
- the candidate must reproduce both frozen rational local-response blocks;
- its strut half must equal the accepted corrected pure-strut carrier;
- summing all scale columns must give the direct homogeneous derivative;
- exact support, rank and staircase-schedule equivariance must be checked
  without reading an action/Hessian target;
- a deliberately corrupted endpoint coefficient must be rejected.

### OPEN

- whether an independent construction reproduces the four generic
  endpoint coefficients;
- whether the generic formulas remain nonsingular and correctly conditioned
  at the accepted curved slab;
- the exact rank and intrinsic singular-value census of the resulting
  `1560 x 240` numerical carrier;
- whether pulling it through the frozen Regge action derivative leaves any
  canonical direction;
- gauge/constraint/dynamical classification and continuum mode labels;
- external novelty of this explicit nonhomogeneous 600-cell formula.

The pre-computation search found no primary source giving this exact
endpoint-supported response on a tetrahedral 600-cell slab.  Search absence
does not establish novelty; external novelty remains **OPEN**.

## Framing attack

The expression “full 240-dimensional dynamic carrier” would prejudge the
physics.  Before an action is applied it is only the **complete kinematic
boundary-data carrier** in the frozen coordinate choice.  Rank 240 would
prove independent parametrization, not 240 physical degrees of freedom.

Likewise, large coefficients near `lambda=1` may indicate a singular
coordinate chart rather than large physical response.  The verifier must
report conditioning and homogeneous identities; no physical conclusion may
be inferred from coefficient magnitude alone.

## Gate

Only after this file and the exploratory disclosure are committed may a
target-disclosed protocol be written.  Only after that protocol is committed
may the independent geometric verifier be implemented and run.  The action,
Hessian, strong equations, representation-sector targets and continuum
labels remain forbidden until the geometric artifact is frozen.

