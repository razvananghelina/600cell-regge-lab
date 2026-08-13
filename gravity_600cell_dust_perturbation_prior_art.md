# Prior-art gate: local perturbations of the 600-cell dust sandwich

Date: 2026-08-13

Status: **targeted primary-literature map completed before preregistration and
before evaluating the dust-solution Jacobian**.

This is not a novelty proof.  The search was targeted at the same carrier,
action, solution and differential operator; absence from these searches leaves
external novelty **OPEN**.

## 1. Exact proposed object

Use the already reproduced De Felice--Fabri time-symmetric dust sandwich on
the fixed 600-cell staircase slab.  For either of the two ordered-schedule
parities, let

```text
x = (30 positive staircase-diagonal squares,
      5 positive pole squared magnitudes),       dim x = 35
y = (30 positive final-boundary edge squares),  dim y = 30.
```

The old boundary is fixed and regular.  The dust mass `M` is fixed to the
published value.  The total one-slab action is the already certified
Lorentzian Regge curvature action plus the published dust world-line term.
Define the 35 local internal equations

```text
R(x,y) = (1/24) * partial S_total / partial x.
```

The proposed calculation is the logarithmic internal Jacobian

```text
J_x = partial R / partial log(x/x0)
```

at the published stationary point `(x0,y0)`, together with the final-boundary
response `J_y`.  If `J_x` is nonsingular, the ordinary implicit-function
theorem gives a locally unique `x(y)` within this 65-coordinate invariant
sector.

This formulation corrects an earlier loose phrase.  It does **not** hold the
lapse fixed: the five pole magnitudes are among the solved internal variables.
Holding all five poles fixed would leave 30 internal unknowns for 35 equations
and would require a separately justified equation/gauge selection.  No such
selection is made here.

## 2. What is already known

### KNOWN: implicit tent evolution and the 600-cell control

Barrett, Galassi, Miller, Sorkin, Tuckey and Williams formulated local
implicit tent evolution, related the underdetermination to approximate
Bianchi identities, and illustrated it on a homogeneous dust-filled
600-cell:

[A Parallelizable Implicit Evolution Scheme for Regge Calculus](https://arxiv.org/abs/gr-qc/9411008)

De Felice and Fabri corrected the schedule to five classes of 24 and allowed
more length variables, but retained extra equalities/shift choices and solved
five successive `4 x 4` systems:

[The Friedmann universe of dust by Regge Calculus](https://arxiv.org/abs/gr-qc/0009093)

[Singularities of the closed RW metric in Regge Calculus: a generalized evolution of the 600-cell](https://arxiv.org/abs/gr-qc/0106077)

Those papers establish the solution used here as a **CONTROL**.  They do not,
in the formulations located, print or diagonalize the complete order-24
`35 x 35` one-slab Jacobian proposed above.

### KNOWN: Hessian, gauge and pseudo-constraint interpretation

For linearized Regge calculus about a flat background, exact vertex
displacements give null directions and first-class constraints.  Beyond that
linearized flat regime the symmetry is generally broken and the canonical
conditions become background-dependent pseudo-constraints:

[From covariant to canonical formulations of discrete gravity](https://arxiv.org/abs/0912.1817)

On curved Regge solutions, Bahr and Dittrich explicitly find no exact generic
gauge symmetry; small Hessian eigenvalues can remain and cause numerical
ill-conditioning:

[(Broken) Gauge Symmetries and Constraints in Regge Calculus](https://arxiv.org/abs/0905.1670)

Hamilton's principal function, pre/post momenta, a priori free data and later
constraints are standard in canonical simplicial gravity:

[Canonical simplicial gravity](https://arxiv.org/abs/1108.1974)

The identification and counting of vertex-displacement and lattice-graviton
modes for linearized four-dimensional Regge calculus is also known:

[Canonical linearized Regge Calculus: counting lattice gravitons with Pachner moves](https://arxiv.org/abs/1411.5672)

Therefore a numerical full rank is not, by itself, a count of physical
degrees of freedom.  It only addresses local solvability for fixed old
boundary and source on this selected finite carrier.

### KNOWN but different: inhomogeneous closed lattice cosmology

Liu and Williams perturb one mass in a Collins--Williams closed lattice
universe and find a well-behaved perturbed evolution:

[Regge calculus models of closed lattice universes](https://arxiv.org/abs/1502.03000)

This is direct prior art for inhomogeneous Regge cosmology.  Its global
continuous-time symmetry reduction, matter placement and perturbation are not
the 35/65-variable Sorkin staircase boundary-value operator above.  Thus it
prevents any broad claim that "inhomogeneous polytopal Regge evolution" is
new, while leaving the exact present matrix comparison **OPEN**.

## 3. KNOWN / CONTROL / OPEN split

### KNOWN

- local implicit Sorkin/tent evolution;
- approximate Bianchi underdetermination and gauge fixing;
- the corrected five-stage 600-cell dust evolution;
- action Hessians and their relation to gauge modes, pseudo-constraints and
  canonical evolution;
- inhomogeneous closed-lattice Regge cosmologies in other reductions.

### CONTROL

- the unrounded published dust mass and time-symmetric sandwich;
- all 35 stationary internal orbit equations, already independently
  reproduced for both schedule parities;
- the already certified 65-coordinate action and boundary response.

### OPEN

- the rank, spectrum and conditioning of `J_x` at this exact dust solution;
- whether both schedule parities give the same local linear response;
- whether all 30 permitted final-boundary directions have locally solvable
  internal continuations in this order-24 invariant sector;
- separation of those directions into gauge and physical modes;
- extension from one slab to multiple ticks;
- any theorem on the unreduced 840-internal-edge carrier;
- external novelty of the exact `35 x 35` computation.

## 4. Framing attack and decision

The proposed result is deliberately narrower than "nonhomogeneous gravity
exists".  The 35/65 coordinates are invariant under the order-24 pointwise
stabilizer of one chosen five-phase schedule; they are not the full 840/720
edge spaces.  Dust remains confined to the 120 prescribed pole world-lines,
its mass is an input, and only one slab is tested.  A full-rank `J_x` would
neither select a perturbation nor derive a clock, mass, Lorentz invariance or
continuum graviton.

Nevertheless, `J_x` is the next clean falsifiable object.  A robustly
nonsingular matrix would justify a local family by the implicit-function
theorem at the analytic level, subject to obtaining a trustworthy numerical
nonsingularity certificate.  A null or badly unresolved direction would be
equally informative and must be compared with the known Regge gauge/Bianchi
structure rather than discarded as a solver inconvenience.

**Decision:** preregister the matrix, coordinate scaling, differentiation
scheme, convergence checks and outcome labels before its first evaluation.
Do not choose a boundary perturbation or solve a displaced point until this
linear audit has been committed.

## 5. Post-result search after the unresolved first Jacobian

The first frozen audit exposed thirty regular directions and five nearly null
pole/lapse directions.  A second primary-source search used those technical
terms rather than the earlier broad query.

Dittrich and Hoehn discuss a symmetry-reduced five-valent tent move for which
the action Hessian has a nonzero eigenvalue associated with the lapse when
curvature is present; that eigenvalue tends to zero toward the flat
configuration.  In the flat case Hessian degeneracy violates the hypotheses
of the implicit-function theorem:

[Canonical simplicial gravity](https://arxiv.org/abs/1108.1974)

This sharpens the **KNOWN** boundary.  Small pole/lapse eigenvalues on a
curved Regge background are not, by themselves, new and are precisely the
published pseudo-constraint phenomenon.  The single-zero-mode interpretation
of a Regge Hessian as a diffeomorphism remnant also has older explicit
examples:

[Linearized dynamics from the 4-simplex Regge action](https://arxiv.org/abs/0707.4513)

What remains **OPEN** is narrower:

- the exact five-dimensional spectrum for this 600-cell dust sandwich;
- whether four relative phase-lapse modes have genuine curvature-induced
  stiffness while the collective mode is exactly null;
- whether that split survives an arbitrary-precision action-only audit;
- the relation, if any, between this split and the five-stage schedule rather
  than a generic tent-move gauge choice.

No located source prints the present `30+5` block decomposition or its
600-cell spectrum.  That absence is not a novelty proof.

## 6. Prior-art gate for the gauge-fixed quotient response

The exact collective lapse path changes the correct linear problem.  The
internal Hessian must be restricted to the complement of its gauge tangent,
and the mixed internal--boundary block must satisfy the corresponding left
null compatibility condition before boundary data can be propagated.

This quotient/gauge-fixing procedure is **KNOWN STRUCTURE**, not a new method.
For flat linearized Regge calculus, vertex-displacement vectors are null
vectors of the Hessian and generate first-class constraints; gauge fixing is
needed before inversion:

[From covariant to canonical formulations of discrete gravity](https://arxiv.org/abs/0912.1817)

Canonical simplicial gravity formulates the same issue through effective
actions, pre/post constraints and the failure of the ordinary
implicit-function theorem when the Hessian is degenerate:

[Canonical simplicial gravity](https://arxiv.org/abs/1108.1974)

Boundary effective Regge actions are obtained by eliminating regular bulk
directions after treating internal null vectors by gauge fixing; their mixed
bulk--boundary coupling must annihilate genuine internal gauge vectors for
the result to be gauge independent.  This is used explicitly in linearized
four-dimensional boundary-graviton calculations:

[Holographic description of boundary gravitons in (3+1) dimensions](https://arxiv.org/abs/1811.11744)

The **OPEN** object is not this linear algebra.  It is the explicit result for
the present curved, dust-filled, five-phase 600-cell slab:

- whether the exact collective lapse tangent annihilates the complete
  `35 x 30` internal--final-boundary mixed block;
- if it does not, which one-dimensional boundary compatibility condition is
  induced;
- whether the `34 x 34` gauge-fixed internal quotient is regular using the
  already certified thirty large and four relative small modes;
- the dimension and spectrum of the resulting local boundary-to-bulk
  response in this order-24 invariant sector.

No located primary source computes those matrices for this carrier.  External
novelty remains **OPEN**.  The next protocol must distinguish all-thirty
boundary compatibility, one boundary constraint, and numerical
non-resolution before inspecting the mixed block.

## 7. Post-result gate for the nearly homogeneous boundary row

The first gauge-quotient run found one nonzero compatibility row whose
normalized direction has cosine `0.999999999770358` with the homogeneous
all-ones final-boundary scale direction.  This numerical alignment was not a
preregistered target and is therefore **PATTERN**, not yet an exact result.

The broad interpretation is **KNOWN**.  In closed Regge--FLRW models, global
variation of an entire orbit of equal spatial edges supplies the homogeneous
scale equation, while variation of the struts supplies the Regge Hamiltonian
constraint/initial-value equation.  Liu and Williams explicitly compare
global with local variation, study when the Hamiltonian constraint is a first
integral, and derive the initial-value equation at the moment of time
symmetry:

[Regge calculus models of the closed vacuum Lambda-FLRW universe](https://arxiv.org/abs/1501.07614)

Consequently, finding that the obstructed boundary direction is homogeneous
scale would not establish a new constraint mechanism.  It would identify the
known cosmological constraint inside this dust-filled, locally varied
600-cell staircase.

The following remain **OPEN** before the precision correction:

- whether the thirty components are exactly equal, rather than merely equal
  within the current fourth-order mixed-difference error;
- whether the orthogonal twenty-nine-dimensional zero-sum shape space is the
  exact compatible boundary subspace;
- whether both schedule parities have a regular 34-dimensional quotient once
  the independently certified high-precision relative Schur form replaces
  the unresolved binary64 weak cluster.

No located source prints this mixed internal--boundary row, its exact
uniformity, or the corresponding 29-direction response for the present
dust-filled 600-cell carrier.  External novelty remains **OPEN**.

## 8. Post-result gate for nonlinear continuation and base refinement

The precision quotient exposed a load-bearing scale separation, but the first
post-result audit of that separation contained a coordinate error.  The four
relative curvatures are about `4.605e-8` in logarithmic coordinates.  The
published residual table records `(partial S/partial x)/24`, not the
logarithmic equation `x*(partial S/partial x)/24` differentiated by the
recorded Hessian.  Projecting the former directly onto the latter falsely
predicted an odd-schedule correction of `8.675e-3`.

After inserting the missing factors `x`, the correctly conditioned predicted
base corrections are

```text
even: 2.521e-10,
odd : 9.019e-7.
```

Both are below the preregistered `1e-5` weak-scale base tolerance.  Thus the
specific claim that the printed odd base is displaced by order `1e-2` is
**RETRACTED / REFUTED BY COORDINATE CONSISTENCY**.  The matrix result at the
printed base was not undermined by that audit.

This danger is **KNOWN** in Regge calculus.  Dittrich and Hoehn show that
higher-order dynamics can break a linearized gauge symmetry: the lowest
nonlinear equations impose consistency conditions on background gauge
parameters, and the quadratic constraints become lapse-dependent
pseudo-constraints:

[From covariant to canonical formulations of discrete gravity](https://arxiv.org/abs/0912.1817)

The more general canonical framework likewise emphasizes that data which are
free in an early move can be fixed a posteriori by later constraints:

[Canonical simplicial gravity](https://arxiv.org/abs/1108.1974)

Hence a nonlinear calculation still cannot assume that the collective path
remains a gauge orbit away from the exactly regular boundary.  The correct
local question remains a Lyapunov--Schmidt reduction.  However, base
stationarity must be audited with logarithmic equations and an
arbitrary-precision error smaller than the weak scale.  The following remain
**OPEN**:

- whether the entire printed collective path is stationary to a precision
  which resolves the four soft modes, rather than merely compatible with zero
  inside a larger finite-difference error;
- whether the collective parameter is free or fixed by the scalar
  pseudo-constraint once the four soft directions are resolved;
- whether the `34+1` tangent response survives recomputation at that refined
  point;
- only after those pass, whether zero-sum boundary perturbations possess
  nearby nonlinear continuations.

These are standard consistency questions, not a novel method.  No located
source performs this weak-scale audit for the present five-stage 600-cell
dust sandwich.
