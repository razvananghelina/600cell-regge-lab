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
