# Prior-art and framing gate: finite-height carrier quadratic canonicity

Date: 2026-08-22.

Status: **PRE-CALCULATION GATE.**  No Hessian pullback, parity comparison,
spectrum or target value has been evaluated for the background below.

## Exact object and complete hypotheses

Fix the first positive-height slab generated from the preregistered incoming
state `v=3/2` by the homogeneous cellular 600-cell dust action.  In the
normalization `L_minus=1`, write

```text
h=h1,
q=q1,
lambda=L_plus/L_minus=1+h*q,
rho=h^2,
M=mu(3/2).
```

The already certified representative is

```text
h      = 0.2040549716108237129281133802550001672...,
q      = 9.618002653341898097389452520652463553...,
lambda = 2.962601258380508161014504264675487873....
```

The numerical values are frozen inputs, not targets.  The calculation must
also reconstruct the root from the exact homogeneous equations and compare
it with the committed interval-certified representative.

For each of the two global staircase parities, let `S_p(y)` be the same bare
Lorentzian Regge action with conserved total dust mass and zero cosmological
constant.  Its complete logarithmic signed-squared-edge carrier is

```text
y = 720 old boundary + 840 internal + 720 new boundary
  = 2280 variables.
```

Freeze the 720 old-boundary variables.  Let `H_p` be the `1560 x 1560`
internal-plus-new block of the complete action Hessian at the finite-height
background.  Let

```text
G_p : R^240 -> R^1560
```

be the already derived complete rank-240 scale-plus-strut tangent carrier,
rebuilt at this background from its generic geometric coefficients.  Its
columns are fixed by logical 600-cell vertex labels:

```text
120 upper scale data + 120 strut data.
```

No eigenspace, continuum mode, desired speed or fitted linear combination may
be used to choose or modify `G_p`.

The target-free object is the pair of real bilinear forms

```text
Q_p = G_p^T H_p G_p.
```

They are compared only after the two schedule-specific row carriers have
been identified by physical edge labels and the 240 columns by the same
logical vertex labels.

## Coordinate-covariance caveat

A Hessian away from a critical point is not a tensor.  For a nonlinear
embedding `Phi_p` with first and second jets `(G_p,K_p)`, the actual second
derivative is

```text
d^2(S_p o Phi_p) = G_p^T H_p G_p + grad(S_p) dot K_p.
```

Therefore a bare comparison of `Q_even` and `Q_odd` would be meaningless if
the omitted second-jet terms could differ.  The proposed comparison is
admissible only after the verifier establishes all of the following at the
same background:

1. every artificial/internal diagonal gradient vanishes for both schedules;
2. the old and new boundary gradients agree edge by edge under the physical
   edge identification;
3. the physical old/new boundary embedding is common to both schedules.

Under these conditions the schedule-dependent part of
`grad(S_p) dot K_p` vanishes: schedule-specific second jets occur only in
artificial internal diagonals, whose gradients are zero, while the common
physical-boundary term cancels in the difference.  Thus `Q_even-Q_odd` is a
well-defined necessary test of quadratic schedule independence.  Failure of
any first-gradient condition gives `CONTROL_FAILED`, not a physical result.

The 240-column carrier is presently an exact infinitesimal compatibility
space, not a proved finite nonlinear configuration manifold.  Consequently,
even equality of the two forms would establish only a canonical quadratic
tangent candidate.  It would not prove nonlinear integrability of the
carrier.

## Repository controls and non-duplication

The current repository already establishes the following distinct results.

- The homogeneous cellular action equals both staircase actions throughout
  the regular homothetic domain, including vanishing artificial-diagonal
  equations.  This does not cover nonhomogeneous quadratic variations.
- The generic scale-plus-strut construction gives an exact rank-240
  kinematic carrier at an older near-static background.  It was not pulled
  through the action Hessian in the comparison defined above.
- At that older background, the intersection of the carrier with a
  fixed-input weak canonical graph and all full equations is zero.  This is a
  local fixed-input closure statement, not a parity comparison of quadratic
  forms and not a statement at the new finite-height branch.
- On the older order-24 quotient, even and odd canonical boundary maps agree
  to first order but differ at second nonlinear order on 32 frozen rays.
  That quotient and background are different from the complete 240-column
  carrier and the finite-height state fixed here.
- A refined static barycentric carrier has its own all-schedule Hessian
  programme.  It is not a refinement of this expanding finite-height slab.

An exhaustive repository text search found no existing computation of
`G_p^T H_p G_p` for both parities on this background.  This establishes only
non-duplication inside this repository.

## Primary prior art

### KNOWN

- Dittrich and Steinhaus derive quadratic Regge actions under Pachner moves
  and determine classical triangulation dependence in three and four
  dimensions; see the abstract and the 4D analysis of
  [arXiv:1110.6866](https://arxiv.org/abs/1110.6866), especially the
  quadratic expansion around a background solution.
- Dittrich, Kaminski and Steinhaus explain why 4D discretization independence
  is restricted and connect it directly with the Regge Hessian; see Sections
  I--II and equations (5)--(10) of
  [arXiv:1404.5288](https://arxiv.org/abs/1404.5288).  In particular, a 4D
  Regge action is not invariant under every Pachner move, and the Hessian is
  the relevant quadratic object.
- Bahr and Dittrich show that curved Regge backgrounds generically break
  exact discrete gauge symmetry and yield pseudo-constraints in
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670).
- Improved or perfect actions are a known route to recovering continuum
  dynamics and gauge symmetry on a discretization; see
  [arXiv:0907.4323](https://arxiv.org/abs/0907.4323).
- Hoehn's Sections 6, 10 and 11 distinguish Hessian degeneracies, gauge
  variables and propagating lattice-graviton variables for flat-background
  linearized canonical Regge calculus in
  [arXiv:1411.5672](https://arxiv.org/abs/1411.5672).
- Homogeneous and selected inhomogeneous closed lattice universes are already
  studied in the Collins--Williams framework; see
  [arXiv:1502.03000](https://arxiv.org/abs/1502.03000).  The symmetric
  600-cell dust evolution itself is also established prior art in
  [arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093).

### CONTROL

The known flat-cell subdivision identity predicts exact equality on the
two-dimensional homogeneous tangent inside the 240-column carrier.  Both
parities must reproduce that sub-block.  Deliberately corrupting one
incidence coefficient or one schedule Hessian entry must produce a resolved
nonzero difference.

### OPEN

The literature search did not locate this exact Lorentzian 600-cell dust
background, its two staircase parity representatives, or the complete
rank-240 geometric carrier comparison.  Search absence is not a novelty
proof; external novelty remains **OPEN**.

It is also **OPEN** whether the two forms agree, whether the tangent carrier
integrates nonlinearly, whether a common quadratic form has stable local
modes, and whether any later boundary evolution is schedule independent.

## Falsifiable meaning of the next calculation

- `QUADRATIC_PARITY_DEPENDENT`: a resolved nonzero `Q_even-Q_odd` closes the
  current bare length-only route to a canonical nonhomogeneous linear theory
  on this carrier.  Choosing one parity afterward would be an additional
  theory input.
- `QUADRATIC_PARITY_INDEPENDENT`: equality within a calibrated and hostilely
  tested error bound supplies a necessary schedule-independent quadratic
  tangent candidate.  It does not yet supply a boundary evolution map,
  gravitons, a continuum limit, `c`, `G` or a physical tick.
- `CONTROL_FAILED` or unresolved numerical separation leaves the question
  **OPEN**.

No outcome of this calculation alone is evidence for a desired continuum
dispersion law.
