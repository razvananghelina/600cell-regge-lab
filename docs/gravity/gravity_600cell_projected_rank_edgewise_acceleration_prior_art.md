# Prior-art gate: homogeneous Regge calibration on the canonical projected carrier

Date: 2026-08-19

## Exact object, operator, carrier and hypotheses

This mission asks a deliberately narrow calibration question before any
inhomogeneous Hessian is attempted.

Let

```text
K_0 = P(sd K_600),
K_1 = P(Esd_2(sd K_600)),                 P(x)=x/||x||,
```

be the two complete carriers certified in
`gravity_600cell_projected_rank_edgewise_carrier_result.md`.  Their f-vectors
are respectively

```text
(2640,17040,28800,14400),
(19680,134880,230400,115200).
```

On each carrier, use the same direct cellular Lorentzian Regge action already
audited in the projected-red calculation.  Every spatial vertex is
homothetically scaled by one common variable `s`; every strut has common
squared proper length `rho`; every tetrahedron sweeps out its own flat
Lorentzian frustum; all lateral hinges and both spatial-boundary terms are
summed without angular averaging.

The matter term is the same global minisuperspace dust term

```text
-8*pi*M*sqrt(rho).
```

At each level separately, the one conserved total mass is selected by the
exact static global lapse constraint at unit volume radius.  This is a global
dust control only.  It does not define local dust weights and cannot support
an inhomogeneous matter Hessian.

With

```text
R = s * (V_bar/(2*pi^2))^(1/3),
eta = tau/R,
log(R_1/R_0) = a*eta^2 + O(eta^4),
```

the blind stage will compute the two coefficients `a_0,a_1` and freeze them
before loading the continuum value.  A later comparison stage may test them
against the closed-dust FLRW half-step coefficient `-1/2`.

## Primary prior art

Regular 5-, 16- and 600-cell closed Friedmann evolutions with frustum-like
world tubes are established Regge-calculus constructions:

- P. A. Collins and R. M. Williams, *Dynamics of the Friedmann Universe Using
  Regge Calculus*, DOI `10.1103/PhysRevD.7.965`.

Subdivided closed-FLRW Regge carriers, the distinction between global and
local variation and their finite-resolution behaviour are also established:

- R. G. Liu and R. M. Williams, *Regge calculus models of the closed vacuum
  Lambda-FLRW universe*, DOI `10.1103/PhysRevD.93.024032`.

Projected subdivisions of the 600-cell and their pseudo-regular continuum
model are established geodesic-4-dome prior art.  The authors explicitly
replace the cumbersome direct irregular Regge system by angular averages:

- R. Tsuda and T. Fujiwara, *Oscillating 4-Polytopal Universe in Regge
  Calculus*, arXiv:`2011.04120`, DOI `10.1093/ptep/ptab079`, Section 6.

The later inhomogeneous problem is not just a larger homogeneous calculation.
On curved Regge backgrounds, exact diffeomorphism gauge symmetry is generally
broken and canonical constraints become background-dependent
pseudo-constraints:

- B. Bahr and B. Dittrich, *(Broken) Gauge Symmetries and Constraints in
  Regge Calculus*, arXiv:`0905.1670`, DOI
  `10.1088/0264-9381/26/22/225011`.

Regge Hessians are a standard tool for linearized dynamics and gauge-mode
analysis:

- B. Dittrich, L. Freidel and S. Speziale, *Linearized dynamics from the
  4-simplex Regge action*, arXiv:`0707.4513`, DOI
  `10.1103/PhysRevD.76.104020`.

The search found no primary source evaluating the complete direct irregular
action on the exact projected rank-edgewise carrier above.  Search absence is
not a novelty proof; external novelty remains **OPEN**.

## KNOWN / CONTROL / OPEN

### KNOWN

- Regge approximations to closed FLRW and projected 600-cell geodesic domes
  are established prior art.
- The repository already derives a two-level improvement toward FLRW for four
  disclosed projected-red regulators.
- The repository now derives that none of those projected-red diagonal rules
  is `H4`-equivariant, and separately certifies the rank-edgewise replacement.
- The continuum target and the previous projected-red coefficients are not
  needed to compute the new blind coefficients.

### CONTROL

- Both carriers must reproduce their frozen topology and geometry hashes or
  frozen scalar certificates before the action is evaluated.
- The static cellular action must equal `tau*sum_e l_e delta_e`; the boundary
  angles must be `pi/2`; the selected dust mass must cancel the static global
  lapse equation.
- Every sampled frustum metric must have Lorentzian inertia `(3,1)` and the
  action must be real on the chosen branch.
- The seam and lapse routes must independently reconstruct the same dynamic
  coefficient within frozen numerical tolerances.
- A tetrahedron relabelling must leave the action invariant.

### OPEN

- The two blind coefficients on the canonical carriers.
- Whether the fine coefficient is closer to `-1/2` than the base coefficient.
- Any asymptotic convergence order; one refinement cannot establish one.
- A local dust discretization, local lapse equations, a constraint-reduced
  inhomogeneous Hessian, tensor modes and a propagation speed.

## Framing attack and decision

This calibration cannot discover new gravitational dynamics.  It varies only
one global scale and one global lapse, so the Friedmann sector is built into
the ansatz.  A positive result merely shows that the new canonical carrier is
compatible with the already known homogeneous Regge continuum trend.

Nevertheless it is a necessary gate.  It prevents us from investing in a
large anisotropic Hessian on a carrier whose direct action already fails the
known homogeneous control.

The outcome hierarchy is therefore:

1. **CALIBRATION FAILURE:** an internal action/control gate fails, or the fine
   carrier does not improve over the base carrier after the blind commit.
   Stop this carrier's dynamics route and diagnose before any Hessian.
2. **CALIBRATION PASS:** all action gates pass and the fine coefficient moves
   toward the continuum coefficient.  Proceed to the separate local-dust and
   inhomogeneous-Hessian mission; do not call the calibration new physics.

No particle, Planck-scale, limiting-speed or continuum-dispersion target is
part of this mission.
