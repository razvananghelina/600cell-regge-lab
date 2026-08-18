# Prior-art gate: analytic weak-lapse dynamics of the cellular 600-cell action

Date: 2026-08-17

Status: written after the action-invariance result `87faad5` and before any
new differentiation, series expansion or evaluation of the closed cellular
action on committed tick data.

## 1. Exact object

Use only the already certified homogeneous cellular action

```text
Delta = L_plus-L_minus,
h = sqrt(rho+Delta^2/4),
c = (Delta^2+2 rho)/(2(Delta^2+3 rho)),
b = Delta/sqrt(8(Delta^2+3 rho)),

S_grav = 360(L_minus+L_plus)h[2 pi-5 acos(c)]
       + 600 sqrt(3)(L_minus^2-L_plus^2)asinh(b),

S_total = S_grav-8 pi M sqrt(rho).
```

Here `L_minus,L_plus,rho>0`, all 600 frusta are on the same certified
Lorentzian branch, and `M` is one conserved total dust mass.  The action
already includes the spatial Regge boundary term.

The inherited static normalization is

```text
epsilon3 = 2 pi-5 acos(1/3),
M = (90/pi) epsilon3 L0.
```

It makes `S_total(L0,L0,rho)=0` for every positive `rho`; therefore a static
slab alone cannot select a lapse.

## 2. New question

The committed weak-lapse calculations found, numerically and before the
closed action was available,

```text
u_n/u_1       -> n,
v_n/v_1       -> 2n-1,
a_n/u_1       -> n(n+1)/2,
r_n/v_1       -> n^2,
p_post,n/k    -> 2n+1
```

for `n<=4` as the inherited lapse is scaled to zero.  The new mission asks:

> Do these integer limits follow analytically from the Euler--Lagrange and
> canonical seam equations of `S_total`, and what discrete acceleration do
> they represent relative to closed Friedmann dust cosmology?

This is not a search for another sequence and not a fit to the existing
integers.  The coefficients must be obtained from a target-independent
series solution of the closed action before the committed numerical ratios
are read for comparison.

## 3. Canonical equations and complete hypotheses

For slab `n`, write

```text
S_n = S_total(L_(n-1),L_n,rho_n;M).
```

The physical internal equation is the strut/lapse equation

```text
E_rho,n = partial S_n/partial log(rho_n) = 0.
```

At an intermediate spatial boundary, additivity gives the seam equation

```text
E_L,n = partial(S_n+S_(n+1))/partial log(L_n^2) = 0,
```

equivalently equality of the post-momentum of slab `n` and the pre-momentum
of slab `n+1` with the already fixed signs.  The initial canonical datum is
the exact momentum of the regular static slab with lapse
`tau=lambda tau0`.

The weak-lapse ansatz must be derived, not assumed after seeing the outputs:

```text
a_n = log(L_n/L0) = A_n lambda^2+O(lambda^4),
r_n = log(rho_n/(lambda^2 rho0))
    = R_n lambda^2+O(lambda^4).
```

Only the branch continuous from the committed contracting root is in scope.
Expansion, time-reversed and other nonlinear roots are separate branches.

## 4. Primary prior art

### 4.1 Closed dust Regge cosmology is established

Collins and Williams, [*Dynamics of the Friedmann Universe Using Regge
Calculus*](https://doi.org/10.1103/PhysRevD.7.965), introduced regular
tetrahedral Cauchy surfaces with 5, 16 or 600 cells as discrete closed
Friedmann models.  Thus interpreting the homogeneous 600-cell variable as a
minisuperspace scale factor is **KNOWN**.

Brewin, [*Friedmann cosmologies via the Regge
calculus*](https://doi.org/10.1088/0264-9381/4/4/023), constructs dust-filled
time-symmetric Regge cosmologies, discusses the full leg/strut/diagonal
equations and warns that imposing symmetry before variation need not retain
the unrestricted local Regge equations.  This directly motivates the
repository's previous checks that all artificial-diagonal equations vanish.

De Felice and Fabri, [*The Friedmann universe of dust by Regge Calculus:
study of its ending point*](https://arxiv.org/abs/gr-qc/0009093) and
[*Singularities of the closed RW metric in Regge
Calculus*](https://arxiv.org/abs/gr-qc/0106077), evolve the 600-cell with a
Sorkin scheme and identify a later causality-breaking endpoint.  Multiple
600-cell dust steps are therefore also **KNOWN**.

### 4.2 Modern Lorentzian shell controls

Dittrich, Gielen and Schander, [*Lorentzian quantum cosmology goes
simplicial*](https://arxiv.org/abs/2109.00875), study closed Lorentzian Regge
shells including dust and compare different discretizations with continuum
cosmology.  Their results make refinement dependence load-bearing: agreement
of one fixed 600-cell with Friedmann is not a continuum theorem.

Tsuda and Fujiwara, [*Oscillating 4-Polytopal Universe in Regge
Calculus*](https://doi.org/10.1093/ptep/ptab074), derive the regular
4-polytopal frustum equations from the cellular action and identify strut
variation with the Hamiltonian constraint and spatial-length variation with
the evolution equation.  The role assigned here to `E_rho` and `E_L` is
therefore standard.

### 4.3 Canonical map and pseudo-constraint

Dittrich and Hoehn, [*Canonical simplicial
gravity*](https://arxiv.org/abs/1108.1974), show that a complete additive
one-step action generates pre/post momenta and a discrete canonical map.
Bahr and Dittrich, [*(Broken) Gauge Symmetries and Constraints in Regge
Calculus*](https://arxiv.org/abs/0905.1670), explain why lapse-like gauge
directions on flat backgrounds can become weak pseudo-constraints on curved
discretizations.

Thus even an analytically nonzero weak eigenvalue would not by itself make
the lapse a fundamental physical clock.

### 4.4 Continuum control

For a closed continuum `S^3` with pressureless dust and zero cosmological
constant, the Friedmann equation gives a time-symmetric maximum radius and a
negative proper-time acceleration.  Expanding a scale factor around that
turning point produces a quadratic position law, linearly changing velocity
and odd half-step momenta.  Therefore the *shape* of the integer recurrence
is expected constant-acceleration kinematics, not prima facie new physics.

The informative comparison is the coefficient and its discretization error,
not the existence of triangular numbers.

## 5. KNOWN / CONTROL / OPEN

### KNOWN

- The homogeneous action is exactly independent of the two staircase
  parities and has the closed cellular form above.
- All artificial diagonal gradients vanish on the homothetic family.
- Four weak-lapse ticks numerically exhibit the quadratic integer law.
- Closed Friedmann dust has a quadratic turning-point expansion.
- Regular 600-cell dust cosmology and its finite-resolution departures are
  established prior art.

### CONTROL

- Differentiate the closed formula independently by symbolic algebra and by
  high-precision finite differences.
- Reproduce the exact static action and static pre/post momentum before any
  nonstatic series is accepted.
- Derive the lowest nonzero weak-lapse system without reading the stored
  integer targets.
- Require the coefficient system to have stated rank and a unique
  contracting solution after the initial canonical datum is fixed.
- Only after committing that coefficient solution compare it with the
  `n<=4` artifacts.
- Re-evaluate the committed finite-lapse states with the cellular equations
  and require their residuals to match the original staircase certificates.
- Compare the derived acceleration coefficient with the continuum closed
  Friedmann value after explicitly stating the radius/edge conversion and
  mass convention.

### OPEN

- An analytic derivation of the integer weak-lapse law.
- The exact leading coefficient and its continuum discrepancy.
- Whether the finite 600-cell lapse is gauge or a resolved
  discretization-dependent pseudo-constraint in the cellular variables.
- Spatial-refinement convergence.
- Anisotropic stability and propagating tensor modes.
- An absolute clock, causal limiting speed, Planck time, Planck mass or
  particle masses.
- External novelty of the exact coefficient audit.

## 6. Framing attack

Three shortcuts are forbidden.

1. Expanding only the already observed integer sequences would be circular.
   The expansion must start from `S_total` and its canonical equations.
2. Equality with continuum constant-acceleration *shape* would carry little
   evidence.  The coefficient and a spatial-refinement law are what could
   support recovery of Einstein dynamics.
3. A selected finite-carrier lapse is not automatically physical time.  If
   its selection vanishes in the weak-lapse or refinement limit, it is gauge
   restoration, not a fundamental tick.

A successful derivation closes an internal analytic gap and makes the four
ticks one theorem-controlled local trajectory.  It does not constitute a new
theory of gravity, because the broad Collins--Williams mechanism is known.
External novelty remains **OPEN** regardless of outcome.
