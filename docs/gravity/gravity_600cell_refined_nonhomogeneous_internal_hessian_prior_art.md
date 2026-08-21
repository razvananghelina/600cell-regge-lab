# Prior-art gate: complete nonhomogeneous internal Hessian of the refined slab

Date: 2026-08-21

Status: completed before assembling or diagonalizing any full local Hessian.

## 1. Exact object and complete hypotheses

Use the already certified projected barycentric carrier

```text
K0=P(sd K_600),  f=(2640,17040,28800,14400),
```

its rank-derived chordal metric, the supplied proper duration
`tau0=0.0102`, all 24 colour-ordered staircase triangulations of `K0 x I`,
the corrected complex Lorentzian Regge action with boundary terms, and the
accepted curvature-matched conserved vertex masses

```text
m_v=K_v/(8*pi).
```

The masses are frozen on the vertical worldlines.  No independent dust-density
perturbation is introduced.  This is a load-bearing **STRUCTURAL** matter
hypothesis, not a derivation that dust has no perturbations.

For each schedule use every actual edge of its slab, in logarithmic absolute
signed-squared-length coordinates.  The common labelled boundary contains

```text
17040 old + 17040 new = 34080
```

variables.  The schedule-specific interior contains

```text
17040 cross diagonals + 2640 vertical edges = 19680
```

variables.  Hence the complete local action Hessian has size
`53760 x 53760`; the present gate constructs only the complete internal block

```text
C_s = d^2 S_s / (d i_s d i_s),  C_s in R^(19680 x 19680),
```

with both spatial boundaries fixed.  It does not restrict perturbations to
rank orbits, a selected Laplacian band, a subgroup sector or a proposed
continuum polarization.

The product geometry and curvature-matched masses form a stationary interior
for every schedule.  Varying the common product duration while preserving the
induced cross-edge lengths gives an analytically determined internal tangent
`n_s`.  The target-free question is

```text
ker(C_s) = span(n_s)  for every one of the 24 schedules?
```

No gauge quotient is assumed in advance.  The question is asked before any
boundary Schur reduction, dispersion comparison, wave-speed extraction or
continuum interpretation.

## 2. KNOWN from primary literature

- Linearized four-dimensional Regge calculus on a flat background has exact
  vertex-displacement symmetry and gauge-invariant curvature degrees of
  freedom ("lattice gravitons").  Their count and propagation depend on the
  Pachner moves used in the evolution: P. A. Hoehn,
  [Canonical linearized Regge Calculus: counting lattice gravitons with
  Pachner moves](https://arxiv.org/abs/1411.5672), especially the abstract and
  the canonical construction.
- On curved Regge solutions, exact discrete diffeomorphism symmetry is
  generically broken and canonical constraints become background-dependent
  pseudo-constraints: B. Bahr and B. Dittrich,
  [(Broken) Gauge Symmetries and Constraints in Regge
  Calculus](https://arxiv.org/abs/0905.1670).
- Linearized Regge actions have special factorisations under Pachner moves,
  but four-dimensional triangulation independence is not generic: B. Dittrich
  and S. Steinhaus,
  [Path integral measure and triangulation independence in discrete
  gravity](https://arxiv.org/abs/1110.6866).
- The action of a simplicial slab generates its canonical boundary relation
  only after the bulk equations are imposed, and singular discrete systems
  require explicit pre/post constraints: B. Dittrich and P. A. Hoehn,
  [Canonical simplicial gravity](https://arxiv.org/abs/1108.1974).
- A Sorkin evolution of the unrefined 600-cell with more than one homogeneous
  variable is known, and it encounters a causal stopping point: A. De Felice
  and E. Fabri,
  [Singularities of the closed RW metric in Regge Calculus: a generalized
  evolution of the 600-cell](https://arxiv.org/abs/gr-qc/0106077).

These sources establish the interpretation and the warnings.  They do not
determine the kernel of the present `19680 x 19680` block on
`P(sd K_600)`, with the selected rank masses and all 24 bare staircases.

## 3. Repository controls already established

- **DERIVED COMPUTATIONAL / STRUCTURAL:** every individual internal equation
  vanishes at the curvature-matched product seed for all 24 schedules.
- **DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED:** on the complete
  `H4`-invariant ten-dimensional internal sector, the Hessian has inertia
  `(9,1,0)` and its unique null line is the analytic product-duration tangent.
- **DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED:** that invariant null
  line couples to the boundary and imposes one compatibility condition; it is
  not a null direction of the complete invariant Hessian.
- **DERIVED EXACT:** the 24 schedules have distinct cross-diagonal edge sets;
  none may be selected after observing a desirable spectrum.

The invariant result does not bound the kernel in nontrivial `H4` sectors.
Inferring the complete kernel from ten orbit sums would be a forbidden
minisuperspace-to-full-system extrapolation.

## 4. Proposed complete test

For each schedule:

1. rebuild all actual slab edges, triangles and pentachora;
2. assemble `C_s` sparsely from local area Hessians and local dihedral-angle
   derivatives, using the Regge--Schlaefli first-variation identity;
3. verify every individual internal gradient, symmetry of `C_s`, Lorentzian
   branch continuity and the analytic residual `C_s n_s`;
4. pull `C_s` back to the ten orbit-constant directions and reproduce the
   accepted invariant block;
5. test nonsingularity of the bordered matrix

   ```text
   K_s = [[C_s,n_s],[n_s^T,0]].
   ```

   Since `C_s` is symmetric and `C_s n_s=0`, `K_s` is nonsingular if and only
   if `ker(C_s)=span(n_s)`;
6. compare every schedule with its exact time reverse and report the complete
   near-zero spectral census before attaching any physical label.

Two independently chosen local differentiation scales and an entrywise
assembly-roundoff envelope must bound the operator error.  A smallest
bordered singular/eigenvalue counts as resolved only if its certified lower
bound remains positive.  A deliberate local-stencil corruption must fail the
orbit pullback or null test.

## 5. Outcome hierarchy fixed before calculation

- If stationarity, the invariant pullback or the analytic null line fails,
  the proposed local extension is inconsistent and the refined propagation
  route stops at this gate.
- If any schedule has additional error-compatible null directions, the
  ordinary internal elimination remains forbidden.  Their complete symmetry
  content and boundary coupling become the next problem; no wave claim is
  allowed.
- If every bordered matrix is certified nonsingular, the complete internal
  kernel is exactly the product-duration line.  This licenses, but does not
  establish, the later constrained nonhomogeneous boundary response.
- Schedule-dependent nonzero spectra are reported as **STRUCTURAL**.  They do
  not by themselves prove schedule dependence of the eliminated boundary
  dynamics.

## 6. Framing attack and scope

This is a nullity/solvability calculation, not yet a graviton calculation.
On a curved discrete background, absence of additional zero modes may mean
that discretization has lifted continuum diffeomorphism symmetry; it is not
evidence for the two polarizations of general relativity.  Conversely,
additional null modes are not automatically gauge modes.

Eigenvalue magnitudes depend on the declared logarithmic edge coordinates and
their Euclidean norm.  Kernel dimension, the bordered nonsingularity test and
the existence of a stationary elimination are the coordinate-invariant
content.  No comparison with a spatial Laplacian, continuum harmonic,
desired wave equation, `c`, `G`, Planck scale or particle spectrum is permitted
in this gate.

The search found the general canonical and 600-cell constructions above but
not this exact full refined kernel census.  Search absence is not proof;
external novelty remains **OPEN**.
