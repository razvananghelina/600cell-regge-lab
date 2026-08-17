# Preregistration: full anisotropic canonical rank of one 600-cell dust slab

Date: 2026-08-17

Prior-art gate commit: `f266cc2`

Status: frozen before evaluating any full-space Hessian entry, any nontrivial
`2T` block, or any new singular value.

## 1. Frozen inputs and scientific question

Require the exact accepted first-tick artifact

```text
reproducible/gravity_600cell_dust_homothetic_canonical_lapse.json
SHA-256 4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9
```

and the exact reduced dynamic-tangent control

```text
reproducible/gravity_600cell_dust_dynamic_tangent.json
SHA-256 1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5.
```

Use both already derived five-stage carriers, `even` and `odd`, with no
schedule search.  The new question is:

> At the accepted non-static homogeneous dust slab, is the complete
> pre-Legendre map locally invertible when all 720 old, 840 internal and 720
> new individual squared-edge magnitudes are allowed to vary?

This is a rank census.  It does not form a propagator, compare a continuum
spectrum or parse a desired speed.

## 2. Complete carrier, variables and matter hypothesis

For each schedule require

```text
vertices in the two layers                  240
Lorentzian four-simplices                  2400
old boundary edges                          720
internal diagonal edges                     720
internal pole edges                         120
new boundary edges                          720
all edge variables                         2280
triangular hinges                          6240.
```

Let

```text
z=(log q_old[720], log x_internal[840], log q_new[720]).
```

All entries are logarithms of positive magnitudes.  The actual squared length
of a pole is `-rho`; all other squared lengths are positive.  At the accepted
background read `(s,r)` from the frozen first-tick artifact and set

```text
q_old = L0^2,
q_new = exp(2s) L0^2,
x_diagonal = exp(s) L0^2-rho,
rho = rho0 exp(r).
```

Keep the same total dust mass and distribute it equally on the 120 pole
world-lines:

```text
S_dust=-(8*pi*M/120) sum_p sqrt(rho_p).
```

No matter perturbation or density transfer between vertices is introduced.
This is the unique uniform localization under the full vertex transitivity,
but it is still a stated discretization hypothesis rather than a derivation
of general inhomogeneous dust.

## 3. Independent local Hessian assembly

For every triangular hinge `h`, use the already certified branch and write

```text
epsilon_h = c_h + sum_(sigma contains h) theta_(sigma,h),
c_h = pi on a boundary triangle and 2*pi otherwise,
S_grav = -i sum_h A_h epsilon_h.
```

The signs and complex-angle branch are inherited unchanged from the complete
one-slab action.  In logarithmic positive-magnitude coordinates, the
Schlaefli identity gives

```text
g_e = -i sum_h epsilon_h A_(h,e) + g_dust,e,

K_(e,f) = -i sum_h [epsilon_h A_(h,ef)
                    + A_(h,e) sum_(sigma contains h) theta_(sigma,h,f)]
            + K_dust,(e,f).
```

Here a derivative vanishes when its edge is absent from the local hinge or
simplex.  The triangle-area first and second derivatives are evaluated
analytically from Heron's polynomial after the signed squared lengths are
composed with `x=sign*magnitude*exp(z)`.  The dust Hessian is diagonal on the
120 poles and is differentiated analytically.

Only first derivatives of four-simplex angles are needed.  Evaluate them at
100 decimal digits by centered differences in each of the ten local
logarithmic edge magnitudes, using exactly

```text
operational primary   1e-20
operational shadow    1e-15
validation primary    3e-20
validation shadow     3e-15.
```

Local derivatives are computed once per geometric four-simplex orbit and
then assembled over all 2400 individual simplices.  No global finite
difference of 2280 coordinates is used as the primary Hessian.

Build four complete sparse Hessians, one for every frozen derivative step.
Require the operational-primary/validation-primary difference entrywise to
be at most ten times the sum of the corresponding primary-shadow proxies plus
an arithmetic floor `1e-70` before conversion to binary64.

## 4. Mandatory action and assembly controls

All of these gates precede a rank interpretation.

1. Every base and displaced representative simplex has one timelike Gram
   direction, positive nonzero leading-minor magnitudes and angle-argument
   modulus above `1e-6`; maximum imaginary contamination in the action
   gradient/Hessian is below `1e-60` before binary64 conversion.
2. The full base gradient is constant on every free 24-edge orbit.  Every
   individual internal equation is below `1e-25`, and the old/new momenta
   reproduce the accepted artifact below `1e-40` per edge.
3. The assembled Hessian is reciprocal.  Its spectral antisymmetry must be
   below ten times the combined derivative proxy plus the measured sparse
   assembly roundoff.
4. Contract the full Hessian with the normalized constant vector on each
   orbit.  The resulting `95 x 95` trivial-sector Hessian must reproduce the
   stored canonical singular values to relative `2e-8` and must reconstruct
   the committed `60 x 60` tangent matrix to relative Frobenius `2e-8`.
   This is a control against the independently implemented, globally
   differentiated reduced calculation.
5. Direct high-precision centered differences of the complete full gradient
   must reproduce `K w` on exactly these target-independent directions:

   ```text
   the normalized global old-boundary scale,
   the normalized global internal-lapse direction,
   the normalized global new-boundary scale,
   one lexicographic old-edge orbit contrast,
   one lexicographic internal-edge orbit contrast,
   one lexicographic new-edge orbit contrast.
   ```

   Use steps `1e-8` and `5e-9`, Richardson extrapolation, and require relative
   agreement below `2e-6`.  These six controls are not selected after seeing
   a singular vector.

Failure of any item gives `FULL_CANONICAL_ASSEMBLY_CONTROL_FAILED`; no rank is
reported.

## 5. Geometry-derived `2T` block decomposition

Enumerate the schedule stabilizer before loading a Hessian.  Require its
order spectrum

```text
1^1, 2^1, 3^8, 4^6, 6^8
```

and seven conjugacy-class sizes `1,1,4,4,4,4,6`.  Require each of the 95 edge
orbits to be free.  Index every orbit by the unique group element taking its
lexicographically first edge to the requested edge.

Construct the 24-dimensional left regular matrices.  Sort conjugacy classes
by `(size, element order, member tuple)` and form the fixed central normal
matrix

```text
Z = sum_(class index r=0..6) 2^r sum_(g in class r) L_g.
```

Its seven eigenspace dimensions must be exactly

```text
1,1,1,4,4,4,9,
```

corresponding to irreducible dimensions `1,1,1,2,2,2,3`.  Within each
isotypic component, enumerate noncentral group elements in the same fixed
order and choose the first Hermitian matrix

```text
Y_g = i (L_g-L_g^*)
```

whose restriction has `d` distinct eigenvalues, each repeated `d` times.
The anti-Hermitian part is load-bearing: for the defining `SU(2)` doublet,
`L_g+L_g^*` can be scalar and would not split the representation coordinate.
The eigenspace of the smallest eigenvalue supplies a deterministic
`d`-column subspace `W_d`.

Lift `W_d` identically across the 65 row and column orbit types.  The seven
minimal canonical blocks therefore have sizes

```text
65,65,65,130,130,130,195.
```

Require their dimensions to sum with representation multiplicity to 1560,
require off-block leakage below ten times the derivative/basis proxy, and
require every singular value of a size-`65d` block to occur `d` times in the
corresponding full isotypic restriction.  This prevents favorable deletion
of a copy or sector.

The invariant block is identified mechanically by overlap above
`1-1e-12` with the constant regular vector.  It must be the block that passes
the reduced-control comparison; no block is chosen by its singular values.

## 6. Calibrated singular-value rule

For each of the seven blocks construct

```text
Jop, Jop_shadow, Jval, Jval_shadow
```

from the four full Hessians.  For the operational SVD triplet `(u_k,s_k,v_k)`
define

```text
epsilon_global = ||Jop-Jop_shadow||_2
               + ||Jval-Jval_shadow||_2
               + ||Jop-Jval||_2
               + epsilon_basis
               + epsilon_svd,

epsilon_svd = max_k (||Jop v_k-s_k u_k||_2
                    + ||Jop^* u_k-s_k v_k||_2).
```

`epsilon_basis` is the measured maximum intertwining/off-block leakage, not a
fixed machine-epsilon guess.  Repeat the SVD with both LAPACK drivers
`gesdd` and `gesvd`; their maximum singular-value difference is added to
`epsilon_svd`.

Classify every singular direction by the global rule:

- resolved nonzero if `s_k > 100 epsilon_global`;
- error-consistent zero if all four estimates are below
  `10 epsilon_global`;
- numerically open otherwise.

Also report the directional matrix-difference proxy along every one of the
ten weakest triplets, but it cannot override the global rule.  Print every
singular value in every block, the block condition numbers, the full rank
obtained with irrep multiplicities, and ranks at relative thresholds
`1e-7,1e-9,1e-11,1e-13,1e-15` only as diagnostics.

If an error-consistent nullspace exists, report for each block the singular
values of its projection onto the last `30d` new-boundary columns.  A
nonzero projected direction is resolved only above the same calibrated
subspace perturbation bound used for the parent block.  No null is called
gauge without a separately verified vertex-displacement generator.

## 7. Frozen predictions and outcomes

The prior-art-based **STRUCTURAL PREDICTION** is that curvature plus dust lifts
the flat-background vertex-displacement constraints and all seven blocks are
full rank, with at least one small but resolved pseudo-constraint sector.
This prediction is not an acceptance condition.

Assign one outcome in this order:

1. `FULL_CANONICAL_ASSEMBLY_CONTROL_FAILED` if any carrier, branch,
   derivative, reciprocity, orbit restriction or symmetry-block control
   fails;
2. `FULL_CANONICAL_RANK_NUMERICALLY_OPEN` if any direction lies between the
   frozen zero and nonzero bands;
3. `FULL_CANONICAL_LEGENDRE_REGULAR` if all seven blocks and hence all 1560
   directions are resolved nonzero;
4. `FULL_CANONICAL_DEGENERACY` if at least one direction is
   error-consistent zero and none is numerically open.

A resolved value with `s_k/s_max<1e-6` is additionally labelled
`PSEUDOCONSTRAINT_CANDIDATE`.  This is a scale diagnostic, not a gauge claim
and not part of the rank decision.

## 8. Acceptance boundary for the next physical calculation

Only `FULL_CANONICAL_LEGENDRE_REGULAR` licenses construction of the complete
`1440 x 1440` canonical tangent map without an added gauge choice.  A regular
map still does not prove 720 physical configuration degrees of freedom; the
small sectors must be compared with geometric vertex displacements and with
refinement.

`FULL_CANONICAL_DEGENERACY` requires an explicit null-generator audit before
any inverse is taken.  `FULL_CANONICAL_RANK_NUMERICALLY_OPEN` requires a
precision correction, not a looser threshold.

No outcome in this mission establishes gravitational waves, Einstein's two
local polarizations, a continuum dispersion relation, a limiting speed, an
absolute tick or a Planck scale.
