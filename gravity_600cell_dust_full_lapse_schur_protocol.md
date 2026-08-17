# Preregistered protocol: full vertex-lapse Schur sector

Date: 2026-08-17  
Prior-art commit: `58f14e1`.  No Schur value, determinant, singular value,
subspace angle, or post-result physical comparison has been evaluated at the
time of this protocol.

## 1. Frozen inputs and claim boundary

Use exactly:

```text
reproducible/gravity_600cell_dust_homothetic_canonical_lapse.json
SHA256 4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9

reproducible/gravity_600cell_dust_full_anisotropic_legendre_rank.json
SHA256 7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226

reproducible/verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py
SHA256 834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5
```

The imported full-rank verifier must retain `18/18` and the input artifact
must report both parities as `FULL_CANONICAL_LEGENDRE_REGULAR`, rank 1,560,
zero error-consistent nullity, and zero open directions.

The new calculation asks only whether the 120 pole/lapse directions survive a
precision Schur audit and whether a separately defined geometric lapse carrier
matches their canonical lifts.  It does not parse a continuum spectrum,
Einstein polarization count, wave speed, Planck scale, or experimental target.

## 2. Frozen row/column partition

Retain the full canonical ordering

```text
rows    = 840 internal equations, then 720 negative old momenta,
columns = 840 internal variables, then 720 new-boundary variables.
```

The five internal edge-orbit types consisting entirely of the 120 pole edges
define five weak row and column orbit types.  The remaining 60 canonical orbit
types are strong.  No singular vector is used to choose this partition.

For every irreducible `2T` representation dimension
`d=1,1,1,2,2,2,3`, form

```text
J_d = [ A_d  B_d ],
      [ C_d  D_d ],

A_d: 60d x 60d,
S_d = D_d-C_d A_d^(-1) B_d: 5d x 5d.
```

The full pole Schur carrier has dimension

```text
sum_irreps d * (5d) = 5 sum_irreps d^2 = 120.
```

## 3. Independent high-precision assembly

Use 100 decimal digits for Lorentzian geometry and derivatives.  Retain, as
arbitrary-precision complex numbers rather than converting to binary64, the
four already frozen centered-angle derivative estimates:

```text
operational primary   1e-20
operational shadow    1e-15
validation primary    3e-20
validation shadow    3e-15.
```

Assemble only the 65 representative canonical rows against all `65 x 24`
columns.  Complete every `24 x 24` orbit block by the exact equivariance rule

```text
J[(a,r),(b,c)] = K[a,b,r^(-1)c].
```

The local Hessian is the same Schlaefli-reduced expression and analytic dust
term used by the full verifier.  All base and displaced simplices must retain
one timelike Gram direction, positive nonzero leading-minor magnitudes, angle
argument modulus above `1e-6`, and imaginary contamination below `1e-70`.

## 4. High-precision `2T` Fourier basis

Construct the left regular matrices from geometry before loading any Hessian.
Use the same class-sorted central matrix `Z` as the full census, but diagonalize
the Hermitian matrix

```text
C = Re(Z) + sqrt(2) Im(Z)
```

at 100 decimal digits.  It must have seven eigenspaces of dimensions
`1,1,1,4,4,4,9`.  Within a `d^2` component, use the first nonidentity group
element for which

```text
Y_g=i(L_g-L_g^*)
```

has `d` eigenvalues of multiplicity `d`; the smallest eigenspace gives the
same target-independent `d`-column representation coordinate as before.

Require orthonormality, central/splitter residuals, and invariance under all
right regular matrices below `1e-70`.  Require every high-precision minimal
block to reproduce the corresponding stored binary64 singular multiset with
maximum normalized discrepancy below `2e-10`.  This comparison is a control,
not a rank criterion.

## 5. Ball solve and calibrated rank rule

Convert at least 85 significant decimal digits of every high-precision block
entry into `python-flint` complex balls and perform the linear algebra at 80
decimal digits.  For all four derivative estimates:

1. `det(A_d)` must exclude zero and the ball solve `A_d X=B_d` must succeed;
2. form `S_d=D_d-C_d X` without a binary64 solve;
3. record whether `det(S_d)` excludes zero and the maximum output-ball radius.

For both `A_d` and `S_d`, convert ball midpoints only for diagnostic SVDs and
define

```text
epsilon_step = ||Mop-Mop_shadow||_2
             + ||Mval-Mval_shadow||_2
             + ||Mop-Mval||_2,

epsilon_ball = Frobenius upper bound from all entry radii,

epsilon_svd  = maximum two-sided SVD residual
             + gesdd/gesvd singular-value discrepancy,

epsilon_global = epsilon_step + epsilon_ball + epsilon_svd.
```

Classify a singular direction mechanically:

- resolved nonzero if all four determinants/solves relevant to its block are
  certified and `s_op > 100 epsilon_global`;
- error-consistent zero if all four ordered singular estimates are below
  `10 epsilon_global`;
- numerically open otherwise.

Ordering of clustered singular values is descending for all four estimates.
The determinant exclusion is an additional block certificate; it cannot
override an open singular direction.

## 6. Frozen geometric vertex-lapse carrier

For each new vertex `v+120`, define a 1,560-component column variation:

```text
delta log(rho_v) = 1,
delta log(q_(u,v+120)) = -rho/(exp(s)L0^2-rho)
                         for every spacelike cross-edge ending at v+120,
all other components = 0.
```

Require the 120 columns to be independent, `2T`-equivariant, and their sum to
equal the preregistered collective internal-lapse direction after
normalization below `1e-14`.

In every minimal irrep sector compare three `5d`-dimensional right subspaces:

1. the canonical Schur lifts `(-A_d^(-1)B_d, I)`;
2. the frozen geometric vertex-lapse columns;
3. the weakest `5d` right singular vectors of the complete `J_d`.

Use orthonormal bases and report projector spectral distances.  Assign each
pair independently:

```text
IDENTIFIED       distance < 1e-8,
SEPARATED        distance > 1e-4,
NUMERICALLY_OPEN otherwise.
```

These labels do not affect Schur rank.  In particular, a regular Schur
operator can coexist with an identified lapse carrier: that is precisely the
pseudo-constraint possibility on a curved discrete background.

## 7. Mechanical outcomes

Assign in this order:

1. `FULL_LAPSE_SCHUR_ASSEMBLY_CONTROL_FAILED` for any provenance, branch,
   high-precision group-basis, representative-kernel, or stored-spectrum
   control failure;
2. `FULL_LAPSE_SCHUR_STRONG_BLOCK_OPEN` if an `A_d` solve/determinant or its
   calibrated nonzero classification fails;
3. `FULL_LAPSE_SCHUR_NUMERICALLY_OPEN` if any Schur singular direction is
   open;
4. `FULL_LAPSE_SCHUR_REGULAR` if all 120 Schur directions are resolved
   nonzero;
5. `FULL_LAPSE_SCHUR_DEGENERATE` if at least one direction is
   error-consistent zero and none is open.

The prior-art-based **STRUCTURAL PREDICTION** is `REGULAR`: curvature and dust
should lift exact lapse gauge directions into pseudo-constraints.  The
separate subspace prediction is that the canonical lifts will be close to,
but need not equal, the simple vertex-lapse carrier because the staircase
schedule breaks full vertex symmetry.  Neither prediction is an acceptance
condition.

## 8. Acceptance boundary

Only `FULL_LAPSE_SCHUR_REGULAR` with all controls passing removes the present
precision objection to construction of the complete `1440 x 1440` tangent
map.  It does not identify 720 physical configuration degrees of freedom.

A degenerate result requires an explicit quotient by certified geometric
generators.  An open result requires more precision or interval derivatives,
not a looser threshold.  No outcome here establishes gravitons, stability,
Lorentz dispersion, a limiting speed, an absolute tick, or a Planck scale.
