# Preregistration: canonical Legendre rank of the 600-cell dust slab

Date: 2026-08-16

Prior-art gate commit: `31b1690`

Corrected two-slab control: `a766740`

Status: **frozen before evaluating any new pre-momentum derivative, complete
`95 x 95` action Hessian or `65 x 65` canonical-inversion spectrum**.

This protocol performs a local rank census at the published stationary slab.
It performs no root search and does not produce a new frame.

## 1. Frozen carrier, action and coordinates

Use only the two already derived ordered-schedule parity representatives

```text
even,
odd.
```

For each, reconstruct the certified complete Lorentzian Regge curvature plus
De Felice--Fabri dust action.  Reconstruct the published constants from their
displayed formulas, with no optimized parameter:

```text
q_old = q_new = l0^2 on all 30 boundary edge orbits,
x_1,...,x_30 = l0^2-tau^2,
x_31,...,x_35 = tau^2,
tau = 0.0102,
M = (90/pi)(2*pi-5*acos(1/3))*l0.
```

Use the 95 dimensionless logarithmic coordinates, in this exact order,

```text
z = (log q_old[30], log x[35], log q_new[30]).
```

For every coordinate, define the per-edge logarithmic action gradient

```text
g_i = (1/24) partial S / partial z_i.
```

The factor `1/24` is common because every retained edge orbit has size 24.
The temporal pole coordinates are positive magnitudes; their actual squared
edge lengths are negative and the existing `-1` edge Jacobian is retained.

## 2. Independent 100-decimal gradient evaluator

At 100 decimal digits, reconstruct each simplex Gram matrix, complex angle,
triangle curvature, area and area derivative.  Use the Schlaefli identity to
evaluate all 95 components of `g` analytically from the area derivatives,
including the exact dust derivative on the five poles.

This evaluator must first reproduce:

1. the published base action within the inherited `5e-8` relative action
   class and all 35 internal stationary equations within the published
   per-edge absolute `1e-7` gate;
2. every already certified old and final per-edge momentum from the corrected
   two-slab artifact within absolute `1e-20` componentwise;
3. the one-timelike-direction Gram signature by Jacobi sign changes of
   `(1,Delta_1,...,Delta_4)` for every evaluated simplex;
4. minimum complex angle-argument modulus above `1e-6`;
5. maximum imaginary contamination below `1e-70`.

The analytic-gradient code must not import a stored Jacobian as its output.
Stored corrected internal blocks are used only for the later independent
cross-check.

## 3. Frozen derivative calibration

Perturb every one of the 95 logarithmic coordinates independently.  For a
step `h`, compute

```text
K(h)[:,j] = [g(z+h e_j)-g(z-h e_j)]/(2h).
```

Use exactly the already calibrated disjoint pairs

```text
operational: (1e-20, 1e-15),
validation : (3e-20, 3e-15).
```

In each pair the smaller step is primary and primary minus shadow is the
signed stability proxy.  Thus there are exactly

```text
95 coordinates * 4 step magnitudes * 2 signs = 760
```

gradient evaluations per parity, plus the base evaluation.  No step may be
changed after seeing a singular value.

Let

```text
K_op, K_op_shadow, K_val, K_val_shadow
```

denote the four matrices.  Require every entry of `K_op-K_val` to lie within
ten times the sum of the corresponding operational and validation proxy
magnitudes plus the arithmetic floor `1e-70`.

The complete action Hessian is `K=K_op`.  Check its reciprocity independently:

```text
||K-K^T||_2
```

must lie below ten times the combined calibrated spectral error.  Report the
Frobenius and spectral antisymmetry, not only a pass flag.

## 4. Canonical-inversion matrix

Use index sets

```text
O = old boundary rows/columns, dimension 30,
X = internal rows/columns,     dimension 35,
N = new boundary rows/columns, dimension 30.
```

At fixed old geometry and pre-momentum, the unknown is

```text
y = (log x[35], log q_new[30]).
```

With

```text
e       = g_X,
p_pre   = -g_O,
F       = (e, p_pre-p_target),
```

extract the frozen `65 x 65` matrix

```text
J_can = [[ K_XX,  K_XN],
         [-K_OX, -K_ON]].
```

This sign and block order are fixed before evaluation.  Also construct the
same matrix from `K_val`; operational/validation agreement is a required
control.

Independently compare `K_XX` and `K_XN` against the already committed
precision-corrected internal Hessian and final-boundary block.  Their relative
Frobenius discrepancies must each be below `1e-6`.  This is only an upstream
coordinate/normalization control; the new `K_OX` and `K_ON` blocks decide the
new result.

## 5. Arbitrary-precision singular spectrum and error rule

Compute the full singular system of `J_can` at 100 decimal digits.  Also
compute a binary64 SVD as an independent leading-spectrum audit, but never use
binary64 alone to classify the weakest modes.

Define

```text
D_op    = J_op - J_op_shadow,
D_val   = J_val - J_val_shadow,
D_cross = J_op - J_val,

epsilon_global = ||D_op||_2 + ||D_val||_2 + ||D_cross||_2 + 1e-70.
```

For every singular triplet `(u_k,s_k,v_k)`, also report

```text
epsilon_k = |u_k^T D_op v_k|
          + |u_k^T D_val v_k|
          + |u_k^T D_cross v_k|
          + 1e-70.
```

The global rule controls the verdict; directional ratios are diagnostics:

- `s_k > 100 epsilon_global`: resolved nonzero;
- `s_k < 10 epsilon_global`: error-consistent zero;
- otherwise: numerically open.

Print all 65 values, all `s_k/epsilon_global`, all `s_k/epsilon_k`, ranks at
relative thresholds `1e-7,1e-9,1e-11,1e-13,1e-15`, the condition number when
defined and the operational/validation principal angles of the
five-dimensional weakest singular subspaces.  A threshold rank is never the
final verdict.

## 6. Frozen scale/lapse/shape decomposition

Use the analytic collective internal tangent from the exact lapse path.  In
the present logarithmic `x` coordinates its unnormalized components are

```text
w_x = (-rho/(l0^2-rho) repeated 30 times,
        1 repeated 5 times),
rho = tau^2.
```

Normalize it and build its complement by the deterministic Householder rule
already certified upstream.  On each boundary use

```text
s = (1,...,1)/sqrt(30)
```

as homogeneous scale and the same deterministic Householder complement as
the 29-dimensional shape basis.

Transform rows and columns of `J_can` into the ordered sectors

```text
rows:    internal lapse [1], internal transverse [34],
         pre-scale [1], pre-shape [29];
columns: internal lapse [1], internal transverse [34],
         final-scale [1], final-shape [29].
```

Report the `4 x 4` table of Frobenius block norms.  For every unresolved or
error-consistent-zero singular vector, and for the five weakest resolved
vectors, report squared overlaps with all four row and column sectors.

For an error-consistent nullspace, report the rank and singular values of its
projection onto the 30 final-boundary columns.  This distinguishes an
internal-only history/gauge degeneracy from non-unique geometry evolution.
No favorable basis rotation inside a degenerate singular cluster is allowed.

## 7. Preregistered prediction and outcome hierarchy

The **STRUCTURAL PREDICTION**, recorded before the new bottom blocks are
evaluated, is

```text
rank(J_can) = 65,
```

with the weakest singular pair dominated by the collective-lapse /
homogeneous-scale sectors.  The prediction is not an acceptance condition.

Assign outcomes mechanically, in this order:

1. `CANONICAL_GRADIENT_CONTROL_FAILED` if base gradients, old/new momenta,
   reciprocity, upstream blocks, branches or derivative calibration fail;
2. `CANONICAL_RANK_NUMERICALLY_OPEN` if any singular direction falls between
   the frozen zero and nonzero bands;
3. `CANONICAL_LEGENDRE_REGULAR` if all 65 are resolved nonzero;
4. `ONE_CANONICAL_LAPSE_NULL` only if exactly one is error-consistent zero,
   its final-boundary projection is error-consistent zero and its internal
   overlap with analytic `w_x` exceeds `0.999999`;
5. `ADDITIONAL_CANONICAL_DEGENERACY` for any other resolved nullity.

Small but resolved values are additionally labelled
`RESOLVED_PSEUDOCONSTRAINTS`; operationally this diagnostic label means
`s_k/s_max < 1e-6` while the frozen global error rule still resolves the mode
nonzero.  This diagnostic cutoff does not affect the rank verdict and no such
value is rounded to gauge.

## 8. Degrees-of-freedom and claim boundary

Before evaluation, the expected reduced local count is:

```text
1 homogeneous scale direction + 29 schedule-invariant shape directions.
```

If `J_can` is regular, all thirty have a locally unique canonical image for
admissible `(q_old,p_pre)` data.  If a nullspace projects with rank `r>0` onto
the final boundary, at least `r` reduced configuration combinations are not
uniquely predicted at linear order.  Internal-only nulls are not subtracted
from boundary degrees of freedom until their gauge generator is established.

This is not a count of continuum gravitons.  It excludes boundary modes
outside the order-24 quotient and does not establish constraints on the full
720-edge phase space.

## 9. Acceptance, kill and next step

**Acceptance boundary for continuing canonical inversion:** both parities
pass every implementation gate and return either
`CANONICAL_LEGENDRE_REGULAR` or `ONE_CANONICAL_LAPSE_NULL`, with the latter
quotiented by its preregistered analytic generator.

**Kill boundary for the immediate solve:** any implementation-clean
`CANONICAL_RANK_NUMERICALLY_OPEN` or `ADDITIONAL_CANONICAL_DEGENERACY` blocks a
nonlinear next-frame solve until the rank is resolved or the additional
canonical data are supplied.

All five scientific outcomes are valid verifier outputs.  The verifier exits
nonzero only for a failed mechanical certificate; an implementation-clean
OPEN or kill-boundary result exits zero while recording that nonlinear
continuation is not accepted.

Even acceptance derives no new frame.  The next separate protocol must first
use `p_target=p_pre(published)` as a reproduction control.  Only after exact
reproduction may it use `p_target=p_post(published)` as a candidate forward
tick.  Expansion, contraction, inflation, a clock, `c`, Planck scales and
particle masses are outside this protocol.
