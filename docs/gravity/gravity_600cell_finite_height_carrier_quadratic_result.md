# Finite-height carrier quadratic canonicity result

Date: 2026-08-22.

## Verdict

**DERIVED COMPUTATIONAL, MECHANICALLY DIFFERENTLY REPLICATED:** under the
complete hypotheses below, the two staircase schedule parities induce no
resolved difference in the one-sided quadratic form

```text
Q_p = G_p^T H_p G_p
```

on the exact rank-240 scale-plus-strut tangent carrier at the first physical
finite-height slab from `v=3/2`.

At 180 decimal digits, the independently reconstructed relative parity
differences are

```text
d_R01 = 1.0946049420189856e-101
d_R12 = 6.8412808876186602e-103
d_R23 = 4.2758005547616626e-104

e_step  = 2.8057204102663833e-101
10 e_total = 2.8057204102663833e-100.
```

Thus every observed difference lies inside the preregistered, target-free
error gate.  This is a positive necessary canonicity result, not a derived
nonhomogeneous evolution law.

**STRUCTURAL:** the 240-column carrier is an exact infinitesimal compatibility
space selected by the 600-cell geometry.  Its nonlinear integrability is not
proved.  The global H4 census has two schedule orbits of size 60, represented
by the even and odd constructions; using their orbit representatives relies
on the already certified H4 covariance.

**OPEN:** nonlinear carrier integrability, a reduced boundary evolution map,
gauge reduction, local propagating degrees of freedom, stability, a continuum
limit and external novelty.

**NOT DERIVED:** a graviton, a dispersion relation, a physical tick, a
limiting speed, `c`, `G`, Planck units or particle masses.

## Complete hypotheses and object

The result fixes all of the following:

1. the homogeneous tetrahedral-frustum 600-cell action with zero cosmological
   constant and conserved global dust;
2. the first positive-height slab reconstructed from the fixed incoming state
   `v=3/2`, with `L_minus=1`;
3. the already certified Lorentzian angle branch and boundary-term convention;
4. the complete `2280` logarithmic signed-squared-edge carrier, split into
   `720` fixed old-boundary, `840` internal and `720` new-boundary variables;
5. the two certified H4 staircase parity representatives, covering schedule
   orbits of size `60+60`;
6. the active `1560 x 1560` internal-plus-new Hessian block `H_p`;
7. the geometry-selected sparse map `G_p : R^240 -> R^1560`, with columns
   `120` upper scale variations and `120` strut variations;
8. no quadratic eigenspace, continuum mode, desired speed or fitted
   coefficient used to define the carrier or the comparison.

The finite-height background was independently reconstructed as

```text
q      = 9.6180026533418980973894525206524635528838...
h      = 0.2040549716108237129281133802550001672207...
lambda = 2.9626012583805081610145042646754878733901...
rho    = 0.04163843143909407142842261006057832848428...
M      = 11.7524154611017941257872622934980734615402...
```

The adversarial reconstruction obtained an exact-equation residual of
`8.70e-80` at the precision carried through the frozen primary artifact.

## Why the Hessian comparison is admissible here

A Hessian away from a critical point is not a tensor.  If a nonlinear carrier
embedding has jets `(G,K)`, then

```text
d^2(S o Phi) = G^T H G + grad(S) dot K.
```

The primary verifier therefore established before interpreting `Q_p` that

```text
maximum artificial/internal gradient       8.46e-114,
maximum old-boundary parity difference      7.75e-121,
maximum new-boundary parity difference      4.65e-120.
```

Schedule-specific second jets can occur only on artificial internal
diagonals, where the gradient vanishes, while the physical boundary
embedding is common to the two schedules.  The schedule-dependent part of
`grad(S) dot K` therefore cancels in this comparison.  This licenses a
necessary quadratic parity test; it does not prove the existence of `Phi`
beyond first order.

## Primary calculation

The registered primary verifier assembled four complete binary64
`2280 x 2280` Hessians for each parity and pulled their active blocks through
dense carrier matrices.  It returned

```text
FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_PRIMARY
22/22 PASS
```

with

```text
normalized parity difference       2.2169149789731003e-14
independent gate                    3.2428479263755228e-7
carrier-corruption effect           2.8390557866167732e-4
synthetic rank-one effect           2.8216370331471513e-1.
```

The large primary error envelope was deliberately conservative.  This result
remained primary-only until the adversarial route below completed.

## Mechanically different adversarial reconstruction

The adversarial verifier never materialized a complete ambient Hessian and
never called the primary binary64 assembler.  At 180 decimal digits it:

1. assembled only the identity-row representative `2T` Hessian kernel
   `(row orbit,column orbit,relative group element)`;
2. reconstructed entries with the independently checked group law
   `H[(r,a),(c,b)]=K[r,c,a^-1*b]`;
3. pulled each kernel entry directly into a `240 x 240` form through sparse
   carrier dictionaries;
4. used four logarithmic derivative steps and three Richardson levels;
5. loaded the primary matrices only after fixing the high-precision parity
   classification in memory;
6. independently differentiated the complete scalar action in two frozen
   data directions for each parity.

It returned

```text
FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_ADVERSARIALLY_REPLICATED
18/18 PASS.
```

The high-precision and primary `R12` matrices agree relatively at

```text
even  3.3576639237426950e-14
odd   3.3175764983958699e-14.
```

The four direct complete-action derivatives reproduce their kernel
quadratic values with relative complex residuals between `1.03e-88` and
`3.94e-85`.  The individual scalar-action evaluations remain on the same
Lorentzian branch and below the frozen action-reality gate.

## Hostile controls and convention sensitivity

The equality is not caused by a verifier that is insensitive to its inputs:

```text
+1/10 carrier-coefficient corruption effect   4.9925094698501382e-5
synthetic rank-one Hessian effect              3.6162624011862327e-2
reversed group-product diagnostic              6.4434995338181637e-1.
```

The last value is diagnostic rather than an acceptance gate, but it confirms
that the `a^-1*b` convention matters strongly on this carrier.

## Preserved first failure

Implementation commit `dcb0d1c` first returned `CONTROL_FAILED` because the
protocol mistakenly applied a reality threshold to raw complex Lorentzian
boost angles instead of to the physical action-derived kernel.  It then
failed JSON serialization on a nested `mpmath.mpf`.

The failure, its byte-identical pre-exception matrix and the pre-rerun protocol
correction are preserved in
[the first-failure note](gravity_600cell_finite_height_carrier_quadratic_adversarial_first_failure.md).
No derivative step, parity threshold or decisive outcome gate changed.  The
accepted matrix is byte-identical to the first-failure matrix:

```text
sha256 8a3ea0c3b8ee720d8ffdf07e7486aefdd0247ca1cfdbeb99f443091376f31729.
```

## Provenance ledger

```text
5fab3f5  prior-art and framing gate
01382fb  primary preregistration
a6c93d4  primary registry entry before implementation
a020978  schedule-orbit provenance amendment
2f9361b  complete input-provenance amendment
810d4bd  primary verifier implementation
d2796de  primary result
0235162  adversarial preregistration
6fde152  adversarial registry entry before implementation
dcb0d1c  first adversarial implementation and failed execution
e0d06b3  preserved failure and corrected physical reality gate
0dcacb4  corrected implementation before rerun
dd57c72  accepted adversarial artifacts
```

The ordering proves that the objects, steps, directions and decisive
thresholds preceded their corresponding evaluations.  It does not prove
external novelty or turn this necessary gate into sufficient physics.

## What changed, and what did not

**DERIVED COMPUTATIONAL:** the complete set of 120 certified staircase orders
does not introduce a resolved ambiguity in this one-sided quadratic tangent
form at the frozen finite-height background.  The earlier fear that choosing
even versus odd staircase parity would already choose different linearized
local physics on this carrier is not realized at this order.

**STRUCTURAL:** the result makes the common quadratic form a defensible next
object to study.  It may ultimately be explained by a boundary-preserving
combinatorial equivalence rather than by a new physical principle; either
way, it is a canonicity check rather than a prediction.

**POST-RESULT FRAMING CORRECTION:** bare local integrability of a linear map
`G` is automatic: `Phi(x)=y0+Gx` is a local embedding into the open
nondegenerate length domain.  The meaningful second-order question would be
uniqueness of a geometry-selected normal second jet after exact nonlinear
boundary scaling and face gluing are fixed.  That question is logically
later than action compatibility.

**OPEN NEXT GATE:** before inspecting any continuum dispersion target,
compute the `840 x 240` internal-equation derivative

```text
R_p=H_p[internal,active] G_p
```

for both parities, first on the 720 diagonal rows and then with all 120 pole
rows.  Full column rank closes this carrier route at first order.  A common
nonzero kernel selects the only directions worth carrying into nonlinear
integrability, canonical reduction and eventual spatial-mode analysis.

No full-suite run was performed; only the registered primary and adversarial
verifiers relevant to this gate were executed.
