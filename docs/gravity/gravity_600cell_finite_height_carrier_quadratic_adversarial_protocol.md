# Adversarial preregistration: high-precision orbit-kernel reconstruction

Date: 2026-08-22.

Primary result commit: `d2796de`.

Status: **FROZEN BEFORE THE FIRST ADVERSARIAL ORBIT-KERNEL ASSEMBLY.**

The primary result is known:

```text
FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_PRIMARY,
22/22 PASS.
```

This protocol is allowed to test that claim but may not change the primary
thresholds or relabel the primary computation.  No quadratic eigenmode,
continuum target or desired physical result may be read.

## 1. Frozen inputs

```text
reproducible/gravity_600cell_finite_height_carrier_quadratic.json
  0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9

reproducible/gravity_600cell_finite_height_carrier_quadratic_matrices.npy
  e01bfb28d4c5313b466118315f8ca22c16c2cdc4e94ab05f30c730a136d81cb2

reproducible/verify_gravity_600cell_finite_height_carrier_quadratic.py
  bbe7112270a7f2bcb2d443fab45ca450598e7234250bd335b14a4ed7869443a5

reproducible/verify_gravity_600cell_dust_full_boundary_tangent.py
  c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf

reproducible/gravity_600cell_full_scale_strut_symbolic_gap_resolution.json
  ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179
```

The finite-height background is read from the frozen primary artifact and
then independently checked against the exact formulas; no primary quadratic
matrix is read until the high-precision forms and their parity classification
have been frozen in memory.

## 2. Mechanically different decisive route

The primary calculation assembled four complete dense `2280 x 2280`
Hessians in binary64 and multiplied them by dense carrier matrices.  The
adversarial route is forbidden to call that assembler.

For each schedule parity, use the free regular `2T` ordering of every
24-element edge orbit.  Assemble at 180 decimal digits only the Hessian rows
whose row-group coordinate is the identity.  Store each representative entry
as

```text
(row orbit, column orbit, relative group element) -> value.
```

The expected number of nonzero representative entries is not frozen and is
not an outcome criterion.  Reconstruct active Hessian entries solely through
the independently checked group law

```text
H[(r,a),(c,b)] = K[r,c,a^{-1}b].
```

Equivalently, a kernel key `(r,c,k)` contributes at the 24 positions

```text
((r,a),(c,a*k)),  a in 2T.
```

Pull each such entry directly into the `240 x 240` data form using sparse
row dictionaries for the scale-plus-strut carrier.  Do not materialize a
full ambient Hessian.  All accumulation remains in `mpmath` arbitrary
precision; conversion to binary64 is allowed only after the adversarial
classification.

This route shares the already audited local Lorentzian angle and area
primitives, but it does not reuse the primary decisive assembly, row order,
dense multiplication or arithmetic precision.

## 3. Frozen derivative hierarchy

At 180 decimal digits use centered logarithmic angle derivatives at

```text
h0 = 1e-25,
h1 = 5e-26,
h2 = 2.5e-26,
h3 = 1.25e-26.
```

Let `Q0,Q1,Q2,Q3` be the resulting directly pulled-back forms.  Construct
three Richardson forms

```text
R01=(4 Q1-Q0)/3,
R12=(4 Q2-Q1)/3,
R23=(4 Q3-Q2)/3.
```

The common area-Hessian and dust terms are retained because their Richardson
coefficient is `(4-1)/3=1`.

Require every base and displaced simplex to remain on the same Lorentzian
branch and agreement of the two finest Richardson local derivative tables
within their observed truncation hierarchy.  Require imaginary contamination
below `1e-140` on the assembled physical orbit kernel and on the complete
scalar action used in Section 6.  A branch or precision failure is
`CONTROL_FAILED`.

### Pre-rerun correction after the first failed execution

The first execution of implementation commit `dcb0d1c` exposed that the
phrase "raw imaginary contamination" was attached to the wrong intermediate
object.  Lorentzian dihedral angles and their derivatives are complex boost
variables and can have a nonzero physical imaginary component; multiplication
by the explicit `-i` in the Lorentzian Regge action is what produces the real
physical kernel.  The original high-precision implementation frozen as an
input to this protocol likewise gates the final kernel, not the raw angle
table.

The failed run found a maximum raw angle/derivative imaginary component of
`105.028...`, while the final even and odd orbit-kernel residues were only
`9.66e-155` and `1.52e-154`.  All four complete-action second derivatives
also had residual imaginary parts below the already frozen `1e-140` gate and
reproduced the kernel quadratic forms between `1.03e-88` and `3.94e-85`
relatively.  Consequently the corrected gate above tests reality of the
physical action-derived objects.  The threshold, derivative steps, parity
criterion and every decisive outcome gate remain unchanged.  The failed run
is preserved separately and must not be relabelled as a scientific pass.

## 4. Carrier reconstruction

Rebuild the carrier independently as sparse dictionaries indexed by physical
edge labels and logical vertex columns.  Use the exact generic coefficients
from the accepted symbolic theorem, not the primary stored matrices.  Require
the same exact support census, pole identity, connected non-bipartite rank
proof and collective identities as the primary protocol.

The sparse dictionary reconstruction must agree with the primary carrier
formula only after it has passed these internal controls.  Changing the
lexicographically first oriented-diagonal source-scale coefficient by
`+1/10` must change the finest adversarial form by more than `1e-12`
relatively.  This threshold is frozen before the adversarial value is known.

## 5. High-precision parity verdict

For each Richardson level `R` define

```text
Delta_R = R_even-R_odd,
N = max(1, ||R12_even||_F, ||R12_odd||_F),
d_R = ||Delta_R||_F/N.
```

The adversarial error floor is

```text
e_step = max_p(
           ||R01_p-R12_p||_F,
           ||R12_p-R23_p||_F
         )/N,
e_total = e_step+1e-145.
```

No observed primary difference enters this bound.

After the high-precision parity classification is fixed in memory, load the
primary operational matrices.  Require each primary parity form to agree
with `R12` to relative Frobenius error below `1e-10`.  Failure is not
automatically a refutation; it invokes the outcome hierarchy below.

## 6. Direct-action controls

Use the complete scalar action, not an assembled Hessian, on exactly two
preregistered unit data directions:

```text
x0: sigma_v=1/sqrt(120) for every v; all c_v=0,

x1: (sigma_0,sigma_1,c_2,c_3)=(1,-1,1,-1)/2;
    every other component is zero.
```

For each parity and direction, map `x` to active logarithmic edge variation
with its independently rebuilt carrier.  At 180 digits evaluate the complete
action at `delta=0,+t,-t,+t/2,-t/2` with `t=1e-20`, and form the Richardson
second derivative

```text
D_R=(4 D(t/2)-D(t))/3.
```

Require every displaced simplex to retain Lorentzian inertia `(3,1)` and
require

```text
|D_R-x^T R12 x|/max(1,|D_R|,|x^T R12 x|) < 1e-55.
```

These action-level controls are the independent check of the group
convolution and local derivative assembly.  No failed direction may be
replaced.

## 7. Hostile arithmetic control

Add to the odd active Hessian the conceptual rank-one perturbation

```text
u*u^T,
u=G_odd[:,0]/||G_odd[:,0]||_2.
```

Compute its induced carrier form algebraically.  Its relative Frobenius norm
must exceed `1e-6`; otherwise the adversarial sensitivity control fails.

Reversing the group product in the convolution is recorded as a convention
diagnostic but is not an acceptance gate, because a symmetric convolution
may make the two conventions accidentally equivalent.

## 8. Outcome hierarchy

### `FINITE_HEIGHT_QUADRATIC_ADVERSARIAL_CONTROL_FAILED`

Use if provenance, branch, group, carrier, action-level or hostile controls
fail.

### `FINITE_HEIGHT_QUADRATIC_PARITY_INDEPENDENT_ADVERSARIALLY_REPLICATED`

Use only if every control passes,

```text
max(d_R01,d_R12,d_R23) <= 10*e_total,
```

and both primary matrices agree with `R12` below `1e-10` relatively.

### `FINITE_HEIGHT_QUADRATIC_PRIMARY_REFUTED`

Use only if every control passes, the three Richardson differences agree
pairwise within `10*e_total`,

```text
min(d_R01,d_R12,d_R23) > 100*e_total,
```

and the direct-action controls agree with the adversarial rather than the
primary forms.  Any weaker disagreement is **OPEN**, not a refutation.

### `FINITE_HEIGHT_QUADRATIC_ADVERSARIAL_OPEN`

Use for the threshold gap, precision instability, or a primary/adversarial
matrix mismatch not resolved by the direct-action controls.

## 9. Interpretation firewall

Even successful replication proves only that the two staircase parities
induce the same one-sided quadratic form on this exact infinitesimal carrier
at this one finite-height background.  It does not prove nonlinear carrier
integrability, schedule-independent nonlinear evolution, gauge reduction,
gravitons, stability, a continuum limit, a physical tick, `c`, `G`, Planck
scales or particle masses.

Only the new adversarial verifier and static registry checks may run.  The
full suite remains forbidden unless the user explicitly requests it.
