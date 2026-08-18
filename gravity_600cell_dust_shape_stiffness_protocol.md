# Blind protocol: action-relative shape stiffness census

Date: 2026-08-18

Status: **PREREGISTERED BEFORE ANY SHAPE-STIFFNESS SPECTRUM WAS COMPUTED.**

The source archive was inspected only for array names and dimensions.  No
restricted `V`, `Omega`, generalized eigenvalue, inertia count, desired
degeneracy, continuum spectrum or sign target was read before this protocol.

## Frozen provenance

```text
prior-art gate commit
  e37d80c

reproducible/gravity_600cell_dust_centered_jacobi.json
  fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56

reproducible/gravity_600cell_dust_centered_jacobi.npz
  1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef

reproducible/gravity_600cell_dust_conformal_shape_dynamics.json
  c5bbeaa2a64d07688061bc5098a33361dc2f5300d637e44a10b6cccbbd1bb162

reproducible/verify_gravity_600cell_dust_conformal_shape_dynamics.py
  52857835b37722db51c03587a9583426b26caaf2cb6b2d55c4fee05419883112

reproducible/verify_gravity_600cell_dust_conformal_supermetric.py
  d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4

reproducible/verify_gravity_600cell_dust_full_boundary_tangent.py
  c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf

commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

## Complete hypotheses

The result will be conditional on all of the following:

1. the fixed regular 600-cell and its `720` logarithmic signed-squared spatial
   edge variables;
2. the first two accepted nonstationary fixed-mass dust-Regge slabs;
3. the literal adjacent-slice edge identification of the centered Jacobi
   equation;
4. the seven frozen minimal binary-tetrahedral sectors of dimensions
   `3,2,2,2,1,1,1`, both staircase schedules and four derivative variants;
5. the unsigned conformal incidence map `(C sigma)_uv=sigma_u+sigma_v`;
6. the action-relative split

   ```text
   K=im C,  S_H=ker(C*H_M),  H_M=(M+M*)/2;
   ```

7. the committed componentwise midpoint/radius enclosures and the same
   floating-point safety multiplier `1000*eps*n`;
8. no dust perturbation variables, constraint quotient, proper-time unit,
   continuum refinement or spatial tensor-harmonic target.

## Enumeration count

There is no choice of carrier or subspace.  The complete census contains

```text
2 schedules * 7 sectors * 4 variants = 56 shape pencils.
```

Their shape dimensions are `25d`; counting multiplicity, this is

```text
2 * 4 * sum(25d) = 4,800 generalized eigenvalue instances.
```

The actual normalized shape blocks contribute the same number of eigenvalue
instances.  All are recorded; none may be discarded as an outlier.

## Reconstructed carrier

For each schedule and sector, reconstruct the conformal basis `U` from the
literal incidence map exactly as in the certified closure verifier.  For each
derivative variant reconstruct an orthonormal basis `W` of

```text
ker(U* H_M)
```

by a full SVD.  Reuse the committed `10/100` singular-value bands and the
Wedin-style subspace bounds.  Any unresolved carrier rank or failed direct
sum is a control failure; no spectrum is interpreted in that cell.

## Two objects that must remain distinct

For each of the 56 cells form

```text
M_S     = W* H_M W,             H_M=(M+M*)/2,
V_S     = W* H_V W,             H_V=(V+V*)/2,
B       = -M_S,
A       = -V_S,
Omega_S = W* Omega W.
```

The overall minus sign makes the already certified shape kinetic convention
positive; multiplying the complete action by `-1` changes neither the pencil
nor its generalized eigenvalues.

### Object 1: Hermitian action pencil

If `B` is positive definite, solve

```text
A x = lambda B x
```

with a Hermitian definite solver.  This is algebraically the same pencil as
`V_S x=lambda M_S x`.

Do not infer signs from condition-sensitive ordered generalized eigenvalues.
By Sylvester inertia under the congruence

```text
B^(-1/2) A B^(-1/2),
```

their positive/zero/negative counts equal the inertia of `A`.  Classify the
ordinary Hermitian eigenvalues of `A` against the complete restricted-form
error.

### Object 2: actual normalized recurrence block

Compute all eigenvalues of `Omega_S`.  Its eigenvector condition number enters
the Bauer--Fike error.  Classify reality and real sign independently from the
Hermitian pencil.

Also compute the action-compatibility residual

```text
R = M_S Omega_S - V_S.
```

This residual, rather than visual agreement of plotted eigenvalues, tests
whether the actual normalized block realizes the symmetrized action pencil on
the certified shape carrier.

When both sign lists in a cell are completely resolved, their inertia counts
must agree.  A resolved count disagreement is an additional action-mismatch
flag; near-zero/open entries are not forced into either sign.

## Frozen error propagation

For every source midpoint `X` with re-enclosed component radius `R_X`, use

```text
epsilon_X = ||R_X||_F
            + 1000 eps n max(1, ||X||_2).
```

For a restriction to a basis with subspace error `eta_S`, use

```text
epsilon_XS = epsilon_X
             + 2 eta_S (||X||_2 + epsilon_X)
             + 1000 eps n max(1, ||X||_2).
```

For `R=M_S Omega_S-V_S`, use

```text
epsilon_R = epsilon_MS (||Omega_S|| + epsilon_OmegaS)
            + ||M_S|| epsilon_OmegaS
            + epsilon_VS
            + 1000 eps s max(1, ||M_S|| ||Omega_S||, ||V_S||),
```

where `s=25d`.

For the actual eigenvalues use one common cell bound

```text
epsilon_eig = kappa(Z) epsilon_OmegaS
              + 1000 eps s max(1, ||Omega_S||),
```

where `Z` is the right-eigenvector matrix.  A nonfinite or singular `Z` makes
the cell `OPEN`.

For the ordered Hermitian generalized eigenvalues, let

```text
b = lambda_min(B),       b_lower = b - epsilon_B.
```

When `b_lower>0`, the min--max/Rayleigh-quotient perturbation envelope is

```text
epsilon_pencil = epsilon_A / b_lower
                 + ||A|| epsilon_B / (b b_lower)
                 + 1000 eps s max(1, ||B^-1 A||).
```

If `b_lower<=0`, the pencil is open.  This bound is used for ordered spectrum
and schedule comparisons; sign multiplicities continue to be taken from the
more direct inertia test on `A`.

## Frozen classifications

Every scalar `x` with complete error `e` uses the existing bands:

```text
|x| <= 10e       ZERO_CONSISTENT
|x| > 100e       NONZERO_RESOLVED
otherwise        OPEN.
```

For a Hermitian eigenvalue of `A`:

```text
x >  100e        POSITIVE_RESOLVED
x < -100e        NEGATIVE_RESOLVED
|x| <= 10e       ZERO_CONSISTENT
otherwise        OPEN.
```

For an eigenvalue `z` of `Omega_S`, classify its imaginary part first:

```text
|Im z| <= 10e    REAL_CONSISTENT
|Im z| > 100e    COMPLEX_RESOLVED
otherwise        REALITY_OPEN.
```

Only `REAL_CONSISTENT` values receive the same real-sign classification.

The kinetic form `B` is `POSITIVE_DEFINITE_RESOLVED` only if its smallest
eigenvalue exceeds `100 epsilon_B`.  A negative eigenvalue below
`-100 epsilon_B` is a control contradiction; the intermediate case is open.

The compatibility residual is zero/nonzero/open by its operator norm and
`epsilon_R` using the same bands.

## Schedule comparisons

For every sector, variant and object, compare the complete ordered even/odd
eigenvalue lists, including multiplicity.  The comparison error is the sum of
the two `epsilon_pencil` values for the action pencil, or the two
`epsilon_eig` values for `Omega_S`, plus the floating-point floor.

```text
distance <= 10 error    SCHEDULE_ROBUST
distance > 100 error    SCHEDULE_DEPENDENT
otherwise               SCHEDULE_OPEN.
```

There are `7*4*2=56` declared schedule comparisons.  A dependent comparison
blocks a schedule-independent physical reading but does not delete either
spectrum.

## Outcome ladder fixed before execution

Apply the first matching outcome:

1. `SHAPE_STIFFNESS_CONTROL_FAILED` if provenance, carrier reconstruction or
   definite-kinetic controls contradict their frozen expectations;
2. `SHAPE_STIFFNESS_CARRIER_OPEN` if any carrier or kinetic control is open;
3. `SHAPE_STIFFNESS_COMPLEX_OR_ACTION_MISMATCH` if any normalized shape
   eigenvalue is resolved complex, any compatibility residual is resolved
   nonzero, or two completely resolved sign lists disagree;
4. `SHAPE_STIFFNESS_SCHEDULE_DEPENDENT` if any schedule comparison is
   dependent;
5. `SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED` if either complete sign census
   has any resolved negative entry;
6. `SHAPE_STIFFNESS_SIGN_OPEN` if any required sign, reality, compatibility or
   schedule classification remains open;
7. `SHAPE_STIFFNESS_NONNEGATIVE_WITH_ZERO_MODES` if all remaining signs are
   resolved and at least one is zero-consistent;
8. `SHAPE_STIFFNESS_POSITIVE_CENSUS` if all `4,800+4,800` entries are resolved
   real positive and every other gate passes.

Outcomes 7 or 8 are **DERIVED COMPUTATIONAL NECESSARY DIAGNOSTICS**, not a
graviton theorem.  Outcome 5 refutes blanket positive shape stiffness but
does not by itself diagnose a ghost.  No outcome licenses fitting to a desired
continuum spectrum.

## Explicit forbidden comparisons

The verifier must not load or compare:

- an `S^3` scalar, vector or tensor Laplacian spectrum;
- desired degeneracies or a two-polarization count;
- a target dispersion relation;
- the speed of light, Planck time, Planck mass or particle data;
- a refinement trend not already constructed.

These flags must be written as `false` in the artifact.
