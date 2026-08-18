# Preregistered protocol: shifted negative-shape persistence

Date: 2026-08-18

Status: **TARGET-DISCLOSED, PREREGISTERED BEFORE THE SHIFTED SHAPE
STIFFNESS OR ANY OLD/NEW PROJECTOR DISTANCE WAS COMPUTED.**

Prior-art gate commit: `33da8dd`.

The old result is disclosed: on the first centered nonstationary recurrence,
sectors `4` and `5` each contain a resolved `15`-dimensional negative
stiffness space, giving a certified `30`-position autonomous subsystem.  The
new evidence is whether that construction survives at the next independently
accepted slab.  Reproducing `30` is not sufficient: the ambient subspaces and
their dynamic invariance must also be tested.

## 1. Frozen inputs

| input | SHA-256 |
|---|---|
| `commons/cell600.py` | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |
| `gravity_600cell_dust_conformal_supermetric.json` | `b38d55f9f575ddffd34edeaa5e835d9e10919e6d96a0c284d73c31a072675025` |
| `verify_gravity_600cell_dust_conformal_supermetric.py` | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| `gravity_600cell_dust_centered_jacobi.json` | `fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56` |
| `gravity_600cell_dust_centered_jacobi.npz` | `1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef` |
| `verify_gravity_600cell_dust_centered_jacobi.py` | `359b8d7642746c2dc22e304353e3b83104874badd86755de4f8f9e6f25e56a20` |
| `gravity_600cell_dust_shifted_centered.json` | `265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47` |
| `gravity_600cell_dust_shifted_centered.npz` | `c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8` |
| `verify_gravity_600cell_dust_shifted_centered.py` | `a3c45e3e636057d83a663d3248dd023f7d04ec6e544c698f9116307822be337a` |
| `gravity_600cell_dust_conformal_shape_dynamics.json` | `c5bbeaa2a64d07688061bc5098a33361dc2f5300d637e44a10b6cccbbd1bb162` |
| `verify_gravity_600cell_dust_conformal_shape_dynamics.py` | `52857835b37722db51c03587a9583426b26caaf2cb6b2d55c4fee05419883112` |
| `gravity_600cell_dust_shape_stiffness.json` | `03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868` |
| `verify_gravity_600cell_dust_shape_stiffness.py` | `d4f0a9a805910de37011ba70f407907daa2d11c650aeea22e571ab867282a44c` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |

Require the upstream outcomes `CENTERED_JACOBI_CERTIFIED`,
`SHIFTED_CENTERED_CERTIFIED`,
`CONFORMAL_SHAPE_DYNAMICS_DECOUPLED_POWER_CERTIFIED` and
`SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED` with their exact pass counts.

## 2. Complete hypotheses

Every result is conditional on:

1. the fixed labelled regular 600-cell and its `720` logarithmic
   signed-squared spatial-edge variables;
2. the first three independently accepted nonstationary fixed-total-mass
   dust-Regge slabs;
3. the literal adjacent-slice edge identification of the Jacobi recurrence;
4. both staircase schedules, seven frozen minimal binary-tetrahedral sectors
   of dimensions `3,2,2,2,1,1,1`, and all four derivative variants;
5. the unsigned conformal incidence `C[e=(u,v),w]=delta_uw+delta_vw`;
6. the action-relative split

   ```text
   K=im C,
   S_j=ker(C* H(M_j)),
   H(X)=(X+X*)/2;
   ```

7. the committed component midpoint/radius enclosures, the floating floor
   `1000 eps n max(1,norm)` and the existing `10/100` bands;
8. no independent dust perturbations, exact constraint quotient, proper-time
   unit, continuum refinement or tensor-harmonic target.

The old `30` is a certified lower count because other old sectors contain
open stiffness signs.  The test therefore concerns persistence of the
specific two resolved old sector carriers, not a theorem that the full
negative spectral space has dimension exactly `30`.

## 3. Stage A: shifted conformal/shape closure

Reconstruct the exact conformal carrier and, for the shifted centered archive,
the action-relative shape basis `W_2` using the same SVD ranks and error
formulae as the committed closure verifier.  For every stored shifted
`Gamma_2` and `Omega_2`, test invariance of both `im C` and
`ker(C*H(M_2))`.

The complete census is

```text
2 schedules * 7 sectors * 4 variants
* 2 operators * 2 carriers = 224 residuals.
```

Use the frozen labels `ZERO_CONSISTENT`, `NONZERO_RESOLVED`, `OPEN`.  A single
resolved cross residual refutes shifted closure.  No stiffness spectrum is
computed in this stage.

## 4. Stage B: blind shifted stiffness census

Only after the Stage-A source and registration commit, reconstruct all

```text
2 schedules * 7 sectors * 4 variants = 56
```

shifted shape pencils, without loading the old stiffness artifact or the old
sector targets.  For each cell form

```text
M_S = W_2* H(M_2) W_2,
V_S = W_2* H(V_2) W_2,
B_2 = -M_S,
A_2 = -V_S,
Omega_S = W_2* Omega_2 W_2.
```

Apply verbatim the definite-kinetic, Hermitian-inertia,
normalized-recurrence, compatibility and schedule rules from the committed
blind stiffness protocol.  Record all `2,400` representative and `4,800`
full-multiplicity signs per object.  Commit this blind census before the
old/new comparison is executed.

## 5. Stage C: target-disclosed persistence comparison

Only after Stage B is committed may a comparison verifier load both old and
shifted centered archives and the two blind stiffness artifacts.

### 5.1 Rank and sector support

For each schedule and variant, reconstruct the old and new stiffness
eigenspaces independently.  Record every sector's complete resolved inertia.
The disclosed old target is

```text
sector 4: 15 negative resolved + 10 positive resolved,
sector 5: 15 negative resolved + 10 positive resolved.
```

For a comparison carrier to exist, the shifted cells in both sectors must
have the same fully resolved `15+10` split and no zero/open entry.  Additional
resolved shifted negative sectors are recorded as an expansion of the
certified carrier and refute an *exactly-30-only* reading, but do not erase a
persisting embedded sector-4/5 carrier.

### 5.2 Ambient projectors, not eigenvector columns

For a qualifying sector/cell, let `E_j` be the orthonormal negative
eigenbasis of `A_j` and define in the common literal edge-sector coordinates

```text
F_j = W_j E_j,
P_j = F_j F_j*.
```

If `g_j` is the negative/positive eigengap and `epsilon_Aj` the complete
restricted-form error, require `g_j>2 epsilon_Aj` and set

```text
eta_Ej = 2 epsilon_Aj/(g_j-2 epsilon_Aj)
         + 1000 eps (25d),
eta_Pj = 2 eta_Sj + 2 eta_Ej
         + 1000 eps (30d).
```

The deliberately conservative factor two covers lifting the internal
projector through the uncertain shape embedding.  Compute

```text
distance = ||P_2-P_1||_2 = sin(theta_max),
error    = eta_P1 + eta_P2
           + 1000 eps (30d) max(1,distance).
```

Classify every one of the `16` disclosed cells:

```text
distance <= 10 error    COMMON_CONSISTENT,
distance > 100 error    ROTATION_RESOLVED,
otherwise               PERSISTENCE_OPEN.
```

Also record all principal angles and `||P_2-P_1||_F`; neither may replace the
operator-norm decision.

### 5.3 Shifted dynamic invariance

For `X=Gamma_{S,2},Omega_{S,2}`, test

```text
||(I-E_2 E_2*) X E_2||_2
```

using the identical inherited subspace/error propagation from the old
negative-dynamics protocol.  There are `16*2=32` classifications.  A reduced
shifted recurrence exists only when both operators are
`INVARIANT_CONSISTENT` in a cell.

### 5.4 Non-autonomous product gate

A product of consecutive reduced companions may be constructed only if all
`16` projector comparisons are `COMMON_CONSISTENT` and all `32` shifted
invariance tests pass.  Use one of the already computed orthonormal bases as
the common coordinate system and propagate projector, leakage and matrix
errors into both steps before multiplying

```text
C_21 = C_2 C_1.
```

Its singular values and complex eigenvalues are diagnostics only.  One
two-step product is neither a Lyapunov exponent nor a physical stability
theorem.  If any projector is `ROTATION_RESOLVED`, do not use a polar,
Procrustes or overlap unitary: the action-selected temporal connection remains
`OPEN`.

## 6. Frozen outcome hierarchy

Stage A:

1. `SHIFTED_CONFORMAL_SHAPE_CONTROL_FAILED`;
2. `SHIFTED_CONFORMAL_SHAPE_CARRIER_OPEN`;
3. `SHIFTED_CONFORMAL_SHAPE_MIXING_REFUTED`;
4. `SHIFTED_CONFORMAL_SHAPE_DYNAMICS_OPEN`;
5. `SHIFTED_CONFORMAL_SHAPE_DYNAMICS_DECOUPLED`.

Stage B uses the old blind stiffness ladder with the prefix
`SHIFTED_SHAPE_STIFFNESS_`.

Stage C applies the first matching branch:

1. `SHIFTED_NEGATIVE_PERSISTENCE_CONTROL_FAILED`;
2. `SHIFTED_NEGATIVE_RANK_OR_SECTOR_CHANGED` if either disclosed shifted
   sector lacks its complete `15+10` split;
3. `SHIFTED_NEGATIVE_CARRIER_OPEN` if a required gap or projector enclosure
   is unresolved;
4. `SHIFTED_NEGATIVE_DYNAMICS_MIXED` for any resolved shifted leakage;
5. `SHIFTED_NEGATIVE_DYNAMICS_OPEN` for unresolved shifted leakage;
6. `SHIFTED_NEGATIVE_BUNDLE_ROTATED` for any resolved old/new projector
   rotation;
7. `SHIFTED_NEGATIVE_PERSISTENCE_OPEN` for an otherwise unresolved projector
   comparison;
8. `SHIFTED_NEGATIVE_COMMON_CARRIER_PRODUCT_FAILED` if common-carrier
   prerequisites pass but the finite product construction contradicts its
   algebraic identities or error regularity;
9. `SHIFTED_NEGATIVE_COMMON_CARRIER_CERTIFIED` only if the same carrier and
   both consecutive reduced recurrences are certified in all `16` cells.

Every correctly resolved scientific branch, including refutation, is a valid
passing verifier outcome.

## 6a. Disclosed post-first-run implementation correction

The first blind Stage-B execution produced the complete numerical census but
assigned `SHIFTED_SHAPE_STIFFNESS_COMPLEX_OR_ACTION_MISMATCH`.  The artifact
was not committed; its SHA-256 was

```text
102aeaced347bb5684da8437cbc139d2e865f5baeb293692562cff3f3b38cb18.
```

Inspection showed `56/56` compatibility residuals `ZERO_CONSISTENT`, no
resolved complex eigenvalue, and `16` count flags caused solely by inherited
implementation logic that treated `ZERO_CONSISTENT` as a resolved sign.
That is incorrect: a zero-consistent interval may contain either sign, so a
positive-resolved classification from the better-conditioned Hermitian form
does not contradict a zero-consistent classification from the normalized
matrix.

Before any old/new target comparison, the implementation is corrected to
compare inertia counts only when *every* entry in both lists is strictly
`POSITIVE_RESOLVED` or `NEGATIVE_RESOLVED`; any `ZERO_CONSISTENT` or `OPEN`
entry disables the count comparison for that cell.  This is the literal
meaning of the already-preregistered phrase "both sign lists ... completely
resolved."  No midpoint, error, carrier, sign census or target is changed.
The initial blind census itself is disclosed unchanged:

```text
Hermitian full multiplicity: 2000 positive, 0 negative,
                              2336 zero-consistent, 464 open;
normalized full multiplicity: 0 positive, 0 negative,
                              4592 zero-consistent, 208 open.
```

Thus the corrected expected branch is selected by the frozen hierarchy from
the unchanged facts; it is not a new fitted numerical test.

## 7. Interpretation boundary

- **DERIVED COMPUTATIONAL** may describe only the finite operators and their
  certified subspaces.
- Same rank/sectors with resolved rotation is **STRUCTURAL**, not a common
  propagating carrier.
- A common carrier would be genuine progress toward a finite nonautonomous
  wave subsystem, but identifying it with gravitational waves remains
  **OPEN**.
- No branch derives a continuum dispersion relation, `c`, a physical tick,
  Planck units, two polarizations or particle masses.
