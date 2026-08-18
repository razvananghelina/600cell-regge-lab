# Preregistration: blind tangent map about the first dynamic dust tick

Date: 2026-08-17

Prior-art gate: `25722d9`.

Status: frozen before evaluating any dynamic-slab Hessian entry, tangent-map
entry, singular value or eigenvalue.

## 1. Frozen inputs and carrier

Require exact SHA-256 values

```text
accepted first dynamic tick artifact:
4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9,

old-to-final gluing artifact:
a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77,

audited canonical-Hessian source:
396c491fe51a9f5e04fa8402e2e5b16884fe23fc5057d8ded325e6064fbd3b9e,

audited one-slab action source:
ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf.
```

The first artifact must have outcome
`HOMOTHETIC_CANONICAL_LAPSE_SELECTED`, pass `7/7`, and contain both `even`
and `odd` schedules.  The gluing artifact must have outcome
`TWO_SLAB_GLUING_CONTROL_PASSED` and provide a permutation of all thirty
boundary orbits for each parity.

Use `DPS=100` and exactly the existing fixed mass, carrier, action, dust
world-lines and complex-angle branch.  No continuum harmonic, desired
degeneracy, speed or experimental number may be loaded.

## 2. Dynamic background

For each parity, read the accepted state `(s,r)` and reconstruct

```text
q_old[i] = L0^2,
q_new[i] = exp(2s)*L0^2,
rho      = rho0*exp(r),
x[0:30]  = exp(s)*L0^2-rho,
x[30:35] = rho.
```

At this point independently evaluate the complete action gradient in the 95
logarithmic orbit coordinates

```text
z=(o[30],x[35],n[30]),
g=(1/24) partial S/partial z.
```

Require all 35 internal equations below `1e-25`, reproduction of every stored
pre/post momentum component below `1e-45`, one timelike Gram direction for
all 2400 simplices, minimum leading minor positive, minimum angle-argument
modulus above `1e-6`, and maximum imaginary contamination below `1e-70`.

## 3. Complete calibrated Hessian

For each of the 95 coordinates, perturb the dynamic base multiplicatively in
its logarithmic coordinate.  Use exactly

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Thus each parity evaluates `95*4*2=760` displaced gradients plus its base.
Every displaced geometry must retain the frozen Lorentzian branch.

For each primary/shadow pair form the centered derivative matrix.  Require
every operational-primary minus validation-primary entry to lie below ten
times the sum of the corresponding operational and validation stability
proxies plus `1e-70`.  Require Hessian reciprocity in spectral norm below ten
times the combined calibrated spectral error.

No stored dynamic Hessian exists and none may be imported.

## 4. Canonical implicit matrix

With index sets

```text
O = old boundary, dimension 30,
X = internal,     dimension 35,
N = new boundary, dimension 30,
```

extract from each calibrated Hessian `K`

```text
J = [[ K_XX,  K_XN],
     [-K_OX, -K_ON]].
```

Define

```text
epsilon_J = ||Jop-Jop_shadow||2
          + ||Jval-Jval_shadow||2
          + ||Jop-Jval||2 + 1e-70.
```

Classify a singular value resolved nonzero only when it exceeds
`100*epsilon_J`, error-consistent zero only below `10*epsilon_J`, and open
otherwise.  The tangent map is constructed only if all 65 singular values
are resolved nonzero.

## 5. Tangent-map formula

Let the input column be `(delta o,delta p_pre)`.  The linearized internal and
pre-momentum equations give

```text
J [delta x; delta n] = R [delta o; delta p_pre],

R = [[-K_XO, 0],
     [ K_OO, I]].
```

After solving `Y=J^-1 R`, form the unpermuted output

```text
delta n      = Y_N [delta o;delta p_pre],
delta p_post = [K_NO,0][delta o;delta p_pre]
             + [K_NX,K_NN]Y[delta o;delta p_pre].
```

Apply the committed final-to-next-old permutation `P` to both output blocks:

```text
T = diag(P,P) [delta n;delta p_post].
```

This fixes one and only one `60 x 60` self-map per schedule.  Report

```text
N_maps=2,
N_distinct_maps,
```

where distinctness is measured only after the calibrated comparison below.

## 6. Canonical controls

In phase order `(delta z,delta p)`, use

```text
Omega = [[0,I],[-I,0]],
Phi(T)=T^T Omega T-Omega.
```

Build `Top`, `Top_shadow`, `Tval`, `Tval_shadow`.  Define

```text
epsilon_T = ||Top-Top_shadow||2
          + ||Tval-Tval_shadow||2
          + ||Top-Tval||2 + 1e-70,

epsilon_sym = ||Phi(Top)-Phi(Top_shadow)||2
            + ||Phi(Tval)-Phi(Tval_shadow)||2
            + ||Phi(Top)-Phi(Tval)||2 + 1e-70.
```

Require

```text
||Phi(Top)||2 <= 10*epsilon_sym.
```

For the descending singular values `sigma_i`, form the reciprocal-pair vector

```text
c_i(T)=sigma_i*sigma_(59-i)-1,  i=0,...,29.
```

Apply the identical operational/shadow/validation difference construction to
`c(T)` and require `||c(Top)||_infinity` below ten times that proxy.  This is
an independent consequence of symplecticity, not a fitted spectrum.

## 7. Frozen scale/shape decomposition

Let

```text
s=(1,...,1)/sqrt(30)
```

and construct its deterministic Householder complement `Y_shape`.  The
orthogonal boundary basis `B=(s,Y_shape)` gives the canonical phase basis
`C=diag(B,B)`.

In `T_sector=C^T T C`, the two scale-phase indices are `(0,30)` and the other
58 indices are the zero-sum shape phase space.  Let `M(T)` concatenate both
off-diagonal scale/shape blocks and define its calibrated proxy from all four
maps exactly as in `epsilon_T`, using the Frobenius norm.

Per parity classify:

- `SCALE_SHAPE_INVARIANT` if `||M(Top)||F <= 10*epsilon_mix`;
- `SCALE_SHAPE_MIXED` if `||M(Top)||F > 100*epsilon_mix`;
- `SCALE_SHAPE_OPEN` otherwise.

Only in the invariant case may the `58 x 58` shape block be called a closed
one-step subsystem.  Its singular values and eigenvalues are still blind
linear-algebra data, not graviton frequencies.

## 8. Blind spectra and schedule comparison

For every complete map record:

- all 60 singular values and reciprocal products;
- all 60 complex eigenvalues and moduli;
- spectral radius, determinant and eigenvector condition number;
- counts of eigenvalues consistent with, open around, or resolved away from
  the unit circle using the calibrated eigenvalue comparison below.

If scale/shape is invariant, record the same data for the 2-dimensional scale
block and 58-dimensional shape block.

Match complex eigenvalue multisets with the Hungarian algorithm.  For each
parity, define `epsilon_eig` as the maximum optimal-matching distance between
the operational spectrum and each of its shadow/validation spectra, plus the
scaled independent binary64 eigensolver floor

```text
1e-15*max(1,max absolute operational eigenvalue).
```

Define `epsilon_svd` analogously from the largest componentwise difference of
the ordered singular spectra, with floor
`1e-15*max(1,sigma_max)`.  Classify an operational eigenvalue as
unit-consistent below `10*epsilon_eig`, resolved off-unit above
`100*epsilon_eig`, and open otherwise.  This is an empirical sensitivity
certificate, not a pseudospectral theorem.

For schedule comparison, optimally match even and odd spectra and separately
compare their ordered singular values.  Let the uncertainty be the sum of
the two within-schedule calibrated uncertainties.  Assign:

- `SCHEDULE_ROBUST` only if both spectral distances are below ten times their
  uncertainty;
- `SCHEDULE_DEPENDENT` if either exceeds one hundred times its uncertainty;
- `SCHEDULE_OPEN` otherwise.

No eigenvalue cluster is compared with continuous `S^3` before this entire
blind artifact is committed.

## 9. Mechanical outcomes

Assign exactly one:

1. `DYNAMIC_TANGENT_CONTROL_FAILED`;
2. `DYNAMIC_CANONICAL_RANK_OPEN`;
3. `DYNAMIC_CANONICAL_DEGENERATE`;
4. `DYNAMIC_TANGENT_SYMPLECTICITY_FAILED`;
5. `DYNAMIC_SCALE_SHAPE_MIXED`;
6. `DYNAMIC_SCALE_SHAPE_OPEN`;
7. `DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT`;
8. `DYNAMIC_SHAPE_TANGENT_SCHEDULE_OPEN`;
9. `DYNAMIC_SHAPE_TANGENT_BLIND_CENSUS_CERTIFIED` only if both maps are
   calibrated symplectic, scale/shape invariant and schedule robust.

The verifier passes when it reconstructs the frozen object and assigns the
outcome mechanically.  Only outcome 9 licenses a later continuum-harmonic
comparison inside this quotient.  A schedule-dependent result remains useful
as a discretization diagnostic but cannot be called universal physics.

## 10. Claim boundary

No outcome identifies gauge-invariant tensor modes, covers the full 720-edge
boundary, proves long-time stability, derives a dispersion relation, selects
a proper-time unit or yields a limiting speed.  A continuous compact `S^3`
already has a discrete spectrum; discreteness, a gap or degeneracy alone is
not evidence for new physics.

Only the new targeted verifier will be run.  The full suite will not be run.
