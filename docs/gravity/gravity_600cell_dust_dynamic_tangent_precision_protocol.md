# Preregistration: high-precision audit of the dynamic tangent spectrum

Date: 2026-08-17

Prior-art gate: `16f4310`.

Status: frozen before evaluating any high-precision eigenvalue or power-trace
comparison between schedules.

## 1. Frozen input and exclusions

Use only
`reproducible/gravity_600cell_dust_dynamic_tangent.json`, whose required
SHA-256 is

```text
1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5.
```

Require its blind provenance commits `25722d9` and `0bceb9b`, outcome
`DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT`, `12/12` tests, and two 60 by 60
decimal tangent matrices.  Preserve the original artifact and classifier.

Do not load a continuum spectrum, desired degeneracy, limiting speed,
experimental number, or full-carrier result.  This audit asks only whether the
two already committed finite matrices have a calibration-resolved spectral
difference.

## 2. Shape matrices

Read the 50-decimal tangent matrices into `mpmath` and reconstruct exactly the
same deterministic scale/shape basis as the blind verifier:

```text
s=(1,...,1)/sqrt(30),
v=s+e_0,
H=I-2 vv^T/(v^T v),
B=(s,H[:,1],...,H[:,29]),
C=diag(B,B).
```

Delete scale-phase indices `(0,30)` from `C^T T C` to obtain each 58 by 58
shape matrix.  Recheck that the omitted mixing blocks remain below the stored
calibration.

The maximum entry-rounding uncertainty of a stored 60 by 60 matrix is bounded
in Frobenius, hence spectral, norm by

```text
epsilon_store = 60 * 0.5e-50.
```

For schedule `q`, set

```text
delta_q = stored epsilon_t_q + epsilon_store.
```

The stored `epsilon_t` is a calibrated proxy, not a theorem about the unknown
exact continuum action.

## 3. Basis-independent power traces

At 100 and 160 decimal digits, use the fixed normalization

```text
R = 2^20
A_q = shape_q/R.
```

For every `k=1,...,58`, compute

```text
p_q(k)=Tr(A_q^k).
```

The first 58 power sums determine the characteristic polynomial of a 58 by 58
matrix through Newton identities.  They are used here directly, without
forming the potentially cancellation-prone polynomial coefficients.

Let `sigma_q` be the stored largest shape singular value plus its stored
`epsilon_svd`.  The calibration propagation bound is

```text
u_q(k) = 58*k*((sigma_q+delta_q)/R)^(k-1)*(delta_q/R).
```

Add the 100-versus-160-digit discrepancy for each schedule and `1e-140` to
obtain `u(k)`.  Compare

```text
d(k)=abs(p_even(k)-p_odd(k)).
```

Classify the trace family as:

- `TRACE_SPECTRUM_CONSISTENT` if every `d(k) <= 10*u(k)`;
- `TRACE_SPECTRUM_DEPENDENT` if at least one `d(k) > 100*u(k)`;
- `TRACE_SPECTRUM_OPEN` otherwise.

Record all 58 differences, bounds and ratios.  No subset may be selected after
inspection.

## 4. High-precision eigendecomposition

Independently compute all shape eigenvalues and right eigenvectors at 100 and
160 decimal digits with `mpmath.eig`.  Match multisets with the Hungarian
algorithm; use binary64 only to select the assignment, then evaluate every
matched distance with the high-precision numbers.

At 160 digits require the normalized right-eigenpair residual

```text
||A V - V diag(lambda)||F /
max(1, ||A||F ||V||F)
```

below `1e-80`.  Let

```text
kappa_F(V)=||V||F*||inverse(V)||F.
```

This is a conservative consistent-norm Bauer--Fike factor.  For each schedule
define

```text
e_q = kappa_F(V_q)*delta_q
      + max matched 100-versus-160-digit eigenvalue distance
      + 1e-100.
```

For the optimally matched even/odd eigenvalue distance `d_eig`, classify:

- `EIG_SPECTRUM_CONSISTENT` if `d_eig <= 10*(e_even+e_odd)`;
- `EIG_SPECTRUM_DEPENDENT` if `d_eig > 100*(e_even+e_odd)`;
- `EIG_SPECTRUM_OPEN` otherwise.

## 5. Frozen combined verdict

- If both independent classifiers are `DEPENDENT`, return
  `SCHEDULE_DEPENDENCE_CONFIRMED_HIGH_PRECISION`.
- If both are `CONSISTENT`, return
  `SCHEDULE_SPECTRUM_NOT_RESOLVED`.
- Otherwise return `SCHEDULE_SPECTRUM_OPEN_NONNORMAL`.
- Any provenance, reconstruction, precision-convergence, inverse, or residual
  failure returns `PRECISION_AUDIT_CONTROL_FAILED` and no physical conclusion.

`SCHEDULE_SPECTRUM_NOT_RESOLVED` means only that the committed finite
calibration cannot distinguish the two characteristic spectra.  It is not a
proof of exact isospectrality and not a claim of continuum universality.

## 6. Outputs

Write a deterministic JSON artifact containing all controls, all trace-family
data, both high-precision eigenvalue multisets, condition factors,
classification labels, and the combined outcome.  Register the verifier before
its first execution.  Run only this targeted verifier, not the full suite.
