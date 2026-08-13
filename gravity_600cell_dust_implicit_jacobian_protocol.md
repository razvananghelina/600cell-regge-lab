# Preregistration: implicit Jacobian at the 600-cell dust sandwich

Date: 2026-08-13

Prior-art gate commit: `31717a8`

Upstream external-control commit: `66a6465`

Status: **protocol frozen before the first evaluation of any derivative of
the dust-solution residual map**.

No boundary direction, eigenvector, singular value, rank or condition number
has been inspected before this protocol.  No measured quantity, `a1=5`,
Planck scale or desired speed is compared.

## 1. Frozen carrier, source and equations

For each of the two certified ordered-schedule parities, reuse without change:

- the 2400-simplex Lorentzian staircase slab and its 100 order-24 simplex
  orbits;
- the exact angle branch and Regge curvature action;
- the unrounded De Felice--Fabri dust mass and time-symmetric sandwich;
- the published dust world-line action on the five pole orbits;
- the old regular boundary, held fixed.

Write the positive internal squared-magnitude variables and positive final
boundary squared lengths as

```text
x = (q_0,...,q_29,rho_0,...,rho_4),  dim x = 35
y = (f_0,...,f_29),                  dim y = 30.
```

Here the actual Lorentzian pole squares are `-rho_k`.  The five `rho_k` are
solved internal variables, not a fixed lapse.  The dust mass is fixed.

Use target-free logarithmic coordinates around the published solution:

```text
u_i = log(x_i/x0_i),    v_j = log(y_j/y0_j).
```

For the common orbit multiplicity 24 define the logarithmic local equations

```text
E_i(u,v) = (x_i/24) * partial S_total / partial x_i.
```

They have the same zero set and the same local rank as the ordinary internal
Regge equations.  At an exact stationary point,

```text
J_x = partial E / partial u
```

is `1/24` times the Hessian in logarithmic internal coordinates and is
symmetric.  Also compute

```text
J_y = partial E / partial v,                         shape 35 x 30
K_x = partial[(y/24) partial S_total/partial y]/partial u,
```

for which Hessian reciprocity requires `J_y = transpose(K_x)` at the
stationary point.

## 2. Frozen differentiation and convergence audit

Differentiate the already certified analytic complete-action gradient by
centered differences in logarithmic coordinates at exactly

```text
h = 1.0e-3, 5.0e-4, 2.5e-4.
```

For each step, evaluate every plus/minus point for all 65 coordinates.  At
every point require all 100 representative simplices to remain Lorentzian,
nondegenerate and on the certified real branch.

Let `J_mid` and `J_fine` denote the internal matrices at the last two steps.
Use the frozen centered-difference Richardson estimate

```text
J_R = (4 J_fine-J_mid)/3,
epsilon_emp = norm_2(J_R-J_fine).
```

`epsilon_emp` is an empirical truncation indicator, not a rigorous error
bound.  Report:

- entrywise and spectral-norm step convergence;
- imaginary contamination relative to the largest matrix scale;
- antisymmetric norm of every `J_x` and `J_R`;
- cross-reciprocity error between `J_y` and `transpose(K_x)`;
- singular values, eigenvalues and inertia of the symmetrized `J_R`;
- ranks at relative thresholds `1e-7`, `1e-9`, `1e-11` for all three steps
  and `J_R`;
- `condition_2(J_R)`, `s_min(J_R)` and `s_min/epsilon_emp`;
- the corresponding ranks of `J_y` and `[J_x | J_y]`;
- the response matrix `X=-J_R^{-1}J_y` when numerically available, its
  linear residual and its singular spectrum;
- phase-parity differences without assuming equality.

Frozen implementation tolerances are

```text
relative imaginary contamination       < 3e-7
relative internal antisymmetry          < 3e-6
relative cross-reciprocity error        < 3e-6
relative fine-versus-mid change         < 3e-4
minimum Gram modulus at all points      > 1e-8
minimum angle-argument modulus          > 1e-6.
```

A tolerance failure is retained and classified; no step may be changed after
seeing the spectrum.

## 3. Frozen independent action-only curvature checks

The matrix above differentiates an analytic gradient.  Independently use the
already certified 60-decimal action-only implementation to compute centered
second differences of the complete action at logarithmic steps `5e-4` and
`2.5e-4`, followed by the same Richardson formula.

Test these four preregistered unit directions in the internal 35-space:

1. the normalized all-ones vector;
2. the normalized alternating vector `(-1)^i`;
3. the normalized vector supported equally on the five pole variables;
4. the unit eigenvector of the symmetrized `J_R` having smallest absolute
   eigenvalue, with sign fixed by making its first nonzero component positive.

The fourth direction is result-dependent but not selected for agreement: it
is deterministically the weakest direction of the frozen matrix and is
specified before evaluation.  Compare each action curvature with

```text
24 * w^T J_R w.
```

Report the relative errors.  The control target is `3e-4`; failure makes the
nonsingularity verdict numerically unresolved rather than licensing a new
step or direction.

## 4. Frozen outcome labels

Assign exactly one label per parity.

### ROBUST_NUMERICAL_REGULARITY

All implementation controls pass, every step and `J_R` has rank 35 at all
three relative thresholds, and

```text
s_min(J_R) > 100 * epsilon_emp,
```

with the weakest-direction 60-decimal curvature control passing.

This is strong computational evidence for nonsingularity.  It is **not** an
interval proof that the exact determinant is nonzero.

### RESOLVED_NUMERICAL_NULLITY

A rank below 35 is stable across all steps and thresholds, the null cluster is
separated from the remaining spectrum by at least `100*epsilon_emp`, and the
60-decimal weakest-direction curvature is compatible with that cluster.

### NUMERICALLY_UNRESOLVED

Any other spectral situation, including a smallest singular value comparable
with the empirical error, unstable rank, failed reciprocity or failed
action-only curvature control.

The verifier passes when it reconstructs the carrier, evaluates the frozen
audit and applies these labels consistently.  Full rank is a scientific
outcome, not a test-suite PASS target.

## 5. Acceptance boundary and claim discipline

If both parities have `ROBUST_NUMERICAL_REGULARITY`, the analytic
implicit-function theorem says that exact nonsingularity would yield a unique
local internal response `x(y)` for every sufficiently small permitted
30-coordinate final-boundary perturbation.  Until a validated error bound or
interval certificate proves exact nonsingularity, label this **DERIVED
COMPUTATIONAL EVIDENCE**, not a theorem.

If that numerical gate passes, the next separately preregistered calculation
may choose deterministic boundary basis vectors and solve displaced roots.
No displaced root is searched here.

Scope remains strict:

- only the order-24 invariant sector of one chosen five-phase schedule;
- not the unreduced 840 internal edges;
- one slab, not a continuum or multi-tick spacetime;
- fixed externally supplied dust mass and world-line placement;
- no physical/gauge mode separation and no graviton claim;
- no derivation of time, light speed, inertia or a Planck scale.
