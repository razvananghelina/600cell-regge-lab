# Precision correction preregistration: five lapse directions

Date: 2026-08-13

Prior-art gate and post-result update: `0882934`

Unresolved first-run record: `7d5e9fc`

Original frozen protocol: `41acf7b`

Status: **correction frozen before any arbitrary-precision reconstruction of
the five-dimensional Schur complement**.

The original `3e-4` convergence threshold is not altered and its `13/15`
result remains a failure.  This correction does not recalculate or search the
already separated thirty-dimensional regular spectrum.

## 1. Frozen input and block reduction

For each parity load the recorded symmetrized logarithmic internal matrix

```text
H = [[A, B],
     [B^T, C]],

A : 30 x 30 staircase-diagonal block,
C :  5 x  5 pole/lapse block.
```

First report the singular spectrum of `A`.  If `A` is not rank 30 at all
relative thresholds `1e-7, 1e-9, 1e-11`, stop this correction as
`BLOCK_REDUCTION_INVALID`.

With no fitting, define the five canonical Schur-lifted directions

```text
L_i = (-A^{-1} B e_i, e_i),    i=0,...,4,
```

where `e_i` is the printed pole-coordinate basis.  For the recorded matrix,
the quadratic form on these lifts is the Schur complement
`C-B^T A^{-1}B`.  The arbitrary-precision action, not the recorded matrix,
will determine its corrected values.

## 2. Frozen 80-decimal action reconstruction

Use the independent complete action-only implementation at 80 decimal
digits.  For a fixed logarithmic direction `w`, evaluate

```text
Q_h(w) = [S(x0 exp(h w))-2 S(x0)+S(x0 exp(-h w))]/h^2
```

at exactly

```text
h = 5.0e-4, 2.5e-4, 1.25e-4.
```

Apply

```text
R_coarse = (4 Q_2-Q_1)/3,
R_fine   = (4 Q_3-Q_2)/3,
Q_6      = (16 R_fine-R_coarse)/15.
```

No step may be changed after evaluating an action.  Use the fifteen frozen
directions

```text
L_i                         (5 directions)
L_i+L_j, 0 <= i < j <= 4  (10 directions).
```

Reconstruct the symmetric action quadratic matrix `Q` by polarization:

```text
Q_ii = Q_6(L_i),
Q_ij = [Q_6(L_i+L_j)-Q_ii-Q_jj]/2.
```

The corrected logarithmic per-edge Schur matrix is

```text
S_5 = Q/24.
```

Use

```text
epsilon_5 = norm_2(S_5-S_5_fine)
```

as an empirical, non-rigorous precision envelope, where `S_5_fine` is
reconstructed from `R_fine`.  Report the full matrices, eigenvalues,
eigenvectors, parity difference and all action imaginary parts.

## 3. Frozen collective/relative decomposition

Do not select a mode after seeing `S_5`.  Use the phase-collective vector and
its orthogonal projector

```text
c = (1,1,1,1,1)/sqrt(5),
P = I-c c^T.
```

Report:

- collective curvature `lambda_c=c^T S_5 c`;
- collective-relative mixing `norm(P S_5 c)`;
- the four eigenvalues of `P S_5 P` on a deterministic QR basis for
  `c^perp`;
- deviation from the permutation-symmetric form `a P+b c c^T`, reported but
  not imposed.

## 4. Frozen lapse-family control

There is a target-independent candidate for the collective tangent.  Keep
`l0` and the dust mass fixed and vary the common positive pole square `rho`
while imposing the published time-symmetric relation on every diagonal:

```text
q = l0^2-rho.
```

Its logarithmic tangent at the control is

```text
w_lapse = (-rho/q on all 30 diagonal coordinates,
            1       on all  5 pole coordinates).
```

Compare the normalized `w_lapse` with the smallest recorded eigenvector and
with the Schur lift `L c`.  Independently evaluate the certified analytic
35-equation residual at

```text
rho = rho0 * exp(eta),   eta in {-1e-3, 0, +1e-3},
q   = l0^2-rho,
y   = l0^2.
```

Use the unchanged mass and old/final boundaries.  The frozen stationarity
gate is maximum absolute per-edge residual `1e-7`.  Passing all three points
is computational evidence for a lapse family; it is not an analytic proof.

## 5. Frozen outcome labels

After implementation and branch controls, assign exactly one label per
parity.

### ONE_COLLECTIVE_NULL_FOUR_STIFF

All three lapse-family points pass, and

```text
abs(lambda_c)              <= 10 epsilon_5,
norm(P S_5 c)              <= 10 epsilon_5,
min(abs(relative eigs))    > 100 epsilon_5.
```

Interpretation: one collective lapse direction is numerically compatible
with gauge freedom, while four relative phase-lapse directions have resolved
curvature/pseudo-constraint stiffness.

### FIVE_NULL

All five absolute eigenvalues are at most `10 epsilon_5` and the lapse-family
control passes.

### FIVE_STIFF

All five absolute eigenvalues exceed `100 epsilon_5`.

### BLOCK_REDUCTION_INVALID

The thirty-dimensional `A` block is not robustly invertible.

### NUMERICALLY_UNRESOLVED

Every other situation.

The physical label is not a verifier PASS target.  PASS means only that the
frozen reconstruction and classification completed consistently.

## 6. Claim boundary

Even `ONE_COLLECTIVE_NULL_FOUR_STIFF` would not prove new physics.  Lapse
null modes and curvature-induced pseudo-constraints are **KNOWN** in Regge
calculus.  The result would establish only their explicit `1+4` realization
in this order-24, five-phase 600-cell dust sandwich.

It would also correct the next dynamical question.  With the five pole/lapse
coordinates prescribed, the regular `30 x 30` block can solve the thirty
diagonal evolution equations locally, while the five pole equations must be
treated as constraints/pseudo-constraints; the collective combination may be
redundant.  This is still one slab, fixed dust and a restricted invariant
sector.  It does not derive a clock, physical perturbations, a continuum
graviton, light speed or Planck units.
