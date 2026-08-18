# Preregistration: precision-corrected gauge quotient

Date: 2026-08-13

Prior-art gate: `44bd4cf`

Exact collective lapse result: `790fc7f`

Five-direction high-precision Schur result: `dc927a5`

First gauge-quotient result: `14a4517`

Status: **frozen before evaluating the new `1.25e-4` mixed-action
rectangles or constructing the corrected quotient**.

The old `OPEN NUMERICALLY` verdict is not changed retroactively.  This is a
precision correction with deterministic inputs and maps.  It does not relax
the failed `0.3%` comparison gate and does not select a new weak eigenvector.

## 1. Complete hypotheses and fixed carrier

Use exactly the published time-symmetric 600-cell dust sandwich already
certified in `verify_gravity_600cell_published_dust_control.py`, separately
for the even and odd five-stage schedules.  The logarithmic internal carrier
is

```text
R^35 = R^30 staircase diagonals + R^5 positive pole magnitudes,
```

and the varied final boundary is the thirty printed squared-edge
coordinates in logarithmic form.  The old boundary, unrounded published
`l0`, fixed external `tau`, and fixed dust mass are unchanged.

The following independently frozen data are the only matrix inputs:

- the recorded binary64 `35 x 35` internal Hessian and `35 x 30` mixed block
  from `7d5e9fc`;
- the action-derived five-dimensional Schur matrix and its empirical error
  from `dc927a5`;
- the exact analytic collective lapse tangent established by the complete
  action family in `790fc7f`.

No Standard-Model, cosmological or dimensional target enters this audit.

## 2. New boundary-row precision level

Repeat the complete-action mixed rectangles of protocol `25d9ee9` at the
three frozen steps

```text
h1 = 5e-4,  h2 = 2.5e-4,  h3 = 1.25e-4
```

using 100 decimal digits.  For each of the thirty boundary coordinates and
each schedule parity, evaluate all four signs on the exact collective path

```text
rho(t)=tau^2 exp(t),   q(t)=l0^2-rho(t),
f_j(v)=l0^2 exp(v).
```

After dividing by `24*norm(w_raw)`, define

```text
D1 = D(h1), D2 = D(h2), D3 = D(h3),
R12 = (4 D2-D1)/3,
R23 = (4 D3-D2)/3,
c6  = (16 R23-R12)/15,
epsilon6 = norm_2(c6-R23).
```

The sixth-order construction and all steps are fixed before the new action
values are inspected.  Every displaced simplex must remain Lorentzian and
off all branch boundaries.  Maximum imaginary contamination must be below
`1e-80`.

Classify the nonzero compatibility equation using

```text
floor6 = max(epsilon6,1e-35).
```

- `ALL_BOUNDARY_DIRECTIONS_COMPATIBLE` if `norm(c6)<=10*floor6`;
- `ONE_BOUNDARY_CONSTRAINT` if `norm(c6)>100*floor6`;
- `BOUNDARY_NORM_UNRESOLVED` otherwise.

## 3. Frozen uniform-scale diagnostic

The all-ones direction is a target suggested by the already inspected first
result, so this is not a blind discovery test.  It is a preregistered
cross-resolution check and retains **PATTERN provenance** even if numerically
confirmed.

Let

```text
e = ones(30)/sqrt(30),
r = (I-e e^T)c6.
```

Assign exactly one label:

- `UNIFORM_WITHIN_FROZEN_ERROR` if `norm(r)<=10*floor6`;
- `RESOLVED_NONUNIFORM` if `norm(r)>100*floor6`;
- `UNIFORMITY_UNRESOLVED` otherwise.

Report the mean, spread, cosine with `e`, `norm(r)`, and their ratios to the
error envelope.  Finite numerical agreement is not an analytic proof that
all components are identical.

If the first label holds, the response boundary basis is fixed to the exact
zero-sum complement of `ones(30)`.  If the second holds, use the deterministic
Householder complement of `c6`.  If uniformity or the norm is unresolved, do
not construct a response.

## 4. Deterministic precision correction of the internal Hessian

Write the recorded symmetrized Hessian in blocks

```text
H = [[A,G],[G^T,C]],  A:30x30, G:30x5.
```

Let `u` be the thirty diagonal components of the exact raw tangent and let
`p=ones(5)` be its pole components.  Correct only the collective component of
`G` by the unique minimum-Frobenius-norm update satisfying the exact top
null equation:

```text
rG = -A u-G p,
G* = G + rG p^T/(p^T p).
```

Let `P5=I-p p^T/(p^T p)` and let `S_hp` be the independently frozen
high-precision five-direction Schur matrix.  Remove only the collective row
and column required by the exact action family:

```text
S* = P5 [(S_hp+S_hp^T)/2] P5,
C* = S* + G*^T A^-1 G*,
H* = [[A,G*],[G*^T,C*]].
```

These formulas are fixed and contain no fitted Schur coefficient.  They are
the minimum-norm collective correction to `G`, the orthogonal projection of
the certified Schur form, and the exact inverse Schur reconstruction.
Require:

- `A` has rank 30 at relative thresholds `1e-7,1e-9,1e-11`;
- all four eigenvalues of `S*` on `p^perp` exceed
  `100*epsilon_S` in magnitude;
- `H* w_raw` is zero to normalized tolerance `1e-12`;
- `H*` has one and only one singular value below absolute `1e-9`;
- each relative correction `norm(G*-G)/max(1,norm(G))` and
  `norm(H*-H)/max(1,norm(H))` is below `1e-6`.

Failure of a correction-size or source gate is
`PRECISION_CORRECTION_UNRESOLVED`; it is not repaired by changing a
threshold.

## 5. Deterministic correction of the mixed block

The new action row controls the exact compatibility target.  With
`w=w_raw/norm(w_raw)` and the recorded mixed block `B`, set

```text
target_raw = norm(w_raw)*c6,
B* = B + w_raw (target_raw-w_raw^T B)/(w_raw^T w_raw).
```

This is the unique minimum-Frobenius-norm update that enforces `w^T B*=c6`.
Require the relative correction `norm(B*-B)/max(1,norm(B)) < 1e-6` and report
it.  The construction is rejected as unresolved if this gate fails.

## 6. Quotient and response

Construct the deterministic Householder complement `Q` of the exact
normalized tangent and set

```text
H_Q = Q^T H* Q.
```

The quotient is `QUOTIENT_REGULAR` only if all gates in sections 4 and 5 pass
and `H_Q` has 34 singular values above absolute `1e-9`.  This rank statement
must also agree with the block-factorization count `rank(A)+rank(S*|p^perp)
=34`.

For a resolved boundary basis `Y`, solve

```text
H_Q X = -Q^T B* Y,
delta_u = Q X.
```

Require relative quotient and corrected unprojected residuals below `1e-7`.
Report the response singular spectrum, ranks at `1e-7,1e-9,1e-11`, norm,
condition, its projection on the four relative pole directions, and the
corresponding residual against the uncorrected recorded matrices.  A large
uncorrected residual is reported, never optimized away.

## 7. Frozen outcomes and claim boundary

Per parity use exactly one of:

- `REGULAR_QUOTIENT_29_ZERO_SUM_RESPONSES`;
- `REGULAR_QUOTIENT_29_NONUNIFORM_RESPONSES`;
- `REGULAR_QUOTIENT_ALL_30_RESPONSES`;
- `PRECISION_CORRECTION_OR_BOUNDARY_UNRESOLVED`.

If both parities give the first outcome, the result is **DERIVED
COMPUTATIONAL at linear order** that the precision-corrected restricted
carrier propagates twenty-nine zero-sum boundary shape perturbations modulo
collective lapse, while the homogeneous scale is constrained.  Because the
uniform target was learned from the first run, its provenance remains
**PATTERN / cross-resolution confirmed**, not a blind selection.

Even the strongest outcome does not prove a nonlinear displaced solution,
identify a graviton, derive a clock, cover the unreduced 840-edge carrier,
establish a continuum constraint algebra, or propagate multiple slabs.  The
known Regge Hamiltonian/initial-value interpretation is **STRUCTURAL**; an
explicit nonlinear continuation is the next possible falsification test.
