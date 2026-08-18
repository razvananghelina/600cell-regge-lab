# Preregistration: additive overlay-transfer rank audit

Date: 2026-08-17

Prior-art commit: `c77a87c`.

Status: frozen before constructing any row of the transfer matrix.

## 1. Frozen input

Read only

```text
reproducible/gravity_600cell_universal_staircase_overlay.json
SHA-256 0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc.
```

Require the certified outcome, `148` feasible sign words, `24` staircase
orders, `14` `S4 x C2` fine chamber orbits and `12/12` passing controls.  Do
not read a gravity action, experimental quantity, desired rank or continuum
target.

## 2. Reconstructed aggregation matrix

Decode each sign word in mask order `1,...,14`.  Recompute, without importing
the generating verifier, the unique staircase position `k` for every chamber
and every vertex order.  Construct the integer matrix

```text
R in {0,1}^{96 x 148},
R[(o,k),C] = 1 iff C is assigned to k for order o.
```

Require for every order:

- every column appears in exactly one of its four rows;
- the row sums are exactly `(19,55,55,19)`;
- the sum of the four rows is the all-ones vector.

Store a SHA-256 digest of the row-major byte matrix.

## 3. Exact full-space rank

Compute `rank_Q(R)` by exact rational row reduction and set

```text
nullity_Q(R) = 148-rank_Q(R).
```

Cross-check the rank using independent modular Gaussian elimination for the
primes

```text
1000003, 1000033, 1000037.
```

All modular ranks must equal the rational rank.  Require the analytic bounds

```text
rank_Q(R) <= 73,
nullity_Q(R) >= 75.
```

Construct an exact rational nullspace basis and require `R N=0`, the basis has
the printed nullity, and its columns are independent.

## 4. Exact symmetry reduction

Reconstruct the `S4 x C2` action on sign words, with

```text
vertex permutation: sign(A) -> sign(g(A)),
time reflection:    sign(A) -> -sign(A^c).
```

Re-enumerate the fine chamber orbits and require exactly the frozen size
distribution.  Form the `148 x 14` orbit-indicator matrix `C` and the invariant
aggregation matrix

```text
R_inv = R C  in Z^{96 x 14}.
```

On coarse labels use

```text
g:(o,k) -> (g o,k),
time reflection:(o,k) -> (reverse(o),3-k).
```

Require exactly two coarse orbits, of sizes 48 each, and require rows of
`R_inv` to be identical within each coarse orbit.  Compute the exact rational
rank and nullity of `R_inv`, cross-check them modulo the same three primes and
require

```text
rank(R_inv) <= 2,
nullity(R_inv) >= 12.
```

## 5. Constructive positive nonuniqueness witness

Choose the first exact null vector of `R_inv` in the deterministic SymPy
nullspace ordering.  Let `v_j` be its 14 orbit coordinates and lift it to 148
chambers by orbit membership.  Set

```text
epsilon = 1/(2 max_j |v_j|),
x_plus  = 1 + epsilon v,
x_minus = 1 - epsilon v.
```

Require exactly:

- `v != 0` and `R v=0`;
- both `x_plus` and `x_minus` are strictly positive;
- both vectors are `S4 x C2` invariant;
- `x_plus != x_minus`;
- `R x_plus = R x_minus`;
- the total fine weights agree.

Store the 14 rational orbit coordinates, `epsilon`, the common 96 coarse
totals and the minimum fine weight.  This is a reproducible infinite-family
witness, since every rational parameter in a sufficiently small interval
around zero gives another positive lift.

## 6. Mechanical outcome

- All controls pass and both the full and invariant nullities are positive,
  with the constructive witness passing:
  `POSITIVE_INVARIANT_ADDITIVE_TRANSFER_NONUNIQUE`.
- Full-space nullity is positive but invariant nullity is zero:
  `SYMMETRY_REMOVES_ADDITIVE_TRANSFER_KERNEL`.
- Full-space nullity is zero:
  `ADDITIVE_TRANSFER_INJECTIVE`.
- Any source, reconstruction, exact-rank, modular-rank, group-action or
  witness check fails:
  `ADDITIVE_TRANSFER_CONTROL_FAILED`.

No rank other than the analytic upper bounds is preregistered.

## 7. Scope and kill boundary

The first outcome kills uniqueness only under the fully stated additive
top-cell hypothesis, even after positivity and finite symmetry.  It proves
that the common carrier alone does not select a fine scalar action density.
It does not kill hinge-local, nonlocal, metric-volume, Galerkin or dynamically
perfect actions.  Any continuation must add one of those principles explicitly
and test whether it is derived rather than fitted.

Run only the registered verifier for this mission.  Do not run the full suite.

