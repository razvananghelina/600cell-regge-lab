# Adversarial protocol: finite-height internal carrier rank

Date: 2026-08-22.

Primary implementation commit: `1472213`.

Primary result commit: `2b9c124`.

Status: **FROZEN AFTER THE PRIMARY RESULT AND BEFORE THE FIRST FULL-REAL-SPACE
REPLICATION RANK, SINGULAR VALUE OR NULL VECTOR IS EVALUATED.**

## 1. Claim under attack

The primary 180-digit orbit-sector calculation reports, for both staircase
parities, that the finite-height internal-equation map

```text
R_p = H_p[internal,active] G_p : R^240 -> R^840
```

has rank `239` and a common one-dimensional right kernel.  It also reports
that the diagonal-only `720 x 240` map has rank `119` and nullity `121`.

This is a result-aware adversarial replication, not a second target-blind
search.  Its purpose is to falsify the claimed ranks, common kernel, row split
or orbit-sector implementation by a mechanically different construction.

## 2. Frozen provenance

The verifier must reject any hash mismatch:

```text
reproducible/gravity_600cell_finite_height_internal_carrier_rank.json
  513fdea33f6b868efa6d6f2b2526bade7ce615ea949f955588916a8d0baee0c8

reproducible/gravity_600cell_finite_height_internal_carrier_rank_matrices.npz
  97f5b8318be2b3ccf843db87e678ac1ac6ce402db262023c6bbc63a7b647321b

reproducible/verify_gravity_600cell_finite_height_internal_carrier_rank.py
  fff2c70dc3685562b4e5f1e7646886c828a5df1aa7cbed792cef8b19afdf8c62

reproducible/verify_gravity_600cell_finite_height_carrier_quadratic.py
  bbe7112270a7f2bcb2d443fab45ca450598e7234250bd335b14a4ed7869443a5

reproducible/verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py
  834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf

reproducible/gravity_600cell_finite_height_carrier_quadratic.json
  0ec142bfc68d04498992a6cdba7437933560b860244573d187cb6e018ece78f9
```

The primary rank JSON and matrices remain unread until the replication has
fixed its ranks, nullities and kernel projectors in memory.  Their hashes may
be checked first, because a cryptographic hash exposes no scientific value.

## 3. Mechanically different reconstruction

Do not use the primary rank verifier's identity-row orbit kernel, seven
minimal representation sectors, 180-digit SVD or determinant minors.

For each parity:

1. reconstruct all `2,400` four-simplices and all `2,280` physical edges;
2. assemble each complete `2,280 x 2,280` Hessian in physical edge order from
   every local triangle and simplex contribution;
3. construct the `1,560 x 240` carrier as a full sparse physical matrix;
4. multiply directly in real space;

```text
T_p = H_p[active,active] G_p,
R_p = T_p[internal rows].
```

The Hessian assembly uses 120-digit local arithmetic followed by explicit
binary64 storage.  It uses the already frozen derivative steps

```text
1e-20, 1e-15, 3e-20, 3e-15,
```

which differ from the primary rank calculation's 180-digit
`1e-25,...,1.25e-26` hierarchy.  Classify the global matrices directly; no
group Fourier transform is allowed.

Identify pole and oriented-diagonal equation rows from physical edge labels.
Do not assume their orbit positions.

## 4. Frozen scaling and error interval

For each parity and each of the diagonal/full maps, let the midpoint response
be the average of the `1e-20` and `3e-20` matrices.  Derive separate positive
scale/strut column factors from that midpoint exactly as in the primary
protocol and apply them to all four matrices.

Set

```text
N = max(1, ||R_mid||_2),
e_step = max_level ||R_level-R_mid||_2 / N.
```

If `b` is the maximum certified binary64 Hessian entry error and `G_j` is a
carrier column, use the conservative response Frobenius bound

```text
e_round_abs = b * sqrt(number_of_rows) * sqrt(sum_j ||G_j||_1^2),
e_round = e_round_abs/N + 500 eps_machine max(matrix dimensions),
e_total = e_step + e_round.
```

At every derivative level, classify singular values by

```text
ZERO     sigma/N <= 10 e_total,
NONZERO  sigma/N >  100 e_total,
OPEN     otherwise.
```

All four levels must have no `OPEN` singular value and must give the same
nullity.  Binary64 is legitimate here only because its explicit roundoff
allowance enters the classifier.

As an independent rank diagnostic, pivoted QR of each midpoint must put every
accepted diagonal pivot above `100 e_total N` and every rejected pivot below
`10 e_total N`.  QR may corroborate or make the result open; it cannot rescue
an SVD failure.

## 5. Global kernel comparison

For a resolved positive full-map nullity, form the right-kernel projector
directly from each full physical SVD.  Its within-parity uncertainty is the
maximum projector difference over the four derivative levels plus

```text
500 eps_machine max(matrix dimensions).
```

The parity projectors agree below ten times the combined uncertainty, differ
above one hundred times it, and are open between those gates.

Only after this classification is fixed may the verifier load the primary
artifact.  Replication requires:

```text
full nullities       even=odd=1,
diagonal nullities   even=odd=121,
primary/replication kernel projectors agree under the same 10/100 rule.
```

## 6. Direct nonlinear-gradient secant

This control bypasses Hessian assembly.  Let `x0` be the replicated full-map
kernel vector and `x1` the right singular vector for the smallest classified
nonzero singular value.  In the complete physical edge coordinates evaluate
the action gradient directly at

```text
delta = +/- h G x,     h in {1e-4, 5e-5, 2.5e-5}.
```

Every displaced simplex must retain Lorentzian inertia `(3,1)`.  The centered
secants must converge with the expected second-order hierarchy, allowing

```text
||D_(h/2)-D_(h/4)|| <= 0.4 ||D_h-D_(h/2)|| + 1e-20.
```

At the finest step the secant for `x0` must be below `1e-6` of that for `x1`,
while the `x1` secant must agree with the assembled linear response to
relative error below `1e-5`.  Failure refutes or leaves open the claimed
linear null direction even if the matrix SVD appears favorable.

## 7. Hostile controls

1. Under the actual `e_total`, a zero matrix must have full nullity and an
   embedded identity must have zero nullity.
2. Adding `+0.1` to the deterministic first carrier coefficient must change
   the complete response by more than `1e-12` relatively.
3. Adding `+0.1` to a deterministic internal Hessian diagonal entry must
   change the complete response by more than `1e-12` relatively.
4. Reversing the physical diagonal/pole row labels must fail the registered
   row-count control.

No failed hostile control may be removed after execution.

## 8. Outcome hierarchy

### `FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_ADVERSARIALLY_REPLICATED`

Use only if every control passes, both global rank censuses reproduce the
primary nullities, parity and primary kernel projectors agree, and the direct
nonlinear-gradient secants validate the null direction.

### `FINITE_HEIGHT_INTERNAL_CARRIER_PRIMARY_REFUTED`

Use if every prerequisite/control is valid and a resolved real-space rank,
projector or direct secant contradicts the primary result.

### `FINITE_HEIGHT_INTERNAL_CARRIER_ADVERSARIAL_OPEN`

Use if a singular value, QR pivot, projector or direct secant remains within
its open band.

### `FINITE_HEIGHT_INTERNAL_CARRIER_ADVERSARIAL_CONTROL_FAILED`

Use for provenance, geometry, branch, stationarity, carrier, row-label or
hostile-control failure.

## 9. Interpretation firewall

Even a successful replication selects only a one-dimensional homogeneous
first-order direction inside this exact carrier at this exact slab.  Whether
it is physical evolution, lapse/gauge, a constraint direction or an artifact
of the homogeneous dust ansatz remains **OPEN**.

The separate invariant-region result gives a unique successor at every later
finite step.  It is not a theorem of convergence, infinite proper duration or
geodesic completeness.

No outcome here derives gravitons, a wave equation, a continuum limit, a
physical tick, `c`, `G`, Planck units, particle masses or Standard-Model
physics.  Run only the adversarial verifier, never the full suite.
