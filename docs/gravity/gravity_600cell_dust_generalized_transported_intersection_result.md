# Result: the transported generalized phase intersection is exactly zero

Date: 2026-08-18

## Status ledger

| Statement | Status |
|---|---|
| The generalized rank-15 configuration fibers rotate between the two middle slices | **DERIVED COMPUTATIONAL** (conditional on the frozen four schedules) |
| Their full rank-30 cotangent lifts are transported by the action tangent | **DERIVED COMPUTATIONAL REFUTATION** |
| A nonzero subset of the old lift is transported into the shifted lift | **DERIVED COMPUTATIONAL REFUTATION** |
| The transported intersection is zero in all 16 parity/sector/schedule cells | **DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED** |
| The generalized-pencil phase-fiber route supplies propagating modes | **CLOSED** |
| The full Regge tangent has no physical perturbations | **NOT TESTED** |
| A direct Hessian/constraint-selected Regge mode sector exists | **OPEN** |
| Dispersion, inertia, mass or limiting speed | **OPEN / NOT COMPUTED** |
| External novelty | **OPEN** |

## 1. Exact question

With the already fixed phase fibers and canonical second-slab tangent,

```text
F0 = E_old direct-sum E_old*,
F1 = E_shifted direct-sum E_shifted*,
K0 = F0 intersection T_2^-1(F1),
```

the mission computed `dim(K0)` without selecting a graph or inspecting a
desired rank.

For

```text
R=(I-Q1) T_2 Q0,
```

the exact right factor has rank 30, so `rank(R)<=30` and

```text
dim(K0)=30-rank(R).
```

Resolving 30 nonzero singular values is therefore sufficient to certify
`K0={0}`; the other 30 singular values are ambient zeros from the right
projector and do not need an exact numerical-zero claim.

## 2. Exact execution chronology

- prior-art/framing commit: `656b1d7`;
- preregistered protocol: `1e6d8ce`;
- initially registered verifier: `6e785d6`;
- execution-history/serialization addendum: `999f8fa`;
- repaired verifier: `895e224`;
- exact artifact commit: `775440c`.

The first execution terminated `6/7` with
`TRANSPORTED_INTERSECTION_CONTROL_FAILED`.  The recomputed SVD norm differed
from the committed phase norm by `2.56e-30...3.46e-29`, while the physical
operator error was approximately `2e-47`.  The cause was not the matrix: the
committed norm had been serialized to only 30 significant digits, so the
protocol asked a text string to support an impossible 47-digit replay check.

The failure was recorded before rerunning.  The repair adds

```text
1e-29 max(1,|committed norm|)
```

only to the norm-string overlap control.  It is not added to the singular-value
errors and changes no rank label or outcome criterion.  The `480/480` pattern
from the failed run was treated as diagnostic only.

Two subsequent executions passed `7/7` and produced byte-identical SHA-256

```text
207cbe61bfaaf2b13d62cc3dbbb2ed5ea4931b7aab13cd47a8dd2802410c55c0.
```

No full suite was run.

## 3. Exact rank result

Across all 16 parity/sector/schedule cells, the complete 960-value census was

```text
480  SINGULAR_NONZERO_RESOLVED
480  SINGULAR_ZERO_CONSISTENT
```

Every cell contained exactly 30 resolved nonzero singular values and 30
zero-consistent ambient values.  The smallest resolved value was

```text
4.188398429739906...e-7,
```

at a minimum separation of

```text
1.820740239880615...e40
```

complete error units.  Thus every exact residual has rank at least 30; the
right-factor theorem gives rank at most 30.  Consequently

```text
rank(R)=30,
dim(K0)=0
```

in all 16 cells.

Accepted exact outcome:

```text
TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL
```

## 4. Independent adversarial replication

Project rule 4 required a separate implementation:

- independence gate: `0d3a46c`;
- adversarial protocol: `08ff16c`;
- registered verifier: `e7daa19`;
- adversarial artifact: `c30135a`.

This audit did not use the 100-digit projectors, exact tangent reconstruction
or exact singular values.  It rebuilt the bases through the earlier float64
null-space/Cholesky path, loaded the earlier tangent archive and formed the
independent square leakage matrix

```text
L=W1_perp^H T_2 W0.
```

Both runs passed `7/7` and produced byte-identical SHA-256

```text
aa6d4cf2ec248ce823d0b22467c1c0337e2341e333fd593cc56e9fae221a5fc0.
```

All 480 singular values were `AUDIT_SINGULAR_NONZERO_RESOLVED`; every one of
the 16 matrices had numerical rank 30.  The smallest singular value range was

```text
4.1883846e-7 ... 4.1884669e-7,
```

or `811.1...811.2` times the preregistered float64 floor.  Synthetic maps with
full intersection and zero intersection were classified correctly in every
cell, and basis reversal/rephasing preserved the spectra.

Accepted adversarial outcome:

```text
ADVERSARIAL_INTERSECTION_ZERO_CORROBORATED
```

The float64 audit is **STRUCTURAL INDEPENDENT CORROBORATION**, not a replacement
for the exact certificate.

## 5. Kill boundary and physical meaning

**DERIVED COMPUTATIONAL NEGATIVE.** No nonzero phase datum in the proposed old
generalized fiber is carried by the canonical tick into the corresponding new
fiber.  Therefore neither the full lift nor a smaller graph contained inside
it can serve as the propagating carrier.  The generalized-pencil phase-fiber
route is closed.

The upstream time-dependent configuration fiber remains a valid algebraic
property of the kinetic--stiffness pencil, but it has no demonstrated status as
a dynamically transported matter, graviton or inertia bundle.

This does not say that the complete 60-dimensional Regge phase space has no
dynamics.  It only rejects one imposed spectral subspace.  Higher-rank or
different hand-selected bundles would reopen the fitting freedom and require
an independent action-selection theorem before calculation.

## 6. Post-result literature and next direction

Bahr and Dittrich show that on curved Regge backgrounds exact vertex gauge
symmetries are generically broken and become pseudo-constraints; small Hessian
eigenvalues diagnose approximate rather than exact gauge symmetry
(<https://doi.org/10.1088/0264-9381/26/22/225011>,
<https://arxiv.org/abs/0905.1670>).

Hoehn's canonical linearized Regge analysis identifies propagating curvature
degrees of freedom and counts lattice gravitons for Pachner evolution around a
flat background (<https://arxiv.org/abs/1411.5672>).  Our background is curved
and contains dust, so those flat-background constraints and counts cannot be
imported as a result here.  Dittrich and Hoehn's covariant-to-canonical
construction likewise derives the constraints from the action rather than an
external spectral fiber (<https://arxiv.org/abs/0912.1817>).

The next honest programme is therefore:

1. audit what the repository has already established about the complete action
   Hessian and pre/post Legendre ranks on the curved dust solution;
2. identify exact null and quantitatively weak Hessian directions without a
   desired count;
3. distinguish gauge/pseudo-constraint directions from physical curvature
   directions using the action's pre/post maps;
4. transport only the quotient or relation selected by that analysis across
   both available ticks;
5. only after that ask for wave spectra or dispersion.

No claim about a graviton, inertia, mass or `c` survives merely from the closed
generalized-fiber route.

