# Primary result: exhaustive Coxeter blocks for internal schedule 0

Date: 2026-08-21

Status: **15/15 primary checks passed; acceptance pending an adversarially
independent replication.**

## Frozen provenance

- prior-art gate: `295e90d`;
- target-free protocol: `acfe795`;
- registered implementation before first execution: `69d2293`;
- primary artifact:

  ```text
  reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_blocks.json
  SHA-256 640f07a3d13ae3692761243cb62ace3ac2fd38f646d03b7750df05883d3f0267
  ```

Only schedule 0 was assembled.  Neither the full verifier suite nor the old
twelve-pair sparse factorization census was run.

## Construction controls

The verifier independently rebuilt the rank-coloured chamber graph and found
the Coxeter relations `(3,3,5)`.  The colour-preserving left action generated
by the frozen word `(0,1,2,3)` has order 30 on chambers, barycentric cells and
internal slab edges.  The forbidden right-product convention does not descend
to rank cells: it produces `48,960` consistency contradictions.

On the `19,680` internal edges the left action consists of exactly

```text
656 cycles x 30 edges/cycle.
```

Consequently all 30 cyclic Fourier sectors have dimension 656.  Including the
border coordinate in the invariant sector gives the exhaustive weighted
dimension

```text
657 + 14*(2*656) + 656 = 19,681.
```

The rebuilt sparse matrix has the exact frozen CSR digest.  Its group-average
distance bound is `4.14715e-16`, compared with the preregistered covariance
gate `2.38686e-10`.  A single `1e-4` diagonal corruption raises that bound to
`9.66667e-5` and is detected.  The analytic duration tangent has no
non-invariant Coxeter component within its multiplication envelope.

## Exhaustive primary spectrum

All 16 independent Hermitian blocks (`k=0,...,15`) were diagonalized in full;
the remaining sectors are their explicitly checked conjugates.  Trace and
Frobenius/Parseval controls agree exactly at the printed precision, and the
explicitly built `k=29` block agrees with the conjugate of `k=1` inside its
roundoff envelope.

Every one of the `19,681` bordered eigenvalues is separated from zero under
the frozen criterion.  The least favourable sector is `k=1`, with

```text
minimum |lambda| = 1.4556490539e-9,
zero-exclusion gate = 2.4245031150e-10,
margin ratio = 6.0039.
```

The post-hoc comparison with the old sparse solver reproduces its first eight
Ritz values with maximum difference `2.56e-16` under a `2.44e-10` gate.

The complete spectrum also explains the failed 32-vector diagnostic.  Using
the already frozen zero-exclusion scale to group its smallest values, the
first cluster multiplicities are

```text
4, 9, 16, 25, 36, ...
```

The old window of 32 necessarily ended three vectors into the fourth cluster,
and its two ARPACK runs returned different incomplete bases near that cut.
This is a diagnostic explanation, not a physical interpretation of the
multiplicities.

## Primary verdict and scope

Conditional on the numerical error model, the primary computation gives

```text
K0 is nonsingular,
therefore ker(C0) = span(n0).
```

Labels before replication:

- **DERIVED COMPUTATIONAL, PRIMARY ONLY:** exhaustive cyclic decomposition,
  covariance, and a positive zero-exclusion margin for schedule 0;
- **PATTERN:** the low cluster multiplicities `4,9,16,25,36,...`;
- **OPEN:** acceptance of the schedule-0 kernel claim until a mechanically
  different replication passes;
- **OPEN / NOT TESTED:** the other eleven schedule pairs, boundary
  propagation, graviton interpretation, a selected tick, `c`, `G`, or a
  Planck scale;
- **OPEN:** external novelty of this exact instance.  The finite-group block
  method itself is known prior art.

The next admissible step is an adversarial verifier with a different group
action construction and a different block assembly/diagonalization route.
It may not merely rerun this script.
