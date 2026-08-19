# Result: the static gradient kernel is exactly the prism-shift carrier

Date: 2026-08-19

## Verdict

**DERIVED EXACT RECONCILIATION.**  The recently derived 119-dimensional
static variable-connection kernel and the previously certified equal-scale
prism-shift family are one embedded global carrier written in two local
coordinate systems.

For every tetrahedron, the Cartesian translation `s_T` and the old covector
evaluations `a_T` obey

```text
a_T=Q_T s_T,
Q_T=D_T^T,
det(Q_T)=-16.
```

Globally,

```text
Q G=B,
C_new=C_old Q
```

hold entrywise over the rationals.  The second identity was checked on all
2,400 ordered face-trace equations.  Consequently

```text
Q ker(C_new)=ker(C_old),
dim=119.
```

This is not another physical sector and not a new discovery of 119 modes.  It
is the exact removal of duplicate project labels.

## Provenance

| stage | commit |
|---|---|
| prior-art/framing gate | `fe9fc18` |
| original protocol | `e2c633e` |
| registered verifier | `d2498e0` |
| preserved inert-control failure | `a442429` |
| disclosed protocol correction | `f5f315a` |
| corrected verifier | `2b36d94` |
| passing artifact | `9cd8938` |

The passing verifier source has SHA-256

```text
daad236db940430536909bfc785e1053292dc7dbee7e6f7c3d1f923581df296b.
```

The passing artifact has SHA-256

```text
2de7c4594ae1dc458fea8abd23461f06b25a038cb0cf7ee290d61d7aed83bb27.
```

Only the targeted verifier was run.  It returned `11/11`; the full suite was
not run.

## Transparent first failure

The first execution returned `10/11` and
`RECONCILIATION_CONTROL_FAILED`.  Every target identity and rank had passed,
but the negative control replaced the transport by the identity on a face
where the legitimate transport was already identity on both tested
tangents.  That artifact is preserved with SHA-256

```text
96126215b507f6fac8e054d1634401400a9f8925e15223194b4c9dcabf0490ff.
```

Before rerunning, the protocol was changed to the lexicographically first
face whose derived tangential transport is non-identity.  This selection uses
only local geometry and did not alter the target equalities, ranks or outcome
tree.  The corrected control then failed as intended.

## Exact audit ledger

The successful run found, at both primes `1000003` and `1000033`,

```text
rank(B)=rank(G)=119,
rank(C_old)=rank(C_new)=1681.
```

It additionally required:

- exact closed 600-cell incidence `f=(120,720,1200,600)`;
- invertibility of every local coordinate block;
- entrywise potential and face intertwining;
- failure after a one-axis frame corruption;
- failure after removing a genuinely nontrivial target transport;
- persistence under the odd canonical relabelling `(0 1)`.

Thus equality is not inferred from dimension, common notation or a fitted
basis.

## Consequences

All established prism-shift results apply to the new static-gradient name:

1. **DERIVED:** on the equal-scale branch, its Regge Hessian is

   ```text
   [2*pi-5*acos(1/3)]/(L*sqrt(rho)) * Delta_600.
   ```

2. **DERIVED:** for unequal homogeneous scale with common struts, the local
   branch equation `(q-1)s=0` removes every nonconstant potential.
3. **DERIVED COMPUTATIONAL:** on the complete accepted dust slab, the
   canonical Schur equations eliminate all 119 relative directions.
4. **STRUCTURAL:** the most defensible reading is an auxiliary or
   pseudo-constraint-like longitudinal variable, not a propagating scalar.

Therefore the new variable-face closure calculation is an independent
geometric derivation and a correction to frozen-connection rigidity, but its
static 119-dimensional kernel does not reopen the already closed free-scalar
route.

## Physical ledger

- duplicate 119-dimensional carriers: **CLOSED / IDENTICAL**;
- free scalar propagation from this carrier: **CLOSED**;
- time or physical tick from this carrier: **NOT DERIVED**;
- limiting speed from this carrier: **NOT DERIVED**;
- particle inertia, mass, `G` or Planck units: **NOT DERIVED**;
- external novelty of the abstract gradient theorem: **KNOWN PRIOR ART**;
- external novelty of the project-specific exact 600-cell reconciliation:
  **OPEN**, with no physics novelty claimed.

## Correct continuation

Do not differentiate or diagonalize this 119-dimensional carrier again.  The
coarse 600-cell has already refuted exact gauge symmetry and every imposed
spectral/York proxy tested so far.  The remaining load-bearing question is
whether action-derived weak/pseudo-constraint directions approach geometric
vertex-displacement directions under the already selected `H4`-equivariant
refinement family.

That refinement mission must report singular-value and subspace scaling
before imposing any quotient.  Only a recovering constraint carrier can
legitimately expose a tensor sector whose temporal/spatial coefficients could
later define an effective propagation speed.

