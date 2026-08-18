# Preregistration: can fibre incidence split the golden `D5` doublets?

Date: 2026-08-11

## Disclosed candidate and scope

This is not blind discovery.  During the post-verdict audit of
`hopf_six_crossed_real_galois_verdict.md`, the following candidate was found
before computation.

Fix one of the six derived Hopf fibrations.  Its effective `A5` stabilizer is
the already certified group `D5`.  Let `r` be the order-five stabilizer
element induced by either of the two 600-cell neighbours of a vertex along
its decagonal fibre.  Reversing the edge replaces `r` by `r^-1`, so the
self-adjoint group-algebra element

```text
u_edge = r+r^-1
```

does not require an orientation.  The candidate observation is that its two
scalars on the golden `D5` doublets should be

```text
phi-1 > 0,  -phi < 0.
```

If fibre incidence selects the inverse pair `{r,r^-1}` globally, rather than
the other inverse pair `{r^2,r^-2}`, then real spectral order would distinguish
the two `M12(R)` crossed-product blocks.  This would refute the stronger
wording that the current geometry has no mechanism at all for distinguishing
the golden pair.  It would not construct a Krajewski support, intersection
form or finite spectral triple.

No Hessian, Standard-Model module, mass, coupling or later target will be
inspected in this audit.

## Frozen canonicity tests

The verifier will independently rebuild the 600-cell, its six left-coset
Hopf fibrations, the 600-cell edge relation and the conjugation action of
`A5` on the six fibration labels.

For each of all six fibrations:

1. recover its identity fibre, hence its order-ten binary subgroup, without
   choosing a generator;
2. enumerate every one of the 120 undirected 600-cell edges lying inside its
   twelve decagonal fibres;
3. for either orientation of every such edge, compute the relative binary
   element and its effective conjugation action;
4. require the resulting actions to be exactly one inverse pair in the
   order-five rotation subgroup of the `D5` stabilizer;
5. require the other order-five inverse pair to occur at fibre graph distance
   two and never at distance one.

The candidate is **not canonical** if the selected inverse pair depends on a
base vertex, a fibre, an edge orientation or an enumeration convention.

## Frozen representation-theory tests

Using only the incidence-selected inverse pair, the verifier will form the
integer symmetric matrix of `u_edge` in the left regular representation of
`D5`.  It must establish exactly:

```text
charpoly(u_edge)=(lambda-2)^2 (lambda^2+lambda-1)^4,
spectrum: 2 [multiplicity 2], phi-1 [4], -phi [4].
```

It will construct the real spectral projectors and check:

- the negative projector is central, self-adjoint, idempotent and has rank
  four;
- it selects exactly one real `M2(R)` group-algebra block, hence exactly one
  `M12(R)` block after six-point Morita amplification;
- golden conjugation exchanges it with the positive-doublet projector;
- their sum is rational, but neither individual projector descends to
  `Q[D5]`;
- the alternative chord element `r^2+r^-2` exchanges the two projectors,
  making the edge-versus-chord incidence distinction load-bearing.

The unique reflection class also gives the central self-adjoint element

```text
v_ref = sum_(s reflection) s.
```

The joint central spectrum of `(u_edge,v_ref)` will be tested against

```text
trivial:       (2,  5)
reflection-sign: (2, -5)
first doublet:   (phi-1, 0)
second doublet:  (-phi,  0),
```

up to naming the two doublets by their displayed eigenvalues.  This tests
whether the geometry labels all four real Wedderburn nodes without arbitrary
Schur coefficients.

## Frozen interpretation boundary

Acceptance has two levels.

- **DERIVED NODE SEPARATION:** all incidence and exact spectral tests pass.
  Then ordered-real functional calculus canonically separates the conjugate
  blocks, and the earlier statement that no current geometric separator
  exists must be corrected.
- **NO FINITE-TRIPLE CLAIM:** even on acceptance, no nondegenerate KO6
  intersection matrix is selected.  The old Galois-compatible no-go remains
  valid under rational descent; what fails is treating that extra descent
  requirement as mandatory after the real geometry has supplied an ordered
  spectral separator.

If incidence does not distinguish the inverse classes uniformly, or if the
regular spectrum/projectors fail any exact test, record a **DERIVED
NEGATIVE** and retain the previous no-separator verdict.

Only the new targeted verifier will be run.  The full suite is deliberately
excluded from this continuation at the user's request.
