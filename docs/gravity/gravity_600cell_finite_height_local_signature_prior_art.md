# Prior-art gate: local stability of the finite-height branch signature

Date: 2026-08-22.

Status: completed before preregistering or implementing the local certificate.

## Exact proposed object

Fix the homogeneous tetrahedral-frustum 600-cell action, zero cosmological
constant, conserved global dust, the committed pre/post momentum convention,
positive proper height and positive endpoint scale. Restrict the initial data
to the already derived analytic incoming curve

```text
(m,pi)=(mu(v),p(v))
```

near the exact representative `v0=3/2`.

The proposed statement is only local: there exists an unspecified
`epsilon>0` such that every physical incoming `v` with
`|v-3/2|<epsilon` has the same ordered physical tree through slab four as the
accepted representative: one nontrivial first slab, two physical second
slabs, one branch dying at the next update and one branch entering the already
proved invariant half-strip `D` at slab four.

No explicit radius, basin measure, generic-state statement or physical
selection law is proposed. The exact diagonal zero-height solution `q=v`
persists identically and is recorded but excluded from the physical tree.

## Known mathematics and physics

1. The real-analytic implicit-function theorem gives a locally unique analytic
   continuation of a simple root. Local constancy of a complete real-root
   census additionally requires control of critical points, the origin and
   both tails; merely checking the roots already seen is insufficient. This
   mathematical mechanism is **KNOWN**, not a new theorem.

2. Fredrik Johansson, *Arb: efficient arbitrary-precision midpoint-radius
   interval arithmetic*, IEEE Transactions on Computers 66 (2017),
   [arXiv:1611.02831](https://arxiv.org/abs/1611.02831), provides the
   outward-rounded ball arithmetic used for the strict sign certificates. Arb
   is a verification method, not evidence for the physical model.

3. Dittrich and Höhn, *Canonical simplicial gravity*,
   [arXiv:1108.1974](https://arxiv.org/abs/1108.1974), establish the broader
   canonical setting in which later moves can impose a-posteriori constraints.
   They do not state the present scalar relation or local branch signature.

4. De Felice and Fabri,
   [arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093) and
   [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077), show that dust
   600-cell evolution can stop at a causal boundary. Their result is a control
   against inferring continuation from a few regular steps, not the local
   theorem proposed here.

## Boundary after the search

### KNOWN

- simple analytic roots and strict inequalities persist locally;
- a nondegenerate stationary tangency that is an exact identity can be
  separated from other roots locally;
- later consistency conditions are structurally admissible in canonical
  simplicial gravity.

### CONTROL

- the complete accepted tree at `v=3/2`;
- the finite discovery census showing the same terminal multiset at adjacent
  diagnostic nodes;
- the invariant theorem on `D`.

The discovery census may motivate this gate but is not its proof.

### OPEN

- whether every load-bearing root at `v=3/2` is simple after quotienting the
  exact diagonal tangency;
- whether every physical, nonphysical, tail and invariant-entry inequality is
  strict under rigorous outward-rounded evaluation;
- the existence of the claimed open neighbourhood;
- any explicit maximal radius or global basin;
- the physical status of complete extendibility;
- external novelty of this exact application.

## Framing warning

A positive local theorem removes the literal single-point objection to
`v=3/2`. It does not make that state derived, show that the neighbourhood is
large, classify the full incoming domain, or turn future extendibility into a
local law. The discovery margins already show that the `D` entry is close to a
boundary, so no genericity language is admissible.

