# Global audit of the chamber one-form minimum

Date: 2026-07-28

## Question

Is

`dim_C Omega_D^1 = 126`

the global minimum over every orientable, first-order, integral
Poincare-dual partition algebra on the 120 oriented chambers?

## Verdict

**REFUTED.  The value 126 is not the global minimum.**

The later exact `C5`-symmetric witness printed by
`reproducible/verify_chamber_symmetry_sat.py` satisfies first order,
orientability and integral Poincare duality, and has

`dim_C Omega_D^1 = 110`.

Its intersection matrix has dimension and rank 30 and determinant 1.  Thus
it is a concrete counterexample to the proposed lower bound 126.  The global
minimum itself remains **OPEN**, with the current certified upper bound 110.

The value 126 is:

- **DERIVED:** attained by an explicit registered integral witness;
- **DERIVED:** the unique minimum in the completely enumerated Hamming
  radius-two neighborhood of the original integral witness;
- **DERIVED:** a local upper bound, superseded globally by the 110 witness;
- **NOT DERIVED:** a global lower bound.

The equality with the number of `E7` roots remains **PATTERN**.

## Exact combinatorial formulation

The 180 chamber edges form 90 pairs `{e,J(e)}`.  For a partition algebra:

- first order requires at least one edge of every pair to be internal;
- orientability forbids both paired edges from becoming internal after
  transitive closure.

Thus every orientable first-order partition is represented by an exact
one-of-two contraction choice on the 90 pairs, modulo choices producing the
same component partition.

For a completed partition, `dim_C Omega_D^1` is the number of directed
nonloop edges in the quotient chamber graph.  Integral Poincare duality is
checked by the exact determinant of its intersection matrix.

## Attempt 1: exact branch-and-bound

`reproducible/verify_chamber_global_minimum.py` uses union-find propagation.
A branch is rejected immediately when both edges in a `J` pair become
internal.

The largest run used limits of 600 seconds and 10,000,000 nodes.  It
reported:

- nodes: `10,000,000`;
- conflicts: `2,779,520`;
- complete leaves: `2,220,453`;
- distinct completed partitions: `2,220,453`;
- status: `INCOMPLETE`.

The tree was not exhausted.  Therefore this is not a global certificate.

## Attempt 2: global SMT

`reproducible/verify_chamber_global_smt.py` encodes an arbitrary set
partition by a restricted-growth string.  It imposes:

- first order on all 180 chamber edges;
- orientability on all `60*60=3600` opposite-sign chamber pairs;
- at most 62 distinct unordered quotient edges, equivalent to
  `dim Omega <=124`.

If this formula were `UNSAT`, the registered 126 witness would prove the
global minimum.  Z3 5.0 returned:

`unknown: timeout`

after 600 seconds.  No UNSAT certificate was produced.

## Honest boundary

Superseded licensed statement (before the `C5` search):

> There exists an integral finite KO6 chamber triple with 126 one-form
> directions, and 126 is the unique minimum in its fully enumerated
> radius-two neighborhood.

Unlicensed statement:

> 126 is the global minimum over all legal chamber partitions.

Current licensed statement:

> 126 is a certified local minimum in the stated radius-two neighborhood,
> but not a global minimum.  An exact legal integral witness attains 110.

The assertion that 126 might be the global minimum is now **REFUTED**; the
exact global minimum remains **OPEN**.

## Next rigorous routes

- a stronger SAT encoding including fixed block count and mod-2
  nondegeneracy of the intersection form;
- isomorph-free generation of component partitions rather than binary edge
  choices;
- a mathematical lower bound relating unimodular antisymmetric Poincare
  pairings to quotient-graph edge count;
- an external proof-producing SAT/CP solver with decomposition by chamber
  edge type.
