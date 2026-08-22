# First direct-bisection adversarial result: OPEN

Date: 2026-08-22.

Verifier commit: `bc30c63`.

Status: **ADVERSARIAL OPEN; PRIMARY NOT YET CONSOLIDATED**.

## Frozen result

The first execution returned

```text
RESULT: 4/11 checks passed
OUTCOME: LOCAL_SIGNATURE_ADVERSARIAL_OPEN
```

Artifact SHA-256:

```text
139dcee2e9ee021c131aae1090433fe16bd70c9f2b10ec52d32b0c5ebd7748a7.
```

The protocol provenance, radical factorization, positive `sqrt(2)` control,
negative polynomial control and exact diagonal nondegeneracy passed. The
gravity recursion stopped at the initial state.

## Exact failure mode

The rational stationary bracket `[5,6]` has strict opposite endpoint signs for
`p(q)-pi`, but direct evaluation of `p'` on the entire unit interval returned

```text
[+/- 3.34e+2],
```

which contains zero. Likewise, the root bracket `[9,10]` has strict opposite
endpoint signs for `E`, but direct evaluation of `E_q=p(q)-pi` over the whole
unit interval returned

```text
[+/- 2.00e+2].
```

The wide derivative balls prevent uniqueness certification. Because the root
intervals were not narrowed, the physical `h` gate also remained unresolved,
so the implementation correctly produced zero accepted physical children and
did not construct later states.

This is interval dependency in the raw formulas, not a contradictory root or
sign. Endpoint signs prove existence in both brackets but not uniqueness. The
frozen protocol required uniqueness on the complete bracket and therefore had
to return `OPEN`.

## Evidential consequence

The primary `10/10` result remains a primary certificate only. It must not be
described as adversarially corroborated while this disagreement is unresolved.

A permissible resolution must preserve all integer brackets and avoid
interval Newton. It may separately preregister a fixed dyadic cover that
certifies derivative signs cell by cell, then use those signs in direct
bisection. Changing brackets after seeing their roots or importing primary
root balls is forbidden.

