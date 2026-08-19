# Preserved first failure: the declared face control was geometrically inert

Date: 2026-08-19

The first targeted execution of
`verify_gravity_600cell_static_gradient_prism_shift_reconciliation.py`
returned

```text
10/11 checks passed
RECONCILIATION_CONTROL_FAILED.
```

The failed artifact is preserved with SHA-256

```text
96126215b507f6fac8e054d1634401400a9f8925e15223194b4c9dcabf0490ff.
```

The executed verifier source has SHA-256

```text
ccb55d7efb98a9753abc9c32e7b705dc54868f78f37e7d2036facd39dab9b6ea.
```

All scientific identities and ranks had already passed:

```text
Q G = B                                      exact,
C_new = C_old Q                              exact on 2400 rows,
rank(B)=rank(G)=119                           at both primes,
rank(C_old)=rank(C_new)=1681                  at both primes,
odd relabelling                               passed.
```

The sole failure was the preregistered negative control that replaced the
target transport by the identity on the lexicographically first face.  On
that face, the source and target canonical labels happen to agree on the two
shared tangential directions.  Consequently its legitimate transport acts
as the identity on precisely the two covectors tested by the face equations.
The replacement is geometrically inert and has no falsification power.

This is a protocol/control failure, not evidence for or against the target
intertwiner.  Before another execution, the protocol must select the
lexicographically first face whose derived target transport is non-identity
on at least one of its two tangents.  That selection uses only the already
constructed geometry and does not inspect the global equality or a desired
rank.  The failed artifact will not be overwritten in history.

