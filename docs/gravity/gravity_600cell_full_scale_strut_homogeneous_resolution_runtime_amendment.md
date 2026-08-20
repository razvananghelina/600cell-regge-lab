# Runtime amendment: exact homogeneous weak-pole resolver

Date: 2026-08-20  
Status: **preregistered implementation-only amendment after a disclosed failed run**

## Observed failed run

The first targeted execution of source commit `5889c2e` was started from the
clean worktree with

```text
/home/razvan/science/.venv/bin/python
reproducible/verify_gravity_600cell_full_scale_strut_homogeneous_resolution.py
```

It printed two PASS lines and then

```text
[FAIL] the exact momentum and all three carrier differentials vanish identically
```

The process subsequently remained CPU-bound in `sympy.factor` while constructing
`wrong_expected = factor(P_S*P_Z*(1-LAMBDA))`.  It was interrupted after about
11 minutes.  The traceback ended inside SymPy polynomial factorization.  No
result artifact was written.

This run is a **METHOD FAILURE**, not an accepted positive or negative result.
The failed Boolean bundled four elementary differential identities after
substitution of the fully expanded action derivatives, so it did not identify
which identity SymPy failed to normalize.  Separately replacing the two action
derivatives by algebraically independent symbols makes all four differences
reduce exactly to zero; this scratch diagnostic is not itself accepted evidence.

## Frozen implementation correction

The scientific formulae, inputs, thresholds, controls and outcome hierarchy in
the original protocol remain unchanged.  Only the normalization strategy is
changed:

1. Construct the actual action derivatives `p_s` and `p_z` without asking SymPy
   to factor the full expressions.  They remain the expressions used by the
   independent P160/P220 numerical bridge.
2. Prove the generator, carrier-differential and corruption identities over the
   rational function field in algebraically independent symbols `P_s,P_z`.
   This is stronger with respect to the action: the identities then hold for
   every values of the two derivatives, rather than because of an accidental
   special relation in this action.
3. Test nontriviality of the corruption controls only after substituting the
   actual automatic derivatives at the frozen background, as already required
   by the original protocol.
4. Do not call `factor`, `simplify`, `trigsimp` or an equivalent global
   normalizer on a product containing the expanded action derivatives.

The amended verifier must pin this file by SHA-256.  The correction is frozen
before editing or rerunning the verifier.  It does not permit changing any
numerical tolerance or interpreting the omitted pole equation.

