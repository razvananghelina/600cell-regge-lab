# Protocol: full Hessian of `Tr(Box^4)`

Date: 2026-08-12

## Status and provenance

This is a **corrective, post-anomaly protocol**, not a blind discovery
protocol.  An exploratory audit has already shown that the matrix called the
Hessian in `reproducible/verify_gravity.py` may omit terms.  Therefore commit
ordering can freeze the correction test, but cannot make the eventual result
target-blind.

No physical target, graviton mode count, Planck scale, or phenomenological
number is used below.

## Frozen hypotheses

Let `B` be the real symmetric `120 x 120` matrix used by
`verify_gravity.py`:

```text
B = 5 A_fiber - A_cross,
```

where `A_fiber` contains the 120 edges internal to the twelve derived Hopf
decagons and `A_cross` contains the other 600 edges of the 600-cell graph.
For every undirected edge `e = {i,j}`, define

```text
E_e = c_e (|i><j| + |j><i|),
c_e = 5 on a fiber edge and -1 on a cross edge.
```

The independently variable weighted operator is

```text
B(w) = B + sum_e w_e E_e.
```

This is the complete hypothesis list.  In particular, the calculation does
not quotient gauge directions and does not assume that the variables `w_e`
are metric variables, gravitons, or physical degrees of freedom.

## Exact object to compute

For

```text
F(w) = Tr(B(w)^4),
```

the full mixed Hessian at `w = 0` is fixed before implementation as

```text
H_ef = 4 Tr(E_f B^2 E_e + B E_f B E_e + B^2 E_f E_e).
```

Equivalently, for real symmetric `B`, `E_e`, and `E_f`,

```text
H_ef = 8 Tr(B^2 E_f E_e) + 4 Tr(B E_f B E_e).
```

The notation in the first display will not be used in code; the second
formula is the frozen implementation formula.  It follows by differentiating
the noncommutative polynomial, not by differentiating eigenvalues alone.

The legacy matrix to audit is

```text
G_ef = 12 sum_k lambda_k^2 S_ke S_kf,
S_ke = <psi_k, E_e psi_k>.
```

It retains only diagonal eigenbasis sensitivities.  It is the full Hessian
only if exhaustive equality with `H` is established.

## Required checks

The registered verifier must perform all of the following.

1. Reconstruct the 600-cell, the twelve Hopf fibers, `B`, and all 720
   matrices `E_e`; certify the `120 + 600` edge split.
2. Certify the exact `120 x 120` rank of integer `B` over the rationals by a
   full-rank minor or modular rank, and report its nullity.
3. Assemble every entry of the integer `720 x 720` matrix `H` from the frozen
   trace formula and certify symmetry exactly.
4. Compare selected entries and deterministic directions against a direct
   matrix-polynomial expansion of `Tr((B+tE)^4)`.  No finite-difference fit is
   admissible.
5. Certify the exact rational rank of `H` by modular elimination.  Numerical
   eigenvalues may report its inertia only after the analytic positivity
   statement below is checked.
6. Use the standard spectral Hessian identity for real symmetric `B`:

   ```text
   d^2 Tr(B^4)[E,E]
     = sum_i 12 lambda_i^2 E_ii^2
       + sum_{i<j} 8(lambda_i^2 + lambda_i lambda_j + lambda_j^2) E_ij^2
   ```

   in a `B` eigenbasis.  Thus the full Hessian is positive semidefinite on all
   symmetric perturbations.  Combined with exact full rank on the 720 edge
   directions, this implies positive definiteness on that restricted space.
7. Reproduce the legacy `G` construction, compare `G` and `H`, and rotate the
   numerical eigenbasis inside degenerate eigenspaces.  If `G` changes, label
   it basis-dependent and therefore noncanonical.  A numerical eigensolver's
   chosen basis is not geometry.

## Decision rule

- If `G = H` exhaustively and its rank/inertia claims survive, retain the old
  result.
- If `G != H`, or its reported rank is not the full Hessian rank, the claims
  "101 positive modes", "619 zero modes", and "graviton propagates at order
  4" are withdrawn.  `G` must be relabelled as a diagonal-eigenvalue-
  sensitivity Gram matrix, and every current summary that uses it as a
  graviton Hessian must be corrected.

Even if `H` is positive definite, this alone is only an algebraic stiffness
result.  A physical graviton interpretation additionally requires a derived
field dictionary, gauge quotient, dynamics, continuum limit, and coupling to
sources; none is assumed here.

## Labels fixed in advance

- Exact polynomial derivative and exact modular ranks: **DERIVED**.
- Numerical eigenvalue ranges used only as controls: **STRUCTURAL**.
- Identification of any edge-weight mode with a graviton: **OPEN** unless a
  separate derivation supplies the missing physical dictionary.
