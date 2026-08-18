# The golden-labelled algebra has a robust orientable KO6 carrier

Date: 2026-08-11

Protocol commit: `12c45eb`.

Registered verifier:
`reproducible/verify_hopf_six_spectral_krajewski.py`.
Targeted exact result: `13/13`.

No Hessian, particle module, mass, coupling or Standard-Model target was used.

## Headline

For the incidence-labelled real algebra

```text
B_R=M6(R)+M6(R)+M12(R)+M12(R),
```

take the off-diagonal full enveloping bimodule

```text
H_off=direct_sum_(i != j) C^(n_i) tensor C^(n_j)*.
```

Let `J` transpose the two central indices and conjugate coefficients.  Orient
each unordered pair using any of the eight signed lexicographic readings of
the exact joint node spectrum `(u_edge,v_ref)`, and give the transposed cell
the opposite grading.

Every one of the eight readings gives:

- a faithful left and opposite action;
- order zero;
- `J^2=+1` and `J gamma=-gamma J`;
- an explicit metric-dimension-zero orientation cycle;
- an exactly unimodular antisymmetric intersection form;
- nonempty first-order-compatible odd block positions;
- a connected graph of possible central links.

This is a **STRUCTURAL ROBUST CARRIER EXISTENCE** result.  It removes
orientability and Poincare duality as obstructions for this particular
four-node real arena.  It is not yet a finite spectral triple because the
geometry has not selected a Dirac operator.

## Complete ambiguity census

The exact joint labels are

```text
node                 (u_edge,v_ref)
trivial              (2,5)
reflection sign      (2,-5)
positive doublet     (phi-1,0)
negative doublet     (-phi,0).
```

They do not specify whether `u_edge` or `v_ref` has lexicographic priority,
or the direction of either ordering.  The preregistered census therefore
contains exactly

```text
2 priorities * 2 u-directions * 2 v-directions = 8 readings.
```

All eight produce distinct node orders.  Simultaneously reversing both
directions sends `gamma` and the intersection form to their negatives, so the
eight readings form four grading-reversal pairs.

No favorable ordering was selected after calculation.  The complete gate hit
fraction is `8/8`.

## Carrier and real structure

With node sizes `(6,6,12,12)`, the complex Hilbert dimension is

```text
dim H_off = 2 sum_(i<j) n_i n_j = 936,
dim H_+ = dim H_- = 468.
```

On a cell `H_ij`, the algebra acts on the first tensor coordinate and the
opposite algebra on the second.  The verifier exhausts every potentially
nonzero pair of matrix units on every basis vector:

```text
sum_(i!=j) n_i^3 n_j^3 = 9,051,264 exact action cases.
```

Every commutator vanishes.  This is a genuine exhaustive order-zero check,
not an assertion from block labels alone.

The antiunitary

```text
J: H_ij -> H_ji
```

swaps tensor coordinates and conjugates coefficients.  Therefore `J^2=1`.
The grading is `+1` on the six spectrally oriented cells and `-1` on their six
transposes, giving `J gamma=-gamma J`.  Since gamma is scalar on each central
cell, `[gamma,B_R]=0`.

## Orientability and Poincare duality

Let `z_i` be the central unit of the `i`-th simple block.  The explicit
metric-zero Hochschild cycle is

```text
gamma = sum_(i!=j) gamma_ij pi(z_i) J pi(z_j) J^-1.
```

On `H_kl`, only the term `(i,j)=(k,l)` survives, so the equality is exact.

Using one minimal projection per real matrix block, one oriented copy of cell
`(i,j)` contributes `+1` and its `J` transpose contributes the antisymmetric
entry.  Hence

```text
Q=mu-mu^T.
```

Every spectral order is a transitive tournament on four nodes.  The eight
exact results are

```text
Pf(Q)=+1 four times,
Pf(Q)=-1 four times,
det(Q)=1 eight times.
```

Thus the ordering ambiguity never threatens Poincare duality.  **DERIVED
conditional on the full off-diagonal carrier.**

## What the first-order calculation does and does not say

For a block from positive cell `(i,j)` to the `J` transpose of positive cell
`(k,l)`, the central first-order rule permits it only when

```text
i=l or j=k.
```

Each reading has exactly eight permitted ordered cell-block positions.  Their
induced possible central-link graph is connected in all eight readings.

This proves only that first order and connectedness are not combinatorially
impossible.  It does not choose a nonzero rectangular matrix in any permitted
position, prove `JD=DJ`, or prove that the commutant of an actual `D` is only
the scalars.

## Hostile framing audit

1. The carrier is canonical as the off-diagonal part of the full enveloping
   bimodule, but the theory has not proved that this maximal correspondence is
   the physical Hilbert space.  Therefore the result is **STRUCTURAL**, not a
   selected matter sector.
2. Robustness under all eight spectral orders removes a look-elsewhere
   concern for Poincare duality, but it does not remove the order ambiguity
   itself.
3. Dimension `936` is a carrier dimension, not a particle count, generation
   count or spacetime dimension.
4. A connected graph of allowed Dirac positions is not a connected spectral
   triple.  Choosing generic coefficients would be an existence proof only
   and would violate the present no-fitting standard.
5. The earlier regular and standard-sheet enveloping doubles remain failed.
   This carrier evades that negative because `J` transposes central cells and
   gamma changes sign between `(i,j)` and `(j,i)`; it is a genuinely different
   KO6 construction, not a reinterpretation of the failed sheet profile.

## Status ledger

- **DERIVED:** eight legitimate signed lexicographic readings exist and give
  eight distinct orders.
- **DERIVED:** every reading has exact carrier dimensions `468+468=936`.
- **DERIVED:** order zero passes all 9,051,264 matrix-unit action cases.
- **DERIVED:** the explicit central 0-cycle represents gamma for all eight.
- **DERIVED:** every intersection form is unimodular, with determinant one.
- **DERIVED:** every support has eight legal odd blocks and a connected
  possible-link graph.
- **STRUCTURAL ROBUST EXISTENCE:** the full off-diagonal carrier passes all
  zero-order, KO-sign, orientability and Poincare gates for `8/8` readings.
- **OPEN:** a geometry-selected `D`, its adjoint/reality signs, actual first
  order, connectedness and nonzero inner one-forms.
- **OPEN:** why this 936-dimensional maximal carrier, or a derived reduction
  of it, is physically selected.
- **NO TARGET COMPARISON:** no Hessian or matter character was inspected.

## Next kill gate

Project every already-defined incidence operator of the six-fibration
transformation groupoid onto the eight permitted odd cell-block positions.
The list must be frozen before testing connectedness.

If every projected canonical operator vanishes, violates `JD=DJ`, or leaves a
nontrivial algebra commutant, this carrier does not produce a finite spectral
triple.  If one survives, its coefficients must be fixed by incidence itself;
a generic fill of the eight legal blocks does not count.

## Subsequent Dirac result

The preregistered continuation is recorded in
`hopf_six_equivariant_dirac_verdict.md`.  Exactly `4/8` spectral readings have
unique equivariant intertwiners on every required link.  All 32 normalized
sign variants pass reality, oddness, first order and nonzero forms, but every
one fails connectedness with algebra commutant dimension 109 or 141.  Thus
the canonical equivariant rook route is closed on this carrier.
