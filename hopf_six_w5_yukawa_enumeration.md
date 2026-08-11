# STEP 1: blind census of `W5`-covariant Dirac couplings

Date: 2026-08-11

Protocol commit: `b9623c4`.

Registered verifier:
`reproducible/verify_hopf_six_w5_yukawa_enumeration.py`.
Targeted exact result: `8/8`.

This commit contains no simplex-vacuum matrix, no Dirac commutant and no
physical target comparison.

## Complete exact census

For

```text
V0=1+5,
V1=3+3',
V2=3+4+5,
V3=3'+4+5,
W5=5,
```

exact character inner products give

```text
m_ij = dim Hom_A5(W5, Hom_R(V_i,V_j))

M_W5 = [[ 4, 2,  6,  6],
        [ 2, 4,  6,  6],
        [ 6, 6, 12, 12],
        [ 6, 6, 12, 12]].
```

The six off-diagonal multiplicities have multiset

```text
2 (x1), 6 (x4), 12 (x1).
```

There are no zero-dimensional coupling spaces, but there are also no
one-dimensional ones.

## All eight legal supports

The verifier independently reconstructs the eight signed lexicographic
gradings.  Every reading again has six positive cells, eight first-order
legal odd cell positions and three central node links.  On every one of the
eight readings, every required `W5` coupling space has dimension greater than
one.

```text
all-unique readings       0/8,
zero-link readings        0/8,
ambiguous readings        8/8.
```

Therefore representation theory plus the existence of the `W5` order
parameter does not select a Yukawa/Dirac tensor.  Picking a Clebsch--Gordan
line after evaluating connectedness would be forbidden fitting.

## Precise interpretation

- **DERIVED:** the full multiplicity matrix and off-diagonal multiset above.
- **DERIVED:** all eight legal supports are coupling-ambiguous.
- **DERIVED CANONICITY NEGATIVE:** `W5` covariance alone selects no tensor
  line on any required central link.
- **NOT TESTED YET:** whether the entire maximal covariant span can connect
  the carrier at a nonzero simplex vacuum.
- **OPEN:** whether the already derived full label Hessian supplies a
  particular functorial contraction inside these large spaces.
- **NO TARGET COMPARISON:** no vacuum, matter character, mass or coupling was
  inspected.

## Phase-2 boundary

The next calculation may impose the entire `W5`-covariant span at a simplex
vacuum.  If even this deliberately overcomplete family leaves a non-scalar
commutant, the linear order-parameter route is dead.  If it connects, that
is only a **STRUCTURAL OPENING**: the current geometry still owes a rule that
selects one point in coupling spaces of dimensions 2, 6 and 12.

The blind enumeration in this file must remain fixed regardless of the
Phase-2 outcome.
