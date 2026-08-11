# The ineffective binary crossed product also fails KO6 Poincare duality

Date: 2026-08-11

Protocol commit: `32d00c7`.

Registered verifier:
`reproducible/verify_hopf_hessian_binary_crossed_real.py`.
Targeted exact result: `15/15`.

No Hessian or selector target is used.

## Complete result

Retain the full binary group, including the centre which acts trivially on
the five points:

```text
B_bin=R(P) crossed_product 2I.
```

The action is transitive with point stabilizer the exact binary tetrahedral
subgroup `2T` of order 24.  Therefore

```text
B_bin ~= M5(R[2T]).
```

Under KO6 and strict Poincare duality, no finite real spectral triple over
this algebra exists, arbitrary multiplicities and Dirac operators included.

## Exact binary stabilizer

The verifier reconstructs all 120 quaternion products.  The maximum squared
matching error back to the indexed group is `3.14e-21`.  Taking the preimage
of the exact effective stabilizer chain gives

```text
Q8 normal in 2T,
|Q8|=8,
|2T|=24,
2T/Q8=C3.
```

This identifies the stabilizer internally from the binary action rather than
importing a presentation of `SL(2,3)`.

## Characters and real type

The three one-dimensional quotient characters, the defining two-dimensional
`SU(2)` character and its two twists, and the three-dimensional rotation
character give seven exact orthonormal irreducibles:

```text
degrees =1,1,1,2,2,2,3.
```

Their squares sum to 24.  The exact Frobenius--Schur indicators are

```text
1,0,0,-1,0,0,1.
```

Thus:

- the trivial and three-dimensional representations are real;
- the untwisted two-dimensional spinor is quaternionic;
- the two nontrivial one-dimensional characters form a complex pair;
- their two-dimensional twists form a second complex pair.

Consequently

```text
R[2T] ~= R + C + H + M2(C) + M3(R)
```

and

```text
B_bin ~= M5(R)+M5(C)+M5(H)+M10(C)+M15(R).
```

Its real dimension is

```text
25+50+100+200+225=600=5*|2I|.
```

Complexification produces the seven blocks

```text
M5^3 + M10^3 + M15,
```

whose squared dimensions also sum to 600.

## KO6 obstruction

There are five real simple summands, hence

```text
rank K0(B_bin)=5.
```

KO6 makes the intersection form antisymmetric.  The generic five-by-five
alternating matrix has identically zero determinant and maximum rank four.
It cannot pair a rank-five `K0` group nondegenerately.

Therefore:

> **DERIVED FULL-ARENA POINCARE NO-GO.**  The canonical ineffective binary
> crossed product admits no strict KO6 Poincare-dual finite real triple,
> regardless of bimodule multiplicities or `D`.

## Comparison with the effective quotient

| groupoid | real algebra type | real `K0` rank | KO6 PD |
|---|---|---:|---:|
| `R(P) crossed A5` | `M5(R)+M5(C)+M15(R)` | 3 | impossible |
| `R(P) crossed 2I` | `M5(R)+M5(C)+M5(H)+M10(C)+M15(R)` | 5 | impossible |

Keeping binary isotropy adds a quaternionic block and another complex block,
but it changes odd rank three into odd rank five rather than curing parity.

## Scope

The no-go assumes:

- the canonical real groupoid algebra;
- KO6, hence antisymmetric intersection form;
- strict Poincare duality.

Changing to the complexification, changing KO dimension, dropping Poincare
duality or selecting a different real subalgebra changes the arena.  The
result does not claim those altered theories are impossible.

## Status ledger

- **DERIVED:** the exact stabilizer is `2T`, with `Q8` normal and quotient
  `C3`.
- **DERIVED:** all seven characters and FS indicators are exact.
- **DERIVED:** `R[2T]=R+C+H+M2(C)+M3(R)`.
- **DERIVED:** the binary crossed product has five real summands and real
  dimension 600.
- **DERIVED FULL-ARENA NO-GO:** strict KO6 Poincare duality fails by odd
  `K0` rank, arbitrary multiplicities and `D`.
- **DERIVED NEGATIVE:** retaining the ineffective binary centre does not
  rescue the five-point matter route.
- **OPEN:** a different geometry-derived real algebra with even `K0` rank.

## Programme consequence

Both canonical groupoid choices are now closed under the same physical
axioms.  Quotienting the centre gives rank three; retaining it gives rank
five.  Neither supplies the required even-rank real finite algebra.  A future
continuation must change the algebraic construction itself rather than
toggle between `A5` and `2I`.
