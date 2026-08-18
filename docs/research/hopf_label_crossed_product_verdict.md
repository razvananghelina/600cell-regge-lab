# The canonical crossed product does not derive label superselection

Date: 2026-08-10

Protocol commit: `d15e7fa` (tests frozen before the algebra type and selector
comparison were computed).

Verifier: `reproducible/verify_hopf_label_crossed_product.py`.  Targeted
result: `14/14`.

## Complete hypotheses and scope

Let `F` be the six unoriented Hopf fibrations already derived from the binary
icosahedral geometry.  The centre of `2I` acts trivially on `F`, so the
effective action is the transitive `A5` action.  This audit concerns exactly
the transformation-groupoid algebra

```text
B = C(F) crossed_product A5.
```

It does not identify this algebra with the larger crossed product by the
ineffective `2I` action, and it does not choose a subalgebra after inspecting
the desired selector.

## Exact algebra type

The stabilizer of one fibration has ten elements, with element-order census

```text
1^1, 2^5, 5^4.
```

It is `D5`.  Its four conjugacy classes have sizes `1,2,2,5`, and its complex
irreducible degrees are `1,1,2,2`.  The transitive-groupoid isomorphism gives

```text
C(A5/D5) crossed_product A5
    ~= M6(C[D5])
    ~= M6(C) + M6(C) + M12(C) + M12(C).
```

Both sides have complex dimension `360`, and the centre of the crossed
product has dimension four.

**DERIVED:** this is a geometrically fixed noncommutative Wedderburn type.  It
is not the Standard-Model algebra: its simple blocks have sizes `6,6,12,12`.

## What happens on the six-label carrier

In the natural covariant representation on `H_F=C^6`,

```text
delta_i -> E_ii,
u_g     -> P_g,
```

the exact span of all `E_ii P_g` is the full `M6(C)`: it contains all 36
matrix units and has scalar commutant.  Consequently this representation is
not faithful; its kernel has dimension `360-36=324` and it selects the `M6`
block induced from the trivial stabilizer representation.

Every label projection is noncentral.  Group unitaries conjugate the six
projections transitively, so they are mutually equivalent and have common
central support one.

**DERIVED NEGATIVE:** the canonical crossed product does not turn the six
fibrations into superselection sectors.  It supplies every transition between
them.  Expecting the transformation-groupoid algebra itself to forbid those
transitions reverses its mathematical role.

## Conditional expectations do not repair the result

A `C(F)`-bimodular map back to `C(F)` kills off-diagonal corners, but each
diagonal corner is a copy of `C[D5]`.  Before equivariance, the relevant
linear-functional space has dimension `6*10=60`.  `A5` equivariance reduces
the choice to a conjugation-invariant functional on `D5`, a four-dimensional
space; imposing the value on the unit leaves a three-dimensional affine
family.

In particular, the regular coefficient trace and the trivial-character state
are distinct positive, unital, conjugation-invariant choices: on a nonidentity
stabilizer element their values are respectively zero and one.  Thus even
positivity and equivariance do not select a unique diagonal expectation.

The familiar coefficient-of-the-identity expectation is canonical only after
the regular trace (equivalently, the regular crossed-product construction) is
included as input.  It is not forced by the abstract algebra and action alone.

**DERIVED NEGATIVE:** projecting dynamics onto the diagonal requires an
additional trace/state choice.

## Consequence for the Hopf selector

The already constructed overlap map `Phi(X)` belongs to the diagonal
subalgebra `C(F)`, but is noncentral in `B` for nonzero `X`.  Therefore the
diagonal form of `D_aux` is not a superselection statement in this algebra;
it is a chosen conditional expectation or an additional locality axiom.

This reaches the preregistered kill boundary for deriving diagonal label
superselection from the canonical crossed product.

## Status ledger

- **DERIVED:** `C(F) crossed_product A5` has type
  `M6 + M6 + M12 + M12`, dimension 360 and centre dimension four.
- **DERIVED:** its natural six-label image is all of `M6`, with a
  324-dimensional kernel.
- **DERIVED NEGATIVE:** the six label projections are equivalent and
  noncentral, not superselection sectors.
- **DERIVED NEGATIVE:** equivariant positive diagonal expectations are not
  unique without selecting a state/trace.
- **DERIVED NEGATIVE:** the crossed product does not license the diagonal
  auxiliary selector.
- **STRUCTURAL ADVANCE:** the noncommutative algebra type itself is selected
  by the derived action, rather than fitted.
- **OPEN:** whether a separately derived real spectral triple over the full
  crossed product satisfies KO6, orientability, first order and the other
  physical gates.  Such a triple would describe transition sectors; it would
  not retroactively prove diagonal superselection.

This open real-triple question is now sharpened in
`hopf_six_crossed_real_galois_verdict.md`.  The split real algebra is
`M6(R)^2+M12(R)^2` and passes even-rank KO6 parity, but exact golden-Galois
descent forces the intersection Pfaffian to vanish.  The minimum faithful
left double, regular/enveloping doubles and natural label double all fail a
cheap gate.  Thus no canonical carrier reaches selector comparison.

## Honest next boundary

The two canonical label algebras now give complementary no-go results:

1. `C^6` can be diagonal, but the full `A5`-equivariant KO6,
   metric-dimension-zero orientability requirements forbid every nonzero
   bimodule in the audited arena.
2. `C(F) crossed_product A5` is intrinsically noncommutative, but its natural
   representation contains all label transitions and has no intrinsic
   diagonal expectation.

Continuing toward the six-point selector therefore requires a new physical
input that actually means locality or superselection.  Merely changing
between the function algebra and its full transformation groupoid cannot
supply it.

The two parameter-free dynamical candidates already present have now also
been audited in `hopf_dynamical_superselection_verdict.md`.  The `A5` twirl
fixes `C+C`, while the complete cubic-Hessian evolution preserves only scalar
diagonal observables at the six `Box_i`.  Neither derives `C(F)`.
