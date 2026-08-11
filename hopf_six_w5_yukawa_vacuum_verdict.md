# The six Hopf vacua retain too much symmetry to connect the carrier

Date: 2026-08-11

Protocol commit: `b9623c4`.

Target-blind coupling enumeration commit: `c148821`.

Registered Phase-2 verifier:
`reproducible/verify_hopf_six_w5_yukawa_vacuum_gate.py`.
Targeted exact result: `16/16`.

No matter character, mass, coupling or Standard-Model target was used.  The
full suite was not run, by explicit user instruction.

## Headline

The derived five-component order parameter does break the static `A5`
symmetry when it chooses one of the six Hopf points, but it breaks it only to
the point stabilizer

```text
A5 -> D5.
```

That residual symmetry is already fatal to connectedness on the fixed
936-state carrier.

The verifier gives the proposal maximal possible freedom: on every legal
central link it includes the **entire** `D5`-equivariant Hom space available
at the vacuum, including both the invariant affine baseline and every linear
`W5` coupling.  Even this overcomplete span has

```text
maximal algebra-commutant dimension  11 (x4 readings), 12 (x4 readings),
maximal commutator-map rank         349 (x4),          348 (x4),
connected readings                   0/8,
connected signed simplex vacua       0/12.
```

This is a **DERIVED LINEAR-FIELD NO-GO** under the preregistered scope, and a
broader **DERIVED RESIDUAL-SYMMETRY NO-GO** for any first-order covariant
operator evaluated at the same six-point orbit.

## Independent reconstruction of the six-point orbit

The verifier does not import rotation eigenvectors or a precomputed
six-fibration permutation table.  It constructs `A5` as all 60 even
permutations of five symbols, finds its 24 order-five elements and their six
Sylow-five subgroups, and lets `A5` act on those six subgroups by conjugation.

The action is faithful and transitive.  Every point stabilizer has exact
order census

```text
identity:       1,
order two:      5,
order five:     4,
```

hence is `D5` of order ten.  Orbit--stabilizer gives `60/10=6`, independently
recovering the abstract orbit type of the six derived fibrations.

The negative signed orbit has the same stabilizers, so changing the unresolved
affine-action sign does not change the obstruction.

## Why the blind multiplicities become full `D5` Hom spaces

The six-point permutation module is

```text
R[A5/D5] = 1 + W5.
```

For any pair of node modules, Frobenius reciprocity gives the evaluation
identity

```text
Hom_A5(R[A5/D5], Hom(V_i,V_j)) = Hom_D5(V_i,V_j).
```

Splitting the left-hand side into its constant and centered parts yields

```text
Hom_D5 = Hom_A5 + Hom_A5(W5,Hom).
```

The exact matrices found in the two committed stages are

```text
Hom_A5 = [[2,0,1,1],
          [0,2,1,1],
          [1,1,3,2],
          [1,1,2,3]],

W5 linear couplings =
         [[4,2, 6, 6],
          [2,4, 6, 6],
          [6,6,12,12],
          [6,6,12,12]],

Hom_D5 =
         [[6,2, 7, 7],
          [2,6, 7, 7],
          [7,7,15,14],
          [7,7,14,15]].
```

Thus Phase 2 does not select a favorable point in any of the preregistered
2-, 6- or 12-dimensional coupling spaces.  It imposes all of them, plus the
constant invariant part, simultaneously.  Any actual affine Dirac family
has a commutant at least as large as this maximal common commutant.

## Exact residual-symmetry obstruction

Restriction of the four node modules to the vacuum stabilizer has
multiplicity rows, in the order
`(1, reflection-sign, positive-doublet, negative-doublet)`,

```text
V0|D5 = (2,0,1,1),
V1|D5 = (0,2,1,1),
V2|D5 = (1,1,3,2),
V3|D5 = (1,1,2,3).
```

For each of the four real `D5` irreducible types, put its isotypic projector
in every simple matrix block where it occurs and zero elsewhere.  These four
elements belong to

```text
B_R=M6(R)+M6(R)+M12(R)+M12(R),
```

are nonzero, mutually orthogonal, sum to the identity, and commute with every
`D5` intertwiner on every legal link.  They are checked exactly as explicit
0/1 matrices in all eight maximal constraint systems.

Therefore the algebra commutant contains at least four independent lines and
cannot be the scalar line required by connectedness.  The complete rank
calculation sharpens the lower bound to dimensions 11 or 12.

There is also a direct general formulation.  If a covariant family obeys

```text
D(gX)=U_g D(X) U_g^-1
```

and `hX=X`, then `D(X)` commutes with the residual `D5` action.  First order
puts its odd blocks in tensor-factor form, so the corresponding residual
representation elements and isotypic projectors lie in the represented
algebra commutant.  This argument does not depend on linearity in `X`.

## Hostile framing audit

1. The result does not say that spontaneous symmetry breaking is useless.
   It says the particular six-vacuum breaking `A5 -> D5` is incomplete for
   this connectedness axiom.
2. The large Phase-1 multiplicities already killed canonicity before this
   gate.  The Phase-2 negative is stronger: even arbitrary choices inside
   all those spaces cannot repair connectedness while `D5` survives.
3. The affine Hessian action remains **STRUCTURAL**, not a licensed inner
   fluctuation.  The no-go grants it anyway, so this provenance weakness
   cannot be blamed for the failure.
4. A generic point of `W5` may have trivial stabilizer, but the certified
   selector does not choose a generic point; it chooses exactly the six
   `D5`-fixed directions.  Replacing its minima after seeing this result
   would change the theory.
5. A second order parameter could break the residual `D5` further.  No such
   field or action is presently selected, so adding one is an open new route,
   not a repair already contained in the calculation.

## Status ledger

- **DERIVED:** the six-point `A5` orbit and all six order-ten `D5`
  stabilizers, reconstructed combinatorially.
- **DERIVED:** the four exact `D5` restriction rows and full `Hom_D5` matrix.
- **DERIVED:** `Hom_D5=Hom_A5+W5-couplings` entry by entry.
- **DERIVED:** four explicit nontrivial residual-isotypic projectors commute
  with every maximal legal vacuum span.
- **DERIVED:** maximal commutant dimensions are 11 or 12 for all eight
  readings.
- **DERIVED LINEAR-FIELD NO-GO:** connected hit fractions are `0/8` readings
  and `0/12` signed vacua.
- **DERIVED RESIDUAL-SYMMETRY NO-GO:** any first-order covariant `D` at these
  `D5`-stabilized vacua is nonconnected on this carrier.
- **OPEN:** a derived second field/vacuum that breaks `D5` completely, or a
  different carrier whose residual action is scalar.
- **NO TARGET COMPARISON:** no desired matter spectrum was used.

## Programme boundary

The current six-Hopf selector cannot simultaneously serve as the vacuum that
makes the finite spectral triple connected.  It selects too symmetric a
point.

The next honest question is no longer whether the five-component field can
be coupled more cleverly.  It is whether the geometry supplies a second,
independently derived order parameter whose common stabilizer with one Hopf
point is trivial.  Before testing any Dirac rank, the complete list of such
geometric fields and their stabilizer intersections must be enumerated.
