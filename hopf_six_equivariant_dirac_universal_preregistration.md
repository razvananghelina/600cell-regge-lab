# Preregistration: universal equivariant-Dirac connectedness bound

Date: 2026-08-11

## Disclosed starting point

The normalized-rook audit in commit `0c7ca39` tested 32 concrete Dirac
operators on the 936-state carrier.  It found algebra-commutant dimensions
109 or 141, but its protocol fixed equal coefficient magnitudes and retained
only spectral readings whose three central links had one-dimensional
equivariant Hom spaces.

The following facts are therefore known before this protocol:

```text
B_R = M6(R) + M6(R) + M12(R) + M12(R),
H_off = direct_sum_(i != j) V_i tensor V_j*,
V0 = 1+5,
V1 = 3+3',
V2 = 3+4+5,
V3 = 3'+4+5,
Hom Gram = [[2,0,1,1],
            [0,2,1,1],
            [1,1,3,2],
            [1,1,2,3]].
```

The candidate theorem was noticed before implementation: naturality under
the full automorphism group of the incidence-labelled six-fibration geometry
forces a geometrically canonical scalar operator to be `A5`-equivariant.
First order then restricts each odd cell block to `T tensor I` or
`I tensor T`, with `T` in the corresponding `Hom_A5(V_i,V_j)`.

This protocol asks whether the previous failure extends from the 32
normalized examples to the complete space of such equivariant first-order
operators.  The expected negative is disclosed; the calculation is not
blind.

## Frozen scope

The verifier will cover all eight signed lexicographic gradings already
derived from `(u_edge,v_ref)`.  For each grading it will:

1. reconstruct all first-order-compatible odd cell positions from the
   central-support rule, without importing the previous candidate list;
2. reconstruct the exact `A5` character table, the four induced modules and
   their complete equivariant Hom Gram matrix;
3. determine the central node links occurring in those positions;
4. place the **entire** equivariant Hom space on every legal link, not one
   chosen operator;
5. impose simultaneously all commutation equations
   `A_j T = T A_i` and their adjoints for a basis of every such Hom space;
6. compute the exact rank of this maximal constraint system on the 360 real
   coefficients of `B_R`.

Using the entire Hom span deliberately gives more constraints than any one
Dirac operator can give.  Therefore, if its kernel has dimension greater
than one, every individual equivariant first-order `D` has a non-scalar
algebra commutant.  Zero coefficients, unequal magnitudes and arbitrary
linear combinations inside higher-dimensional Hom spaces cannot improve the
bound.

The verifier must also exhibit a concrete non-scalar element of the maximal
commutant for every reading and check every defining equation exactly.

## Naturality hypothesis

The universal conclusion has the following complete hypothesis:

> `D` is a scalar operator on the fixed 936-state carrier, is odd for one of
> the eight derived gradings, satisfies first order for the full matrix
> algebra `B_R`, and is natural under every orientation-preserving
> automorphism in the derived `A5` action.

The last clause means `D` commutes with the induced diagonal `A5` action.  It
is not inferred merely from the word "geometric".  A choice of one vertex,
one chamber, one fibre or one component of an `A5`-covariant field violates
this hypothesis unless the theory separately selects that choice.

## Acceptance and kill boundaries

- **REFUTATION / ADVANCE:** at least one reading has maximal equivariant
  commutant dimension one.  Then an equivariant connected `D` is not ruled
  out, and the complete coefficient space on that reading must be analysed.
- **DERIVED UNIVERSAL EQUIVARIANT NO-GO:** all eight maximal commutants have
  dimension greater than one.  Then no `A5`-natural incidence tensor on this
  carrier can satisfy connectedness, irrespective of normalization.

This does **not** rule out a dynamically selected non-invariant field, a
spontaneously chosen vacuum, an enlarged carrier, or a different algebra.  A
single basepoint-dependent tensor will be labelled **STRUCTURAL**, not
canonical, unless an independent selection mechanism is supplied.

No Hessian, particle representation, mass, coupling or Standard-Model target
will be inspected.  Only the new targeted verifier will be run; the full
suite remains excluded by user instruction.
