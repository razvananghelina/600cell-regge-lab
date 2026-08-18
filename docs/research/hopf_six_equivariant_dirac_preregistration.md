# Preregistration: canonical equivariant Dirac gate on the 936-state carrier

Date: 2026-08-11

## Disclosed input and candidate

This protocol follows the structural carrier result committed in `dc48128`.
It is not blind to the following representation-theory observation, obtained
before implementation.

For the four `D5` stabilizer irreps `(1,sgn,rho,rho_sigma)`, Frobenius
reciprocity is expected to give the real `A5` induction decompositions

```text
V0=Ind(1)         = 1+5,       dim 6,
V1=Ind(sgn)       = 3+3',      dim 6,
V2=Ind(rho)       = 3+4+5,     dim 12,
V3=Ind(rho_sigma) = 3'+4+5,    dim 12.
```

Consequently the expected off-diagonal equivariant Hom dimensions are

```text
dim Hom_A5(V0,V1)=0,
dim Hom_A5(V2,V3)=2,
all other off-diagonal dimensions=1.
```

These values are disclosed candidates, not results of the verifier.

For each of the eight preregistered spectral readings, the 936-state
Krajewski carrier has eight first-order-compatible ordered odd cell blocks
and a three-edge possible central-link graph.  The candidate Dirac operator
uses the unique normalized `A5` intertwiner on each required node link,
tensored with the identity on the unchanged bimodule factor, in every legal
cell position.  Adjoint blocks are then forced by self-adjointness.

No arbitrary rectangular matrices, Schur coefficients or target comparison
are allowed.

## Frozen character calculation

The verifier must independently reconstruct the exact `A5` and `D5`
characters in `Q(sqrt(5))`, restrict every `A5` character to the selected
stabilizer classes, and derive all four induction decompositions by exact
inner products.  It must then compute the complete `4 x 4` Hom Gram matrix.

A spectral reading is eligible for a canonical equivariant rook operator
only if all three links required by its connected possible-link graph have
one-dimensional Hom spaces.

- a zero-dimensional required Hom kills the reading;
- a Hom space of dimension greater than one makes that link noncanonical and
  kills the reading for this protocol;
- no favorable line may be chosen from a higher-dimensional Hom space.

The hit fraction among all eight readings must be reported before any Dirac
gate is interpreted.

## Frozen normalized intertwiner construction

Use one fixed real model for each multiplicity-one `A5` irrep and realize the
four induced modules as the displayed orthogonal direct sums.  On a
one-dimensional Hom space, the unique intertwiner is the identity between
the common irreducible summands and zero on their complements.  This fixes it
as a partial isometry; only a sign remains.

For every eligible reading:

1. put this partial isometry into every one of the eight legal odd cell
   positions, tensoring it with the identity on the unchanged factor;
2. use equal coefficient magnitude one in all positions;
3. add the transpose blocks required by self-adjointness;
4. enumerate all independent signs on the three central-link intertwiners;
5. record whether the sign choices are related by cellwise orthogonal gauge
   transformations and whether any gate depends on them.

Equal magnitude is part of the parameter-free rook construction.  This
protocol does not classify arbitrary relative magnitudes.

## Frozen exact gates

For every eligible reading and every inequivalent sign choice, compute on the
full 936-dimensional carrier:

- `D*=D`;
- `gamma D=-D gamma`;
- `JD=DJ` for the already fixed cell-transpose real structure;
- order zero;
- first order on a spanning matrix-unit set of the full real algebra;
- rank of `a -> [D,pi(a)]` on all 360 algebra matrix units;
- dimension of `{a in B_R:[D,pi(a)]=0}`;
- whether that kernel is exactly the one-dimensional scalar line;
- dimension of the represented inner one-form space, or at minimum a
  certified nonzero basis element if the full rank is too costly;
- rank and connected components of the actual, nonzero central-link graph.

The already certified orientability and Poincare form are carrier data and
will not be counted again as Dirac successes.

## Acceptance and kill boundaries

- **ADVANCE:** at least one preregistered spectral reading has all required
  Hom dimensions one, and its incidence-normalized equivariant rook operator
  passes reality, first order, connectedness and nonzero forms for every
  inequivalent sign choice.
- **PATTERN:** only a proper subset of the eight readings is eligible; report
  the exact fraction even if an eligible `D` passes.
- **KILL FOR CURRENT EQUIVARIANT ROOK ROUTE:** every reading has a zero or
  ambiguous required Hom, or every eligible canonical `D` fails an axiom.

A failure is scoped to the normalized `A5`-equivariant rook construction on
the full off-diagonal carrier.  It does not rule out a new non-equivariant
geometric field that explicitly breaks the residual symmetry.

No Hessian, Standard-Model character, mass or coupling is used.  Only the new
targeted verifier will be run; the full suite remains excluded by user
instruction.
