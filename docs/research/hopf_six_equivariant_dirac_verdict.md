# The canonical equivariant rook Dirac is never connected

Date: 2026-08-11

Protocol commit: `705edeb`.

Registered verifier:
`reproducible/verify_hopf_six_equivariant_dirac.py`.
Targeted exact result: `15/15`.

No Hessian, Standard-Model module, mass or coupling target was used.

## Headline

The incidence-labelled algebra admits a robust orientable, Poincare-dual
936-state KO6 carrier.  Its most economical `A5`-equivariant Dirac candidate
can now be constructed exactly on half of the preregistered spectral
readings.

It passes self-adjointness, oddness, reality, first order and nonzero forms,
but it is never connected:

```text
eligible spectral readings                    4/8,
normalized sign variants                       32,
connected variants                              0,
algebra commutant dimensions              109 or 141.
```

Therefore this continuation reaches a **DERIVED SCOPED KILL**: the canonical
normalized `A5`-equivariant rook operator on the full off-diagonal carrier is
not a finite connected spectral triple.

The result does not close a future non-equivariant geometric `D` that carries
additional symmetry-breaking data.

## Exact induced modules

Restrict the exact `A5` character table to the selected `D5` stabilizer, whose
classes are identity, five reflections, two fibre-edge rotations and two
distance-two rotations.  Exact Frobenius reciprocity gives

```text
V0=Ind_D5^A5(1)         = 1+5,
V1=Ind_D5^A5(sgn)       = 3+3',
V2=Ind_D5^A5(rho)       = 3+4+5,
V3=Ind_D5^A5(rho_sigma) = 3'+4+5.
```

Their dimensions are `(6,6,12,12)`, and the complete equivariant Hom Gram
matrix is

```text
G = [[2,0,1,1],
     [0,2,1,1],
     [1,1,3,2],
     [1,1,2,3]].
```

Thus:

- `Hom(V0,V1)=0`;
- `dim Hom(V2,V3)=2`;
- every other off-diagonal Hom is one-dimensional.

These are exact representation-theory statements, not numerical ranks.

## Look-elsewhere result: exactly four of eight readings

For each spectral reading, first order leaves eight ordered odd cell blocks
whose possible central-link graph is a three-edge tree.

The complete result is

```text
u_edge first, v_ref tie-breaker: 4/4 eligible,
v_ref first, u_edge tie-breaker:  0/4 eligible.
```

Every reflection-first tree contains the required edge `(V0,V1)`, whose Hom
space is zero.  Every edge-first tree uses three one-dimensional Hom spaces.
The ambiguous two-dimensional `(V2,V3)` space occurs in none of the eight
trees, so no line was fitted inside it.

The hit fraction is therefore `4/8`.  This is labelled **PATTERN**, because
the joint spectrum did not independently settle which central operator has
lexicographic priority.  The calculation does show precisely what a future
priority principle must select: fibre-edge incidence first.

## Parameter-free operator

On a one-dimensional Hom line, the induced modules share exactly one
multiplicity-one real `A5` irrep.  The normalized intertwiner is the identity
on that common summand and zero on the orthogonal complements.  It is a
partial isometry and is unique up to sign.

For each eligible reading, the verifier puts this map into all eight legal
odd cell positions:

```text
T tensor I  when the left node changes,
I tensor T  when the right node changes.
```

Every block is compared exactly with its Kronecker formula; there are no
additional entries.  Adding its transpose fixes self-adjointness.

There are `2^3=8` signs on the three tree links and therefore

```text
4 readings * 8 signs = 32 operators.
```

All signs are related to the all-positive operator by the explicit cellwise
orthogonal gauge

```text
U|H_ij=(g_i g_j) I,
```

because a tree has no sign holonomy.  Hence the 32 raw cases introduce no
physical coefficient fitting.

## Gates that pass

All 32 exact 936-dimensional sparse matrices satisfy

```text
D*=D,
gamma D=-D gamma,
JD=DJ.
```

Order zero was already exhaustive on the carrier.  First order follows here
from, and is checked by, the exact `T tensor I` / `I tensor T` support of
every block.  Every operator is nonzero and the commutator map has positive
rank, providing a certified nonzero inner-one-form witness.

## Connectedness failure

For each tree link with normalized partial isometry `T:V_i->V_j`, commuting
algebra elements must solve exactly

```text
A_j T = T A_i,
A_i T* = T* A_j.
```

The verifier builds the full integer linear system on all

```text
6^2+6^2+12^2+12^2=360
```

matrix-unit coefficients.  Its two exact outcomes are

```text
commutator-map rank 251, kernel dimension 109: 16 operators,
commutator-map rank 219, kernel dimension 141: 16 operators.
```

Connectedness would require kernel dimension one.  None comes close.

The reason is structural: each canonical intertwiner sees only the common
`A5` irreducible summand.  Large orthogonal complements remain free, leaving
whole matrix algebras in the commutant.  Changing signs cannot repair this,
and changing relative nonzero magnitudes does not change the same linear
commutation equations.

## Hostile framing audit

1. The `4/8` eligibility is not a prediction of edge priority.  Selecting
   edge-first after seeing this result would be a look-elsewhere choice.
2. Passing first order is meaningful: every block has the exact tensor-factor
   form required for all matrix units.  But it does not imply connectedness.
3. Nonzero one-forms are only a minimal gate.  Their existence does not cure
   a 109- or 141-dimensional commutant.
4. Adding arbitrary maps on the invisible complementary subspaces could make
   the operator generic and perhaps connected, but those maps lie outside the
   unique equivariant Hom lines.  Without a geometric construction, this
   would be fitted Schur data and is forbidden.
5. The result is scoped to the full off-diagonal carrier and the normalized
   equivariant rook rule.  A field that explicitly breaks residual `A5`
   covariance, or a different derived carrier, is not covered.

## Status ledger

- **DERIVED:** the four exact induced-module decompositions displayed above.
- **DERIVED:** the complete Hom Gram matrix.
- **PATTERN:** exactly `4/8` spectral readings have unique intertwiners on all
  required links; these are precisely the edge-first readings.
- **DERIVED:** all 32 sign variants are cellwise gauge equivalent.
- **DERIVED:** all 32 pass self-adjointness, oddness, `JD=DJ`, first order and
  possess nonzero one-forms.
- **DERIVED NEGATIVE:** connected count is `0/32`; commutant dimensions are
  109 or 141.
- **DERIVED SCOPED KILL:** the canonical equivariant rook construction does
  not produce a connected finite spectral triple.
- **OPEN:** a geometry-selected non-equivariant operator filling the unseen
  complementary sectors without fitted coefficients.
- **NO TARGET COMPARISON:** no Hessian or matter module was inspected.

## Next legitimate boundary

The central separator and equivariant intertwiners have now given everything
they can: node labels, a KO6 carrier, Poincare duality, first order and nonzero
forms, but not connectedness.

A continuation must exhibit a noncentral incidence tensor already present in
the six-fibration groupoid whose projections cover the complementary
`A5` summands.  Before computing a commutant, the complete list of such
tensors and their normalizations must be frozen.  If the repository contains
no such tensor, adding generic rectangular matrices is not a continuation of
the theory.

## Subsequent provenance result

The preregistered inventory is recorded in
`hopf_six_existing_operator_lift_verdict.md`.  None of the seven existing
operator families has both a faithful canonical lift to the 936-state carrier
and a nonzero odd part.  Five have no defined lift; the two crossed-product
families generated by the algebra preserve every central cell and are even.
The next step therefore requires a genuinely new noncentral incidence tensor.
