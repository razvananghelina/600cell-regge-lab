# Matter module construction: hypercharge no-go

Date: 2026-07-22 (session 2)

## Decision

The already-derived nonabelian seeds construct a Standard-Model-shaped
15- or 16-dimensional vector space, but none of the already-derived integer
gradings defines hypercharge on it.  Consequently a genuine simultaneous
derived `u(1)+su(2)+su(3)` matter module has **not** been constructed.

This is a sharpened no-go, not a claim that the Standard Model module is
mathematically impossible.  The missing datum is an equivariant map from the
derived discrete flavor space to chiral Weyl summands, together with an
everywhere-defined integer operator in the commutant of the nonabelian action.

## Nonabelian construction and dimension enumeration

Let `V2` be the defining complex `SU(2)` module restricted to `2I`, and let
`V3` be the defining complex `SU(3)` module whose restriction uses the derived
`3'` embedding `A5 -> SO(3) -> SU(3)`.  Their external tensor products give
the block types

`(3,2), (3,1), (1,2), (1,1)`

of dimensions `6,3,2,1`.

- **DERIVED:** the two defining actions and the `3'` embedding exist.
- **STRUCTURAL:** placing independently derived actions on an external tensor
  product is canonical representation theory, but the repository has no
  discrete matter space identified with that tensor product.
- **DERIVED:** solving `6q+3c+2l+s=15` and `=16` over nonnegative integers
  gives respectively 42 and 48 dimension decompositions.
- **DERIVED:** imposing the requested SM multiplet inventory selects
  `(q,c,l,s)=(1,2,1,1)` for 15 dimensions and `(1,2,1,2)` for 16 dimensions.
- **OPEN:** the real `A5` restriction does not distinguish the complex color
  fundamental from its conjugate.  The chiral choices `3` versus `bar(3)` are
  not selected by the finite `A5` data used here.

Thus the candidate vector spaces are

`M15=(3,2)+(bar3,1)+(bar3,1)+(1,2)+(1,1)`

and `M16=M15+(1,1)`.  Their existence as abstract nonabelian modules is not
the missing step; the derived abelian action is.

## Exact Standard-Model benchmark

For left-handed fields

`Q=(3,2)_(1/6), u^c=(bar3,1)_(-2/3), d^c=(bar3,1)_(1/3),`
`L=(1,2)_(-1/2), e^c=(1,1)_1 [, nu^c=(1,1)_0]`,

and `T(fundamental)=1/2`, exact computation gives

`(T1,T2,T3)=(10/3,2,2)`, hence ratio `5:3:3`.

The anomaly tuple

`(Tr Y, Tr Y^3, Tr(Y T_2^2), Tr(Y T_3^2))`

is exactly `(0,0,0,0)`.  The weighted number of weak doublets is
`3+1=4`, hence the Witten parity is zero.  A sterile `nu^c` changes none of
these results.

- **DERIVED:** all of the preceding trace and anomaly statements.
- **PATTERN:** the old prefactor target `8:5:2`; a physical generation gives
  `5:3:3`, so the old target is not pursued.

## Audit of the three proposed abelian gradings

The exact charged-flavor data by generation are

| generation | lepton | up | down | exponent triple |
|---|---|---|---|---|
| 0 | `e:(0,0)` | `u:(3,-2)` | `d:(1,0)` | `(0,3,5)` |
| 1 | `mu:(1,1)` | `c:(2,1)` | `s:(1,1)` | `(11,16,11)` |
| 2 | `tau:(1,2)` | `t:(4,1)` | `b:(-1,4)` | `(17,26,19)` |

Here `n=5a+6b`.

### McKay exponent `n`

A `u(1)` action commuting with `su(2)` must be scalar on each irreducible
weak doublet.  Under the only available flavor identification, the quark
doublet pairs have exponent pairs

`(3,5), (16,11), (26,19)`.

- **DERIVED (negative):** none is equal.  Therefore `n` cannot itself be the
  hypercharge generator on the required weak doublets.
- **OPEN:** the nine mass labels also contain no neutrino labels and do not
  split a Dirac species into its left- and right-handed Weyl summands.  Hence
  they do not define an operator on `M15` or `M16`.

### Hopf/C10 grading

The decagonal character supplied by the same exponent has residues

`(3,5), (6,1), (6,9) mod 10`

on the three quark pairs.

- **DERIVED (negative):** this also violates the weak-doublet commutant
  condition in every generation.
- **OPEN:** no separate derived map from C10 winding eigenspaces to all chiral
  Weyl blocks exists.  Treating arbitrary residues as rational charges would
  be an assignment by hand.

### `Z[phi]` unit grading

For `z=a+b phi`, `N(z)=a^2+ab-b^2`.  Across the nine slots the norms are

`0,-1,1,1,5,1,-1,19,-19`.

- **DERIVED (negative):** a unit exponent is defined only when `|N(z)|=1`.
  It is therefore not an everywhere-defined grading of the nine slots, much
  less of the 15/16 Weyl components.

## Sharp obstruction

- **DERIVED:** the desired nonabelian dimensions and SM anomaly target are
  internally consistent.
- **DERIVED (negative):** `n` and `n mod 10` fail the necessary commutant
  condition for every quark doublet; the unit grading is not globally defined.
- **STRUCTURAL:** an abstract SM-shaped `su(2)+su(3)` tensor module can be
  written using the already-derived defining actions.
- **OPEN:** no derived chiral-slot map, color conjugation choice, neutrino
  placement, or abelian generator exists.  Therefore charges cannot be
  computed without inserting the Standard Model hypercharges by hand.

The missing ingredient is specifically **a derived chiral matter functor**
from the McKay/Hopf/flavor data to `M15` or `M16`, carrying an integer
commutant operator.  Anomaly cancellation would then be a genuine test of its
eigenvalues.  At present the anomaly-cancelling numbers can only be verified
as the external SM benchmark, not derived from the discrete theory.

All finite claims are verified by `reproducible/verify_matter_module.py`.
