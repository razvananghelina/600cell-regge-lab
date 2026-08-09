# Free-cell convolution theorem and the failed fluctuation dichotomy

Date: 2026-07-27

## Decision

The requested freeness and convolution statements are **DERIVED**.  The
proposed two-horn inner-fluctuation dichotomy is **REFUTED**.

The exact oriented 600-cell cochain module is

`C = C^0 direct-sum C^1 direct-sum C^2 direct-sum C^3
   = C[2I] tensor C^22`,

with free-orbit counts `(1,6,10,5)`.  In the explicit orbit chart below every
`2I`-equivariant operator has the unique form

`T = sum_(alpha,beta) R_(w_alpha,beta) tensor E_(alpha,beta)`.

In particular the derived Kähler--Dirac operator `D=d+d*` has this form.
Consequently it commutes with the canonical **left** group-algebra action,
and all left inner one-forms `sum a_i[D,b_i]` vanish.

However, the same structure theorem also exposes a missed third possibility.
The **right** regular algebra is itself fully `2I`-equivariant, and the
verifier finds an exact group element `s` and coefficient block for which

`[D,R_s] != 0`.

Thus an algebra with nonzero fluctuations need not break equivariance and
need not act only on `C^22`.  It can act by right convolution on `C[2I]`.
This refutes the proposed dichotomy; it does not construct a physical finite
spectral triple.  For the orbitwise inversion antiunitary, `JD` has neither
sign, so the new order-zero/first-order algebra does not pass the real
spectral-triple gate.

The inverse problem also has a negative correction: order zero and first
order do not define a unique greatest algebra from `D`.  They require a
specified `J`, representation, and ambient algebra.  For orbitwise inversion,
the full left and full right group algebras are two incomparable admissible
algebras.  No algebra containing both can satisfy order zero.  Therefore
there is no object honestly describable as **the** maximal algebra in the
unqualified inverse question.

All finite claims and the complete coefficient certificate are emitted by
`reproducible/verify_inner_fluctuations.py`.

## 1. Free cell layers and explicit regular-module chart

Vertices are the sorted list of the standard 120 unit quaternions used by the
existing Kähler--Dirac verifier.  A simplex is the increasing tuple of its
vertex-list indices.  With that reproducible convention, one representative
from every left `2I` orbit is:

| `alpha` | degree | representative |
|---:|---:|---|
| 0 | 0 | `(0,)` |
| 1 | 1 | `(0,1)` |
| 2 | 1 | `(0,2)` |
| 3 | 1 | `(0,3)` |
| 4 | 1 | `(0,4)` |
| 5 | 1 | `(0,5)` |
| 6 | 1 | `(0,6)` |
| 7 | 2 | `(0,1,2)` |
| 8 | 2 | `(0,1,3)` |
| 9 | 2 | `(0,1,4)` |
| 10 | 2 | `(0,1,5)` |
| 11 | 2 | `(0,1,6)` |
| 12 | 2 | `(0,2,3)` |
| 13 | 2 | `(0,2,4)` |
| 14 | 2 | `(0,2,7)` |
| 15 | 2 | `(0,3,5)` |
| 16 | 2 | `(0,4,8)` |
| 17 | 3 | `(0,1,2,3)` |
| 18 | 3 | `(0,1,2,4)` |
| 19 | 3 | `(0,1,3,5)` |
| 20 | 3 | `(0,1,4,6)` |
| 21 | 3 | `(0,1,5,6)` |

For every representative the verifier applies all 120 group elements,
including the orientation sign needed to return the image simplex to
increasing order.  Every target list has 120 distinct cells.  The resulting
orbits exhaust respectively `120`, `720`, `1200`, and `600` cells, proving

`C^0=Reg`, `C^1=6 Reg`, `C^2=10 Reg`, `C^3=5 Reg`.

Define

`F(delta_g tensor e_alpha) = L_g [representative alpha]`,

where the right side includes its sorted-simplex orientation sign.  The
exhaustive chart collision and coverage checks prove that `F` is an explicit
signed-basis isomorphism, not just a character comparison.  **DERIVED.**

## 2. Wedderburn blocks and their exact scope

In McKay-chain convention the irreducible dimensions are

`(1,2,3,4,5,6,4,2,3)`.

Exact character orthogonality and the square sum give

`C[2I] = M1(C) + M2(C) + M3(C) + M4(C) + M5(C)
          + M6(C) + M4(C) + M2(C) + M3(C)`,

equivalently

`M1 + M2^2 + M3^2 + M4^2 + M5 + M6`,

of complex dimension `120`.  **DERIVED Wedderburn type.**

The exact Frobenius--Schur indicators in the same order are

`(+1,-1,+1,-1,+1,-1,+1,-1,+1)`.

Thus the first two-dimensional spinor block has quaternionic real form, and
the displayed standard `H subset M2(C)` is available.  The `rho0`, `rho1`,
and `rho8` factors give a corner of real type

`C + H + M3(C)`.

In the McKay-chain convention fixed above the simultaneous nontrivial Galois
action is precisely

`rho1 <-> rho7`, `rho2 <-> rho8`,

and fixes the other five labels.  Thus the two triples are exactly
`(rho0,rho1,rho8)` and `(rho0,rho7,rho2)`.  There is one simultaneous
Galois choice of sheet, not four independent weak/color choices.
**DERIVED.**

Two limitations are binding:

- the three-factor corner has its central corner unit, not the identity of
  all nine Wedderburn factors;
- FS `-1` fixes the quaternionic isomorphism class, but a concrete matrix copy
  of `H` is unique only up to allowed unitary conjugacy.

Accordingly “canonical blocks” is justified; a fully allocated unital
Standard-Model representation on all of `C` is still **STRUCTURAL/OPEN**.

## 3. Structure theorem and the extracted `D`

Let `G` act by `L_g tensor I_m` on `C[G] tensor C^m`.  Write an operator as
blocks `T_alpha,beta:C[G]->C[G]`.  Equivariance says

`T_alpha,beta L_g = L_g T_alpha,beta`.

Such a block is determined by `T_alpha,beta(delta_e)`.  If

`T_alpha,beta(delta_e)=sum_h w_alpha,beta(h) delta_h`,

then equivariance gives

`T_alpha,beta(delta_g)=sum_h w_alpha,beta(h) delta_(g h)`.

This is right convolution under the convention
`R_h delta_g=delta_(g h)`.  Conversely every such right convolution commutes
with every `L_g`.  Therefore

`End_G(C[G] tensor C^m)
 = M_m(C[G]^op)`

and every equivariant operator has the unique asserted coefficient form.
**DERIVED structure theorem.**

For `m=22`, the verifier extracts `D` by applying it only to the 22 orbit
representatives.  It obtains 124 nonzero group coefficients in 112 nonzero
`(alpha,beta)` blocks.  It then reconstructs and compares all 14,880
nonzero signed incidences of `D`; equality is exact over the integers.

The complete machine-readable list is printed between
`D_COEFFICIENTS_BEGIN/END`.  In that list

`w[alpha,beta] h:c`

means that `c R_h tensor E_alpha,beta` occurs.  This executable table is the
explicit extraction; no diagonalization is used.

## 4. Inner fluctuations: vanishing result and refutation

### Canonical left placement

For the spanning set `{L_g:g in 2I}`, exact sparse integer matrices give

`[D,L_g]=0` for all 120 elements.

It follows by linearity that `[D,a]=0` for all `a in C[2I]`, and hence

`Omega_D^1(A_L)=0`, `D_A=D`

for `A_L=L(C[2I])` and every inner fluctuation.  This includes every
left-placed `C+H+M3` corner.  **DERIVED negative:** that placement supplies
no gauge fields or Yukawa couplings.

### The missed right-convolution horn

The structure theorem does not say equivariant operators commute with one
another.  It says they form the noncommutative algebra
`M_22(C[2I]^op)`.  The verifier compares the extracted coefficient maps

`D R_s: h -> s h`, `R_s D: h -> h s`

and finds the exact witness `s=1` in the deterministic group-index convention,
already on coefficient block `(alpha,beta)=(2,0)`.  Their integer
coefficient dictionaries differ, so `[D,R_s]!=0`.

`R_s` commutes with the full left `2I` action and acts on the group-algebra
factor, not only on `C^22`.  Therefore the proposed theorem

“nonzero fluctuation implies broken equivariance or multiplicity action”

is **REFUTED**.  A correct trichotomy must at least include right convolution.
More general elements of `M_22(C[G]^op)` mix both right-convolution and
multiplicity indices, so even “trichotomy” should not be promoted without an
additional naturality axiom.

The finite content of the original horns remains true but incomplete:

- breaking to `C3` enlarges the form-relevant candidate arena as already
  mapped, but the adversarial-audit wording remains binding: it is not a
  constructed real spectral triple;
- pure multiplicity actions can have nonzero commutator with the
  `22 by 22` coefficient matrix of `D`, but selecting their matrix carriers
  has `U(22)`-type freedom;
- right convolution is a third exact, equivariant possibility.

## 5. Inverse problem: no greatest algebra

The inverse question cannot be posed from `D` alone.  Order zero uses `J`,
and “maximal” also needs an ambient algebra and must distinguish a greatest
element from merely inclusion-maximal elements.

Use the explicit orbitwise inversion

`J_iota(delta_g tensor e_alpha)
 = delta_(g^-1) tensor e_alpha`

with coefficient conjugation.  It is an antiunitary signed permutation and
`J_iota^2=+1`.  Then both

`A_L=L(C[2I]) tensor I_22`,

`A_R=R(C[2I]) tensor I_22`

satisfy order zero and first order:

- for `A_L`, `J A_L J^-1=A_R` and `[D,A_L]=0`;
- for `A_R`, `J A_R J^-1=A_L`, while
  `[[D,R_a],L_b]=0` because both `D` and `R_a` are left-equivariant.

Both have exact Wedderburn type

`M1 + M2^2 + M3^2 + M4^2 + M5 + M6`.

They are incomparable.  If a greatest admissible algebra contained both,
order zero applied to `L_s` and `R_t` would require

`[L_s,J R_t J^-1]=[L_s,L_(t^-1)]=0`,

which is false for the exact noncommuting group pair exhibited by the
verifier.  **DERIVED:** there is no greatest order-zero/first-order algebra
for this `J` in `End(C)`.  Consequently no exact Wedderburn type can be
assigned to “the maximal algebra”; the premise of Observation A is false.

Relative answers are exact:

| specified ansatz | largest algebra in that ansatz | SM corner? | fluctuations? |
|---|---|---:|---:|
| left group algebra | `A_L`, Wedderburn type above | yes, with stated corner/Galois scope | zero |
| `R(C[G])` times coordinate-diagonal multiplicity algebra | `A_R`, same type | yes, same scope | nonzero |

The second row is not tautological.  The extracted coefficient-support graph
on the 22 orbit labels is connected.  A unital star-subalgebra of the
coordinate-diagonal `C^22` is a partition algebra.  If one of its blocks
separates the endpoints `alpha,beta` of a nonzero coefficient of `D`, its
characteristic projector gives a nonzero first-order double commutator.
First order therefore forbids every cut edge.  Connectivity forces the
one-block partition, namely scalar multiplicity action.  Hence `A_R` is
exactly maximal in this explicit equivariant tensor/diagonal ansatz.
**DERIVED direct ansatz solve.**

This does not classify every non-diagonal star-subalgebra of the
58,080-dimensional `End_G(C)`.  Claiming that larger classification from the
partition solve would repeat the overstatement errors prohibited by the
adversarial audit.

The right result is a genuine algebraic reopening of the order-one screen,
but not yet of the physical matter program: the next real-structure gate
fails.

## 6. Derived/structural `J` inventory

The table distinguishes candidates defined on this exact primal arena from
objects defined only on another arena.

| candidate and arena | KO signs `(J^2,JD,Jgamma_form)` | order zero for left group algebra | fluctuations | status |
|---|---|---|---|---|
| coefficient conjugation `K` on primal `C` | `(+,+,+)` | **fails** for noncommuting blocks; it maps left to left | zero because `[D,A_L]=0` | **DERIVED** |
| central antipode `L_z K` | `(+,+,+)` | **fails** for the same reason | zero | **DERIVED** |
| orbitwise inversion `J_iota` after the 22 origins above | `(+`, neither `+` nor `-`, `+)` | holds; left maps to right | zero for `A_L`; nonzero for `A_R` | **DERIVED conditional on the explicit orbit origins** |
| `L_z J_iota` | `(+`, neither `+` nor `-`, `+)` | holds | same left/right result | **DERIVED conditional** |
| primal Hodge-star times conjugation | undefined: `120!=600`, `720!=1200` | undefined | undefined | **DERIVED nonexistence on primal `C`** |
| primal-to-dual Hodge star | signs not computed without a specified dual Hilbert space, metric, and dual `D` | open | open | **STRUCTURAL/OPEN** |
| Galois sheet swap on `W+W^sigma` | conditional KO6 `(+,+,-)` for doubled McKay adjacency | SM order zero undefined because no derived SM representation exists | undefined | **DERIVED sign algebra on a different arena; STRUCTURAL intertwiners** |
| signed/other-chirality Galois compositions | exact node-level alternatives include `J^2=-`, `JD=-`, or `Jgamma=+` | same missing representation | undefined | **DERIVED conditional on doubled node data** |
| Galois/Hodge compositions on primal `C` | undefined | undefined | undefined | **OPEN; constituent map absent** |

For `J_iota`, exact sparse comparison gives

`nnz(JD-DJ)=11000`, `nnz(JD+DJ)=20380`.

Thus the only defined primal candidate that makes the nonzero right
fluctuation pass order zero and first order has no `JD=epsilon' DJ` sign.
No combination in the derived inventory simultaneously supplies:

1. a defined antiunitary on the 2640 primal arena,
2. the required real-structure sign with this `D`,
3. order zero and first order, and
4. nonzero inner fluctuations.

**DERIVED/OPEN verdict:** the right-convolution result reopens the algebraic
order-one route, but no real spectral triple or physical matter sector has
been constructed.  `Y`, anomaly forcing, color orientation, Yukawa data, and
the multiplet census remain **OPEN**.

## Status ledger

### Strengthened

- **DERIVED:** every cell layer is a free `2I` set, with the 22 explicit
  representatives above.
- **DERIVED:** an explicit signed-basis isomorphism
  `C=C[2I] tensor C^22`.
- **DERIVED:** the general equivariant-operator/right-convolution theorem.
- **DERIVED:** all 124 integer coefficients of `D` and exact reconstruction
  of its 14,880 incidences.
- **DERIVED:** exact Wedderburn and FS tables; the weak/color alternatives
  are one simultaneous Galois flip.
- **DERIVED:** canonical left fluctuations vanish.
- **DERIVED discovery:** equivariant right fluctuations are nonzero and pass
  order zero/first order for orbitwise inversion.

### Downgraded / refuted

- **REFUTED:** nonzero fluctuations force either equivariance breaking or
  multiplicity-only action.
- **REFUTED:** `D` canonically determines a unique greatest order-one
  algebra.  Left and right admissible algebras are incomparable.
- **DELIMITED:** the `C+H+M3` factors form a canonical-type corner up to the
  simultaneous Galois flip and quaternionic basis conjugacy, not a fully
  allocated unital representation on all nine factors.
- **NEGATIVE:** orbitwise inversion has no `JD=+/-DJ` sign for this `D`.

### Open

- a derived antiunitary and arena satisfying all real spectral-triple signs
  while retaining the nonzero right fluctuation;
- whether `D` or primal-dual geometry canonically aligns the 22 right torsor
  actions without chosen orbit origins;
- classification of inclusion-maximal order-zero/first-order algebras after
  an ambient algebra and valid `J` are fixed;
- a derived matter representation, `Y`, Route-C anomaly input, color
  orientation, generations, gauge fields, and Yukawa couplings.
