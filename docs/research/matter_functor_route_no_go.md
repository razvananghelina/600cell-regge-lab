# Derived chiral matter functor: route audit and exact target

Date: 2026-07-22 (session 3)

## Decision

No derived chiral matter functor is presently obtained.  Route A stops before
the first-order condition is defined; Route B admits the Standard Model
residue tuple but does not select it and cannot recover its integer lift.
Route C nevertheless turns the missing functor into an exact equivariance and
Diophantine target.  In particular, the anomaly equations almost determine
the M15 hypercharges once the nonabelian module and a nonzero primitive charge
unit are supplied.

This note never uses anomaly cancellation to relabel an inserted tuple as
derived.  The distinction is:

- anomaly equations can constrain a candidate integer operator;
- the repository still lacks the construction that produces that operator on
  its finite matter space.

All finite claims below are checked by
`reproducible/verify_matter_functor.py`.

## Route A: finite real spectral triple

### Data actually available

Let

`H_McKay = direct_sum_i V_i`, with dimensions
`(1,2,2,3,3,4,4,5,6)`.

- **DERIVED:** `dim H_McKay=30`.
- **DERIVED:** the affine-E8 graph is bipartite and supplies
  `gamma_F=(-1)^(2j)` with `gamma_F^2=1`; its eigenspace dimensions are
  `(dim H_+,dim H_-)=(16,14)`.
- **DERIVED:** every McKay-edge Dirac topology anticommutes with this grading.
- **DERIVED:** Galois conjugation exchanges the paired 2-dimensional irreps
  and paired 3-dimensional irreps and fixes the other nodes.  This node
  permutation preserves dimensions and commutes with `gamma_F`.
- **DERIVED:** the gauge kernel and bracket results give the compact Lie
  algebra skeleton `u(1)+su(2)+su(3)`.

The 30-dimensional sum is numerically compatible with particle/antiparticle
doubled M15, but dimension equality is not an identification.

- **PATTERN:** `30=2*15`.
- **OPEN:** no equivariant isomorphism identifies this McKay sum with doubled
  M15, and its `16+14` grading is not the `15+15` exchange split that such a
  numerical reading suggests.

### Exact KO6 obstruction on the present grading

For the usual finite Standard-Model KO-dimension-6 sign,
`J gamma_F=-gamma_F J`.  Any invertible `J` with this relation maps `H_+`
bijectively to `H_-`.  Their dimensions must therefore agree.

- **DERIVED (negative):** `16 != 14`, so no invertible antiunitary—or even
  invertible linear map—can anticommute with the present `gamma_F` on
  `H_McKay`.
- **OPEN:** another KO sign or an enlarged Hilbert space can evade this
  dimension obstruction, but neither is selected by the derived data.

### Why Galois conjugation is not yet `J`

The arithmetic automorphism `sigma(phi)=phi'` acts on coefficient fields and
permutes some irreducible labels.  A real structure is instead an antiunitary
endomorphism of one specified complex Hilbert space, with a specified positive
inner product, square, KO signs, and opposite-algebra action

`rho^o(a)=J rho(a*) J^(-1)`.

- **DERIVED:** the node permutation is an involution commuting with the
  bipartite grading.
- **STRUCTURAL:** composing an extension of `sigma` with ordinary complex
  conjugation is a plausible candidate mechanism.
- **OPEN:** the extension to all representation matrices and CG blocks, its
  antiunitarity for a derived Hilbert metric, its KO signs, and its action on a
  particle/antiparticle doubling have not been constructed.

### The first-order system is undefined, not inconsistent

The Connes conditions require matrices for

`rho(A_F)`, `rho^o(A_F)`, `D_F`, `J`, and `gamma_F`,

including

`[rho(a),rho^o(b)]=0`,
`[[D_F,rho(a)],rho^o(b)]=0`,

and unimodularity on the resulting representation.

The Lie-algebra dimensions `1+3+8` do not by themselves construct the real
associative algebra

`A_F=C+H+M_3(C)`

or its representation.  Distinct associative representations can have the
same infinitesimal gauge algebra and different central charges.

- **DERIVED (negative inventory statement):** `gamma_F` and the McKay Dirac
  topology exist, but `rho(A_F)`, the opposite action, and an antiunitary
  endomorphism `J` do not exist as derived matrices on the same space.
- **OPEN:** consequently there is no finite first-order equation to solve for
  `Y`; writing the standard NCG representation would insert the desired
  matter content.

**Minimal missing axiom for Route A:** a faithful real even representation of
`C+H+M_3(C)` on a derived particle/antiparticle Hilbert space, including its
positive metric and an antiunitary `J` satisfying specified KO signs.  Once
this is given, the already-available Dirac topology and grading can be tested
against order-zero, first-order, and unimodularity conditions.

## Route B: section winding on the Hopf C10 fiber

This route is distinct from applying `n mod 10` to mass slots.  Assign one
vertical character to each nonabelian irreducible matter block.  In the order

`(Q,u^c,d^c,L,e^c,nu^c)`,

the external SM benchmark has

`6Y=(1,-4,2,-3,6,0)`

and hence C10 residues

`(1,6,2,7,6,0)`.

- **DERIVED:** these reductions modulo 10 are exact.
- **DERIVED:** one character per irreducible block is automatically scalar on
  each weak doublet and hence passes the commutant gate.
- **DERIVED:** with no block-selection rule, six C10 characters leave
  `10^6` commutant-compatible assignments.
- **PATTERN:** the displayed SM residue tuple is consistent but is not selected
  by fiber geometry, McKay position, or a derived matter-section functor.

There is also an information-theoretic obstruction.  C10 sees only a residue:

`y` and `y+10m` define the same character.

In particular, `u^c` and `e^c` both have residue 6 although their integer
sixth-hypercharges are `-4` and `6`.

- **DERIVED (negative):** C10 winding alone cannot recover integer
  hypercharge or its normalization.
- **OPEN:** a derived lift from C10 characters to primitive integers, and a
  derived assignment of section sectors to the six chiral blocks.

**Minimal missing axiom for Route B:** an equivariant matter-section
construction with a canonical integral lift of vertical characters.  Merely
choosing the six benchmark residues is assignment by hand.

## Route C: exact necessary specification

Write primitive integer charges `y=6Y` on the left-handed blocks as

`(q,u,d,l,e,n)` for `(Q,u^c,d^c,L,e^c,nu^c)`.

Any successful functor must satisfy all of the following.

### Equivariance and representation conditions

1. **DERIVED necessary:** dimension 15 or 16 with SM nonabelian multiplicities
   `(1,2,1,1)` or `(1,2,1,2)` for block dimensions `(6,3,2,1)`.
2. **DERIVED necessary:** `Y` lies in the commutant of
   `su(2)+su(3)`, hence is scalar on every irreducible block and especially on
   both components of each weak doublet.
3. **DERIVED necessary:** three generations carry the same six eigenvalues.
4. **OPEN selection:** the functor must distinguish `3` from `bar(3)` and
   select which color singlet is `u^c` versus `d^c`; real A5 restrictions do
   neither.
5. **DERIVED necessary:** the weighted weak-doublet count is even.  For the
   SM inventory it is `3+1=4`, so Witten parity vanishes.

### Exact anomaly system

The four local conditions are

`6q+3u+3d+2l+e+n=0`,

`6q^3+3u^3+3d^3+2l^3+e^3+n^3=0`,

`3q+l=0`, `2q+u+d=0`.

- **DERIVED:** `(1,-4,2,-3,6,0)` is primitive and solves all four.

For M15, omit `n`.  Eliminating the three linear equations gives

`l=-3q`, `d=-2q-u`, `e=6q`,

and the cubic anomaly factors exactly as

`18 q (2q-u)(4q+u)=0`.

- **DERIVED partial theorem:** for nonzero `q`, anomaly cancellation forces
  `u=2q` or `u=-4q`.  After primitive normalization and sign orientation,
  the two solutions are

  `(1,-4,2,-3,6)` and `(1,2,-4,-3,6)`.

  Thus anomalies recover the SM tuple only up to exchange of the two colored
  singlets.  The exchange is not resolved by their identical nonabelian
  representations.

For M16 the linear equations give `e=6q-n`, while the cubic factors as

`18 q (-n+2q-u)(-n+4q+u)=0`.

- **DERIVED:** anomaly cancellation leaves an integer `n` parameter.  Setting
  `n=0` recovers the M15 branches.
- **OPEN:** neutrality of the extra singlet is not forced by anomalies alone.

### Trace-index target

6. **DERIVED necessary physical check:** after a successful construction,
   ordinary hypercharge must give `(T1,T2,T3)=(10/3,2,2)`, ratio `5:3:3`.
7. **PATTERN only:** `8:5:2` is recorded as the old gauge-prefactor pattern
   and is not a matter-functor target.

## Final obstruction

## Follow-up finite test: is the McKay 16 equal to diagonal M16?

### Exact halves and flavor slots

In the one-based irrep convention of `verify_mckay_chirality.py`, the derived
bipartition is

`H_even = rho_1 + rho_4 + rho_5 + rho_6 + rho_8`

with dimensions `1+3+3+4+5=16`, and

`H_odd = rho_2 + rho_3 + rho_7 + rho_9`

with dimensions `2+2+4+6=14`.  More explicitly:

| half | 2I irrep | construction | dim | assigned charged slot |
|---|---|---|---:|---|
| even | `rho_1` | `Sym^0(std)` | 1 | `e` |
| even | `rho_4` | `Sym^2(std)` | 3 | `d` |
| even | `rho_5` | `Sym^2(std')` | 3 | `b` |
| even | `rho_6` | `std tensor std'` | 4 | `tau` |
| even | `rho_8` | `Sym^4(std)` | 5 | `mu` |
| odd | `rho_2` | `std` | 2 | `u` |
| odd | `rho_3` | `std'` | 2 | `t` |
| odd | `rho_7` | `Sym^3(std)` | 4 | `s` |
| odd | `rho_9` | `Sym^5(std)` | 6 | `c` |

- **DERIVED:** the node split, dimensions, and charged-slot placement under
  the repository's existing Wilson-line assignment.
- **DERIVED (negative):** neither half is generation-local: the even half
  contains five charged species from all three generations, while the odd
  half contains the remaining four.  No neutrino slot occurs in this node
  assignment.

For class order

`(1A,2A,4A,6A,3A,10A,5A,5B,10B)`,

exact character addition gives

`chi_even=(16,16,0,1,1,1,1,1,1)`,

`chi_odd=(14,-14,0,1,-1,1,-1,-1,1)`.

In particular, the central element `-1` acts as `+I_16` on the entire even
half and `-I_14` on the entire odd half.

### Correct diagonal restriction of M16

The derived finite diagonal homomorphism relevant to the proposed test is

`g in 2I -> (rho_5(g),rho_2(g)) in SU(3) x SU(2)`.

Here `rho_5=3'` factors through `2I/{+/-1}=A5` in the repository convention,
and `rho_2=std` is the defining weak spinor.  Therefore an external product
`(R_3,R_2)` restricts by the tensor product of these two pulled-back actions.
This justifies, rather than assumes, the rule in the question.  The real A5
module `3'` is also the restriction of `bar(3)`, so this finite action still
cannot select color orientation.

The McKay edge `rho_5--rho_9` states exactly

`rho_5 tensor rho_2 = rho_9`.

Consequently

`M16|_2I = rho_9 + 2 rho_5 + rho_2 + 2 rho_1`.

Its exact character is

`chi_M16=(16,0,0,3,1,`
`(5-sqrt(5))/2,(7-3sqrt(5))/2,`
`(7+3sqrt(5))/2,(5+sqrt(5))/2)`.

This row is obtained independently both from its irreducible decomposition
and from

`chi_3' chi_2 + 2 chi_3' + chi_2 + 2`.

- **DERIVED (negative):** `H_even` is not isomorphic to diagonal M16.  The
  shortest certificate is already the central character:
  `chi_even(-1)=16` but `chi_M16(-1)=0`.
- **DERIVED:** the equality `16=dim(M16)` is therefore only a dimension
  coincidence for this derived diagonal action, not a generation module.

### The odd 14

Two natural 14-dimensional comparisons also fail:

1. diagonal M15 with its remaining singlet removed has
   `rho_9+2rho_5+rho_2` and central character `-2`;
2. derived gauge content `rho_1+rho_4+rho_5+rho_8` plus the defining spinor
   `rho_2` has central character `10`.

Neither equals `H_odd`, whose central character is `-14` and whose exact
decomposition is `rho_2+rho_3+rho_7+rho_9`.

- **DERIVED (negative, scoped):** the odd half matches neither tested natural
  candidate.
- **OPEN:** no physical matter or gauge interpretation of the odd half is
  derived.  This is not an exhaustive classification of all abstract
  14-dimensional `2I` representations.

### Interpretation of the finite test

- **Strengthened / DERIVED:** the character comparison closes the tempting
  `16=M16` route negatively.
- **PATTERN rejected:** dimension equality alone is not evidence for an SO(10)
  spinor or one-generation module here.
- **OPEN unchanged:** another functor could use a different derived space or
  action, but it must still construct color orientation and the abelian
  generator before the Route C anomaly theorem applies.

## Final obstruction

- **Strengthened / DERIVED:** the present McKay grading has a KO6 dimension
  obstruction; C10 winding loses the integral lift; the M15 anomaly variety
  factorizes and fixes charges up to scale and the `u/d` exchange.
- **Downgraded / PATTERN:** the SM C10 winding tuple is a consistent pattern,
  not a derived assignment.
- **OPEN:** the central object remains a derived associative-algebra
  representation or matter-section functor that selects chiral blocks,
  color conjugation, integer lift, orientation, and generation-blind `Y`.

The sharpest next target is not another charge search.  It is construction of
one of the two missing maps:

`(McKay/Hopf data) -> Rep_real_even(C+H+M3(C))`

or

`(McKay/Hopf data) -> C10-equivariant chiral matter sections with an integral lift`.

Either output can be inserted into the exact specification above and accepted
or rejected without normalization fitting.
