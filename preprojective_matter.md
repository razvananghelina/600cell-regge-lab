# Preprojective smooth fiber: exact boundary no-go

Date: 2026-07-22

## Decision

The regular 120-dimensional fiber is genuinely the first matter candidate in
this program whose equivariant multiplicity algebra is large enough to contain
weak and color matrix blocks.  The hoped-for canonical finite triple does not,
however, follow from it.

Two independent exact constraints kill the canonical construction:

1. inversion on `C[2I]` has `J^2=+1` but **commutes** with the proposed
   integer/spinor isotypic grading.  Galois composition does not change this.
   Thus its sign is `J gamma=+gamma J`, not the KO6 target
   `J gamma=-gamma J`;
2. for the canonical maximal multiplicity algebra
   `B=End_2I(C[2I])`, the regular `B-B` bimodule has only diagonal Krajewski
   vertices `(i,i)`.  Distinct vertices share neither index, so first order
   permits no off-diagonal block.  Oddness then forces `D=0`.

There is a second scope correction.  The factors `M2(C)` and `M3(C)` prove
that quaternionic and color actions are *available*, but do not give a
canonical unital inclusion of `C+H+M3(C)` on the entire regular module.
Putting `H` in one `M2` factor and `M3(C)` in one `M3` factor leaves seven
central factors unassigned.  Extending the unit requires an allocation rule
not supplied by the geometry.

**DERIVED boundary theorem (for the three canonical routes audited):**
Standard-Model matter is not derivable from the McKay node module, the
bidirected edge module, or the canonical smooth preprojective fiber module of
this geometry.  The three certificates are respectively missing `2I` fusion
support / no `M3` in the node commutant; first order plus orientation oddness;
and the wrong canonical KO grading sign together with a zero first-order odd
Dirac for the maximal multiplicity action.

This does not classify every abstract module over the affine-E8
preprojective algebra.  It closes the natural smooth-fiber construction with
the algebra, grading, real structure, and operators specified in the mission.
All finite claims are checked by
`reproducible/verify_preprojective_matter.py`.

**2026-07-27 scope correction (DERIVED).**  “Three canonical routes” and
“real structures specified in the mission” are binding restrictions.  The
later Q8 factor-swap example shows that inversion/Galois candidates do not
classify arbitrary multiplicity-mixing antiunitaries on a free regular
arena.  No such universal no-go is claimed here.

## 1. What the smooth fiber actually supplies

For the affine McKay quiver, the standard preprojective/skew-group
correspondence identifies the relevant affine-ADE category with modules over
`C[x,y]#2I` (up to the usual Morita choice).  Only its following finite
consequence is used here.  At a point with trivial stabilizer, the orbit has
120 points.  Functions on that orbit have delta basis `{delta_g}`, and the
group permutes it freely and transitively.  Therefore the fiber is the regular
module `C[2I]`, of dimension 120.  **DERIVED finite consequence.**

Peter-Weyl gives

`C[2I] = direct_sum_i V_i tensor C^(n_i)`,

where, in McKay-chain order `1,2,3,4s,5,6,4,2',3'`,

`(n_i)=(1,2,3,4,5,6,4,2,3)`.

Thus

`End_2I(C[2I]) = M1+M2+M3+M4+M5+M6+M4+M2+M3`.

The isotypic dimensions are `n_i^2` and sum to 120.  Integer and spinor
isotypic dimensions each sum to 60.  **DERIVED.**

### Frobenius--Schur table

The exact formula `nu(V)=|G|^-1 sum_g chi_V(g^2)` gives:

| irrep | dimension | parity | FS indicator | real Schur type |
|---|---:|---|---:|---|
| `1` | 1 | integer | `+1` | real |
| `2` | 2 | spinor | `-1` | quaternionic |
| `3` | 3 | integer | `+1` | real |
| `4s` | 4 | spinor | `-1` | quaternionic |
| `5` | 5 | integer | `+1` | real |
| `6` | 6 | spinor | `-1` | quaternionic |
| `4` | 4 | integer | `+1` | real |
| `2'` | 2 | spinor | `-1` | quaternionic |
| `3'` | 3 | integer | `+1` | real |

Hence the real commutant of either irreducible weak spinor is isomorphic to
`H`.  **DERIVED isomorphism class.**  An invariant quaternionic structure is
unique only up to its allowed scalar/unitary conjugacy, so a concrete copy
`H subset M2(C)` is not a distinguished matrix embedding without a basis
choice.  **STRUCTURAL choice.**

There are precisely two dimension-two factors (`2,2'`) and two
dimension-three factors (`3,3'`).  Galois exchanges the two weak candidates
and exchanges the two color candidates.  Thus the residual *factor choices*
are elegantly two single Galois flips.  **DERIVED.**  This elegance does not
remove the additional unital allocation of the other seven factors, nor does
`3 <-> 3'` become color conjugation; both three-dimensional representations
are real.  **OPEN selector / DERIVED limitation.**

## 2. Exact signs of canonical inversion

In the group basis define

`J(sum c_g delta_g)=sum conjugate(c_g) delta_(g^-1)`.

It is antiunitary, exchanges the left and right regular actions, and satisfies
`J^2=+1`.  In a Peter-Weyl matrix-coefficient basis it sends `(i,a,b)` to
`(i,b,a)` anti-linearly.  All nine irreps are self-dual, so it leaves the
irrep label—and therefore integer/spinor parity—fixed.  Consequently

`J gamma=+gamma J`.

The exact canonical sign inventory is therefore

`J^2=+1`, `Jgamma=+gamma J`,

with the `JD` sign depending on the chosen operator.  It is not a KO6 table.
The outer/Galois involution preserves the bipartition, so composing it with
`J` still gives the plus grading sign.  **DERIVED decisive negative.**

Obtaining KO6 would require a separately chosen parity-reversing doubling or
antiunitary.  That is the same kind of extra sheet already classified as
STRUCTURAL/OPEN in `galois_doubling_triple.md`; it is not canonical inversion
on the 120 vertices.

## 3. Dirac candidates

The two proposed operator families behave differently, but neither completes
the triple.

### Preprojective coordinate arrows

Multiplication by the two coordinate functions on a free orbit realizes the
doubled-quiver incidence: tensoring by the defining spinor connects precisely
the affine-E8 adjacent isotypic sectors.  Every such block reverses
integer/spinor parity.  **DERIVED:** the edge-space orientation obstruction
does not recur verbatim; these are maps between regular-module isotypic
multiplicity spaces, not maps between distinct oriented Hom blocks.

Nevertheless, the canonical `J` already has the wrong KO6 grading sign.  For
the maximal canonical algebra `B`, the Krajewski support is diagonal `(i,i)`.
A first-order block between `(i,i)` and `(j,j)` with `i!=j` shares neither a
left nor a right index.  Hence all adjacency-shaped odd blocks are illegal and
first order plus oddness gives `D=0`.  **DERIVED for the maximal canonical
action.**

A smaller `C+H+M3(C)` restriction can change Krajewski legality, but no
canonical unital restriction exists: one must decide how all nine central
factors are assigned to the three algebra summands and choose concrete real
embeddings.  Consequently its first-order equations are not uniquely defined,
rather than exhaustively failed.  Testing one hand-picked allocation would
insert the missing matter bimodule.  **OPEN outside the canonical route.**

### Existing 120-vertex operators

The 600-cell adjacency/Laplacian and every left-equivariant convolution
operator preserve left isotypic components.  Since `gamma` is the action of
the central element `-1`, these operators commute with `gamma`; the existing
`Box` constructions likewise preserve the antipodal parity.  They are
gamma-even, not finite Dirac candidates for this grading.  **DERIVED scope
negative.**  The coordinate arrows are odd but fail the canonical first-order
gate above; the old vertex operators pass parity preservation in the wrong
direction.

Thus the edge no-go is *evaded as an orientation argument* but replaced by
two fiber-specific constraints: the inversion KO sign and diagonal-bimodule
first order.

## 4. Why Route C and the generation census do not start

No canonical KO6 real even SM bimodule with a nonzero first-order Dirac
survives.  Therefore there is no derived commutant generator `Y` on which to
impose doublet constancy and generation blindness, and Route C cannot be
applied.  **DERIVED scope limit.**

The numbers 2 and 3 in the regular decomposition are multiplicities of an
irrep in the regular representation—equivalently dimensions of its dual
right-regular multiplicity space.  After choosing an algebra action they are
the carrier dimensions of candidate weak/color matrix factors.  They do not
count copies of a 15- or 16-dimensional chiral matter module and therefore do
not count generations.  **DERIVED interpretation.**

The 120 dimensions admit numerical factorizations `8*15` and do not equal an
integer multiple of 16, but the actual decomposition is the nine isotypic
blocks of dimensions

`1,4,9,16,25,36,16,4,9`.

It is not a sum of uniform M15 or M16 blocks.  Reading eight generations from
`120=8*15`, or three generations from an `M3` multiplicity, would be a forced
pattern.  **PATTERN rejected.**

## Status ledger

### Strengthened

- **DERIVED:** the smooth free-orbit fiber is the 120-dimensional regular
  representation and has the exact matrix-factor commutant above.
- **DERIVED:** all nine FS indicators; every spinor is quaternionic.
- **DERIVED:** both weak and color residual factor choices are single Galois
  flips.
- **DERIVED negative:** canonical inversion commutes with isotypic gamma.
- **DERIVED negative:** maximal multiplicity first order plus oddness forces
  `D=0`.
- **DERIVED:** preprojective arrows evade edge orientation parity but not the
  fiber-specific constraints.

### Downgraded

- `H` and `M3(C)` occur as available real/complex factor types, not as a
  canonical unital Standard-Model algebra action on all of `C[2I]`.
- Multiplicity dimensions 2 and 3 are carrier dimensions, not generation
  numbers.
- The canonical sign table is KO0-like in its `J-gamma` sign, not KO6.

### Open

- a geometrically selected unital allocation of all nine factors to
  `C+H+M3(C)`;
- a parity-reversing antiunitary not inserted as an extra sheet;
- a nonzero first-order odd Dirac for such a selected smaller bimodule;
- color orientation, generation blocks, and a generation-blind `Y`.

Those are new axioms/structures, not information hidden in the canonical
smooth fiber.
