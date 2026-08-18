# The final canonical bimodule arena

Date: 2026-07-27 (fifth session)

## Decision

The authorized final enlargement does not open the matter gate.

There is first a necessary algebra correction.  For the non-simple
semisimple algebra

`A=C[2I]=direct_sum_i M_(d_i)(C)`,

the abstract enveloping algebra `A tensor A^op` has dimension 14,400, but it
is **not**

`End_C(A)=M_120(C)`.

The natural two-sided multiplication map

`a tensor b^op : x -> axb`

has image dimension

`sum_i d_i^4=2628`

and kernel dimension `14400-2628=11772`.  The identification asserted in the
mission is therefore **REFUTED**.  It would be true for a central simple
matrix algebra, not for a nine-block group algebra.

The abstract `A tensor A^op` Hilbert bimodule is nevertheless canonical and
is the arena audited below.  Calling it “universal” is justified in the
enveloping-algebra sense; calling it “the most economical Hilbert bimodule”
is not.  The regular bimodule `A` is smaller, and finite spectral triples
permit many bimodule multiplicities.  Thus the motivation is
**STRUCTURAL**, not a uniqueness theorem.

For every derived Dirac candidate tested here, order zero and first order
hold.  The gate fails because:

1. no derived grading anticommutes with the candidate;
2. every candidate obtained from a fully `2I`-equivariant vertex operator
   commutes with the represented left algebra, so all inner one-forms vanish.

The stopping rule is now binding: no larger arena is proposed.

## 1. Bimodule, real structure, and adjoint action

Use the Hilbert-space basis

`delta_x tensor delta_y^op`, `x,y in 2I`,

with the product inner product.  The represented algebra acts on the first
factor, and its opposite acts on the second:

`pi(a)(u tensor v^op)=au tensor v^op`,

`pi^o(b)(u tensor v^op)=u tensor (vb)^op`.

Associativity gives

`[pi(a),pi^o(b)]=0`

exactly.  Define

`J(u tensor v^op)=v* tensor (u*)^op`

antilinearly.  In the group basis this is flip plus inversion:

`J(delta_x tensor delta_y^op)
 =delta_(y^-1) tensor delta_(x^-1)^op`.

It is antiunitary and

`J^2=+1`.

It exchanges the two represented factors, giving the stated opposite
action.  **DERIVED.**

### Diagonal adjoint decomposition

The diagonal adjoint action used here is conjugation on both group-algebra
factors:

`g:(u tensor v^op) -> gug^-1 tensor (gvg^-1)^op`.

For one group-algebra factor its character is `|C_G(g)|`.  Exact character
inner products give, in McKay order,

`A_Ad = (9,0,7,0,9,0,6,0,7)`,

`(A tensor A^op)_Ad
 = (296,0,736,0,1192,0,932,0,736)`.

The latter dimension closes as

`296*1+736*3+1192*5+932*4+736*3=14400`.

The central element `z=-1` satisfies `zgz^-1=g`; hence it acts trivially
under this adjoint action.  Only the five central-even irreps occur.

The central-parity segregation theorem assumes an odd operator between
source and target sectors with opposite `z` characters.  Here every sector
has character `+1`, so there is no integer/spinor split to segregate.
Precisely that hypothesis fails.  This does not produce an odd Dirac; it
only makes the segregation theorem inapplicable.

## 2. Derived Dirac candidates

### The 600-cell adjacency element

Let `e` be the identity vertex.  Its twelve adjacent vertices are exactly
the group elements with quaternionic trace `phi`.  They form one
inverse-closed conjugacy class.  Therefore

`c=sum_(g adjacent to e) g`

is self-adjoint and central in `C[2I]`.  Right convolution by `c` is exactly
the certified 600-cell vertex adjacency.  The vertex Laplacian is

`Delta_0=12-c`.

No rediagonalization is used.

On the tensor arena define

`D_c^- = L_c tensor 1 - 1 tensor R_c`,

`D_c^+ = L_c tensor 1 + 1 tensor R_c`.

Both are self-adjoint and nonzero: on `1 tensor 1` their two supports occupy
different tensor coordinates.  Flip-star gives

`J D_c^- J^-1=-D_c^-`,

`J D_c^+ J^-1=+D_c^+`.

Thus the plus candidate has the KO6 target `JD=DJ`; the minus candidate does
not.

Because `c` is central,

`[D_c^±,pi(a)]=0`

for every `a`.  First order follows, but

`Omega_D^1(A)=0`.

The “commutator/adjoint derivation” language must be qualified: as a
two-factor operator `D_c^-` is nonzero, but the internal derivation
`[c,x]` on `A` itself is identically zero.

### Class sums and Laplacian

Every derived class sum lies in `Z(C[2I])`.  Its inverse-symmetrization is
self-adjoint.  The same plus/minus conclusions hold:

- plus: `JD=+DJ`;
- minus: `JD=-DJ`;
- both: order zero and first order;
- both: zero inner one-forms.

Polynomials in the adjacency and the vertex Laplacian remain central, so
they add no new axiom row.  **DERIVED family reduction.**

### Hopf/Box vertex operator

The certified Hopf partition consists of twelve right `C10` cosets.  The
vertex `Box_0` is built from fiber and cross-edge Laplacians.  Left `2I`
multiplication preserves the right-coset partition, so `Box_0` is a
`2I`-equivariant vertex operator.  By the regular-module structure theorem
it is right convolution, although unlike `c` its coefficient need not be
central.

Consequently its canonical two-factor plus/minus lifts are self-adjoint,
nonzero, and have the same flip-star `JD` signs.  But they still commute with
the represented left algebra:

`[Box_0 tensor 1,pi(a)]=0`.

Terms on the second factor commute automatically.  Hence their inner
one-forms also vanish.  This is a **DERIVED structural negative**, not a
spectral fit.

The edge/face/cell Box operators act on dimensions 720, 1200, and 600 rather
than on the 120-dimensional vertex factor.  Moving them to this tensor arena
requires an additional carrier map.  No such map is derived, so they are not
candidate operators here.

No further Hopf or McKay operator in the repository supplies a
self-adjoint, odd endomorphism of this exact bimodule without adding a
carrier choice.  **OPEN outside the derived candidate list.**

## 3. Grading search

### Adjoint central parity

The adjoint action of `-1` is the identity.  It is not a nontrivial grading,
has `J gamma=+gamma J`, and commutes with every candidate `D`.

### McKay parity

On a regular group-algebra factor, McKay integer/spinor parity is represented
by the central involution.  Since it is central, it commutes with all
convolution candidates.  Its natural same-factor and two-factor extensions
therefore satisfy `[gamma,D]=0`, not `{gamma,D}=0`.  Flip-star either
exchanges the one-factor choices or commutes with the symmetric two-factor
choice; it does not supply the required derived KO6 row.

### Vertex/Hopf sign grading

The 600-cell graph contains triangles; the verifier prints an explicit one.
It is not bipartite.  Hence no diagonal vertex-sign involution makes the
adjacency odd.  The Hopf fiber cycles alone are even, but `Box_0` has
nonzero weights on the full fiber-plus-cross edge support, which contains
the same odd cycles.  The fiber-alternating sign therefore does not make
`Box_0` odd.

This is the expected perfect-group obstruction in concrete form.  A
generator-sign grading would define a nontrivial homomorphism
`2I -> Z2`, but `2I` is perfect and has none.

### Form parity

Form parity is defined on the 2640-dimensional primal cochain complex and
its 5280-dimensional primal--dual double.  The repository supplies no
canonical map transporting it to the 14,400-dimensional tensor arena.
Choosing one would be a fitted carrier.  It is therefore **not a derived
candidate**.

A spectral sign chosen after seeing `D` remains excluded by the mission and
the stopping rule.  The ban on an extra trivial even doubling is stronger
than policy: it is a **DERIVED lemma**.  If

`H=H0 direct-sum H0`, `pi=diag(pi0,pi0)`,
`gamma=diag(1,-1)`, and `J=[[0,J0],[J0,0]]`,

then every product `pi(a)Jpi(b)J^-1` is `diag(S,S)`.  Hence every
intersection-form entry is

`Tr(gamma diag(S,S))=Tr(S)-Tr(S)=0`.

The intersection form is identically zero for every algebra and every
choice of projections, so Poincare duality fails maximally.  The same
sheet-identical form also prevents a metric-dimension-zero orientability
cycle from representing `gamma`.  A viable grading must therefore come
from genuinely different sheet data, not a trivial copy inserted solely
to manufacture oddness.

## 4. Complete axiom table

`odd` means `{D,gamma}=0`.  A dash means that no derived nontrivial grading
exists, not that an untested grading was assigned a negative value.

| derived operator | self-adjoint/nonzero | `J^2` | `JD` | derived odd `gamma` | `Jgamma=-gamma J` | order zero | first order | inner forms |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `D_c^-` | yes/yes | `+` | `-` | none | no | yes | yes | zero |
| `D_c^+` | yes/yes | `+` | `+` | none | no | yes | yes | zero |
| self-adjoint class-sum `D^-` | yes/nonzero unless scalar cancellation | `+` | `-` | none | no | yes | yes | zero |
| self-adjoint class-sum `D^+` | yes/nonzero | `+` | `+` | none | no | yes | yes | zero |
| adjacency/Laplacian polynomial lifts | yes, except the zero polynomial | `+` | `-/+` by lift | none | no | yes | yes | zero |
| vertex `Box_0^-` lift | yes/yes | `+` | `-` | none | no | yes | yes | zero |
| vertex `Box_0^+` lift | yes/yes | `+` | `+` | none | no | yes | yes | zero |

No row has the KO6 signs `(+,+,-)`, oddness, and nonzero fluctuations
simultaneously.  Therefore no real-even spectral triple is claimed.

## 5. SM-type corners

The two block choices are

`(rho0,rho1,rho8)` and `(rho0,rho7,rho2)`.

Each has complex regular support

`1^2+2^2+3^2=14`.

On the 14,400-dimensional tensor arena, its left corner unit has rank

`14*120=1680`,

not 14,400.  Thus the literal `C+H+M3(C)` action is nonunital on `H`.
The Frobenius--Schur `-1` on the selected two-dimensional block supplies its
quaternionic real form but does not repair the missing unit.

Restricting the algebra does not improve any derived row:

- order zero and first order still hold;
- every full-algebra-zero one-form remains zero;
- no derived grading becomes odd;
- the corner representation is nonunital.

Both Galois choices therefore have the same axiom table and the same failure.
Allocating the other six sectors would be an additional noncanonical
representation choice and is not performed.

## 6. Physics gate

The gate does not open.  There are no represented inner gauge fields or
Yukawa fluctuations, so no `Y`, Route C anomaly forcing, `M15/M16`
multiplet census, generation count, or mass-ratio comparison is licensed.
The registered `Z[phi]` target list receives zero new trials and the
look-elsewhere count is zero.

## 7. Closing boundary of the matter program

Under the binding stopping rule, the audited arena sequence is now closed:

| arena | exact positive content | failed gate |
|---|---|---|
| primal `22 Reg` | derived nonzero form-odd `d+d*`; right convolution has nonzero forms | all geometric order-zero `J` candidates fail KO6/reality; arbitrary multiplicity-mixing `J` remains open |
| primal plus dual `44 Reg` | cellular star supplies KO6 candidates; one-form space is nonzero before the gate | pure star fails order zero/first order; star--inversion loses the `JD` sign; arbitrary mixing remains open |
| icosahedral/orbifold base | exact stabilizers `C10,C4,C6`; nonregular module; scalar cell functions have 240 one-form directions | scalar algebra fails first order and KO6; twisted stabilizer fibers are unspecified |
| final `A tensor A^op` | canonical order zero and flip-star `J`; exact adjoint decomposition | no derived odd grading and all derived equivariant vertex lifts have zero one-forms |

The Q8 counterexample remains important: it proves that free arenas and the
listed axioms are not abstractly incompatible.  It does not canonically
select a 2I solution.

What remains genuinely open is sharply delimited:

1. the continuous multiplicity-mixing `J` varieties on the already existing
   `22 Reg` and `44 Reg` arenas;
2. a specified twisted stabilizer fiber and Dirac on the already identified
   orbifold arena;
3. a derived, not fitted, algebra allocation and grading.

A future matter attempt would have to supply one of those missing structures
from independent geometry and then pass every axiom exactly.  It may not
enlarge the arena merely to create room, insert a fitted `D` or `gamma`, or
reinterpret an available matrix block as a selected physical carrier.

No further arena enlargement is proposed or authorized.

## Status ledger

### Strengthened

- **DERIVED:** exact flip-star real structure and order-zero bimodule.
- **DERIVED:** diagonal-adjoint decomposition and trivial central action.
- **DERIVED:** adjacency is the inverse-closed central 12-class sum.
- **DERIVED:** complete reduction of class-sum/Laplacian candidates.
- **DERIVED:** every audited equivariant vertex lift has zero inner forms.
- **DERIVED:** no audited derived grading makes a derived candidate odd.

### Downgraded / refuted

- **REFUTED:** `C[2I] tensor C[2I]^op=End_C(C[2I])`.
- **DOWNGRADED:** the tensor bimodule is canonical, but not uniquely “most
  economical.”
- **DOWNGRADED:** escaping central-parity segregation does not itself produce
  a grading or matter Dirac.

### Open

- arbitrary multiplicity-mixing `J` on the already authorized free arenas;
- explicitly derived twisted orbifold fibers;
- any future independently derived matter structure satisfying the closed
  boundary conditions above.
