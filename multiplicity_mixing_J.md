# Multiplicity-mixing real structures on the free 2I arenas

Date: 2026-07-27 (fourth session)

## Decision

The order-zero problem is classified exactly, and it corrects a second
over-simplification in the proposed search.  An embedding is not classified
by one nine-entry multiplicity vector.  It is classified by a nonnegative
`9 x 9` bimodule multiplicity matrix with two weighted margin conditions.
Each such discrete type carries a positive-dimensional continuous unitary
orbit.  The `J^2`, grading, reality, and first-order equations depend on the
point of that orbit, not only on its integer label.

The repository-derived geometric strata have been exhausted and remain
negative:

- on `22 Reg`, coefficient conjugation and central/Galois compositions have
  the wrong grading sign; orbitwise inversion has order zero and first order
  but has no `JD=+/-DJ` sign;
- on `44 Reg`, pure cellular star has KO6 signs but fails order zero and first
  order, while star--inversion repairs those two axioms but has no `JD` sign.

The arbitrary multiplicity-mixing strata have **not** been exhaustively
solved.  Therefore neither existence nor a 2I-specific no-go is derived.
The verdict for a non-geometric `J` on `m=22` or `m=44` is
**STRUCTURAL/OPEN**.  This is a boundary result, not outcome (iii).

The Q8 counterexample remains the exact positive control.  Its multiplicity
space contains two complete regular `Q8` modules (`16=2*8`), permitting an
explicit factor swap.  Neither `C^22` nor `C^44` contains a regular `2I`
module because both dimensions are below 120.  This explains why the Q8
formula does not transfer literally, but it is not a theorem excluding
distributed 2I bimodule embeddings.

No claim that a real spectral triple has been constructed is made.

## 1. Exact order-zero classification

Put

`A=C[2I]=direct_sum_i M_(d_i)(C)`,

with McKay-chain dimensions

`d=(1,2,3,4,5,6,4,2,3)`.

On `H=C[2I] tensor C^m`, the right action has irreducible `i` with
multiplicity `d_i m`.  Its commutant is

`A'=L(C[2I]) tensor M_m(C)
   =direct_sum_i M_(d_i m)(C)`.

For `J=UK`, order zero says that the opposite representation

`pi^o(b)=U overline(pi(b)) U*`

is a unital star-representation of `A` in `A'`.  Simultaneously
decomposing the commuting left and right actions gives

`H = direct_sum_(i,j) V_i tensor V_j^o tensor C^(k_ij)`.

Thus the complete discrete invariant is a matrix

`K=(k_ij) in M_9(Z_{\geq 0})`

satisfying

`sum_j d_j k_ij=d_i m` for every `i`,

`sum_i d_i k_ij=d_j m` for every `j`.                 (1)

Conversely every nonnegative integer matrix satisfying (1) constructs a
commuting representation with the required left and opposite marginal
multiplicities.  This is the exact finite classification up to unitaries in
`A'`.  It is a labelled finite set for each fixed `m`, but the margin
equations, rather than a nine-vector or the single total dimension `120m`,
are the honest enumeration.

Examples for both `m=22` and `m=44` are

`k_ij=m delta_ij`,

and `k_ij=m delta_(j,sigma(i))` for every dimension-preserving permutation
`sigma`.  In particular the simultaneous Galois flip

`rho1<->rho7`, `rho2<->rho8`

and the independent interchange of the two four-dimensional labels already
give distinct types having the same global regular character.  Hence a
single multiplicity vector provably loses embedding data.

### Continuous dimensions

For fixed `K`, its conjugacy orbit inside `A'` has real dimension

`dim O_K=sum_i [(d_i m)^2-sum_j k_ij^2]`.              (2)

For a fixed represented image, the unitaries implementing it form a torsor
for the commutant of the original `A` action, of real dimension

`sum_i (d_i m)^2=120m^2`.                              (3)

For the diagonal type, (2) is `111m^2`: 53,724 for `m=22` and 214,896 for
`m=44`.  The fixed-image torsor dimensions are respectively 58,080 and
232,320.  These dimensions are before imposing `J^2`, grading, `JD`, or
first order.  They show why listing integer types is not an exhaustive
search over `J`.

Because every complex 2I irrep is self-conjugate, `J` sends the `(i,j)`
bimodule sector to `(j,i)`.  A necessary condition for a real structure is
therefore `k_ij=k_ji`, after incorporating the real/quaternionic
Frobenius--Schur forms.  It is not sufficient: the FS-negative blocks impose
symplectic rather than real bilinear forms, and the grading split and `D`
must also be respected.  **DERIVED necessary condition; sufficiency OPEN.**

## 2. Remaining equations on each orbit

For a representative embedding of type `K`, write its implementing unitary
as `U`.  In the real signed cell basis, `D` and `gamma` are real.  The
remaining equations are exactly

`U overline(U)=epsilon`,

`U gamma=-gamma U`,

`U D=epsilon' D U`,

`[[D,R_g], U R_h U*]=0` for all `g,h in 2I`.           (4)

Unitarity is `U*U=1`.  These are polynomial equations of degree at most
three in the real and imaginary parts of `U` (quadratic after retaining the
opposite representation as an auxiliary variable).  Nonzero fluctuations
are already certified independently by `[D,R_g]!=0`; they survive whenever
the other axioms admit a `J`.

Equations (1)--(4) are a finite real-algebraic decision problem.  But the
integer type alone does not determine the answer, and no exact
Groebner/real-quantifier-elimination certificate or interval-certified
solution of (4) on the 2640- or 5280-dimensional arena is present in the
repository.  Reporting a zero-dimensional survivor table for all `K` would
therefore be fabricated.  **RULE ZERO boundary.**

### Certified per-stratum table

| stratum | `J^2` | `J gamma=-gamma J` | `JD=+/-DJ` | order zero | first order | one-forms |
|---|---:|---:|---:|---:|---:|---:|
| primal coefficient/Galois conjugation | `+` | no | `+` | no for full `A` | not a passed triple | nonzero algebraically |
| primal orbitwise inversion | `+` | no | neither | yes | yes | nonzero |
| doubled pure cellular star, KO6 variant | `+` | yes | `+` | no | no | dimension 1191 before gate |
| doubled star--inversion | `+` | yes | neither | yes | yes | nonzero candidate space |
| arbitrary symmetric type satisfying (1) | unresolved | unresolved | unresolved | yes by construction | unresolved | nonzero if gate passes |
| Q8 factor-swap control | `+` | yes | `+` | yes | yes | nonzero |

The first four rows are **DERIVED exact negatives** from the registered
integer/signed-permutation verifiers.  The arbitrary row is
**STRUCTURAL/OPEN**, not a numerical negative.

## 3. Full algebra versus the SM corner

The full `C[2I]` representation is unital and the classification above
applies directly.

The proposed standard-model-type blocks are different.  The
`(rho0,rho1,rho8)` corner, or its simultaneous Galois mate
`(rho0,rho7,rho2)`, has complex regular rank

`1^2+2^2+3^2=14`,

not 120.  Its corner unit has rank `14m`, whereas the arena identity has rank
`120m`.  The quaternionic real form in the two-dimensional block does not
alter this support calculation.

Therefore the displayed `C+H+M3(C)` is not by itself a unital algebra action
on the full arena.  A full finite-triple test needs an allocation of its unit
and representations on the other six Wedderburn sectors.  The two Galois
choices do not provide that allocation.  Testing the literal corner gives a
scoped nonunital representation only.  Its arbitrary multiplicity-mixing
`J` problem remains **OPEN** for the same reason as the full-algebra problem,
with additional carrier freedom.

## 4. Verdict and physics gate

The exact alternatives currently supported are:

1. **DERIVED negative:** no repository-derived geometric `J` passes all the
   listed axioms for the derived `D`.
2. **DERIVED scoped positive control:** Q8 proves the listed KO/order axioms
   plus nonzero fluctuations are satisfiable on a free arena.  Its trivial
   chirality doubling has an identically zero intersection form, so it does
   not prove existence of a Poincare-dual, manifold-like finite triple.
3. **OPEN:** whether a non-geometric point in one of the 2I unitary orbits
   solves (4).

Consequently there is no licensed gauge-field, `Y`, anomaly-forcing,
generation, multiplet, or Yukawa extraction.  The registered
`Z[phi]` mass-exponent targets were not tested and the new look-elsewhere
count is zero.

Even if a point of (4) is eventually found, equations (2)--(3) show that
existence alone would not make it **DERIVED**.  A canonical selection from
the surviving moduli, using the Hopf map, McKay grading, Galois involution,
form degree, or orbit chart, would still be required.  Otherwise the result
would be an **EXISTS/STRUCTURAL** carrier choice.

## 5. Scope corrections

The Q8 counterexample does not invalidate the route-specific node, edge, or
preprojective calculations.  It invalidates extrapolating them to every free
arena or every antiunitary:

- the node result concerns its explicit diagonal restriction and node
  commutant;
- the edge theorem concerns the canonical endpoint bimodule and orientation
  grading;
- the preprojective boundary concerns its specified inversion/Galois real
  structures and maximal diagonal multiplicity action;
- the Kähler--Dirac geometric-J negatives do not cover arbitrary solutions
  of (1)--(4).

These boundaries are now appended explicitly to the affected notes.

## Status ledger

### Strengthened

- **DERIVED:** complete weighted-bimodule classification (1) of order-zero
  embeddings.
- **DERIVED:** exact continuous orbit and intertwiner dimensions (2)--(3).
- **DERIVED:** the SM blocks are rank-14 corners, not unital arena actions.
- **DERIVED:** all currently geometric 2I `J` candidates fail at least one
  required axiom.

### Downgraded

- **DOWNGRADED:** “all 2I multiplicity-mixing real structures fail” is not
  proved.
- **DOWNGRADED:** an embedding type cannot be represented by one global
  multiplicity vector.
- **DOWNGRADED:** the SM corner cannot be treated as an already specified
  unital finite algebra.

### Open

- the exact or interval-certified solution of (4) across all symmetric
  matrices satisfying (1), for `m=22` and `m=44`;
- the dimension and connected components of every surviving `U` variety;
- a geometrically selected point if a variety is nonempty;
- every physics extraction contingent on that gate.
