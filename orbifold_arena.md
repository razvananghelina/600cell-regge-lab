# Icosahedral Hopf-base orbifold arena

Date: 2026-07-27 (third session)

## Decision

The non-free icosahedral cochain arena and its stabilizers are **DERIVED**.
The scalar-cochain real-triple gate is **DERIVED negative**.

The Hopf base has:

| layer | cells | `2I` stabilizer | odd-order part |
|---|---:|---|---|
| vertices | 12 | `C10` | `C5` |
| edges | 30 | `C4` | `C1` |
| faces | 20 | `C6` | `C3` |

Thus it is not a multiple of the regular representation.  Its exact
oriented-cochain decomposition contains only integer-spin 2I sectors.

The canonical cell-function algebra `C^62` with coefficient conjugation
satisfies order zero and has 240 nonzero one-form directions, but fails first
order.  Its real-structure signs are `(+,+,+)`, not KO6.  Global scalars pass
first order but have zero fluctuations.

The stabilizer group algebras are abstractly available, but the scalar
cochain space does not carry faithful `C[C10]`, `C[C4]`, or `C[C6]` fibers.
Adding their characters/twisted sectors is additional Hilbert-space data.
Accordingly no real spectral triple, `Y`, anomaly input, multiplet census, or
Yukawa sector is constructed.

All finite claims are checked in
`reproducible/verify_free_orbifold_arenas.py`.

## 1. Full oriented cochain complex

Use the regular icosahedron with standard coordinates given by signed
permutations of `(0,1,phi)`.  Increasing vertex tuples orient its 30 edges
and 20 triangular faces.  The integer coboundaries are

`d0:C^12 -> C^30`,

`d1:C^30 -> C^20`.

Exact checks give:

`d1 d0=0`,

`rank(d0,d1)=(11,19)`,

`b=(1,0,1)`.

The ranks have GF2 lower certificates `(11,19)`.  Connectedness bounds the
first rank by 11, and `d^2=0` plus Euler closure gives the matching second
upper bound.  **DERIVED exact.**

The derived Dirac

`D_ico=d+d*`

is a nonzero, real, self-adjoint `62 by 62` integer matrix with 240 nonzero
entries and anticommutes with form parity.

## 2. Exact stabilizers and non-regular decomposition

The verifier constructs all 60 orientation-preserving rotations of the
icosahedron and doubles them through

`2I -> A5`

because the central element acts trivially on the base.  Direct fixed-cell
counts give A5 stabilizer orders

`(5,2,3)`

on vertices, edges, and faces.  Their full preimages in `2I` have orders

`(10,4,6)=(120/12,120/30,120/20)`.

The unique-involution result already proved for `2I` implies these cyclic
lifts are exactly

`C10`, `C4`, `C6`.

Signed fixed-cell characters on the five A5 classes give:

| layer | `1` | `3` | `3'` | `4` | `5` |
|---|---:|---:|---:|---:|---:|
| `C0` | 1 | 1 | 1 | 0 | 1 |
| `C1` | 0 | 2 | 2 | 2 | 2 |
| `C2` | 1 | 1 | 1 | 2 | 1 |

Dimensions close at 12, 30, and 20.  The central `-1` acts trivially, so no
spinorial/quaternionic irrep occurs.  In particular this 62-dimensional
module is not `n Reg(2I)` for any `n`.  **DERIVED.**

## 3. Canonical algebra candidates

### Cell functions

The completely specified algebra on the existing scalar arena is

`A_cell=C^(12+30+20)=C^62`,

acting diagonally on oriented cell basis vectors.  Coefficient conjugation
`K` has

`(J^2,JD,Jgamma)=(+,+,+)`.

Because `A_cell` is commutative, order zero holds.  Independent endpoint
projectors on any nonzero incidence give an exact nonzero

`[[D,a],KbK^-1]`.

Thus first order fails.

Each directed nonzero incidence matrix unit can be isolated by two cell
projectors, so

`dim_C Omega_D^1(A_cell)=nnz(D)=240`.

The calculus is nonzero but not licensed as physical gauge fields because
first order and KO6 fail.  Restricting to global constants makes first order
hold and forces `Omega_D^1=0`.

### Stabilizer algebras

The full stabilizer algebras are

`C[C10]`, `C[C4]`, `C[C6]`.

The initially suggested `C[C5]` is the algebra of the vertex stabilizer's
odd-order subgroup, not of its full `2I` stabilizer.  Over `C`, all three
cyclic group algebras split into one-dimensional character sectors.

On the scalar cochain module a stabilizer fixes its cell by a scalar
orientation character.  This is not a faithful regular stabilizer fiber and
does not realize the full group algebra.  A twisted-sector construction
would have to add selected stabilizer characters or induced bundles and
specify:

- their Hilbert multiplicities;
- the algebra action and opposite action;
- a grading and antiunitary;
- incidence/Dirac maps between different stabilizer types.

None is selected by the scalar complex alone.  These are
**STRUCTURAL/OPEN**, not silently inserted sectors.

Because cyclic groups have only one-dimensional complex irreducibles,
ordinary projective phases do not by themselves manufacture a canonical
`H` or `M3(C)` factor.  Noncommutative matrix blocks require multiplicities
or induced/twisted bundles.

## 4. Relation to the segregation theorem

Removing the unique central involution from the three cyclic stabilizers
leaves their maximal odd-order subgroups:

`C10 -> C5`,

`C4 -> C1`,

`C6 -> C3`.

The set is exactly

`{C1,C3,C5}`,

the complete odd-order subgroup list independently derived by the central
parity segregation theorem.  The maximal nontrivial cases are again `C3`
and `C5`.

This equality is **DERIVED group theory**, not a numerical coincidence.
However, the statement

“physical matter is forced to live in the vertex/face twisted sectors”

is only a **PATTERN/OPEN mechanism**.  The scalar cochain module has no
spinorial sectors and no valid matter algebra.  A functor from central-parity
escape subgroups to stabilizer-localized matter modules remains to be
constructed.

The structural convergence is therefore real but limited:

- **DERIVED:** identical subgroup types arise from two independent exact
  classifications;
- **PATTERN:** interpreting that equality as localization of generations or
  SM matter.

## 5. Axiom verdict and physics gate

| candidate | KO signs | order zero | first order | fluctuations |
|---|---|---:|---:|---:|
| `C^62`, coefficient `K` | `(+,+,+)` | holds | fails | nonzero, dim 240 |
| global `C` | `(+,+,+)` | holds | holds | zero |
| stabilizer group algebras | representation not defined | open | open | open |

A two-dimensional primal--dual Hodge star cannot repair the KO grading sign:
in dimension two, `k` and `2-k` have the same parity, so Hodge duality
commutes with form parity.  A separately doubled or twisted grading would be
additional data.

The gate does not open.  Consequently:

- no `Y` is selected;
- Route C is not activated;
- no `M15/M16` or generation census is licensed;
- no Yukawa block is defined;
- the frozen `Z[phi]` mass-target comparison is skipped with zero new trials.

## Status ledger

### Strengthened

- **DERIVED:** exact icosahedral cochain complex and S2 homology.
- **DERIVED:** stabilizers `C10,C4,C6` and odd parts `C5,C1,C3`.
- **DERIVED:** exact layerwise isotypic decomposition.
- **DERIVED:** non-regularity and absence of spinorial sectors.
- **DERIVED:** 240-dimensional cell-function one-form calculus.
- **DERIVED:** exact equality with the segregation escape subgroup list.

### Negative / delimited

- **DERIVED negative:** full cell functions fail first order and KO6.
- **DERIVED negative:** global scalars fluctuate trivially.
- **DELIMITED:** scalar fixed cells do not automatically carry their full
  stabilizer group algebras.
- **PATTERN only:** subgroup equality as physical matter localization.

### Open

- a derived twisted-sector Hilbert space and induced algebra;
- stabilizer-changing Dirac maps and a valid real structure;
- noncommutative weak/color blocks, `Y`, anomalies, multiplets, generations,
  Yukawa data, and frozen mass comparisons;
- comparison with a specifically constructed `C^2/2I` orbifold spectral
  triple rather than McKay exceptional-data analogy alone.

