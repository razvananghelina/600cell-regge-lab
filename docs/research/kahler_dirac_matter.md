# Kaehler--Dirac matter on the full 600-cell complex

Date: 2026-07-24

## Decision

**DERIVED positive:** the full primal cochain space

`C=C^0+C^1+C^2+C^3`, with dimensions `(120,720,1200,600)`,

has a nonzero, real, self-adjoint, fully `2I`-equivariant operator
`D=d+d*`.  With `gamma_form=(-1)^k` on `C^k`,
`{D,gamma_form}=0` exactly.

**DERIVED independence:** central/spin parity is the signed cell action of
`z=-1 in 2I`.  It commutes with form parity, is neither form parity nor its
negative, and `[D,gamma_spin]=0`.  Thus the central-parity segregation theorem
is not violated: that theorem kills operators odd for `gamma_spin`; this
operator is even for `gamma_spin` and odd for the independent
`gamma_form`.

**DERIVED algebra content:** as exact `2I` modules,

`C^0=Reg`, `C^1=6 Reg`, `C^2=10 Reg`, `C^3=5 Reg`,

and hence `C=22 Reg`.  The old multiplicity obstruction is completely gone.
The total multiplicities in McKay-chain order are

`(22,44,66,88,110,132,88,44,66)`.

**DERIVED negative at the decisive matter gate:** large matrix factors make
abstract `C+H+M3(C)` actions available, but do not canonically select one.
On an isotypic multiplicity space `C^(22 dim rho)`, its natural `U(22 dim
rho)` basis freedom moves every proposed one-, two-, or three-dimensional
seed carrier.  Requiring the construction to be invariant under this
unbroken freedom leaves only the center on each isotypic block.  That
commutative algebra cannot contain `H` or `M3(C)`.

Consequently order zero and first order can be studied after choosing an
embedding, but there is no uniquely defined seed bimodule on which to perform
the requested Standard-Model test.  Choosing carrier subspaces or allocating
their complements would insert precisely the missing matter data.  The
result is not `D=0`; it is that the canonical `D` does not canonically produce
the requested finite algebra.  Hypercharge, Route C, and the multiplet census
therefore do not start.

All finite checks are in
`reproducible/verify_kahler_dirac.py`.

## 1. Complex, orientations, and Hodge anchors

Every simplex is stored as an increasing vertex tuple.  Its orientation is
that order, and

`partial[v_0,...,v_k]=sum_i (-1)^i[v_0,...,omit(v_i),...,v_k]`.

The matrices called `d_k` have shape `f_(k+1) by f_k`; they are coboundaries,
the transposes of the corresponding boundary matrices in the orthonormal
oriented-cell basis.  The verifier constructs them as sparse integer
matrices.  It proves exactly over `Z`

`d_1 d_0=0`, `d_2 d_1=0`.

Numerical ranks (double precision, singular-value threshold `10^-9`) are
`(119,601,599)`.  They give

`b=(1,0,0,1)`

and the layerwise Hodge counts

| degree | exact | harmonic | coexact | total |
|---:|---:|---:|---:|---:|
| 0 | 0 | 1 | 119 | 120 |
| 1 | 119 | 0 | 601 | 720 |
| 2 | 601 | 0 | 599 | 1200 |
| 3 | 599 | 1 | 0 | 600 |

The incidence identities are **exact**; ranks and Betti numbers are
**high-precision verified numerics**, independently consistent with the
known `S^3` topology and Euler characteristic.

In block form `D` has `d_k` below and `d_k^T` above the diagonal.  It has
`14880` nonzero entries.  Self-adjointness, nonzeroness, and form oddness are
exact sparse-matrix identities.

### Why `2640=c0`

This is the same object, not a coincidence.  In
`verify_spectral_action.py`, the spectral triple Hilbert space is precisely
this oriented-cell cochain space and the zeroth finite spectral coefficient
is defined by

`c0=Tr(I_C)=dim C=sum_k f_k`.

Therefore

`c0=120+720+1200+600=2640`

is a **DERIVED definitional identity for the same complex**.  The separate
factorization `2640=240*11` is an additional numerical identity, not the
reason `c0` equals the cell count.

## 2. Full symmetry and the two gradings

The vertices are the 120 unit quaternions of `2I`.  Left quaternion
multiplication acts freely and transitively on vertices and induces signed
permutations on every oriented cell layer.  The verifier constructs all 120
actions.  It checks coboundary equivariance for a faithful generator exactly;
the construction by functorial oriented boundary makes the same identity
hold for the whole action.  Thus `D` is `2I`-equivariant.  **DERIVED.**

Let `gamma_spin=rho(-1)`.  On cell cochains this is an antipodal signed
permutation, not a degree sign.  Exact sparse identities give

`gamma_spin^2=1`, `[gamma_spin,gamma_form]=0`,

`gamma_spin != +/- gamma_form`,

`[D,gamma_spin]=0`, `{D,gamma_form}=0`.

The old theorem assumed that its Dirac was odd for `gamma_spin`.  The present
operator is not.  It belongs to the theorem's allowed even commutant while
being odd for a different grading.  **DERIVED: no contradiction and no
symmetry breaking.**

## 3. Exact character decompositions

Signed fixed-cell traces on all nine conjugacy classes, paired with the exact
`Q(sqrt(5))` character table, give:

| irrep | dim | `C0` | `C1` | `C2` | `C3` | total |
|---|---:|---:|---:|---:|---:|---:|
| `rho0` | 1 | 1 | 6 | 10 | 5 | 22 |
| `rho1` | 2 | 2 | 12 | 20 | 10 | 44 |
| `rho2` | 3 | 3 | 18 | 30 | 15 | 66 |
| `rho3` | 4 | 4 | 24 | 40 | 20 | 88 |
| `rho4` | 5 | 5 | 30 | 50 | 25 | 110 |
| `rho5` | 6 | 6 | 36 | 60 | 30 | 132 |
| `rho6` | 4 | 4 | 24 | 40 | 20 | 88 |
| `rho7` | 2 | 2 | 12 | 20 | 10 | 44 |
| `rho8` | 3 | 3 | 18 | 30 | 15 | 66 |

Every row is `dim(rho)*(1,6,10,5)`.  Integrality and dimension closure are
verified to residual below `2e-8`; the characters themselves are exact
signed integers and the character table is algebraic exact.  The table is
therefore labeled **exact character computation with numerically evaluated
algebraic inner products**.

The scalar seed can occur on `rho0`, a quaternionic weak seed is available
on either FS-negative two-dimensional spinor (`rho1` or its Galois mate), and
a color-sized complex factor is available on a three-dimensional
multiplicity carrier associated with `rho2` or `rho8`.  Availability is
**DERIVED**; a concrete embedding and complement allocation are
**STRUCTURAL/OPEN**.

## 4. Real structures and the decisive order conditions

There is an important correction to the proposed real-structure route.
Poincare duality on this triangulated `S^3` does not give an invertible
primal-cochain map `C^k -> C^(3-k)`: already

`dim C^0=120 != 600=dim C^3`,

`dim C^1=720 != 1200=dim C^2`.

A combinatorial Hodge star maps primal cells to complementary-dimensional
cells of the **dual cellulation**, or realizes duality on cohomology.  It is
not an antiunitary endomorphism of the stated 2640-dimensional primal arena.
This is a **DERIVED dimension obstruction** to that candidate, not a failure
of Poincare duality.

The endomorphism candidates actually defined on this arena have signs:

| candidate | `J^2` | `JD` | `J gamma_form` | status |
|---|---:|---:|---:|---|
| coefficient conjugation `K` | `+` | `+` | `+` | DERIVED |
| central antipode times `K` | `+` | `+` | `+` | DERIVED |
| primal Hodge star times `K` | -- | -- | -- | nonexistent on this arena |
| orbitwise `g -> g^-1` | `+` | not fixed | `+` | requires choices of origins among 22 free orbits |

The last variant is canonical on `C^0=Reg` after declaring the identity
vertex, but extensions to the other 21 regular orbits require choosing an
origin in each orbit, and commutation with the incidence-coupling `D` is not
automatic.  It is **OPEN/STRUCTURAL**, not a derived real structure.

For the abstract full multiplicity algebra, a standard-form opposite action
can make order zero and first order compatible with equivariant
multiplicity maps.  But reducing that large algebra to
`C+H+M3(C)` demands non-invariant carrier choices.  The exact killing
constraint on a *canonical* seed action is:

`U(m)`-naturality on every multiplicity space implies scalar action by
Schur's commutant theorem.

The surviving canonical algebra is therefore only block-center data.  It
cannot realize the noncommutative weak and color summands.  Any particular
smaller embedding may pass or fail first order depending on its allocation;
testing one is hand insertion forbidden by Rule Zero.

## 5. Taste audit

The exact harmonic kernel has dimension

`dim ker D=sum_k b_k=2`,

one constant 0-form and one volume-class 3-form.  Form parity pairs every
nonzero eigenvalue `lambda` with `-lambda`.  Hence the remaining spectrum is
exactly 1319 positive/negative pairs, counting multiplicity.

There is **no derived uniform taste multiplicity**.  The familiar continuum
`2^d` exterior-algebra count describes the local Clifford module of a smooth
orthonormal frame.  This irregular finite simplicial complex has no verified
global Clifford generators commuting with its Laplacian and enforcing a
uniform degeneracy.  The exact natural algebra found here is the form-parity
`Z2`, the `2I` action, and spectral commutants that vary by eigenvalue.

- `N_gen=3`: **no match is derived**.  The kernel is 2, while nonzero
  degeneracies are representation- and eigenvalue-dependent.
- Distinguished `C3`: **none**.  Full `2I` remains unbroken, and no canonical
  taste algebra selects one of its ten conjugate `C3` subgroups.
- Any identification of a threefold spectral degeneracy with generations is
  **PATTERN** unless a canonical invariant projector and a repeated matter
  module are produced.

## Status ledger

### Strengthened

- **DERIVED:** an unbroken, nonzero, self-adjoint, form-odd Kähler--Dirac
  operator exists on the exact spectral-action Hilbert space.
- **DERIVED:** `2640=c0` is identity of the same cochain object.
- **DERIVED:** spin and form gradings are independent; the segregation theorem
  does not apply to form oddness.
- **DERIVED:** `C=22 Reg(2I)` and all exact layer multiplicities above.
- **DERIVED:** the old small-multiplicity obstruction is removed.
- **DERIVED:** kernel dimension 2 and 1319 nonzero `+/-` pairs.

### Downgraded / killed

- **DERIVED negative:** primal Hodge star is not an endomorphism of this
  cochain arena; a dual cellulation is required.
- **DERIVED negative:** availability of large `M2` and `M3` carriers does not
  canonically select the Standard-Model algebra.
- **DERIVED negative:** invariance under full multiplicity-basis freedom
  reduces canonical seed actions to scalars.
- **PATTERN rejected:** Kähler tastes do not derive three generations or a
  distinguished `C3`.

### Open

- a geometrically invariant `C+H+M3(C)` subalgebra/allocation;
- a real structure on a justified primal-plus-dual arena and its KO signs;
- order-zero/first-order tests after, and only after, such an algebra is
  derived;
- a selected `Y`, Route C anomaly input, color orientation, and multiplet
  census;
- the optional critical-circle CP-phase diagnostic was not needed for this
  decisive audit and remains open.

**2026-07-27 scope correction (DERIVED).**  The subsequent Q8 counterexample
proves that coefficient conjugation, Hodge, inversion, and Galois
compositions do not exhaust antiunitaries on a free arena: an antiunitary may
mix the regular coordinate with multiplicity space.  The geometric-candidate
negatives in this note remain exact, but they are not a no-go for every
`J` on `22 Reg` or `44 Reg`.
