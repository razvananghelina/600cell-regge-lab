# Inflation towers: exact finite Bratteli audit

Date: 2026-07-22

## Decision

The two proposed inclusion systems are **DERIVED**, and their multiplicities
do grow.  Inflation therefore removes the narrow matrix-size obstruction of
the one-cell node module: the rooted McKay floors contain arbitrarily large
matrix blocks.

It does not supply a Standard-Model finite bimodule.  The rooted McKay tower
has a stronger parity obstruction.  Every even floor contains only
integer-spin nodes, including the scalar and color seeds, while every odd
floor contains only spinor nodes, including the quaternionic weak seed.  No
single floor can contain the three requested canonical seed factors.  This is
an all-level consequence of bipartiteness, not merely a failure through the
tested depth.

If one replaces a floor by a cumulative sum of floors, the three factor types
can be placed in different summands, but that replacement is **STRUCTURAL**.
For the canonical independent endpoint actions, every nonzero shift between
successive floors is gamma-odd and fails first order.  Thus the old
orientation argument does not recur word for word, but its endpoint-projector
mechanism does: a source projector and a target projector give a nonzero
double commutator on every nonzero shift block.

Consequently no derived `Y` exists and Route C cannot start.  The inflation
tower strengthens the boundary theorem for the precise rooted-floor and
canonical cumulative-endpoint constructions audited here.  It does not rule
out a different, independently derived representation of the AF algebra with
an opposite action and odd Dirac.

All floor vectors, Smith forms, stable matrices, traces, parity statements,
and the finite first-order witness are checked by
`reproducible/verify_bratteli_tower.py`.

## 1. The two rooted towers

### A. Fibonacci fusion

Use column vectors and

`F = [[1,1],[1,0]]`,  `m(0)=(1,0)`,  `m(n)=F^n m(0)`.

The floors through 12 are

| n | block sizes |
|---:|---|
| 0 | `(1,0)` |
| 1 | `(1,1)` |
| 2 | `(2,1)` |
| 3 | `(3,2)` |
| 4 | `(5,3)` |
| 5 | `(8,5)` |
| 6 | `(13,8)` |
| 7 | `(21,13)` |
| 8 | `(34,21)` |
| 9 | `(55,34)` |
| 10 | `(89,55)` |
| 11 | `(144,89)` |
| 12 | `(233,144)` |

Thus the level algebra is the sum of the two nonzero blocks `M_m(C)` shown in
the row.  The inclusion is **DERIVED** from `tau tensor tau=1+tau`.  No weak,
color, real structure, opposite action, or Dirac is supplied by this fusion
rule, so reading its large blocks as matter is **PATTERN**, not a matter
construction.

### B. Rooted McKay tensor tower

In chain order

`(rho0,...,rho8)=(1,2,3,4s,5,6,4,2',3')`,

the affine-E8 edges are

`01,12,23,34,45,56,67,58`.

For their symmetric adjacency matrix `A`,

`m(0)=e_rho0`,  `m(n)=A^n e_rho0`.

The exact floors are

| n | `(m0,m1,m2,m3,m4,m5,m6,m7,m8)` |
|---:|---|
| 0 | `(1,0,0,0,0,0,0,0,0)` |
| 1 | `(0,1,0,0,0,0,0,0,0)` |
| 2 | `(1,0,1,0,0,0,0,0,0)` |
| 3 | `(0,2,0,1,0,0,0,0,0)` |
| 4 | `(2,0,3,0,1,0,0,0,0)` |
| 5 | `(0,5,0,4,0,1,0,0,0)` |
| 6 | `(5,0,9,0,5,0,1,0,1)` |
| 7 | `(0,14,0,14,0,7,0,1,0)` |
| 8 | `(14,0,28,0,21,0,8,0,7)` |
| 9 | `(0,42,0,49,0,36,0,8,0)` |
| 10 | `(42,0,91,0,85,0,44,0,36)` |
| 11 | `(0,133,0,176,0,165,0,44,0)` |
| 12 | `(133,0,309,0,341,0,209,0,165)` |

The level algebra is `direct sum_i M_{m_i(n)}(C)`, omitting zero blocks.
The dimension vector `d=(1,2,3,4,5,6,4,2,3)` obeys `A d=2d`, and the finite
certificate `d dot m(n)=2^n` holds through level 12.  Both the inclusion and
the rooted initial condition are **DERIVED**.

## 2. Exact dimension groups and states

### Fibonacci: the golden dimension group

`det(F)=-1`, so every bonding map is an automorphism of `Z^2`.  Hence

`K0(A_F)=lim(Z^2,F) = Z^2`

as an abelian group.  With `phi^2=phi+1`, the canonical state is

`tau_F([x,n]) = (phi*x_0+x_1)/phi^(n+1)`.

It sends the rooted order unit `[e0,0]` to 1.  Its range is exactly
`Z[phi]`, because powers of `phi` are units of that ring.  The ordered group
is therefore

`(Z[phi], Z[phi] intersect R_{>0} union {0}, 1)`.

This is **DERIVED**, not imported from the literature: unimodularity, the
characteristic polynomial `x^2-x-1`, all displayed floors, and the PF
functional are finite certificates in the verifier.

### Rooted McKay: dyadic trace plus infinitesimals

It would be wrong to use the full stationary system `lim(Z^9,A)`: it contains
a second parity component not reached by the chosen root.  On even floors the
active vertices are `(rho0,rho2,rho4,rho6,rho8)` and on odd floors they are
`(rho1,rho3,rho5,rho7)`.  Let `C` be the resulting `5 x 4` incidence matrix.
After telescoping two floors,

`B=C C^T =`

```text
1 1 0 0 0
1 2 1 0 0
0 1 2 1 1
0 0 1 2 1
0 0 1 1 1
```

The exact certificates are

`char_B(x)=x(x-4)(x-1)(x^2-3x+1)`,

`SNF(B)=diag(1,1,1,1,0)`,

`ker(B)=Z*(-1,1,-1,0,1)`.

The image is saturated.  In the basis given by the first four columns of
`B`, the stable injective map is

```text
M = 1 1 1  1
    1 2 0 -1
    0 1 3  2
    0 0 1  2
```

with `det(M)=4` and `SNF(M)=diag(1,1,1,4)`.  Thus the honest abstract answer
is the exact stationary presentation

`K0(A_M)=lim(Z^4,M)=union_(k>=0) M^(-k) Z^4 subset Q^4`.

The primitive positive row `l=(1,3,5,4)` satisfies `lM=4l`.  With the stable
rooted order unit represented by `e0`,

`tau_M([z,k]) = l z / 4^k`.

Its range is exactly `Z[1/2]`.  The equality follows already from finite
floor data: even active irrep dimensions have gcd 1, odd active dimensions
have gcd 2, and the denominator at floor `n` is `2^n`.

The trace kernel is rank 3.  An explicit unimodular change of basis makes the
matrix block lower triangular with quotient block `[4]` and kernel block of
determinant 1 and characteristic polynomial
`(x-1)(x^2-3x+1)`.  Therefore there is an exact sequence

`0 -> Z^3 -> K0(A_M) -> Z[1/2] -> 0`.

No splitting is asserted.  This restraint matters: the finite matrix has a
nontrivial lower-left extension term, so writing `Z[1/2]+Z^3` without an
explicit integral splitting would overclaim.

The verifier has `B^4>0`; its other eigenvalues have absolute value below 4.
Consequently the rooted stationary order is the strict trace order:

`K0(A_M)^+ = {0} union {g: tau_M(g)>0}`.

The normalized `tau_M` is the canonical state.  The rank-3 kernel consists of
infinitesimals invisible to it.  These statements concern this exact rooted
tower; no physical scale interpretation is made.

### The `Z[phi]` pattern check

- Fibonacci: the mass-exponent lattice `Z[phi]` agrees canonically with the
  trace lattice because both come from the same fusion polynomial and rooted
  state.  This algebraic embedding is **DERIVED**.  Identifying its order or
  levels with physical masses or energy scales remains **STRUCTURAL**.
- McKay: the trace range is the rational dyadic group `Z[1/2]`.  It cannot
  contain `phi`, so there is no trace-preserving embedding of the golden
  lattice.  Abstract rank-two subgroup embeddings into the rank-four group
  would require choices and do not respect the canonical state.  The result
  is **DERIVED negative**; any claimed golden match here is rejected as a
  **PATTERN**.

## 3. Finite-level matter test

### The algebra-size gate

Matrix sizes grow, so the old statement “the commutant has no block of size
3” is no longer true.  For example `rho2` has multiplicities `3,9,28,91,309`
at levels `4,6,8,10,12`.  Inflation genuinely passes that size gate.

The requested canonical seed test is stricter.  A unital real representation
of `H` in `M_m(C)` requires even `m`, and a unital representation of
`M3(C)` in `M_m(C)` requires `3|m`.  The verifier tests these divisibilities
at levels 2 through 12.  More decisively,

- `rho0`, `rho2`, and `rho8` are even nodes;
- `rho1` is an odd node;
- `A^n e0` has support on only one parity.

Therefore `rho0/C` and either color seed can never coexist on a floor with
`rho1/H`.  **DERIVED all-level no-go for the stated floor algebra.**

There is also a unital-allocation issue familiar from the smooth-fiber audit:
the three selected factors do not by themselves specify how the algebra unit
acts on every other block.  The parity no-go occurs before that choice is
needed.

### Consecutive-floor shift and first order

Let a cumulative path model carry independent endpoint actions on its floor
summands, and let `S:H_n -> H_(n+1)` be any nonzero incidence/shift block.
Choose `a` to be the source-floor left projector and `b^op` the target-floor
right projector.  Direct block multiplication gives

`[[S+S^dagger,L(a)],R(b)^op] = +/- S`.

The verifier contains the minimal exact `1 x 1` shift witness; tensoring it
with any nonzero incidence block preserves nonvanishing.  Meanwhile the level
grading gives `gamma D=-D gamma`.  Hence

`gamma-odd + canonical independent endpoint first order ==> S=0`.

This is **DERIVED** for that cumulative endpoint representation.  The richer
matrix sizes do not help because projectors distinguishing consecutive floors
remain in the relative endpoint algebra.

The phrase “Christensen--Ivan style” must be scoped carefully.  Standard AF
spectral-triple constructions motivate using the filtration, but they do not
derive this shift.  Their Dirac operators are normally assembled from
filtration-difference projections and chosen coefficients.  Here the
sequence `c_n`, a Hilbert-space representation, an opposite action, and a
real structure are all **NOT DERIVED**.  No claim about an AF-limit Dirac is
made from the finite witness alone.

### Hypercharge and Route C

There is no floor carrying the requested derived algebra and no nonzero
first-order shift for the canonical cumulative endpoint action.  Thus no
commutant generator `Y` is licensed.  Generation blindness, doublet
constancy, and Route C anomaly forcing are undefined, not failed anomaly
equations.  **DERIVED scope limit.**

## 4. Status ledger

### Strengthened

- **DERIVED:** both rooted inclusion systems and every floor through 12.
- **DERIVED:** Fibonacci ordered `K0` is the golden dimension group.
- **DERIVED:** rooted McKay `K0` has the explicit stable matrix `M`, dyadic
  canonical trace, and rank-3 infinitesimal kernel.
- **DERIVED:** McKay multiplicities grow without bound and remove the old
  matrix-size obstruction.
- **DERIVED negative:** the golden trace lattice does not embed in the McKay
  trace range.
- **DERIVED all-level negative:** the three canonical matter seeds never
  coexist on one rooted McKay floor.
- **DERIVED scoped negative:** a nonzero consecutive-floor shift fails first
  order for canonical independent endpoint actions.

### Downgraded or rejected

- Large `M2`/`M3`-capable blocks are availability, not a selected unital
  Standard-Model algebra.
- A cumulative-floor replacement is **STRUCTURAL**, not the requested level
  algebra.
- The coefficients `c_n` and any level-to-physical-scale identification are
  **NOT DERIVED / STRUCTURAL**.
- A Christensen--Ivan analogy does not by itself provide a gamma-odd shift.

### Open

- a separately derived AF representation with a compatible opposite action
  in which a nonzero odd first-order operator survives;
- a geometric reason to combine parity-separated floors into one matter
  bimodule;
- a real structure, color orientation, generation blocks, and a canonical
  integral `Y` for such a representation.

The finite house gains a real staircase, but the staircase does not yet carry
matter.
