# Element-local Whitney dynamics leaks completely at assembly

Date: 2026-08-11  
Preregistration commit: `162ce61`

## Result

Duplicating every cochain degree of freedom into its incident tetrahedra makes
the Whitney mass block-local, but the resulting elementwise codifferential
does not preserve the conforming assembled cochains.

> **DERIVED ASSEMBLY-LEAKAGE NO-GO:** all three downward maps leak out of the
> conforming subspace.  The invariance hit fraction is `0/3`.

The leakage maps have exact ranks

\[
(720,1200,600),
\]

equal to their complete input dimensions.  Thus every nonzero edge, triangle
or tetrahedron direction has a nonconforming component after the raw local
codifferential; this is not a small exceptional sector.

There is also a constructive result:

> **DERIVED PROJECTED FACTORIZATION:** metric-orthogonal assembly repairs all
> three degrees exactly.  The global consistent Whitney adjoint is precisely
> element-local evolution followed by the `M_loc`-orthogonal projection onto
> assembled cochains.

Simple equality/Grover averaging of copies fails the Whitney metric in all
three degrees.  The glue must know the mass matrix, not merely which copies
represent the same simplex.

The targeted verifier passes `11/11` in about one second.  No full suite was
run.

## Carrier and exact local data

Each of the 600 tetrahedra carries its 15 nonempty faces.  The duplicated
degree dimensions are

\[
(2400,3600,2400,600),
\]

for a total carrier dimension of 9,000.  The global conforming dimensions
remain

\[
(120,720,1200,600).
\]

The exact occurrence multiplicities are:

| global simplex | copies in incident tetrahedra |
|---|---:|
| vertex | 20 |
| edge | 5 |
| triangle | 2 |
| tetrahedron | 1 |

The local mass matrices are independently integrated from the defining
Whitney forms.  On each congruent regular tetrahedron, all three local metric
adjoints are exactly incidence-proportional:

\[
(\delta_0^{\rm loc},\delta_1^{\rm loc},\delta_2^{\rm loc})
=
\left(rac54d_0^T,\frac52d_1^T,\frac{15}{4}d_2^T\right).
\]

Therefore the failure does not originate in an exotic dense operator inside
one regular element.  It appears when separately computed element
contributions must agree across shared simplices.

## Upward versus downward maps

Restriction to local copies intertwines the exterior derivative exactly:

\[
d_p^{\rm loc}J_p=J_{p+1}d_p,
\qquad p=0,1,2.
\]

All three integer residuals are identically zero.  Exterior differentiation
therefore preserves conformity without a solve.

The adjoint direction behaves differently.  Exact difference operators
compare every copy with a fixed copy of the same global simplex.  Applied to

\[
\delta_p^{\rm loc}J_{p+1},
\]

they give:

| downward degree | exact leakage nonzeros | exact rank | input dimension |
|---:|---:|---:|---:|
| 1-form to 0-form | 10,800 | 720 | 720 |
| 2-form to 1-form | 8,640 | 1,200 | 1,200 |
| 3-form to 2-form | 2,400 | 600 | 600 |

The ranks are certified exactly: sparse structural-rank upper bounds equal
ranks computed over the finite field with prime `1000003`, sandwiching the
rational rank at the same value.

Representative exact witnesses are:

\[
\frac54,\qquad -\frac52,\qquad \frac{15}{4}
\]

for degrees zero, one and two respectively.  The full matrices—not merely
these witnesses—are included in the rank and support census.

## Why metric projection repairs it

Assemble the global masses by

\[
M_p=J_p^TM_p^{\rm loc}J_p.
\]

The verifier checks exactly, after clearing rational denominators,

\[
J_p^TM_p^{\rm loc}\delta_p^{\rm loc}J_{p+1}
=d_p^TM_{p+1}
\]

for all three degrees.  There are zero integer residual entries in every
case.

Consequently, if `P_p^M` denotes the `M_loc`-orthogonal projection onto
`im(J_p)`, then

\[
P_p^M\delta_p^{\rm loc}J_{p+1}=J_p\delta_p,
\]

where

\[
\delta_p=M_p^{-1}d_p^TM_{p+1}
\]

is the global consistent Whitney adjoint.

This is an exact factorization of the missing operation.  It also identifies
the problem: applying a projection is not reversible, while reflecting about
its range requires implementing the metric projector whose formula contains
the assembled inverse mass.

## The cheap Grover glue is insufficient

The Euclidean equality projector

\[
P_p^{\rm eq}=J_p(J_p^TJ_p)^{-1}J_p^T
\]

is attractive because it decomposes into small uniform averages over the
20, 5, 2 or 1 copies of each simplex.  Its reflection is a canonical local
Grover operation.

However, compressing the element-local adjoint with this projector fails the
defining global Whitney identity.  The exact cleared residual nonzero counts
are

\[
(10080,18000,9600).
\]

Thus equality of copies is not enough.  The correct glue must encode the
off-diagonal local mass correlations as well.

## Physical interpretation

The result sharpens the locality problem substantially:

```text
local Whitney evolution
        |
        v
nonconforming element copies
        |
        v
mass-weighted glue/projection  --->  global Whitney Kähler--Dirac
```

The first arrow is exactly local but leaks.  The second restores the accepted
metric exactly but is not yet a reversible local tick.

A promising continuation is to replace instantaneous projection with a
local **constraint/consensus Hamiltonian** whose kernel is precisely the
conforming subspace.  Nonconforming modes would then acquire a gap rather
than being deleted.  In a large-gap or long-consensus limit, the low-energy
sector could approach the metric projection while the microscopic evolution
remains unitary.

That continuation introduces a crucial selection question: whether the
constraint strength and local mass weighting are fixed by the geometry or
are new fitted parameters.  Exact projection may require an infinite-time or
infinite-penalty limit; this must be proved rather than assumed.

## Status ledger

- **DERIVED:** exact 9,000-dimensional duplicated element carrier.
- **DERIVED:** all upward Whitney differentials preserve assembly.
- **DERIVED NEGATIVE:** downward element adjoints preserve `0/3` conforming
  subspaces.
- **DERIVED:** all three leakage maps have full input rank.
- **DERIVED:** mass-orthogonal projection reproduces all three global Whitney
  adjoints exactly.
- **DERIVED NEGATIVE:** simple equality/Grover glue fails in `3/3` degrees.
- **STRUCTURAL:** interpreting a local constraint gap as the physical
  mechanism selecting conforming modes.
- **OPEN:** a coefficient-free reversible realization of mass-weighted glue.
- **OPEN:** whether exact assembly needs an infinite limiting process.
- **NOT CLAIMED:** a refined Dirac cone, mass, inertia or physical `c`.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_element_local_assembly.py
```

Expected result: `11/11`.
