# Whitney trace stiffness is a separator, not the missing dynamics

Date: 2026-08-11

Preregistration commit: `7528f97`

Targeted verifier:
`reproducible/verify_whitney_trace_penalty_role.py`

Targeted result: **8/8 PASS**.  The verifier is registered exactly once.  The
full suite was not run by explicit user request.

## Headline

The exact trace-jump term has two sharply different roles.

> **DERIVED NEGATIVE FOR DYNAMICS:** its compression to the conforming
> Whitney space is identically zero at every refinement level.  Its positive
> spectrum is not a dispersion spectrum of physical assembled fields.

> **DERIVED STRUCTURAL POSITIVE FOR SEPARATION:** on the canonical
> shape-regular tower, its positive gap has the same `1/h` engineering scale
> as the local Kähler--Dirac operator.  Therefore a finite dimensionless
> stiffness can, in principle, keep the conforming and mismatch sectors
> uniformly separated under refinement.

The term is a constraint stabilizer.  It is not the missing law of motion and
does not select the value of its coefficient.

## Exact zero-compression theorem

For one form degree let

\[
 J_h:V_h^{\rm conf}\longrightarrow V_h^{\rm disc}
\]

copy a global Whitney coefficient into every incident tetrahedron, and let
`R_h` subtract the two copies across every shared face.  By construction,

\[
 R_hJ_h=0.
\]

With the positive face-trace Gram matrix `H_h`, define

\[
 B_h=R_h^*H_hR_h.
\]

It follows exactly that

\[
 B_hJ_h=0,
 \qquad
 J_h^*B_hJ_h=0.
\]

Moreover,

\[
 x^*B_hx=\lVert H_h^{1/2}R_hx\rVert^2,
\]

so positivity of `H_h` gives `ker B_h=ker R_h`.  Every occurrence graph is
connected, hence

\[
 \ker R_h=\operatorname{im}J_h.
\]

The verifier constructs `R_h` and `J_h` independently as integer sparse
matrices for all nine combinations of

\[
 k=1,2,4,qquad p=0,1,2,
\]

and obtains exactly zero nonzeros in every product `R_h J_h`.

## This is zero on nonconstant physical fields too

At each level, take a global scalar cochain equal to one at a chosen vertex
and zero elsewhere.  Its assembled copy has:

\[
 R_hJ_hx=0,
\]

but its simplicial derivative is nonzero on 14 incident edges:

\[
 d_0x\ne0.
\]

Thus the penalty does not merely annihilate constants or harmonic modes.  It
annihilates the entire conforming space, including fields with nonzero spatial
variation.

Consequently it cannot be identified with the conforming exterior derivative,
Kähler--Dirac operator, Hodge Laplacian, or a physical gradient energy.  Its
principal/tangent action **restricted to the physical conforming sector** is
zero.

This scope qualification matters.  The penalty has a nonzero symbol on the
discontinuous mismatch branch.  Standard discontinuous Galerkin methods use
jump penalties together with volume derivatives and consistency fluxes, not
as a replacement for the differential operator; see Arnold, Brezzi, Cockburn
and Marini, [Unified Analysis of Discontinuous Galerkin Methods for Elliptic
Problems](https://doi.org/10.1137/S0036142901384162).  The compatible
conforming differential complex is the FEEC/Whitney structure described by
Arnold, Falk and Winther, [Finite element exterior calculus, homological
techniques, and applications](https://doi.org/10.1017/S0962492906210018).

## Complete kernel counts

| `k` | degree | local dimension | conforming dimension | exact jump rank | maximum occurrences |
|---:|---:|---:|---:|---:|---:|
| 1 | 0 | 480 | 30 | 450 | 24 |
| 1 | 1 | 720 | 150 | 570 | 6 |
| 1 | 2 | 480 | 240 | 240 | 2 |
| 2 | 0 | 3,840 | 180 | 3,660 | 24 |
| 2 | 1 | 5,760 | 1,140 | 4,620 | 6 |
| 2 | 2 | 3,840 | 1,920 | 1,920 | 2 |
| 4 | 0 | 30,720 | 1,320 | 29,400 | 24 |
| 4 | 1 | 46,080 | 9,000 | 37,080 | 6 |
| 4 | 2 | 30,720 | 15,360 | 15,360 | 2 |

Every rank is exactly

\[
 \dim V_h^{\rm disc}-\dim V_h^{\rm conf}.
\]

All exact principal minors of all encountered face masses are positive.  The
reported smallest numerical face-mass eigenvalues scale as

\[
 0.0196419,\quad0.00491046,\quad0.00122762,
\]

which is the expected factor four for the degree contributing the global
minimum; the positivity conclusion itself is exact.

## All-level separator theorem

This part is not inferred from the three numerical levels.  It follows from
four explicit hypotheses already established for the rank-edgewise tower.

### 1. Quasi-uniform shape-regular geometry

After the single barycentric rank step, every `Esd_k` child has the same
volume scale and belongs to one of at most three normalized congruence
classes.  Therefore its diameter is comparable to `h=1/k` with constants
independent of `k`.

### 2. Bounded occurrence graphs

Edgewise subdivision has only finitely many local neighbourhood types away
from the boundary of a parent chamber.  Boundary neighbourhoods are induced
from lower-dimensional edgewise subdivisions.  Since the original
barycentric complex is finite, gluing across its faces multiplies these
bounds by a fixed finite incidence number.  Thus the number of tetrahedra
incident to any refined simplex is bounded independently of `k`.

The complete control gives the stable finite counts

\[
 q_0=24,qquad q_1=6,qquad q_2=2
\]

at all three tested resolutions.  These are controls, not the proof and not
the corresponding numerical bound for the full 600-cell.

### 3. Exact affine powers

For Whitney `p`-forms in three dimensions,

\[
 M_{T,p}\sim h^{3-2p},
 \qquad
 H_{F,p}\sim h^{2-2p}.
\]

The powers are certified exactly for every `p=0,1,2`; uniform
shape-regularity turns the proportionalities into level-independent two-sided
matrix bounds.

### 4. Uniform graph singular values

There are only finitely many connected occurrence graphs of bounded order.
Their nonzero incidence singular values therefore have a positive common
minimum.

Combining the four points gives, for the smallest positive generalized
penalty eigenvalue,

\[
 c_p h^{-1}\le g_{h,p}\le C_p h^{-1}.
\]

The local first-order Kähler--Dirac norm satisfies

\[
 a_h\le C_Dh^{-1}.
\]

Hence

\[
 \sup_h\frac{a_h}{g_{h,p}}<\infty.
\]

This proves boundedness, not a limiting numerical value.

## What the computed ratios actually mean

The previous finite-stiffness separation theorem uses the sufficient
condition

\[
 \kappa g_h>2a_h.
\]

The new values therefore measure a sufficient stiffness threshold, not a
physical propagation constant:

| `k` | `h a_h` | `h g_0` | `h g_1` | `h g_2` | worst `2a/g` |
|---:|---:|---:|---:|---:|---:|
| 1 | 11.9369 | 2.29366 | 2.35146 | 2.80283 | 10.4086 |
| 2 | 13.1168 | 1.79877 | 1.91984 | 2.78614 | 14.5842 |
| 4 | 13.1168 | 1.59494 | 1.68318 | 2.80647 | 16.4481 |

The finite thresholds increase over the tested range.  It would be fitting to
declare a limit near 16 or choose `kappa=16.5` from this table.  No such claim
is made.  The theorem supplies only an unspecified finite all-level bound.

Therefore:

- refinement does not force `kappa_h` to diverge merely to preserve a
  separated mismatch sector;
- exact conformity at a fixed finite level still requires
  `kappa -> infinity`;
- geometry still does not select a particular finite `kappa`.

## Consequence for the coupled pencil

For

\[
 W_h+\kappa B_h,
\]

the exact conforming compression is independent of `kappa`:

\[
 J_h^*(W_h+\kappa B_h)J_h=J_h^*W_hJ_h.
\]

At finite stiffness, off-diagonal mixing with the mismatch sector remains.
The penalty can suppress that mixing; it cannot generate new tangent
dynamics because its tangent block is zero.

This reconciles the two earlier results:

1. finite stiffness is a legitimate local mechanism for isolating an
   approximately conforming slow sector;
2. the physical evolution within the exact conforming sector still comes
   entirely from the Whitney/Kähler--Dirac term.

## Physical status

- **DERIVED:** `R_hJ_h=0` and `ker B_h=im J_h` at every level.
- **DERIVED NEGATIVE:** trace stiffness is not a physical spatial derivative
  on assembled fields.
- **DERIVED:** the positive spectra measured so far belong to the mismatch
  sector.
- **DERIVED STRUCTURAL POSITIVE:** `a_h/g_h` is uniformly bounded under all
  stated rank-edgewise hypotheses.
- **STRUCTURAL:** a finite dimensionless stiffness can maintain sector
  separation under refinement.
- **OPEN:** a geometry-selected value of `kappa`.
- **OPEN:** convergence of the complete finite-`kappa` pencil.
- **OPEN:** a local chiral dynamics, Lorentzian time, causal speed, inertia,
  mass, and Planck units.

The next question is therefore no longer “does the penalty flow to a physical
operator?”  It is:

> Does the complete local pencil at fixed dimensionless stiffness converge to
> the assembled Kähler--Dirac dynamics on its low branch, or are additional
> consistency/flux terms required?

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_trace_penalty_role.py
```

Expected result: `8/8`.

