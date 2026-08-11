# Exact Whitney trace stiffness improves, but does not solve, refinement balance

Date: 2026-08-11

Preregistration commit: `b9a4104`  
Paired-control protocol correction: `a92c911`

Targeted verifier:
`reproducible/verify_whitney_trace_stiffness.py`

Targeted result: **13/13 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Result

The exact Whitney trace-jump energy is geometrically and dimensionally better
than the unweighted copy Laplacian, but it is not an exact refinement fixed
point.

> **DERIVED:** under uniform dilation, the mass-orthonormal exact trace
> stiffness and the Kähler--Dirac operator both scale as (h^{-1}) in every
> constrained form degree.

> **DERIVED NUMERICAL NEGATIVE:** on the complete refined control, the three
> degreewise first-step ratios are not one.  Exact first-step refinement
> covariance fails.

> **PATTERN:** relative to the unweighted penalty rebuilt on the same control,
> exact trace weighting reduces the degree-ratio spread from (47.65) to
> (3.05).

This is enough to justify a complete 600-cell calculation.  It is not enough
to claim a refinement law or physical scale.

## Why this operator is canonical

For every shared triangle (F) and form degree (p=0,1,2), the jump of the
two element traces is measured with the exact Whitney Gram matrix on that
face:

\[
 B_p=\sum_F R_{F,p}^*M_{F,p}R_{F,p}.
\]

There is no diagonal mass lumping, fitted degree factor, inverse face-size
power, or selected constraint basis.  The candidate count is exactly one,
up to the same unresolved overall scalar (kappa).

All shared-face metrics agree exactly when computed from either parent.  The
base control has one face metric type; the barycentric refinement has four.
Every exact face mass is positive definite.

Because the occurrence graphs are connected and every face mass is positive,

\[
 \ker B=\operatorname{im}J
\]

exactly.  Complete dense generalized diagonalization confirms the expected
nullities in all six degree/level cases.

## Exact dimensional result

For a dilation by (h), exact symbolic integration gives

\[
 M_{T,p}\mapsto h^{3-2p}M_{T,p},
 \qquad
 M_{F,p}\mapsto h^{2-2p}M_{F,p}.
\]

Therefore

\[
 M_{T,p}^{-1/2}B_pM_{T,p}^{-1/2}\sim h^{-1}
\]

for every (p).  The local mass-orthonormal Kähler--Dirac norm independently
returns the exact numerical ratio (1/2) when all lengths double.

This repairs the dimensional defect of using one raw unit weight across all
form degrees.  It does not guarantee invariance under barycentric refinement,
which changes both element shape and occurrence-graph topology.

## Complete paired control

The two complete carriers are:

| level | f-vector | tetrahedra | duplicated dimension |
|---|---|---:|---:|
| base | `(5,10,10,5)` | 5 | 75 |
| first barycentric | `(30,150,240,120)` | 120 | 1,800 |

The exact copy-constraint ranks are respectively

\[
 (15,20,10)
 \quad\text{and}\quad
 (450,570,240).
\]

The maximum generalized eigenpair residual over both the exact-trace and
unweighted spectra is (1.30\times10^{-14}), well below the preregistered
(10^{-9}) gate.

## Positive gaps

| degree | base exact trace | refined exact trace | base unweighted | refined unweighted |
|---:|---:|---:|---:|---:|
| 0 | 8.660254 | 2.293663 | 7.500000 | 32.233679 |
| 1 | 5.196152 | 2.351464 | 18.000000 | 6.752486 |
| 2 | 3.464102 | 2.802834 | 12.000000 | 1.082291 |

All gaps remain positive.  The stiff approximation mechanism therefore
survives the finite refinement on this control.

## First-step ratios

For each degree define

\[
 s_{r,p}=\frac{a_r}{g_{r,p}},
 \qquad
 R_p=\frac{s_{1,p}}{s_{0,p}},
\]

where (a_r) is the local Kähler--Dirac norm and (g_{r,p}) the positive
constraint gap.

The exact-trace ratios are

\[
 (R_0,R_1,R_2)
 =
 (11.6371,,6.81065,,3.80924).
\]

Their spread is

\[
 \frac{\max R_p}{\min R_p}=3.05498.
\]

The unweighted penalty, reconstructed on the same carriers, gives

\[
 (0.717127,,8.21587,,34.1729)
\]

with spread

\[
 47.6525.
\]

The comparison is paired and therefore valid.  The protocol correction
`a92c911` was committed before computation after noticing that comparing with
the earlier complete-600-cell spread would have mixed two different
geometries.

## Honest interpretation

The trace metric removes most of the catastrophic degree mismatch, but all
three ratios remain substantially above one and differ from one another.
Thus:

- exact uniform dimensional scaling is **DERIVED**;
- exact barycentric refinement covariance is **REFUTED on the control**;
- improved degree balance is a **PATTERN**, because only one refinement and
  one closed control have been tested.

The remaining factor cannot yet be called a renormalization-group flow.  A
single step provides no exponent, fixed point, or continuum limit.

The likely mathematical source is now isolated: barycentric subdivision
changes occurrence-graph topology in addition to local length.  Exact face
mass repairs local dimensional scaling but cannot by itself cancel that
global combinatorial change.

## Physical limitations

The term (kappa B) remains positive and degree-preserving.  Consequently it
breaks exact Kähler--Dirac oddness at finite stiffness.  It also leaves one
overall dimensionless (kappa) unselected.

Nothing here derives time, inertia, a causal cone, (c), or a Planck scale.
The result only identifies a much better microscopic conformity energy.

## Status ledger

- **DERIVED:** unique exact Whitney (L^2) trace-jump bilinear form.
- **DERIVED:** exact face metric agreement and positivity.
- **DERIVED:** exact conforming kernel at both control levels.
- **DERIVED:** uniform (h^{-1}) scaling in all degrees.
- **DERIVED NUMERICAL:** complete base/refined control spectra.
- **DERIVED NUMERICAL NEGATIVE:** exact first-step covariance fails.
- **PATTERN:** paired degree spread improves from `47.6525` to `3.05498`.
- **OPEN:** complete 600-cell trace-stiffness spectra.
- **OPEN:** repeated-refinement behavior.
- **OPEN:** the overall dimensionless stiffness.
- **OPEN:** a chiral local realization of the same constraint energy.
- **NOT CLAIMED:** physical time, mass, inertia, (c), (hbar), Newton's
  (G), or a Planck scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_trace_stiffness.py
```

Expected result: `13/13`.
