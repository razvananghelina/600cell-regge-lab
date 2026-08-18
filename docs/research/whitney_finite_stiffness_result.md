# Finite local stiffness recovers Whitney only as a singular low-energy limit

Date: 2026-08-11

Preregistration commit: `68f839b`

Targeted verifier:
`reproducible/verify_whitney_finite_stiffness.py`

Targeted result: **14/14 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Result

There is a genuine constructive result:

> **DERIVED:** the exact assembled Whitney spectrum is the bounded
> low-energy limit of a canonical local finite-stiffness pencil on duplicated
> element cochains.

But it does not yet provide the missing physical tick:

> **DERIVED SINGULAR-LIMIT NEGATIVE:** exact Whitney recovery requires
> (kappa\rightarrow\infty), while the microscopic generator norm grows at
> least linearly with (kappa).  Geometry has not selected a finite value or
> a refinement law for (kappa).

Thus finite stiffness supplies a mathematically honest approximation
mechanism.  It does not derive a finite causal speed, mass scale, or Planck
scale.

## Canonical local construction

On the duplicated element carrier, let (M>0) be the block-local Whitney
metric, (W=W^*) the block-local weak Kähler--Dirac matrix, (J) the
assembly injection, and (C) the complete canonical neighbour-difference
matrix.

No independent constraint basis or spanning tree is selected.  Every row of
(C) compares two copies of the same simplex in tetrahedra that share a
triangle containing it.  The tested pencil is

\[
 (W+\kappa C^*C)v=zMv,
 \qquad \kappa>0.
\]

The all-row Laplacian (C^*C) is invariant under row order and sign.  On the
boundary-of-4-simplex control it has:

- local dimension (75);
- assembled dimension (30);
- 70 canonical neighbour rows;
- exact rank 45 and nullity 30;
- exactly the conforming space as its kernel.

Every cross-tetrahedron nonzero in the penalty is one of the 70 declared
neighbour pairs.  The weak operator itself is tetrahedron-block local.

## Why the Whitney limit is exact

In mass-orthonormal coordinates define

\[
 A=M^{-1/2}WM^{-1/2},
 \qquad
 L=M^{-1/2}C^*CM^{-1/2}.
\]

Then (ker L=M^{1/2}\operatorname{im}J).  Compression to this kernel gives
exactly the assembled generalized eigenproblem

\[
 J^*WJ\,x=zJ^*MJ\,x.
\]

The verifier reconstructs both sides independently.  The weak matrices agree
exactly, and the independent orthonormal spectral comparison has maximum
residual

\[
 4.89\times10^{-15}.
\]

Let

\[
 a=\lVert A\rVert_2=\sqrt{15}\approx3.872983,
 \qquad
 g=\lambda_{+}^{\min}(L)=7.5.
\]

For (kappa g>2a), exactly 30 eigenvalues remain below the stiff sector.
The min--max principle and completing the off-diagonal square give, in
ordered form,

\[
 0\leq \mu_j-\lambda_j(H_\kappa)
 \leq \frac{a^2}{\kappa g-2a},
 \qquad j=1,\ldots,30,
\]

where (mu_j) are the assembled Whitney eigenvalues.  Every dyadic value
from (kappa=2) through (2^{14}) satisfies the preregistered bound.

## Frozen-grid results

No exponent or preferred (kappa) was fitted.  Representative entries from
the preregistered dyadic grid are:

| (kappa) | max spectral error | max angle sine | max constraint residual | microscopic norm |
|---:|---:|---:|---:|---:|
| 2 | 0.371549 | 0.095495 | 0.427067 | 90.1512 |
| 8 | 0.093695 | 0.024185 | 0.108158 | 360.038 |
| 32 | 0.023437 | 0.006051 | 0.027062 | 1,440.01 |
| 128 | 0.005859 | 0.001513 | 0.006766 | 5,760.00 |
| 1,024 | 0.000732 | 0.000189 | 0.000846 | 46,080.0 |
| 16,384 | 0.0000458 | 0.0000118 | 0.0000529 | 737,280.0 |

The complete multiset is stored in
`reproducible/whitney_finite_stiffness.json`.

The convergence is not inferred from the numerical pattern: it follows from
the analytic bound.  The table is an independent finite control.

Every tested finite separated sector remains nonconforming.  Therefore the
construction does not contradict the earlier exact pure-penalty no-go.  It
approaches conformity but never makes the physical subspace invariant at a
finite (kappa).

## The scale problem is now quantitative

The same spectrum satisfies

\[
 \lVert H_\kappa\rVert_2
 \geq \kappa\lambda_{\max}(L)-a,
 \qquad \lambda_{\max}(L)=45.
\]

At the final frozen point the measured norm is approximately 737,280, while
the rigorous lower bound is approximately 737,276.  Better conformity is
therefore bought by an increasingly fast microscopic sector.

This does not prove a relativistic Lieb--Robinson theorem, because a tensor
product of local quantum systems and a physical time unit have not been
derived.  The justified statement is narrower: the family is not uniformly
bounded in the existing normalization, so the exact limit cannot itself be
declared a finite-speed tick.

## Framing attacks

### The parameter is not derived

The coefficients (+1,-1) fix the graph Laplacian, but they do not fix the
relative coefficient (kappa) between the Whitney weak term and the
constraint term.  Selecting (kappa) to reproduce a desired hierarchy would
be forbidden fitting.

### Refinement units are unresolved

Whitney masses of different form degree scale differently with simplex size.
A single algebraic (kappa) has no established common physical dimension or
canonical law under barycentric refinement.  The small-control result cannot
be promoted to a continuum statement until that law is derived.

### Chiral oddness is lost

The Kähler--Dirac term is odd under form parity, whereas (C^*C) preserves
degree.  Their sum is Hermitian but not an odd Kähler--Dirac operator.  Calling
it the theory's fundamental Dirac operator would be incorrect.  It is a
stiff local Hamiltonian whose low-energy compression reproduces Whitney.

## Physical verdict

This is modestly promising as an **effective mechanism**:

```text
finite local duplicated dynamics
           |
           | increasing stiffness
           v
approximately conforming slow sector
           |
           | singular limit
           v
exact assembled Whitney spectrum
```

It does not solve the central problem.  The geometry now supplies the local
coupling pattern and the limiting slow theory, but still does not supply:

1. a finite stiffness scale;
2. a refinement law keeping the physical and stiff sectors controlled;
3. a uniformly bounded causal tick;
4. a Lorentzian interpretation of time.

## Status ledger

- **DERIVED:** canonical all-neighbour positive stiffness.
- **DERIVED:** exact conformity as its kernel.
- **DERIVED:** exact Whitney compression.
- **DERIVED:** target-free spectral convergence bound.
- **DERIVED CONTROL:** all 14 separated dyadic spectra obey the bound.
- **DERIVED NEGATIVE:** every finite tested sector remains nonconforming.
- **DERIVED NEGATIVE:** exact recovery requires
  (kappa\rightarrow\infty).
- **STRUCTURAL NEGATIVE:** the microscopic norm diverges in that limit.
- **STRUCTURAL:** finite stiffness as an effective low-energy mechanism.
- **OPEN:** geometry-selected (kappa) and its units.
- **OPEN:** uniform refinement scaling.
- **OPEN:** a chiral or enlarged local realization avoiding the positive
  degree-preserving term.
- **NOT CLAIMED:** physical time, mass, inertia, (c), (hbar), Newton's
  (G), or a Planck scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_finite_stiffness.py
```

Expected result: `14/14`.
