# Prior-art and framing gate: kinetic--stiffness generalized modes

Date: 2026-08-18

## Why the previous phase lift failed

The Euclidean negative projector of the stiffness form `A=-V_S` is common
across the two centered slices, but its naive phase lift `E+E*` is not closed:
the `B,D` tangent blocks leak in all cells.  This is not surprising once the
kinetic form is included.  The dynamical mode problem is the Hermitian definite
pencil

```text
A v = lambda B v,     B=-M_S > 0,
```

not the Euclidean eigensystem of `A` alone.

## Primary numerical literature checked

1. The LAPACK Users' Guide section [*Generalized Symmetric Definite
   Eigenproblems*](https://www.netlib.org/lapack/lug/node54.html) gives the
   canonical Cholesky reduction of `A v=lambda B v` with Hermitian `A,B` and
   positive `B` to an ordinary Hermitian problem.  This fixes the generalized
   spectral subspace without an optimized basis alignment.
2. Li, Nakatsukasa, Truhar and Xu, [*Perturbation of Partitioned Hermitian
   Definite Generalized Eigenvalue
   Problems*](https://doi.org/10.1137/100808356), study perturbations of
   Hermitian definite pencils and the special role of off-diagonal coupling.
   It supplies error-analysis context, not a result about this 600-cell
   carrier.
3. Davis--Kahan spectral-subspace control remains applicable after the
   Cholesky reduction to a standard Hermitian problem.

No source found applies this exact pencil to the accepted dust-Regge slabs.
The search is not proof of novelty.

## Action-selected closure criterion

The centered recurrence already provides the normalized operators

```text
Gamma = M^-1 N,
Omega = M^-1 V.
```

The generalized pencil subspace is an eigenspace of the `Omega` part when the
action-selected shape carrier reduces the coefficients.  A genuine
configuration-mode subsystem additionally requires `Gamma` to preserve it.
Thus, for its Hermitian projector `P`, test both fixed leakages

```text
(I-P) Gamma P,
(I-P) Omega P
```

on both centered slices.  This asks directly whether the second-order
recurrence closes; it does not invent a momentum lift.

## Status before execution

- **DERIVED COMPUTATIONAL:** the kinetic form is positive definite on all
  target shape carriers and the stiffness inertia is `15 negative + 10
  positive`.
- **KNOWN:** the Hermitian definite pencil has real eigenvalues and a canonical
  Cholesky reduction.
- **OPEN:** whether its negative generalized fiber is common across slices and
  invariant under both normalized recurrence operators.
- **FORBIDDEN:** choosing a non-Cholesky whitening, rotating the old/new fibers
  to improve overlap, or dropping `Gamma` if it mixes the modes.

Even closure would certify only a finite linear recurrence subsystem.  It
would not yet establish particle inertia, a mass shell, a graviton, dispersion
or a continuum limit.
