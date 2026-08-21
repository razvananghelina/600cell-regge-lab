# Prior-art gate: constrained refined H4 boundary Hessian

Date: 2026-08-21

Status: completed before constructing any constrained effective matrix.

## Exact object and complete hypotheses

Use the internally stationary Lorentzian product over

```text
K0=P(sd K_600),  f=(2640,17040,28800,14400),  tau0=0.0102,
```

with the frozen projected rank geometry, all 24 colour-ordered staircase
triangulations, the corrected Regge action and boundary terms, and the
curvature-selected rank masses `m_r=K_r/(8*pi)`.  In the fixed total-orbit
log-squared-edge coordinates, write the symmetric action Hessian as

```text
H = [[A,B],
     [B^T,C]],

b in R^12,  i in R^10.
```

The accepted computations establish, inside certified numerical envelopes,

```text
inertia(C)=(9 positive,1 zero-compatible,0 negative),
ker(C)=span(n),
c=B n != 0,
c=(1/2)g_boundary,
```

where `n` is the analytically fixed product-duration tangent.  The new object
is the stationary quadratic boundary form only on

```text
S=ker(c^T) subset R^12.                            (1)
```

No unconstrained extension of that form, schedule average, colour-order
selection, continuum target, nonhomogeneous mode, speed or physical constant
is part of the object.

## Finite-dimensional derivation

For a boundary perturbation `b`, the linear internal equation is

```text
C i = -B^T b.                                     (2)
```

Because `C` is symmetric and `ker(C)=span(n)`, its range is `ker(n^T)`.
Equation (2) is soluble exactly when

```text
n^T B^T b = c^T b = 0,
```

which is (1).  If `i_1` and `i_2` are two solutions, their difference is a
multiple of `n`.  On (1), the stationary quadratic value is independent of
that difference because

```text
b^T B(i_1-i_2) = alpha*b^T B n = alpha*c^T b = 0.
```

Thus the bilinear form

```text
q_eff(b_1,b_2)
 = b_1^T A b_2 + b_1^T B i(b_2),  b_1,b_2 in S,   (3)
```

is well defined without selecting a physical gauge.  A Moore--Penrose formula
`A-B C^+ B^T` may represent (3) in the frozen Euclidean coordinates, but its
values away from `S` are not part of the result.  The computation will avoid
that extension and solve a bordered/restricted nonsingular system instead.

Because the base point is not stationary in the boundary variables, (3) is
most safely interpreted as the **linearized boundary-momentum response,
restricted to admissible boundary tangents and modulo the conormal `c`**.  It
is not automatically the intrinsic Hessian of an unknown nonlinear boundary
constraint surface: that stronger interpretation would also require a finite
admissible surface (or its second fundamental form).  No such surface is
assumed here.

## KNOWN from primary literature

- Carlson, Haynsworth and Markham define generalized Schur complements using
  the Moore--Penrose inverse and derive their rank relations:
  *A Generalization of the Schur Complement by Means of the Moore--Penrose
  Inverse*, SIAM J. Appl. Math. 26 (1974), 169--175,
  [DOI 10.1137/0126013](https://doi.org/10.1137/0126013).
- Maddocks treats quadratic forms restricted to a subspace and relates their
  inertia to generalized Schur complements: *Restricted quadratic forms,
  inertia theorems, and the Schur complement*, Linear Algebra Appl. 108
  (1988), 1--36,
  [DOI 10.1016/0024-3795(88)90177-2](https://doi.org/10.1016/0024-3795(88)90177-2).
- Dittrich and Hoehn formulate singular discrete dynamics on constraint
  surfaces and stress that preliminary null directions need not extend to
  Hessian gauge directions; Hamilton's principal function is defined by
  extremizing the bulk action: *Canonical simplicial gravity*, especially
  sections 2 and 6.2, [arXiv:1108.1974](https://arxiv.org/abs/1108.1974),
  [DOI 10.1088/0264-9381/29/11/115009](https://doi.org/10.1088/0264-9381/29/11/115009).
- Dittrich and Hoehn derive canonical constraints and consistency conditions
  from singular Regge dynamics: *From covariant to canonical formulations of
  discrete gravity*, [arXiv:0912.1817](https://arxiv.org/abs/0912.1817),
  [DOI 10.1088/0264-9381/27/15/155001](https://doi.org/10.1088/0264-9381/27/15/155001).

These sources justify constrained reduction as standard mathematics and
discrete canonical mechanics.  They do not evaluate (3) on this projected
600-cell carrier or compare its 24 staircase triangulations.  Search absence
does not establish novelty; external novelty remains **OPEN**.

## CONTROL already known in the repository

- **DERIVED COMPUTATIONAL:** all 24 internal blocks have one isolated
  zero-compatible eigenvalue and nine positive eigenvalues.
- **DERIVED COMPUTATIONAL / STRUCTURAL:** the analytic product tangent spans
  the numerical null line.
- **DERIVED COMPUTATIONAL, adversarially corroborated:** the nonzero common
  compatibility row is `c=B n=tau*C_spatial/8=(1/2)g_boundary`.
- **DERIVED NEGATIVE:** an ordinary unconstrained Schur complement is not
  defined.

These facts license (3), but do not determine its schedule class.

## Proposed representation and invariance test

Freeze a full-rank boundary matrix `P` whose columns span `S`, and an internal
matrix `Q` whose columns span `ker(n^T)`.  Solve

```text
(Q^T C Q) Y = -Q^T B^T P,
K_S = P^T(A P+B Q Y).                             (4)
```

The `9 x 9` matrix in (4) is nonsingular under the frozen inertia hypothesis.
The `11 x 11` matrix `K_S` represents (3).  Changing `P` gives a congruent
matrix; changing `Q` changes the internal representative only by the null
line and must leave (3) unchanged.  Both changes will be tested explicitly.

The primary implementation will use algebraic pivot bases, not fitted or
eigenvector-selected bases:

```text
boundary pivot: old_12 (index 3),
internal pivot: rho_3  (index 9).
```

Alternative pivots `new_12` and `rho_0` are frozen as basis-independence
controls.  The boundary pivots are already known nonzero from the accepted
coupling artifact; the internal pivot values are exactly one.

## Framing attack

1. Calling `K_S` a Hamiltonian constraint would outrun the calculation.  The
   constraint here is a finite invariant-sector solvability condition.
2. A single schedule class would show only that the bare staircase choice is
   invisible to this constrained `H4` quadratic form.  It would not imply
   triangulation independence outside this 11-dimensional sector.
3. A multiple-class result is a clean negative for canonical bare staircase
   evolution on this carrier.  It does not refute Regge calculus or an
   improved/perfect action.
4. The fixed dot products used to write pivot complements are coordinate
   devices.  The evidential object is equality of the restricted bilinear
   forms, verified under explicit changes of both bases.
5. Even a positive result cannot yield `c`: no nonhomogeneous spatial mode or
   dispersion relation appears in the carrier.
6. The linear compatibility hyperplane need not integrate to an exact
   nonlinear eleven-dimensional family.  The present gate tests linearized
   canonical response only; nonlinear integrability remains separate and
   **OPEN**.

## OPEN and next gate

Preregister a target-free, high-precision construction of (4) for all 24
schedules, including complete matrices, error envelopes, time reversal,
class census, direct-action directional checks, corruption controls and the
two basis changes.  Only after the primary result may an adversarial route be
designed.
