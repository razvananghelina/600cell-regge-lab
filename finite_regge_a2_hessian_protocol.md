# Preregistered protocol: full finite-Regge `A2` Hessian at the 600-cell

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- NO FINITE HESSIAN RESULT**

## 1. Framing correction

The smooth round metric is not a point of the finite configuration space of
Euclidean Regge tetrahedra. The canonical point in that space is instead the
**equilateral 600-cell Regge metric**. Consequently this audit cannot answer
whether a finite Hessian “at the smooth round point” agrees with the smooth
conformal calculation. Its well-posed question is:

> Is the equilateral 600-cell a minimum, maximum or saddle of the exact
> equal-volume ordinary full-de Rham conical `A2` when all 720 edge lengths
> are varied?

This is a distinct transverse gate, not a discretization proof for the
smooth `l=2` result.

## 2. Complete configuration space and functional

Let `E` be the 720 edges and use squared edge lengths

```text
x=(x_e)_(e in E),  x_e>0.
```

The admissible neighborhood consists of assignments for which the Gram
matrix of every one of the 600 tetrahedra is positive definite. No embedding
of the deformed complex in `R4` is required and no edge orbit is frozen.

For a tetrahedron `t`, let `Vol_t(x)` be its Euclidean volume and
`theta_(t,e)(x)` its interior dihedral angle at `e`. Define

```text
V(x)=sum_t Vol_t(x),
L_e(x)=sqrt(x_e),
beta_e(x)=sum_(t contains e) theta_(t,e)(x),
C(beta)=16*pi^2/(3*beta)+8*beta/3-8*pi,
F(x)=sum_e L_e(x) C(beta_e(x)),
Ahat(x)=V(x)^(-1/3) F(x).
```

The omitted factor `(2*pi^2)^(1/3)` is positive and cannot change gradient,
nullity or inertia. This is exactly the previously licensed open-edge
ordinary full-de Rham coefficient. Vertices first enter the next heat order,
so no vertex term is silently omitted from `A2`.

The audit is local at the equilateral point `x_e=1`. This choice fixes only a
unit: `Ahat` is invariant under `x -> s*x`, and its Hessian inertia is scale
independent.

## 3. Canonical derivatives

For local vertices `i,j,k,l`, put `e=v_j-v_i`, `a=v_k-v_i`,
`b=v_l-v_i`. Dot products are recovered solely from the six squared lengths.
The interior angle is frozen as

```text
cos(theta_ij)=
 [a.b-(a.e)(b.e)/(e.e)]
 / sqrt([a.a-(a.e)^2/(e.e)] [b.b-(b.e)^2/(e.e)]).
```

It must return `cos(theta)=1/3` at the regular tetrahedron. Volume is
`sqrt(det Gram(v1-v0,v2-v0,v3-v0))/6`.

Second-order automatic differentiation must act on these formulas before
global assembly. The global Hessian must include every term from both
`F V^(-1/3)` factors, including the `C'' d beta tensor d beta` term and the
second derivatives of the dihedral angles. Finite differences alone are
forbidden as the primary Hessian.

## 4. Stationarity, quotient and preregistered probes

The full `H4` symmetry is transitive on edges, so the gradient of `Ahat` at
the equilateral point is proportional to the all-ones vector. Exact scale
invariance makes its contraction with that vector zero. Hence stationarity
is predicted structurally and must also be recovered computationally.

Only the global scale vector is removed. No vertex-displacement or continuum
diffeomorphism direction may be discarded without an independent exact null
identity.

In addition to the complete 719-dimensional quotient spectrum, freeze the
following target-independent directions before calculation:

1. every edge-stabilizer orbit contrast `e_0-e_j`, one representative for
   every orbit of edges relative to the lexicographically first edge;
2. the discrete conformal direction

   ```text
   h_(ij)=f(v_i)+f(v_j),
   f(v)=v_1^2-v_2^2,
   ```

   with its mean removed if exact symmetry does not already make it
   scale-free;
3. the full vertex-quadratic conformal span obtained from all trace-free
   quadratic forms on `R4`, reported as a subspace rather than selecting its
   best vector after seeing the result.

The ambient coordinates label canonical functions on the original
equilateral orbit; they do not constrain the deformed Regge metric to remain
embedded.

## 5. Mandatory hostile controls

The verifier must check:

- f-vector `(120,720,1200,600)` and five tetrahedra around every edge;
- every regular local Gram matrix is positive and every dihedral cosine is
  `1/3`;
- the two independent local angle/volume implementations agree at the
  regular point;
- local first and second derivatives respect all tetrahedral relabellings;
- assembled gradient is zero by the symmetry-plus-homogeneity identity;
- assembled Hessian is symmetric and annihilates the scale vector;
- the quadratic form computed from the dense Hessian agrees with direct
  second-order evaluation on all frozen probe families;
- inertia from symmetric eigendecomposition agrees with an independent
  congruence/LDL count;
- the nonzero sign count is stable when all analytic constants are evaluated
  at 50, 80 and 120 decimal digits.

Numerical inertia may be labelled **DERIVED COMPUTATIONAL** only if the
smallest claimed nonzero magnitude exceeds `10^-7`, cross-precision drift is
below `10^-10`, the scale residual is below `10^-10`, and the independent
counts agree. Otherwise the full inertia is **PATTERN/OPEN**. A strict sign on
an explicit frozen direction may still be reported separately if its error
bound excludes zero.

## 6. Preregistered decision boundary

- **DERIVED FINITE LOCAL MINIMUM:** stationarity holds and the quotient
  Hessian has 719 positive, zero negative and zero additional null modes.
- **DERIVED FINITE LOCAL MAXIMUM:** it has 719 negative directions and no
  additional null modes.
- **DERIVED FINITE SADDLE:** at least one strictly positive and one strictly
  negative non-scale direction are certified.
- **DERIVED EXTRA DEGENERACY:** additional exact null vectors are found and
  identified.
- **OPEN/FAILED CERTIFICATE:** stationarity, derivative validation or robust
  sign separation fails.

Whatever the outcome, it concerns the equilateral **singular Regge** point.
Agreement with the smooth conformal sign is structural evidence for a stable
discretization; disagreement is not a contradiction because the cone
correction is finite and the analytic categories differ.

## 7. Physical boundary

Even a strict finite minimum would concern one asymptotic heat coefficient,
not the complete cutoff action. It would not derive Lorentzian dynamics,
Newton's constant, a Planck scale or a propagating graviton. Conversely a
saddle would close the claim that the positive `A2` term alone supplies a
stable finite gravitational vacuum.

Only the eventual targeted verifier and a non-executing static registry audit
will be run. No full-suite run and no PDF build are authorized.
