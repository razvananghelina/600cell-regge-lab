# Protocol: the ordinary de Rham coefficient on the full round--Regge path

Date: 2026-08-12

## Provenance

This is a post-recognition hostile protocol.  The following facts were known
before it was written:

- on the smooth fixed-volume homogeneous branch, ordinary full-de Rham `A2`
  has the round metric as its unique global minimum;
- after exact conical correction and equal-volume normalization, the round
  endpoint is lower than the fixed 600-cell Regge endpoint by only
  `0.0848366160...`;
- replacing the exact cone term by the linear Regge deficit gives the wrong
  endpoint ordering;
- for a generic interior parameter `0<u<1`, a representative face has a
  nonzero second fundamental form with a factor `u(1-u)`.  Thus a cone-only
  interpolation is already known to be invalid.

No interior value of the complete coefficient has been evaluated.  This is
not a blind discovery protocol, and a small numerical gap will not be called
evidence without a whole-interval certificate.

## Complete hypotheses and operator

1. The carrier is the radially identified boundary of the unit-circumradius
   regular 600-cell.

2. On each open tetrahedral facet, in its Euclidean tangent coordinates, the
   metric family is exactly

   ```text
   g_u=(1-u) I + u (I/r^2 - y y^T/r^4),   0<=u<=1,
   r^2=a^2+|y|^2,
   a^2=(7+3*sqrt(5))/16.
   ```

   It is the already registered affine family, not a geodesic in a newly
   selected metric-space norm.

3. The operator is the self-adjoint Hodge--de Rham Laplacian of the closed
   Hilbert complex on the complete exterior algebra.  Across every open face
   it uses the de Rham transmittal domain

   ```text
   phi_+=phi_-,
   B_0((d+d*)phi)=0,
   ```

   derived by Gilkey--Kirsten--Vassilevich.  Independent Robin, delta-shell or
   point-interaction data are forbidden.

4. The coefficient is defined by

   ```text
   Tr exp(-t Delta_u)
     ~ (4*pi*t)^(-3/2) [A0(u)+t*A2(u)+...].
   ```

   The trace is ordinary and summed over degrees `0,1,2,3`.

5. Every metric is rescaled only after assembling the unnormalized
   coefficient, to the frozen volume `V0=2*pi^2`.  Since `A2` has length
   dimension one,

   ```text
   A2_equal_volume(u)=(V0/V(u))^(1/3)*A2_raw(u).
   ```

## Why the face term is mandatory

For a continuous leading metric with a jump in its normal derivative, the
primary transmittal theorem gives the bracket coefficient

```text
A2_face=(1/6) integral_face Tr[2(L+ + L-)I - 6U].
```

For the full exterior algebra in dimension three, `Tr I=8`.  The de Rham
transmittal endomorphism has

```text
Tr U=4 tr(L+ + L-).
```

Therefore the frozen ordinary all-form face contribution is

```text
A2_face=-(4/3) integral_face tr(L+ + L-).
```

This term vanishes at `u=0` and `u=1`, but not generically between them.  A
calculation omitting it fails the protocol.

Primary source: P. B. Gilkey, K. Kirsten and D. V. Vassilevich, *Heat trace
asymptotics with transmittal boundary conditions and quantum brane-world
scenario*, <https://arxiv.org/abs/hep-th/0101105>, Theorem 2.3 and Lemma 6.1.

## Frozen local decomposition

The raw coefficient must be assembled from all strata contributing at this
order:

```text
A2_raw(u) = B(u)+F(u)+E(u),

B(u)=-(2/3)*600 * integral_tetra R(g_u) dV_u,

F(u)=-(4/3)*1200
     * integral_triangle tr(L+ + L-) dA_u,

E(u)=720 * integral_edge
     [16*pi^2/(3*beta_u)+8*beta_u/3-8*pi] dl_u.
```

Here `beta_u=5*theta_u` is the angle of the tangent cone at the edge point.
The five sectors are forced by the exact `C5` edge link.  The edge formula is
the already verified exact full-de Rham cone coefficient.  Terms containing
extrinsic curvature have the wrong length dimension to enter this edge
coefficient; they first occur at the next heat order.

The following endpoint controls are compulsory:

```text
u=0: B=F=0 and E is the fixed-Regge cone result;
u=1: F=E=0, V=2*pi^2 and B=-8*pi^2.
```

If either fails under refinement, no interior conclusion is admissible.

## Frozen coordinates and quadrature

Let `rho=sqrt(1-a^2)` and use the centered regular tetrahedron

```text
(rho/sqrt(3))*(+,+,+),
(rho/sqrt(3))*(+,-,-),
(rho/sqrt(3))*(-,+,-),
(rho/sqrt(3))*(-,-,+).
```

The bulk metric is radial in these coordinates.  Write

```text
q=1-u+u/r^2,
p=1-u+u*a^2/r^4,
g_u=p ds^2+q s^2 dOmega^2,
f=s*sqrt(q),  dxi=sqrt(p) ds,
R=-4 f_xixi/f + 2(1-f_xi^2)/f^2.
```

The face term is evaluated on the triangle opposite the first vertex, with
the inward normal fixed toward that vertex.  The two incident contributions
are related by the exact face reflection and are included separately through
`L+ + L-`; their sign may not be chosen after inspecting the result.

For the edge between the first two vertices, parameterize

```text
y=(rho/sqrt(3))*(1,t,t),  -1<=t<=1.
```

The two inward face covectors are those pointing toward the two opposite
vertices.  Freeze

```text
cos(theta_u)
 = -(n1^T g_u^{-1} n2)
   /sqrt((n1^T g_u^{-1}n1)(n2^T g_u^{-1}n2)),
beta_u=5*acos(cos(theta_u)).
```

All integrations use tensor Gauss--Legendre quadrature after the standard
Duffy map from `[0,1]^3` to the tetrahedron and its two-dimensional analogue
to the face.  Orders `16,24,32,40,48` are frozen before evaluation.  Edge
orders are twice these values.  The reported digit interval must contain all
of the last three orders.  Polynomial volume/Jacobian controls and both exact
endpoints must pass first.

## Whole-interval decision boundary

The first attack is the preregistered rational grid

```text
u=j/200,  j=0,...,200.
```

This grid can **refute** global round selection if a converged value lies
strictly below `A2_equal_volume(1)`.  It cannot prove the positive claim.

- **REFUTED:** one certified interior interval lies below the round value.
- **PATTERN ONLY:** the frozen grid prefers the round endpoint but no
  derivative/interval certificate covers all `u`.
- **DERIVED PATH SELECTION:** in addition to the converged grid, an interval
  enclosure or analytic inequality proves
  `A2_equal_volume(u)>A2_equal_volume(1)` for every `0<=u<1`.
- **ILL-POSED/INCOMPLETE:** the transmittal, edge-corner or convergence gates
  fail.  No preferred `u` is then reported.

Even the strongest outcome is conditional on the sign and use of this one
asymptotic coefficient.  It does not select a cutoff, scale, Lorentzian time,
Newton constant or source coupling.

