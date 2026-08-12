# Protocol: full-de Rham curvature coefficient at the fixed Regge endpoint

Date: 2026-08-12

## Provenance and hostile scope

This is a **post-recognition** protocol.  Before it was written, the smooth
homogeneous calculation had already shown that the ordinary full-de Rham
coefficient called `A2` in this repository uniquely minimizes at the round
metric at fixed volume.  It was also already known that the fixed Regge
endpoint is locally a product of an edge with a two-dimensional cone, and a
small-deficit estimate suggested that the volume-normalized Regge value might
be close to, or even pass, the round value.  Therefore this is not a blind
endpoint comparison and no evidential claim will be made from numerical
surprise alone.

The purpose is narrower and falsifiable: determine whether the *same* ordinary
full-exterior Hodge heat coefficient is mathematically defined at both
endpoints, and, if it is, compare the two exact values without replacing the
conical coefficient by the smooth Regge-curvature limit.

No heat time, cutoff function, Newton constant, phenomenological target or
fitted normalization may enter the test.

## Frozen operator and coefficient

1. The smooth endpoint is the unit round metric `g_0` on `S3`.

2. The singular endpoint is the boundary of the regular Euclidean 600-cell,
   with unit circumradius and its exact facetwise-flat metric `g_R`, transported
   to `S3` by the already certified radial homeomorphism.

3. The operator is the ordinary Hodge--de Rham Laplacian on the complete
   exterior algebra,

   ```text
   Delta = (d+d*)^2,       p=0,1,2,3,
   ```

   with the closed-Hilbert-complex/Friedrichs domain selected by the existing
   fixed-Regge Whitney limit.  A different point-interaction or ideal-boundary
   extension is outside this protocol.

4. The notation `A2` means the coefficient with length dimension one in

   ```text
   Tr exp(-t Delta)
     ~ (4 pi t)^(-3/2) [A0 + t A2 + ...].
   ```

   It is the curvature coefficient, not the finite moment historically also
   named `A2=Tr(D^4)/2` elsewhere in the repository.

5. The trace is ordinary, summed over all form degrees.  The graded
   supertrace is excluded because its nonconstant terms cancel by the index
   theorem.

## Frozen singular formula to be independently checked

Let a codimension-two stratum have cone angle `beta`, length `L`, and
`gamma=2*pi/beta`.  For the scalar Friedrichs Laplacian, Fursaev and Miele give

```text
S(beta,L) = beta/6 * (gamma^2-1) * L.
```

For the Hodge Laplacian on one-forms in ambient dimension three they give

```text
V(beta,L) = 3 S(beta,L) + 2 (beta-2*pi) L.
```

Hodge duality identifies degrees `3` with `0` and `2` with `1`.  Hence the
candidate exact full-exterior contribution of one open edge is frozen as

```text
C_full(beta,L)
  = 2 S(beta,L) + 2 V(beta,L)
  = [8 * beta/6 * ((2*pi/beta)^2-1)
     + 4*(beta-2*pi)] * L.
```

Before using this expression, the result stage must verify all of the
following:

- the cited formulas are for the Hodge--de Rham operator, not a rough vector
  Laplacian;
- their domain agrees with the regular/Friedrichs choice licensed above;
- the actual 600-cell angle lies in the directly treated range, rather than
  being reached only by an uncontrolled analytic continuation;
- vertex links cannot contribute to the same `t^(-1/2)` order.  Cheeger's
  piecewise-flat heat expansion must place that coefficient on the
  one-skeleton; otherwise the comparison stops as incomplete.

The small-deficit expansion must be checked as a convention control:

```text
C_full(2*pi-delta,L)
  = -(4/3) delta L + O(delta^2),
```

which agrees with the smooth identity

```text
A2 = -(2/3) integral R,
integral R = 2 sum_edges delta_e L_e.
```

Agreement only to first order is not permission to discard the exact conical
correction.

## Frozen 600-cell geometry and normalization

The verifier must derive rather than assume:

```text
N_tetrahedra = 600,
N_edges       = 720,
edge length   = 1/phi,
beta          = 5 arccos(1/3),
V_R           = 600 * edge_length^3/(6 sqrt(2)),
V_0           = 2 pi^2.
```

There are five regular tetrahedra at every edge.  The cone angle must satisfy
`0 < beta < 2*pi`.  Because `A2` scales as length in three dimensions, the
fixed-volume Regge value is

```text
A2_R_equal_volume = (V_0/V_R)^(1/3) * 720 * C_full(beta,edge_length).
```

The smooth round value is frozen independently as

```text
A2_round = -(2/3) * 6 * V_0 = -8*pi^2.
```

No volume, radius or coefficient may be rescaled after inspecting the
comparison.

## Decision boundary

- **ROUND BEATS FIXED REGGE:** `A2_round < A2_R_equal_volume` under the
  already used minimization convention.
- **FIXED REGGE BEATS ROUND:** `A2_R_equal_volume < A2_round`.  This refutes
  any global round-selection reading of the smooth theorem; it does not by
  itself make the fixed Regge endpoint a dynamical vacuum.
- **TIE:** the exact expressions agree.
- **ANALYTIC GATE FAILS:** the domain, conical formula, or skeleton
  localization does not apply.  No endpoint number is then admissible.

Every numerical sign must be backed by an interval or higher-precision
certificate stable under increased precision.  The result will be labelled
**DERIVED CONDITIONAL** at most: one heat coefficient is not a complete
spectral action, and it still supplies neither a cutoff function nor an
absolute scale.

## Primary references frozen before comparison

- D. V. Fursaev and G. Miele, *Cones, Spins and Heat Kernels*,
  <https://arxiv.org/abs/hep-th/9605153>, especially equations (i6a), (i10)
  and the two-dimensional Hodge decomposition.
- J. Cheeger, *Spectral Geometry of Singular Riemannian Spaces*,
  <https://webhomes.maths.ed.ac.uk/~v1ranick/papers/cheeger.pdf>, especially
  Sections 7.1 and 7.5 on piecewise-flat pseudomanifolds and skeleton-local
  heat coefficients.

