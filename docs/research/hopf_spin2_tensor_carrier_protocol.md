# Protocol: Hopf projector fields as a symmetric spin-2 carrier

Date: 2026-08-12

## Provenance

This is a **post-recognition structural audit**, not a blind discovery.  The
repository already certifies the six unoriented fivefold-axis projectors,
their two handed quaternionic lifts, and

```text
span{P_i-I/3} = Sym^2_0(R^3).
```

Before this protocol, those facts had been interpreted as a five-component
order parameter, not tested as transverse-traceless tensor fields on the
round three-sphere.  The expected TT result and connection-Laplacian
eigenvalue were recognized before implementation.  Commit ordering freezes
the formulas and falsifiers but supplies no target-blind discovery claim.

No particle mass, Newton constant, Planck scale, measured speed, or
phenomenological target is used.

## Complete hypotheses

1. The carrier is the **unit round** `S^3`, identified with unit quaternions
   and sectional curvature `+1`.  This is the smooth Hopf carrier, not the
   separately certified fixed piecewise-flat Regge metric.
2. At the identity, `Im(H)=R^3` has its Euclidean metric and orientation.
   For `q in S^3`, use the two already-derived handed tangent frames

   ```text
   L_q(v)=qv,       R_q(v)=vq,       v in Im(H).
   ```

3. Use exactly the six already-certified unoriented fivefold-axis
   projectors `P_i`; no alternative axes, weights, signs, or tensor basis may
   be searched.  Put

   ```text
   T_i=P_i-I/3 in Sym^2_0(R^3).
   ```

4. Lift them as covariant symmetric tangent tensors

   ```text
   T_i^L(q)=L_q T_i L_q^T,
   T_i^R(q)=R_q T_i R_q^T.
   ```

5. In a left-invariant orthonormal frame on the unit sphere,

   ```text
   [e_a,e_b]=2 epsilon_abc e_c,
   nabla_(e_a)e_b=epsilon_abc e_c.
   ```

   If `Gamma_a` is the corresponding skew connection matrix and `H` has
   constant components, freeze

   ```text
   nabla_a H = [Gamma_a,H],
   div(H)_b = sum_a (nabla_a H)_(ab),
   nabla* nabla H = -sum_a [Gamma_a,[Gamma_a,H]].
   ```

   The right-invariant frame has the opposite connection sign; divergence
   and the squared operator are unchanged.

These hypotheses do not posit a Lorentzian time, an Einstein action, a metric
fluctuation, a diffeomorphism quotient, or a stress-energy source.

## Frozen questions

The registered verifier must determine exactly:

1. whether left and right quaternion multiplication give orthonormal tangent
   frames for every unit `q`;
2. whether the six `T_i` have rank-five span and the already-certified tight
   frame `(4/5)I_5`;
3. whether every symmetric tracefree constant-frame tensor has zero
   divergence;
4. whether the connection Laplacian acts by one exact scalar on this
   five-dimensional space and, if so, what that scalar is;
5. whether the two handed constant-coefficient spaces intersect trivially as
   global fields, using the exact `A5` adjoint action rather than numerical
   sampling;
6. whether the geometric exterior-algebra fibre of the existing
   Kaehler--Dirac continuum contains an intrinsic spin-two component.  Under
   proper tangent rotations this must be checked from

   ```text
   Lambda^*(R^3)=Lambda^0 + Lambda^1 + Lambda^2 + Lambda^3,
   ```

   while distinguishing intrinsic tensor type from orbital `l=2` harmonics;
7. whether a bilinear/symmetric-square construction can contain spin two,
   and whether the repository currently supplies a selected nonlinear map or
   action that realizes it.

## Decision rule

- **Kinematic advance:** if the centered Hopf projector fields are TT,
  furnish the full five-dimensional symmetric-tracefree fibre, and the two
  handed homogeneous spaces are distinct exact eigenspaces, record a
  canonical round-`S^3` spin-two **carrier**.
- **Kill:** if tracefreeness, transversality, span, handedness, or the common
  eigenvalue fails, the proposed Hopf spin-two carrier is closed.

Even a positive result is not a graviton derivation.  It supplies the correct
linear tensor type but not:

- a selected variable metric or gravitational action;
- Lorentzian propagation or a massless constraint algebra;
- diffeomorphism gauge symmetry and quotient;
- universal coupling to a conserved stress tensor;
- compatibility with the fixed-Regge rather than round continuum;
- Newton's constant or a Planck scale.

## Labels fixed in advance

- Quaternion-frame, tight-frame, divergence, representation and Casimir
  identities proved exactly: **DERIVED** under the round-`S^3` hypotheses.
- Reading the two handed five-spaces as a candidate spin-two kinematic seed:
  **STRUCTURAL ADVANCE**.
- Calling the seed a physical graviton or emergent gravity: **OPEN**.
- Importing the Einstein--Hilbert/Regge/Lichnerowicz action merely because it
  is the desired gravitational answer: **FITTING/forbidden as selection
  evidence**.
