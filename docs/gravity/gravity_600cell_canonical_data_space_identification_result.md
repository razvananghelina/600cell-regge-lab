# Identification of the modular canonical-data space

Date: 2026-08-19

## Complete hypotheses

Use the frozen complete variable-face flat-frustum compatibility equations
`F f + E e + S s = 0` on the fixed lower regular 600-cell.  Here `f` has
3600 cell-flex coordinates, `e` has 720 upper spatial squared-edge
coordinates, and `s` has 120 strut squared-length coordinates.  The proposed
spatial data are not fitted: for a vertex scalar `sigma`, the edge datum on
`{u,v}` is the derived infinitesimal squared-length variation

```text
e_{uv} = 8 lambda (sigma_u + sigma_v).
```

The test uses both rational representatives `(lambda,tau)=(2,5),(3,11)`,
both frozen primes 1000003 and 1000033, two exact local right-inverse graphs,
reversed face orientation, an odd canonical relabelling, and reversed metric
sign.  It is a linear kinematic compatibility test.  It is not an action,
Hessian, symplectic, constraint, or evolution calculation.

## Frozen provenance and result

The prior-art gate is commit `7fd131d`, the target-disclosed protocol is
commit `8576c84`, and the registered implementation was committed as
`c198268` before its first execution.  The first artifact was frozen in
commit `e9a941b`; its SHA-256 is
`3db0b9ce8c90cba9de3fbbff818129388d79a98e0483a0ca3ae53b2e4d271434`.

The targeted verifier passed 11/11 checks.  Every one of the fourteen
construction/prime combinations gives

```text
rank(F)                 = 3600
rank([F E])             = 4200
rank([F S])             = 3600
rank([F E U])           = 3600
rank([F E U_bad])       = 3601
```

The unsigned incidence map `U` and the deliberately corrupted `U_bad` have
exact rational and modular ranks

```text
rank(U) = 120, rank(U_bad) = 120, rank([U U_bad]) = 121.
```

Thus the negative map is genuinely a different 120-dimensional image, and
the equations distinguish it by one additional rank.

## What follows, and what does not

For either frozen finite field, the preceding target-blind census proved that
the compatible edge-only space has dimension 120 and the entire strut ambient
space is compatible.  The present inclusion, together with `rank(U)=120`,
therefore proves

```text
compatible data = im(U) direct-sum arbitrary strut data
```

over both fields.

**DERIVED (modular).** The 120 spatial compatible data are exactly the
unsigned vertex-scale edge variations, and the other 120 compatible data are
all strut variations.  This is an image equality, not a dimension match: the
one-row-corrupted image is rejected.

**DERIVED NEGATIVE, retained.** The previously tested *local* cell-flex lift
of these same data fails all 3600 face equations.  The present result
identifies the boundary-data image only; it does not resurrect that lift.

**STRUCTURAL.** The result is stable under every frozen graph and convention
attack, but it is still modular and it inherits the recorded blindness
deviation of the preceding projection census.  A material rational carrier
claim requires a mechanically independent exact implementation.

**OPEN.** Equality over `Q`, an explicit globally solved cell-flex lift, the
canonical symplectic phase space, action dynamics, a selected lapse/tick,
tensor propagation, `c`, `G`, and Planck units are not established.

## Physical interpretation

The result supplies a natural coordinate description of the infinitesimal
boundary-data domain:

```text
120 vertex conformal/scale coordinates + 120 arbitrary strut coordinates.
```

It is useful because the next action/Hessian calculation no longer has to be
performed on an unidentified 840-dimensional boundary space.  It is not yet
physical evolution.  In particular, arbitrary compatible struts are closer
to flat-background lapse freedom than to a dynamically selected clock.

## Next falsifiable calculation

Construct the rational lifts without reusing the modular elimination path:

1. solve the exact face equations for the unique 3600 cell-flex response to
   the 120 vertex-scale and 120 strut basis data;
2. verify the full rational residual directly, including a corrupted-image
   attack;
3. expose the lift as a sparse/local formula if one exists, or record its
   unavoidable global support;
4. only after that, restrict the Regge boundary Hessian to the accepted data
   coordinates and separate gauge/null lapse directions from dynamical
   directions.

