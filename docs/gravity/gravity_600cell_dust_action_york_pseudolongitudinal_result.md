# Consolidated result: negative shape modes are pseudo-longitudinal, not exact gauge

Date: 2026-08-19

## Verdict

The proposed exact identity

```text
negative generalized shape carrier
= action-weighted image of tangential vertex displacements
```

is REFUTED at the repository's preregistered high-precision computational
standard on the fixed curved two-slab background.

The correct label is:

```text
DERIVED COMPUTATIONAL / STRUCTURAL:
the thirty negative shape modes are pseudo-longitudinal, not exact gauge.
```

This is not a symbolic or formal interval theorem and makes no claim about
continuum diffeomorphism symmetry.

## Provenance ledger

| stage | commit | outcome |
|---|---|---|
| prior-art gate | `1b076f2` | target and literature disclosed |
| direct protocol | `852cd54`, corrected before execution in `167077b` | method and `10/100` hierarchy frozen |
| registered direct verifier | `5e3c331` | no target result yet |
| direct artifact | `5a450f1` | `17/17`, numerical refutation |
| first adversarial protocol | `c65ed2a` | independent archive/QR route frozen |
| adversarial verifier | `5705d67` | mechanically independent construction |
| disclosed harness repairs | `af9fceb`, `5b38980`; code `ab963ac`, `48ba565` | no residual existed before either repair |
| contradictory adversarial artifact | `33e0fdc` | literal preregistered outcome preserved |
| framing-failure analysis | `0cd3128` | common absolute error shown scale dependent |
| corrected relative protocol | `199696f` | dimensionless classifier frozen |
| corrected verifier | `490b38f` | no corrected result yet |
| corrected artifact | `47172df` | `9/9`, independent confirmation |

Artifact hashes:

```text
direct:
d57351e852ab40eb7809397c84e5f57ff58e5ae0bd31f9dcaf87efdc84be76b5

first adversarial:
e39203741513f128a208f22896abef53daa12db089ee7e43abf9c90643fc579b

scale-invariant adversarial:
aee9088c1b0bb1cd4ebec12707014ed187f5c6cbe7ad84cc6cd74db705b8d20c
```

## Primary direct calculation

The direct verifier rebuilt both physical slabs from all local 4-simplex
Hessians at 120 decimal digits with four newly frozen derivative scales.  It
used exact golden-ratio 600-cell coordinates and no centered midpoint or
radius until the final formula control.

Controls:

```text
direct M,V versus frozen centered archive:  maximum distance 2.946e-11
carrier dimensions:                         5+25, then 15+10
selected cells:                             2 schedules x 2 sectors x 4 variants = 16
```

Every cell returned:

```text
||L* A T||_2                         NONZERO_RESOLVED
||A L-B L(L*BL)^-1(L*AL)||_2        NONZERO_RESOLVED
||P_L-P_negative||_2                 NONZERO_RESOLVED
longitudinal stiffness               NEGATIVE_RESOLVED
transverse stiffness                 POSITIVE_RESOLVED
rotated control                      NONZERO_RESOLVED
```

Representative midpoint/error values are

```text
cross residual     2.96331536e-5 / 9.79842267e-9
image residual     2.97372448e-5 / 9.79842267e-9
projector distance 4.07209597e-4
```

Thus the non-invariance margin is about three thousand empirical error units
before applying the frozen factor 100.

## Independent adversarial calculation

The independent route did not rebuild the local Hessians.  It loaded the
older centered archive, rebuilt normalized binary rigidity, used
column-pivoted QR instead of the primary SVD and tested two different
observables:

```text
distance of A L from im(B L),
leakage of B^-1 A L from L.
```

The first adversarial protocol incorrectly applied one dimensionful absolute
error after `B^-1` had rescaled the second observable.  It therefore returned
the preserved but internally contradictory labels

```text
span residual       NONZERO_RESOLVED 16/16
inverse residual    ZERO_CONSISTENT  16/16
augmented rank      24               16/16.
```

That protocol's literal `ADVERSARIAL_DIRECT_REFUTATION_REFUTED` outcome is a
first-class negative about the classifier, not evidence that the exact
identity holds.  Exact arithmetic makes the two invariance conditions
equivalent.

The preregistered correction normalized each residual by its own operator
norm.  It passed both conditioning inequalities in all sixteen cells and
returned

```text
relative span residual         0.0777626758  NONZERO_RESOLVED 16/16
relative B^-1 A leakage        0.0090188380  NONZERO_RESOLVED 16/16
empirical relative error       about 1.255e-9
relative augmented rank        19 (>15)      16/16
outcome                         SCALE_INVARIANT_DIRECT_REFUTATION_CONFIRMED
```

The exact number assigned to the extra augmented rank depends on the disclosed
relative threshold (`24` under the earlier absolute threshold and `19` under
the corrected one).  Only the threshold-robust statement `rank > 15` is
claimed.

## Physical meaning

DERIVED COMPUTATIONAL / STRUCTURAL:

- the negative carrier is very close in subspace angle to the tangential
  vertex-displacement image (`||P_L-P_-|| about 4.07e-4`);
- it is nevertheless not invariant under the generalized stiffness pencil;
- the negative directions therefore cannot be discarded as an exact gauge
  quotient on this curved finite background.

PATTERN:

- the dimension match `15+10` and the small projector angle identify the
  modes as pseudo-longitudinal in the sense expected when discrete curvature
  breaks vertex-displacement symmetry;
- the relative leakage is not microscopically small (`7.78%` in the span
  normalization and `0.902%` after `B^-1`), so calling the equality "almost
  exact" without stating the norm would be misleading.

OPEN:

- whether these directions are a genuine physical instability, a
  pseudo-constraint artifact or a coarse-lattice effect;
- whether the splitting vanishes in a flat-background or refinement limit;
- the physical tensor quotient, dispersion, propagation speed and continuum
  polarization content.

The result agrees qualitatively with the established distinction between
exact flat-background vertex-displacement symmetry and curved Regge
pseudo-constraints, but that literature does not prove this 600-cell result.

## Next decisive test

The next useful quantity is not another equality test at the same resolution.
It is the scaling of the dimensionless leakage under a controlled curvature or
refinement parameter:

```text
rho_span(curvature, resolution),
rho_comm(curvature, resolution).
```

If both tend to zero toward a flat/refined limit with a preregistered rate,
the pseudo-longitudinal interpretation advances.  If they remain finite, the
negative carrier is not a recovering gauge direction and must be treated as a
physical/coarse instability candidate.  No quotient is authorized before
that test.
