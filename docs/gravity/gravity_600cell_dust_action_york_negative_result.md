# Result: the coarse negative shape carrier is only longitudinal-consistent

Date: 2026-08-19

## Status ledger

| Item | Status | Evidence |
|---|---|---|
| Prior-art gate | DERIVED DOCUMENTARY | commit `9bd990a` |
| Target-disclosed protocol | DERIVED DOCUMENTARY | commit `7165752` |
| Registered verifier before execution | DERIVED COMPUTATIONAL | commit `33f0bf7` |
| First-run harness failure disclosed before repair | DERIVED COMPUTATIONAL | commit `52eb83e` |
| Frozen rank-harness correction | DERIVED COMPUTATIONAL | commit `3bfab80` |
| First complete artifact | DERIVED COMPUTATIONAL | commit `4d7b65e`, SHA-256 `fd0763af779cb02d96f7e1d7a8856b117dd4bf2c9413f01de6246c597743df27` |
| Exact action-weighted longitudinal identity | OPEN | `13/15`; outcome `NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_OPEN` |
| Longitudinal interpretation of the negative carrier | PATTERN | stable in all 16 disclosed cells, but not separated from the propagated enclosure |
| Exact gauge quotient | OPEN / NOT AUTHORIZED | curved finite Regge backgrounds need not have exact vertex-displacement gauge symmetry |

The first execution was stopped at the preregistered global geometry control.
It had omitted the unit normalization used by the frozen rigidity theorem and
had used `numpy.linalg.matrix_rank` instead of the frozen `10/100` singular
bands.  No target verdict or JSON was produced.  The correction changed no
action matrix, carrier definition, selected sector, target threshold or
outcome hierarchy.

## Complete targeted result

The corrected run recovered the frozen carrier controls in both schedules:

```text
rank C = 120
rank R = 470
rank D = 354
dim(im C intersection im D) = 4
open singular values = 0
```

All 56 schedule/sector/derivative cells gave the same action-weighted
longitudinal/transverse dimensions:

| sector | irrep dimension | longitudinal | transverse |
|---:|---:|---:|---:|
| 0 | 3 | 44 | 31 |
| 1 | 2 | 30 | 20 |
| 2 | 2 | 30 | 20 |
| 3 | 2 | 28 | 22 |
| 4 | 1 | 15 | 10 |
| 5 | 1 | 15 | 10 |
| 6 | 1 | 12 | 13 |

Restoring representation multiplicities gives exactly

```text
dim L_H = 350,       dim T_H = 250.
```

This census is DERIVED COMPUTATIONAL.  In particular, the two sectors chosen
by the preceding blind negative-stiffness result acquire dimensions `15+10`
from geometry alone; the match is not a dimension fitted after loading the
negative eigenspaces.

For all 16 selected cells the preregistered classifiers returned:

```text
negative-projector comparison: EQUALITY_CONSISTENT  16/16
A longitudinal/transverse cross block: ZERO_CONSISTENT 16/16
B longitudinal/transverse cross block: ZERO_CONSISTENT 16/16
Gamma leakage: ZERO_CONSISTENT 16/16
Omega leakage: ZERO_CONSISTENT 16/16
transverse stiffness: POSITIVE_RESOLVED 16/16
longitudinal stiffness: OPEN 16/16
same-dimensional rotated control: EQUALITY_OPEN 16/16
```

The midpoint ranges show both why the result is interesting and why it is not
an identity theorem:

| diagnostic over 16 cells | midpoint range | propagated error |
|---|---:|---:|
| `||P_L-P_-||_2` | `4.07209628968e-4` to `4.07209629002e-4` | `3.65988107084e-2` to `3.65988107085e-2` |
| `||L* A T||_2` | `2.96331559511e-5` to `2.96331559537e-5` | `1.58805795544e-5` |
| `||L* B T||_2` | `6.65e-14` to `1.02e-13` | `4.47893333258e-3` |
| largest longitudinal `A` eigenvalue | `-1.95825711706e-4` | sign error `1.57374046319e-5` |
| smallest transverse `A` eigenvalue | `1.74039353911e-2` | same sign error |
| rotated-control projector distance | `0.999934` to `0.999981` | projector error `3.65988107085e-2` |

The longitudinal sign lies between the frozen `10*error` and `100*error`
bands.  The rotated control is genuinely far away at its midpoint, but the
preregistered generic separation rule asks for more than `100*error`; the
current propagated projector enclosure is too broad for that gate.  Neither
failure may be repaired by relabelling or by changing the bands after seeing
the answer.

## Interpretation

DERIVED COMPUTATIONAL: the geometry-selected tangential image has the right
all-sector dimensions, the correct `15+10` split in both disclosed sectors,
and is preserved by the centered first recurrence within the frozen
enclosures.

PATTERN: the negative carrier is extremely close to that longitudinal image
and the result is schedule/variant stable.  This is meaningful evidence that
the thirty coarse negative directions are longitudinal or pseudo-longitudinal
rather than an independently selected tensor carrier.

OPEN: exact equality.  The stable nonzero midpoint values
`||P_L-P_-|| about 4.07e-4` and `||L*AT|| about 2.96e-5` could be a small
curvature/discretization splitting rather than numerical noise.  The present
entrywise ball propagation is too wide to distinguish those cases.

This caution agrees with, but is not proved by, the established Regge
literature.  Exact vertex-displacement gauge symmetry is obtained for the
linearized theory on flat backgrounds, whereas curved discrete solutions
generically replace it by broken symmetries and pseudo-constraints:

- Bahr and Dittrich, [(Broken) Gauge Symmetries and Constraints in Regge
  Calculus](https://arxiv.org/abs/0905.1670).
- Dittrich and Hoehn, [From covariant to canonical formulations of discrete
  gravity](https://arxiv.org/abs/0912.1817).
- Hoehn, [Canonical linearized Regge Calculus: counting lattice gravitons with
  Pachner moves](https://arxiv.org/abs/1411.5672).

Therefore even a later exact carrier equality would be STRUCTURAL on this
fixed slab; it would not by itself prove an exact continuum gauge symmetry.
No quotient, graviton polarization, propagation speed or mass is derived here.

## Next falsifier

The next calculation should not rerun the same broad interval pipeline.  It
must reconstruct only the selected cells directly at high precision from the
local Regge Hessian and the exact golden-ratio 600-cell carrier, before any
binary/entrywise interval wrapping.  Its cheapest decisive observable is

```text
L* A T.
```

If its stable `2.96e-5` midpoint receives a certified error below
`2.96e-7`, the exact longitudinal identity is REFUTED under the frozen
`100*error` rule.  If it instead converges to zero and also resolves the
longitudinal sign and rotated control, the identity can advance.  Until that
independent calculation is preregistered and run, the honest verdict remains
OPEN rather than blocked.
