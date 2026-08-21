# Precision adjudication first result: direct action wins, orbit action disagrees

Date: 2026-08-21

## Provenance

- precision protocol: `093055e`;
- pre-evaluation correction: `287dae9`;
- registered implementation before evaluation: `ab7e491`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_tick_scale_covariance_precision.py`;
- artifact:
  `reproducible/gravity_600cell_tick_scale_covariance_precision.json`;
- artifact SHA-256:
  `7579656a59be5ae8d808bf560b9d231d264b0d14610fb220120433ee639c76f0`.

Only the targeted verifier was run.  It returned **7/10** and exited nonzero.
The full suite and deferred censuses were not run.

## Mechanical outcome

```text
TICK_SCALE_COVARIANCE_HIGH_PRECISION_DISAGREES
```

## Result that survives

The locally reimplemented 80-decimal action loops over all 2400 simplices and
all individual hinges, without orbit multiplicities.  For both parities it
obeys the scale identity at both preregistered nontrivial factors:

| parity | `alpha` | relative scale error |
|---|---:|---:|
| even | `0.6` | `6.67e-75` |
| even | `1.75` | `2.90e-73` |
| odd | `0.6` | `2.52e-74` |
| odd | `1.75` | `5.08e-73` |

All six states contain exactly 2400 Lorentzian simplices on the certified
branch.  Thus direct arbitrary precision independently confirms the action
homogeneity.  This does not rescue the preregistered adjudication outcome,
because its cross-implementation prediction failed.

## Unexpected discriminator

The high-precision direct base action agrees with the binary64 direct action,
not with the 100-decimal orbit-reduced action:

| parity | high precision vs orbit | high precision vs binary64 direct |
|---|---:|---:|
| even | `1.64545e-6` | `2.76622e-9` |
| odd | `1.33937e-6` | `1.80430e-9` |

The protocol predicted that arbitrary precision would select the orbit action.
It selected the opposite implementation.  Therefore the earlier discrepancy
is **not** explained by binary64 precision loss.

## Honest interpretation

- **DERIVED EXACT / independently reproduced:** the direct complete Regge-dust
  action has length degree two under simultaneous geometry and mass scaling.
- **DERIVED COMPUTATIONAL NEGATIVE:** the frozen orbit and direct action
  implementations do not agree on the order-24-invariant off-shell perturbation.
- **OPEN:** whether the fault is in the orbit incidence reduction, in the
  independent direct construction, or in the assumption that the coordinate-
  indexed perturbation is invariant under the action used by both evaluators.
- **NOT ACCEPTED:** no consolidated absolute-tick theorem is issued while the
  complete-action implementations disagree.

The leading hypothesis is an orbit-incidence multiplicity error.  The orbit
action accumulates angles from one representative simplex per simplex orbit;
that shortcut is valid only if representative incidence realizes every
triangle--simplex double incidence with the correct multiplicity.  This has not
yet been proved and must not be assumed from regular-state agreement.

## Next falsification test

Preregister an exact incidence audit before changing any evaluator:

1. enumerate the triangle and simplex orbits under the actual stabilizer;
2. form the complete integer incidence matrix counting, for one triangle in
   each triangle orbit, incident simplices in every simplex orbit;
3. compare it with the multiplicities implicitly used by the orbit evaluator;
4. reconstruct the reduced action using the exact incidence matrix;
5. require agreement with the direct action on a known regular control and the
   frozen off-shell discriminator.

If the matrices differ, this is a real bug in the nonhomogeneous orbit action
and every downstream result using that shortcut must be scope-audited.  If they
agree, the direct arbitrary-precision implementation must be attacked instead.

