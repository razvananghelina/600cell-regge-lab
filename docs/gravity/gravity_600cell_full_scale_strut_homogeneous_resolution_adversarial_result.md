# Adversarial result: unique homogeneous weak-pole line replicated after control repair

Date: 2026-08-20  
Status: **DERIVED COMPUTATIONAL, adversarially replicated after disclosed control repair**

## Outcome and provenance

The fresh P200G run of source commit `7ccc619` reports

```text
HOMOGENEOUS_WEAK_POLE_LINE_REPLICATED_AFTER_CONTROL_REPAIR
7/7 tests passed
```

Artifact SHA-256:

```text
fab74a26ae940cf0e65f26a4f6f167285cc269e282c40d7a630f37d65ba7ab07
```

Verifier source SHA-256:

```text
0d6d50e2c757463eb86e341f8520d5a97d190949aa6a002eac57096a9d0f03d3
```

Only the targeted verifier and a static registry audit ran.  The full suite did
not run.

## Result

The mechanically different verifier re-executed the full-action builder,
rebuilt the complete homogeneous carrier and canonical-lift matrices at 200
decimal digits with Arb at 180 digits, and extracted D/K candidates by fixed
normal-equation solves.  It did not use the primary symbolic generator or the
stored P100 SVD candidate to construct the line.

In both parities:

- all ten D and all fifteen K single-column-deleted Gram determinants exclude
  zero: 50/50 fresh rank certificates;
- normalized D/K residuals lie between `5.57e-83` and `6.99e-83`;
- the independently extracted D and K projectors differ by about `1.20e-85`;
- even and odd D projectors differ by `6.23e-82`;
- the direct physical ratio agrees with the primary exact-action ratio to
  `1.50e-76`;
- the full-action imaginary residue falls to the `3.03e-159`--`9.27e-159`
  range.

Together with the primary exact generator and frozen rank lower bounds, this
establishes exactly one homogeneous weak-pole canonical intersection line.  The
already replicated nonhomogeneous intersection is zero in every sector.

## Hostile controls and disclosed repair

The first P160 adversarial run remains a formal `CONTROL_FAILED` result in
commit `5d43620`; it is not overwritten or relabelled.  Its arbitrary absolute
corruption thresholds were incompatible with the near-pole Frobenius scale even
though the line and rank tests passed.

The repair was preregistered after this disclosure in commit `fa798f5`.  At the
fresh P200G precision:

```text
missing-lambda residual / correct residual = 6.66e69 -- 7.28e69,
wrong-sign residual     / correct residual = 4.28e75 -- 4.68e75.
```

Both corruptions therefore destroy the line by far more than the frozen
scale-free `1e40` requirement.  Because the repaired criterion was chosen after
seeing P160, the qualification “after disclosed control repair” is permanent.

## Interpretation boundary

- **DERIVED COMPUTATIONAL, ADVERSARIALLY REPLICATED:** the complete frozen
  carrier/action intersection is zero in all nonhomogeneous sectors and exactly
  one-dimensional in the homogeneous sector.
- **STRUCTURAL:** the surviving line is the weak-lapse/homothetic canonical
  response preserving the old collective momentum.
- **OPEN:** whether the omitted pole/lapse equation annihilates this line or
  extends it to a tangent of the full equations.
- **NOT DERIVED:** gauge status, propagation, a physical clock/tick, `c`, `G`,
  Planck scales, inertia or particle masses.

This is the narrowest possible candidate direction for homogeneous evolution;
it is not yet evolution.

