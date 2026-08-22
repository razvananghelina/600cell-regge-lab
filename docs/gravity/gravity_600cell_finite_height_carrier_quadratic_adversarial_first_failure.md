# First adversarial execution: misplaced reality gate and serializer failure

Date: 2026-08-22.

Implementation commit: `dcb0d1c`.

Status: **DERIVED INFRASTRUCTURE FAILURE / SCIENTIFIC VERDICT NOT ACCEPTED.**

This note preserves the first execution before any rerun.  The process ran
the complete preregistered numerical calculation, then returned

```text
FINITE_HEIGHT_QUADRATIC_ADVERSARIAL_CONTROL_FAILED
```

for one failed control and finally raised a JSON serialization exception.
It is not counted as a successful adversarial replication.

## 1. Failed control

The implementation incorrectly required the raw Lorentzian dihedral-angle
and derivative entries to have imaginary magnitude below `1e-140`.  Both
parities reported

```text
maximum raw angle/derivative imaginary component = 105.028...
```

This is not contamination: the Lorentzian boost-angle contribution is
complex before the explicit `-i` in the Regge action.  The physical assembled
kernel did pass the frozen reality scale:

```text
even maximum kernel imaginary residue = 9.66126e-155
odd  maximum kernel imaginary residue = 1.51905e-154
```

The input high-precision source
`verify_gravity_600cell_dust_full_boundary_tangent.py` also gates the final
kernel rather than the raw boost-angle table.  The adversarial protocol has
therefore been amended before any rerun so that its `1e-140` reality gate is
applied to the physical orbit kernel and complete scalar action.  No numerical
threshold or decisive parity criterion was changed.

## 2. Scientific diagnostics observed but not yet accepted

The failed run printed

```text
d_R01 = 1.0946049e-101
d_R12 = 6.8412809e-103
d_R23 = 4.2758006e-104
e_step = e_total = 2.8057204e-101 to the displayed precision

primary/adversarial relative matrix error:
  even = 3.357663923742695e-14
  odd  = 3.317576498395870e-14

carrier corruption effect = 4.9925095e-5
synthetic rank-one effect  = 3.6162624e-2
```

All four direct complete-action controls reproduced the orbit-kernel
quadratic values:

```text
even uniform scale       1.0339853e-88
even local scale/strut   3.9318572e-85
odd  uniform scale       1.0339853e-88
odd  local scale/strut   3.9319025e-85
```

These values strongly favour parity independence, but under the frozen
hierarchy they remain **STRUCTURAL / UNACCEPTED** until a corrected registered
execution passes every control.

## 3. Serializer failure and preserved binary artifact

After writing the binary matrix stack, JSON serialization encountered an
unconverted `mpmath.mpf` in a nested pull-control record:

```text
TypeError: Object of type mpf is not JSON serializable
```

The pre-exception matrix stack is preserved as

```text
reproducible/gravity_600cell_finite_height_carrier_quadratic_adversarial_first_failure_matrices.npy
sha256 8a3ea0c3b8ee720d8ffdf07e7486aefdd0247ca1cfdbeb99f443091376f31729
shape  (3,240,240)
order  even_R12, odd_R12, even_minus_odd
```

The serializer fix converts that diagnostic scalar to text only.  It changes
no arithmetic or classification path.
