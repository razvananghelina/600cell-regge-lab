# First multiprecision resolver failure: Flint integer conversion

Date: 2026-08-20

## Frozen execution

- corrected protocol commit: `f011db5`;
- registered source commit: `55ef589`;
- source SHA-256:
  `8fae55534637d972f547f961905d6854e0b407a49c28fb68e51837277148e69e`.

The run passed provenance, the byte-identical action reconstruction, the P100
synthetic determinant controls and P100/even geometry.  It then stopped on
the first actual sector before any actual Gram determinant or singular
spectrum was recorded.

Python-flint returned the exact binary midpoint mantissa and exponent as
`fmpz` objects.  `mpmath.mpf` does not implicitly accept an `fmpz`, producing

```text
TypeError: cannot create mpf from <large fmpz>.
```

No output artifact exists.  **DERIVED SOFTWARE NEGATIVE:** this is an adapter
type failure and carries no rank information.

## Frozen correction

Change only

```text
mp.mpf(mantissa) * 2^exponent
```

to

```text
mp.mpf(int(mantissa)) * 2^int(exponent).
```

The conversion is exact: both quantities are integers and no binary float is
introduced.  No matrix, precision, derivative step, determinant, threshold,
candidate or outcome rule may change.  Commit the repair before rerunning
only the targeted verifier.

