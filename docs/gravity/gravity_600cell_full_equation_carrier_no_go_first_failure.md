# First full-equation carrier no-go run: coverage serialization mismatch

Date: 2026-08-20  
Status: **preserved 5/6 coverage failure; correction frozen before source edit**

## Failed outcome

Source commit `a9e93ff` returned

```text
FULL_EQUATION_CARRIER_COVERAGE_OPEN
5/6 tests passed
```

Failed artifact SHA-256:

```text
ab001d13aeeb1d34bd6005d621748b7f404c5b782a49191689805d4c4bb1329a
```

The fourteen-cell dimension composition, both homogeneous transversality
checks and the pole-null negative control passed.  The failed assertion was
the nonhomogeneous certificate census.

## Exact cause

The upstream adversarial artifact reports `direct_minor_certificate_count=48`
but serializes them as `24` aggregate records:

```text
2 parities x 6 sectors x 2 matrices (D,K) = 24 records.
```

Each record contains two independently cross-precision direct minors:

```text
P100 rows applied to the P160 matrix,
P160 rows applied to the P100 matrix.
```

Thus `24 records x 2 minors = 48 certificates`.  The first verifier incorrectly
required `len(records)==48`.  It did not inspect either nested certificate.
This is a coverage-parser error, not a rank disagreement.

## Frozen correction

Before editing or rerunning the verifier:

1. preserve the failed artifact and pin it by SHA-256;
2. require exactly 24 unique `(parity,sector,matrix)` records covering
   `2 x 6 x 2`;
3. require every record's `pass` flag;
4. require `contains_zero=false` in both nested cross-precision minors;
5. count the resulting 48 nested certificates and compare that count with the
   frozen top-level field;
6. write a new `_corrected.json` artifact without overwriting the failure.

No dimension, outcome hierarchy or physical interpretation may change.

