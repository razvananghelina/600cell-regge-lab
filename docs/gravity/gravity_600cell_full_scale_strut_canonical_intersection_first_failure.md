# First complete-carrier intersection run: classifier failure

Date: 2026-08-20

## Frozen execution

- protocol commit: `b621736`;
- registered source commit: `8e20920`;
- first artifact commit: `61fc3cc`;
- source SHA-256:
  `990af0c836f78e648e96fd6a916a48a4389ccad2dc8ca896c5b3383fcf36904d`;
- artifact SHA-256:
  `6423c3efc03ba6107a82c1b0d813e0226ccf757d242cc3ecc0522003095e97d5`;
- result: `11/13`, `FULL_SCALE_STRUT_CANONICAL_CONTROL_FAILED`.

## What actually happened

**DERIVED SOFTWARE NEGATIVE.**  All 14 actual sectors were numerically open.
The operational smallest scaled singular values included values around
`0.8e-7 ... 1.3e-7` against `epsilon ~= 0.81e-8`, and the homogeneous sector
also contained values around `2e-16 ... 4e-16`.  Thus the preregistered
`10 epsilon / 100 epsilon` rule deliberately refuses a rank decision.

The implementation then made an invalid classifier implication: the
change-of-basis and joined-image controls were coded as Boolean failures
whenever `nullity is None`.  The structural `C=G_strut` control likewise
called an `OPEN` classification a hard failure.  Therefore numerical
uncertainty was incorrectly promoted to `CONTROL_FAILED`, contrary to the
protocol's separate `NUMERICALLY_OPEN` outcome.

The run did **not** find a failed pole identity, wrong input rank, provenance
change or undetected corruption.  It passed byte-identical action rebuild,
two independent carrier reconstructions, all geometry controls and the
source/target corruption.  No intersection dimension is accepted from this
run.

