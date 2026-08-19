# First-run correction: normalize the quaternion vertices and use the frozen rank bands

Date: 2026-08-19

The first execution of
`verify_gravity_600cell_dust_action_york_negative.py` was interrupted after a
preregistered geometry control failed, before accepting any target verdict.
It printed

```text
rank C=120, rank R=470, rank D=450, dim(im C intersection im D)=96
```

instead of the frozen theorem control `120/470/354/4`.

The cause is a harness mismatch.  The audited rigidity construction first
renormalizes every approximate quaternion coordinate to unit length and then
classifies singular values with the repository's `10/100` rank bands.  The new
wrapper omitted the normalization and called `numpy.linalg.matrix_rank`,
whose implicit threshold treated coordinate-level radial leakage as rank.
The printed radial identity residual `6.52e-10` and edge-length spread
`5.60e-11` exposed the mismatch directly.

The preregistered repair is:

1. normalize each `build_600cell()` vertex before forming tangent projectors;
2. replace implicit `matrix_rank` by the same explicit singular census used
   by the frozen rigidity verifier;
3. change no geometry, action matrix, selected sector, comparison threshold
   or outcome hierarchy.

The failed run produced no JSON artifact and was stopped while computing the
first all-sector carrier.  This correction is committed before rerunning.

