# Protocol correction: two-map archive cardinality

Date: 2026-08-18

The first targeted execution is preserved in commit `b3cc397` with outcome

```text
TWO_STEP_FULL_TANGENT_CONTROL_FAILED
15/16 PASS
```

Every geometric, determinant, canonical and schedule test passed.  The sole
failure was the preregistered archive-cardinality assertion `224`.

That count omitted the factor of two for storing both the second-slab tangent
and the two-step product.  The complete declared archive contains

```text
2 schedule parities
* 7 minimal sectors
* 4 derivative variants
* 2 maps (T_2 and T_2 T_1)
* 4 arrays (midpoint, radius, defect midpoint, defect radius)
= 448 arrays.
```

The failed archive already contained exactly 448 uniquely named arrays and
had SHA-256

```text
ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d.
```

This correction replaces `224` by `448` in the protocol and verifier.  It
does not alter or recompute a geometry, Hessian, tangent, product, ball,
threshold, spectral quantity, schedule label or outcome rule.  The passing
rerun must reproduce the same numeric-archive SHA-256; otherwise this
correction is invalid and the result remains control-failed.
