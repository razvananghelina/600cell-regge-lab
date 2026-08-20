# Second multiprecision resolver failure: Hermitian eigensolver dispatch

Date: 2026-08-20

- exact-integer conversion repair: `9d4847f`;
- failing source SHA-256:
  `416933d7a91cba7891bf05eeb7fb51bcbd09413a7a88d5c19189765d33437b8c`.

The repaired run again passed provenance, the action rebuild, P100 synthetic
determinants and P100/even geometry.  On the first actual sector the Flint
midpoint conversion completed, then `mp.eigsy` rejected the complex
Hermitian Gram matrix because `eigsy` is the real-symmetric routine.

No actual spectrum or output artifact was produced.  **DERIVED SOFTWARE
NEGATIVE:** dispatch was wrong; the matrix was not.

The frozen correction is exactly

```text
mp.eigsy(gram) -> mp.eighe(gram).
```

`eighe` is mpmath's Hermitian eigensolver and returns ordered real
eigenvalues plus complex eigenvectors.  No other source line, precision,
matrix, determinant or criterion may change before the next targeted run.

