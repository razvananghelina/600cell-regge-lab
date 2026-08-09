# S08 Closure: Edge Space and Fiber Permutation Module

## Exact claim

For the full six-fibration left-coset Hopf class of the 600-cell, determine
the exact status of:

1. the edge-space kernel of
   \[
   \Box_1(F)=L_{\mathrm{cross}}(F)-a_1 L_{\mathrm{fiber}}(F),
   \qquad a_1=5,
   \]
2. the factorization of the fiber action through `A_5`,
3. the 12-dimensional fiber permutation module.

## Formal setting

- finite line-graph operator on the 720 edges of the 600-cell;
- discrete Hopf fibrations coming from left cosets of order-10 subgroups of
  `2I`;
- exact character decomposition under the left action of `2I`.

## Inputs used

- [s06_hopf_closure.md](D:\infinity\ToE\science\s06_hopf_closure.md)
- [s07_wave_coefficient_closure.md](D:\infinity\ToE\science\s07_wave_coefficient_closure.md)
- [reproducible/verify_edge_gauge_spectrum.py](D:\infinity\ToE\science\reproducible\verify_edge_gauge_spectrum.py)
- [reproducible/verify_s08_edge_fibration_uniformity.py](D:\infinity\ToE\science\reproducible\verify_s08_edge_fibration_uniformity.py)

## Output status

- `Computational fact`
- `Theorem` for the quotient-factorization part
- `Open` for any identification between the two 12-dimensional spaces

## Proof / derivation

The focused verifier `verify_s08_edge_fibration_uniformity.py` checks all six
distinct left-coset Hopf fibrations and confirms:

1. for every fibration `F`,
   \[
   \dim \ker(\Box_1(F)) = 13;
   \]
2. for every `F`, the full `2I`-character decomposition is
   \[
   \ker(\Box_1(F))=\rho_0\oplus 2\rho_5;
   \]
3. for every `F`, the central element `-1` acts trivially on the 12 fiber
   labels, so the permutation action factors through
   \[
   A_5 \cong 2I/\{\pm 1\};
   \]
4. for every `F`, the resulting 12-dimensional fiber permutation module
   decomposes as
   \[
   \mathbf{1}\oplus \mathbf{3}\oplus \mathbf{3'}\oplus \mathbf{5}.
   \]

## What is still not proved

The 12-dimensional nontrivial edge-kernel sector
\[
2\rho_5
\]
and the 12-dimensional fiber permutation module
\[
\mathbf{3}\oplus \mathbf{3'}\oplus \mathbf{5}
\]
live in different ambient spaces:

- edge functions on `R^{720}`;
- fiber-label functions on `R^{12}`.

This closure does **not** prove any canonical `A_5`-equivariant map between
them, and does **not** identify them as the same module.

## Decision

`S08` is closed in the weak uniform form:

- `Computational fact`:
  `ker(Box_1(F))` has dimension 13 and decomposes as `rho_0 + 2 rho_5` for
  all 6 fibrations;
- `Theorem`:
  the fiber action factors through `A_5`;
- `Computational fact`:
  the fiber permutation module is `1 + 3 + 3' + 5` for all 6 fibrations;
- `Open`:
  no identification is asserted between the 12-dimensional edge sector and the
  12-dimensional fiber permutation module.
