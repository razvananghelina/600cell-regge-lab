# S09 Edge-Endomorphism Test: No-Go

## Exact claim

Test whether the exact nontrivial edge sector
\[
G_F := \ker(\Box_1)\cap \rho_0^\perp \cong 2\rho_5
\]
admits an intrinsic gauge-algebra route through the canonical algebra of
`2I`-equivariant endomorphisms.

## Formal setting

- object: the exact 12-dimensional edge sector `G_F = 2 rho_5`;
- symmetry: `2I`;
- canonical algebra:
  \[
  \mathrm{End}_{2I}(G_F).
  \]

## Output status

- `No-go theorem`

## Proof / derivation

Let
\[
D := \mathrm{End}_{2I}(\rho_5).
\]
By Schur's lemma over `\mathbb{R}`, one has
\[
D \in \{\mathbb{R},\mathbb{C},\mathbb{H}\}.
\]
The Frobenius-Schur indicator of the 6-dimensional irrep `\rho_5` is
\[
\nu(\rho_5)=-1,
\]
so `\rho_5` is quaternionic and therefore
\[
D \cong \mathbb{H}.
\]
Hence
\[
\mathrm{End}_{2I}(G_F)\cong \mathrm{End}_{2I}(2\rho_5)\cong M_2(\mathbb{H}).
\]

The canonical compact Lie algebra attached to this multiplicity space is then
\[
\mathfrak{sp}(2)\cong \mathfrak{usp}(4),
\]
which has dimension `10`.

This is incompatible with
\[
\mathfrak{u}(1)\oplus \mathfrak{su}(2)\oplus \mathfrak{su}(3),
\]
which has dimension `12`.

Therefore the intrinsic endomorphism route from the exact edge sector does not
recover the target gauge algebra.

## Computational confirmation

The focused verifier
[reproducible/verify_edge_endomorphism_type.py](D:\infinity\ToE\science\reproducible\verify_edge_endomorphism_type.py)
checks:

1. `ker(Box_1)=rho_0 + 2 rho_5`;
2. `dim(rho_5)=6`;
3. the Frobenius-Schur indicator is exactly `-1`;
4. therefore `End_{2I}(rho_5)=H`;
5. hence the canonical compact Lie algebra is `sp(2)`, not the SM one.

## Decision

This route is closed as `No-go theorem`.

Consequence:

- the direct edge-endomorphism test fails;
- by the agreed decision rule, no further gauge-recovery attempts should be
  made in the current cycle;
- the program should now move to option `A`: exact discrete precursor plus
  non-gauge structural sectors only.
