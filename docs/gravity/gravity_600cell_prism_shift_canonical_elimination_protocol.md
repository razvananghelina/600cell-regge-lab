# Protocol: canonical elimination of the 119 relative prism shifts

Date: 2026-08-19

Prior-art gate commit: `d90a44b`.

This is a target-disclosed composition protocol.  The two input theorems are
already frozen; no new Schur spectrum is being searched for.

## 1. Frozen inputs

Require byte-exact provenance for:

```text
reproducible/gravity_600cell_prism_shift_dynamic_extension.json
SHA-256 32d5269b27756a4c6fec4603855db643106e571007d3f3dd1a0a6c69d33a0095

reproducible/gravity_600cell_dust_full_lapse_schur.json
SHA-256 4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349

reproducible/gravity_600cell_dust_full_anisotropic_legendre_rank.json
SHA-256 7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226

reproducible/verify_gravity_600cell_dust_full_lapse_schur.py
SHA-256 7258899ba96a127515956fa2ea5fb17ad480373765b3f7c88fed40845adc82a6
```

The dynamic-extension artifact must report `13/13` and
`DYNAMIC_SHIFT_EXTENSION_OBSTRUCTED`.  The Schur artifact must report
`18/18`, both parity outcomes `FULL_LAPSE_SCHUR_REGULAR`, and rank `120` with
zero open and zero error-consistent null directions.

No continuum, speed, mass or experimental target may appear in the new
verifier.

## 2. Exact relative-pole carrier

Reconstruct the literal 600-cell vertex graph.  Use the exact map

```text
R : Q^119 -> Q^120,
R e_i = e_i-e_119,  i=0,...,118.
```

Require exact rank `119` over `F_101` and `F_1000003`.  Its image must be the
zero-sum relative-pole hyperplane, and the all-ones collective lapse must be
an independent complement.

The differential between a pole log-magnitude and its squared Lorentzian
length is multiplication by the nonzero background pole magnitude.  It
therefore preserves this rank and relative/collective split.

Verify the 600-cell incidence rank `119` over the same two fields.  This
identifies, dimension and kernel included, the relative pole coordinates
with the global vertex-potential differences from the exact frustum theorem.

## 3. Complete carrier and sector census

From the frozen artifacts require, for both parities,

```text
old boundary       720,
internal           840 = 720 non-pole strong + 120 poles,
new boundary       720,
canonical strong  1440,
canonical weak     120.
```

The five weak orbit positions must be exactly `[30,31,32,33,34]`.  The seven
minimal sector dimensions must be a permutation of

```text
1,1,1,2,2,2,3.
```

For a sector of irrep dimension `d`, independently parse the stored
operational Schur midpoint as a `5d x 5d` complex matrix.  Recompute its
singular values and require:

1. all four stored strong and Schur determinant balls exclude zero;
2. the strong midpoint census has `60d` resolved entries, zero open/zero;
3. the Schur midpoint census has `5d` resolved entries, zero open/zero;
4. every recomputed Schur singular value exceeds the stored
   `100*epsilon_global` boundary;
5. the recomputed singular multiset agrees with the stored one to relative
   `2e-12` with unit absolute floor.

Restore representation multiplicities and require exact counts

```text
sum d*(60d)=1440,
sum d*(5d)=120.
```

## 4. Graph-embedding invariance

Let a geometric polytopal realization place its cross-diagonal variations in
the strong graph `x=Gz+y`.  Verify exactly on a rational nonsymmetric control
that

```text
S_G=(C*G+D)-C*A^-1*(A*G+B)
   =D-C*A^-1*B=S.
```

The note supplies the dimension-independent algebraic proof.  The control
must use nonzero `B,C,G` and a non-diagonal invertible `A`, so the cancellation
cannot pass because all couplings vanished.

Also require from the carrier census that all 720 cross diagonals lie in the
strong internal coordinates, while all 120 poles lie in the weak coordinates.

## 5. Composition and negative control

The decisive implication is

```text
A invertible and S invertible:

A*x+B*z=0,
C*x+D*z=0

=> x=-A^-1*B*z,
=> S*z=0,
=> z=0.
```

Since `R` is injective, `S*R` is injective on all 119 relative coordinates.
Use the smallest certified sector singular value as a positive lower bound;
do not fit a coefficient.

Negative control: replace one frozen Schur sector by an exact zero matrix.
The reconstructed full rank must drop by `d*(5d)` for that sector and the
composition verdict must fail.  This guards against assigning elimination
from dimensions alone.

## 6. Verdicts

Return

```text
RELATIVE_SHIFT_CANONICALLY_ELIMINATED
```

only if all provenance, exact-rank, carrier, Schur, graph-invariance and
negative-control checks pass.

Return `RELATIVE_SHIFT_CANONICAL_STATUS_OPEN` if any Schur direction or
carrier identity is open or inconsistent.  Return
`RELATIVE_SHIFT_FREE_DIRECTION_SURVIVES` if a certified null survives on the
119-dimensional relative carrier.

## 7. Interpretation boundary

A passing result proves only a local homogeneous-equation statement at fixed
old geometry and incoming momentum:

- **DERIVED:** the 119 relative shift/strut modes are fixed, not free
  propagated canonical data;
- **DERIVED:** this statement is independent of the strong diagonal graph
  used to realize their polytopal geometry;
- **STRUCTURAL:** auxiliary/pseudo-constraint sector;
- **OPEN:** sourced boundary response, physical constraint quotient,
  gravitons, dispersion and limiting speed.

Only the new verifier and static guards may run; the full suite is excluded.
