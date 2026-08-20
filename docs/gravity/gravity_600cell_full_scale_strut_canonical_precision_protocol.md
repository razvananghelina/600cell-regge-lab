# Protocol: multiprecision resolution of the complete-carrier intersection

Date: 2026-08-20

Status: the `1e-7` and homogeneous `1e-16` target patterns are disclosed,
but no multiprecision projected carrier, interval Gram determinant or
higher-precision candidate residual has been computed.

## 1. Frozen provenance

Require these SHA-256 values literally:

```text
precision prior-art gate
fd2e230fdc0c0f7aaa771a4781973d0b476b5758d8d82f5b84ba471b391a722c

first OPEN result note
971b1eddcb09e4a35a72ba5d1c359b710ba74673815824c1029ddb68e897bfc5

corrected primary source
a2d5390d39c725a5fb586fefce9da34cede3a1fb84bbe36791f8b0599b3eae42

frozen primary OPEN artifact
b29cc33a9effeb2087fb6133359ee747d100d203778586372a7ceeebc2e4f070

frozen first classifier failure
6423c3efc03ba6107a82c1b0d813e0226ccf757d242cc3ecc0522003095e97d5

resolved symbolic carrier
ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179

action-response source
e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047

action-response artifact
a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff

accepted homogeneous canonical root
4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9
```

The primary artifact must retain `13/13`, the numerically OPEN outcome and
14 unresolved sectors.  No target value may be read from a later artifact.

## 2. Two frozen precision levels

Rebuild all geometry, deterministic minimal-sector bases, action kernels and
carrier projections at both levels:

```text
level P100:
    mp.dps=100, Flint dps=80,
    derivative steps=(1e-20,1e-15,3e-20,3e-15);

level P160:
    mp.dps=160, Flint dps=140,
    derivative steps=(1e-40,1e-30,3e-40,3e-30).
```

At each level rebuild the projected carrier directly from sparse edge roles
and multiprecision coefficients.  Do not create a binary64 global matrix and
do not use the old global spectral-discrepancy proxy.  Project each sparse
entry with the current multiprecision `2T` basis.

Reconstruct the canonical lift as an `acb_mat` solve and retain the full ball
matrix.  Also extract its midpoint only for diagnostic multiprecision Gram
eigenvalues.  Decimal conversion into a ball must use at least `dps-15`
digits.  The two staircase parities remain separate.

Use the already-frozen positive block scalings from the primary artifact;
changing coordinates after seeing the high-precision result is forbidden.

## 3. Full-rank certificates

For both

```text
D = [G_scale,G_strut-C],
K = [G_scale,G_strut,-C],
```

form the Hermitian interval Gram matrix and its determinant.  A determinant
certifies full column rank only if its `acb` ball excludes zero.  Record
midpoint, radius, `abs_lower`, `abs_upper` and zero containment at both
precision levels.

Independently compute midpoint Gram eigenvalues with `mp.eigsy`.  For a
claimed nonzero singular value, P100/P160 relative disagreement must be less
than `1e-20` and the P160 value must exceed the upper propagated ball/midpoint
uncertainty by at least `1e20`.  Failure is OPEN, not zero.

All six non-homogeneous sectors are full rank only if both `D` and `K`
determinant balls exclude zero at both levels and the smallest midpoint
singular values pass this stability rule.

## 4. Homogeneous exhaustive rank and no-refit candidate test

The trivial sector is identified mechanically by constant overlap one, not
by its position in the list.  In addition to the full determinants, compute
every single-column-deleted Gram determinant:

```text
10 minors of the 65 x 10 reduced D,
15 minors of the 65 x 15 joined K.
```

At least one P160 minor in each family must exclude zero to certify ranks at
least 9 and 14.

At P100, choose the normalized smallest right Gram eigenvector of `D` using
the deterministic phase convention that its largest-magnitude component is
positive real.  Round each component to 70 decimal digits and freeze it in
memory before constructing P160.  Do not refit at P160.  Validate:

```text
||D_P160 v_P100||_2 < 1e-50,
||K_P160 (sigma,strut,-strut)||_2 < 1e-50,
both residuals < 1e-40 times the next P160 singular value.
```

For a calibrated one-dimensional homogeneous kernel additionally require:

- full `D` and `K` Gram determinants contain zero at both levels;
- the exhaustive minors certify ranks at least 9 and 14;
- each midpoint spectrum has exactly one value below `1e-50` at P160;
- each next singular value exceeds `1e-8`;
- the smallest singular values decrease by at least `1e20` from P100 to
  P160 while the next values agree relatively within `1e-20`;
- the no-refit candidate tests above pass for both parities;
- the two parity candidate projectors agree within `1e-30` after the frozen
  coefficient ordering is checked.

If full determinants instead exclude zero stably, the homogeneous nullity is
zero.  If neither branch meets every gate, it remains OPEN.  A determinant
merely containing zero never proves a kernel.

## 5. Hostile controls

1. Synthetic full-column-rank and planted-nullity-one matrices at both
   scales must be classified correctly by determinants/minors and midpoint
   spectra.
2. Deleting one pole identity entry from `G_strut` must change the
   homogeneous candidate residual by more than `1e20` relative to its
   uncorrupted value and must change at least one determinant ball.
3. Reversing the first diagonal's source/target roles must change a
   non-homogeneous smallest singular value by more than its complete P160
   uncertainty.
4. `D` and `K` must give the same resolved intersection nullity.
5. Complex conjugation preserves each midpoint singular spectrum within the
   P160 arithmetic bound.

## 6. Outcome hierarchy

1. `FULL_SCALE_STRUT_CANONICAL_PRECISION_CONTROL_FAILED` for provenance,
   reconstruction, interval, parity-ordering or hostile-control failure.
2. `FULL_SCALE_STRUT_CANONICAL_PRECISION_DISAGREEMENT` if a stable result
   contradicts the disclosed binary pattern.
3. `FULL_SCALE_STRUT_CANONICAL_NONHOMOGENEOUS_OPEN` if any non-homogeneous
   sector remains unresolved.
4. `FULL_SCALE_STRUT_CANONICAL_HOMOGENEOUS_OPEN` if all non-homogeneous
   sectors are full rank but the trivial sector is unresolved.
5. `FULL_SCALE_STRUT_CANONICAL_ZERO_INTERSECTION_RESOLVED` if every sector,
   including the homogeneous sector, is certified full rank.
6. `FULL_SCALE_STRUT_CANONICAL_ONE_HOMOGENEOUS_RESOLVED` if every
   non-homogeneous sector is certified full rank and the complete
   homogeneous one-dimensional gate passes.

Outcomes 3 and 4 are honest passing OPEN outcomes.  Outcomes 5 and 6 are
primary material results and require a mechanically different adversarial
verifier before consolidation.

## 7. Claim boundary

Outcome 6 would select one homogeneous first-order canonical direction and
exclude non-homogeneous directions only inside this fixed 240-dimensional
carrier.  It would not classify that direction as gauge or physical, prove
propagation, or derive time, `c`, `G`, Planck units or particle masses.

Run only the new verifier and static registry guards.  Do not run the full
suite.

