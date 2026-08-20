# Protocol: complete scale--strut/canonical intersection census

Date: 2026-08-20

Status: no singular value, rank or nullity of an actual complete-carrier
intersection block has been evaluated.

## 1. Frozen provenance

Require the following literal SHA-256 values:

```text
prior-art gate
134ece5926b429011a1b74428a30454924c1458492f2aed3244bab9258b345c3

accepted carrier result
774dae7fbe3d3becf505867c3272f41f800ab9c917766fcfd347395e36c34ece

finite full-carrier source
e68105df4058f7d2ed39a6913f29e88cd9fe88e123ff52260acf698a2bd7da49

finite full-carrier artifact
6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d

resolved carrier precision artifact
2a2a79271a92fc2ddde343a9d0651402df6eeb4a90efa2697e26f54cafcdf60f

resolved symbolic-gap artifact
ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179

frozen action-response source
e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047

frozen action-response artifact
a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff

frozen tangent archive
816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b

pure-strut zero-intersection artifact
422d8d8cb0fc0d72d842e3bf79609d4d985da6237c58e7c699b5f9cc21b65cec
```

The carrier inputs must retain their exact-rank, precision-resolved and
symbolic-real-gap-resolved outcomes.  The pure-strut artifact must retain 14
resolved zero nullities.  Re-executing the frozen response source must leave
its artifact byte-identical.

## 2. Independent carrier reconstruction

Do not import the carrier builder.  Reconstruct `G` independently from the
artifact's 840 internal and 720 final edge orders, background decimals and
the disclosed endpoint formula.  Also reconstruct it from the artifact's
stored coefficient decimals.  Require equality within 100-decimal
arithmetic before conversion to binary64.

Map its 120 vertex labels into the same five regular-representation pole
orbits used by the action response.  Require complete edge coverage and the
literal blocks

```text
G_scale|pole = 0,
G_strut|pole = I.
```

Rebuild the seven deterministic minimal `2T` sectors and four frozen
action-derivative variants exactly as in the mechanically accepted response
source.  In every sector require full column ranks `rank(G)=10d` and
`rank(C)=5d`, and `C|pole=I` below `1e-13`.

## 3. Target-blind rank objects

For every parity and sector form both

```text
D = [G_scale, G_strut-C]        shape (65d,10d),
K = [G_scale, G_strut,-C]       shape (65d,15d).
```

The nullities of `D` and `K` must agree; that common number is the
intersection dimension in the minimal sector.  The global real/complex
carrier count is reported separately as the sum of each minimal nullity
times the regular multiplicity `d`.  No expected nullity or continuum label
is loaded.

Because raw scale coefficients are large, apply the frozen symmetry-
preserving block equilibration

```text
s_scale = 1/max(1, ||G_scale||_2),
s_strut = 1/max(1, ||G_strut||_2, max_variant ||C_variant||_2).
```

Use `diag(s_scale I_(5d),s_strut I_(5d))` for `D`, and the corresponding
three blocks for `K`.  These scalars are positive and hence cannot change a
nullity.  Record unscaled and scaled condition diagnostics; only scaled
matrices decide the operational rank.

## 4. Frozen numerical calibration

For each sector use both `gesdd` and `gesvd`.  The absolute perturbation
bound is the sum, after the frozen scaling, of:

1. the maximum Frobenius radius of the four Flint lift balls;
2. the largest spectral difference among the four derivative variants;
3. the discrepancy between the two independently reconstructed carriers;
4. the resolved full-carrier high-precision spectral discrepancy multiplied
   by the scaled carrier norm;
5. `200*eps_machine*max(rows,columns)*matrix_norm` for conversion and SVD
   backward error.

Classify each singular value as

```text
ZERO       sigma <= 10 epsilon,
NONZERO    sigma > 100 epsilon,
OPEN       otherwise.
```

A sector is resolved only if neither driver nor any derivative variant has
an `OPEN` value, all variants and both drivers give the same nullity, and
the reduced `D` and joined `K` nullities agree.  Store any resolved kernel as
an orthogonal projector in the equilibrated coefficient coordinates, not as
an arbitrary basis.

## 5. Hostile controls

1. A zero matrix under the same threshold must have nullity `10d`; an
   embedded identity must have nullity zero.
2. Replacing `C` by `G_strut` must give `nullity(D)=5d`, provided the checked
   scale half has rank `5d`.
3. A frozen nonsingular upper-triangular change of all source coefficients
   must preserve the actual nullity.
4. The independent image formula
   `rank(G)+rank(C)-rank([G,C])` must reproduce it.
5. Reverse the source/target roles of the first staircase diagonal; the
   scaled matrix and at least one singular value must change above the
   reconstruction error.
6. Complex conjugation and the two parities are recorded independently;
   equality may be observed but never assumed.

## 6. Outcome hierarchy

1. `FULL_SCALE_STRUT_CANONICAL_CONTROL_FAILED`: provenance, geometry,
   identities, ranks or hostile controls fail.
2. `FULL_SCALE_STRUT_CANONICAL_NUMERICALLY_OPEN`: at least one actual block
   remains in the open band or variants/drivers disagree without a hard
   control failure.
3. `FULL_SCALE_STRUT_CANONICAL_INTERSECTION_RESOLVED`: all 14 actual blocks
   and both independent intersection formulas are resolved.

The verifier passes for outcome 2 or 3 and fails for outcome 1.  Any material
resolved result requires a mechanically different adversarial verifier
before consolidation.

## 7. Claim boundary

A resolved nonzero intersection selects first-order stationary candidates
inside this carrier.  It does not establish gauge invariance, propagation or
a graviton.  A resolved zero intersection closes this particular accepted
carrier on this fixed slab, but does not refute Regge calculus or every
higher-rank/nonlocal carrier.

No outcome supplies a multi-tick solution, a clock, `c`, `G`, Planck units,
particle masses or Standard-Model physics.  Run only the new verifier and
static registry guards; do not run the full suite.

