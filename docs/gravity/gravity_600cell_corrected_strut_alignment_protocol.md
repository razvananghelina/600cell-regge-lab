# Confirmatory protocol: corrected strut carrier versus the weak and hyperbolic sectors

Date: 2026-08-19

Status: written only after commit `dab941b` froze the complete target-blind
carrier artifact.  No corrected-carrier/target angle or projector distance
has been evaluated.

## 1. Provenance and epistemic limit

Require exact SHA-256 values:

```text
corrected carrier source
80f0a17960adee496fe7d51678ea99849280ecd3fca6254efc8acd3753aad348

target-blind corrected carrier artifact
e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7

prior hyperbolic-alignment source and artifact
e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047
a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff

full-boundary tangent JSON, numeric archive and source
4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5
816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b
c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571

full pole-Schur artifact
4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349

accepted tick
4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9
```

The corrected carrier verifier must report `13/13`, outcome
`CORRECTED_STRUT_CARRIER_FROZEN`, exactly one candidate per parity and
`target_artifacts_loaded=false`.  The old alignment audit must retain
`14/14` and outcome `HYPERBOLIC_EXTREME_SUBSPACE_OPEN`.

This is not a blind prediction: the `119+1=120` count motivated the test.
The separate carrier commit proves absence of coefficient fitting, not
absence of hypothesis selection.  An exact subspace identity would be a
confirmatory finite theorem.  A merely small angle is **PATTERN**, regardless
of its visual appeal.

## 2. Reconstruct the corrected target coordinates

For each schedule parity, reconstruct the same 65 orbit types in the
pre-Legendre domain

```text
U = (840 internal,720 new boundary).
```

Map every lexicographically stored carrier row to the exact edge label in the
action ordering; do not match rows by orbit number alone.  Order the 120
columns by the five geometry-selected pole orbits and their frozen group
coordinates.  Require:

1. exact edge-label coverage of all 1,560 domain rows;
2. the pole block is the literal identity in the selected column ordering;
3. projection to every minimal `2T` sector has rank `5d`;
4. summing the corrected columns reproduces the old analytic collective
   lapse column;
5. rebuilding the matrix directly from
   `rho/((lambda-1)q_diag)*(c_u-lambda c_v)` agrees with the committed rows
   below `1e-70` before binary64 projection.

No column rescaling, mixing, graph choice or sector-dependent convention is
allowed.

## 3. Frozen dynamic operators

Use the existing four high-precision derivative variants, seven deterministic
minimal `2T` sectors and both staircase parities.  Reconstruct in each sector

```text
J = [[K_XX,K_XN],[-K_OX,-K_ON]],
R = [[-K_XO,0],[K_OO,I]],
Y = J^-1 R.
```

Require all Flint determinant balls and the old response/branch/reality
controls to pass.  Select the five pole positions from edge geometry before
loading a tangent spectrum.

Form the canonical weak lift

```text
C_weak = reorder_back((-A^-1 B, I)^T)
```

and the corrected geometric carrier `G_corr`.  Separately select, without
loading either candidate, the `k=5d` largest- and smallest-modulus invariant
subspaces of each frozen tangent block and transport them into `U`:

```text
U_plus  = colspace(Y E_plus),
U_minus = colspace(Y E_minus).
```

Retain the old fixed-count Schur selection, direct-eigenvector control and
gap requirement `>2`.  The known homogeneous fifth pair may keep the global
extreme result OPEN; it may not be dropped after inspection.

## 4. Fixed comparisons and calibration

For every parity-sector pair evaluate exactly three projector distances:

```text
G_corr versus C_weak,
G_corr versus U_plus,
G_corr versus U_minus.
```

Thus the look-elsewhere ledger is exactly

```text
2 parities x 7 sectors x 3 comparisons = 42.
```

Do not compare arbitrary mixtures of plus/minus, old/new geometric carriers,
scale columns or sectors and then retain the smallest angle.

Use the same principal-angle definition as the frozen alignment audit:

```text
distance = ||Q1 Q1* - Q2 Q2*||_2 = sin(theta_max).
```

For each named distance combine:

- all four derivative-step discrepancies;
- direct-eigenvector versus reordered-Schur discrepancy after transport;
- first-order Flint radius bounds divided by minimum column singular values;
- the committed carrier's 100-digit/binary64 discrepancy;
- `10 eps_machine` times the largest reported condition number.

Call the sum `epsilon_distance` and assign without alteration:

```text
IDENTIFIED       distance <= 10 epsilon_distance,
SEPARATED        distance > 100 epsilon_distance,
NUMERICALLY_OPEN otherwise.
```

Also report the uncalibrated angle, minimum overlap and all four raw
distances.  A distance of order `1e-5` is not a hit merely because it is
small.

## 5. Negative and convention controls

1. Reproduce the old committed canonical-versus-old-geometric distances
   within `2e-8`; this verifies target ordering without accepting that old
   carrier.
2. The corrected and old geometric carriers must be distinct in every
   non-uniform sector while agreeing on the collective column.
3. Reverse source/target coefficient roles in the first corrected diagonal.
   At least one projected sector and one comparison must change above its
   calibrated carrier roundoff.
4. Even and odd parities are classified separately before any combined
   statement.

## 6. Outcome hierarchy

Assign in order:

1. `CORRECTED_STRUT_ALIGNMENT_CONTROL_FAILED` for any provenance, edge-map,
   rank, response, determinant, reproduction or corruption failure;
2. `CORRECTED_STRUT_EXTREME_SELECTION_OPEN` if any fixed-count/gap gate is
   open, while still reporting every resolved sector comparison;
3. `CORRECTED_STRUT_CANONICAL_AND_PLUS_IDENTIFIED` if both corrected/canonical
   and corrected/plus are `IDENTIFIED` in all 14 parity-sector cases;
4. `CORRECTED_STRUT_CANONICAL_AND_MINUS_IDENTIFIED` for the analogous fixed
   minus branch;
5. `CORRECTED_STRUT_CANONICAL_IDENTIFIED_ONLY` if corrected/canonical is
   identified in all 14 but neither fixed extreme branch is;
6. `CORRECTED_STRUT_ALIGNMENT_REFUTED` if all 42 comparisons are
   `SEPARATED`;
7. `CORRECTED_STRUT_ALIGNMENT_MIXED_OR_OPEN` otherwise.

The verifier passes when the frozen object is reconstructed and classified,
including honest OPEN or negative outcomes.

## 7. Claim boundary and next gate

Even a global subspace identity does not prove gauge freedom: the corrected
carrier is kinematic and the curved dust Schur block is regular, not null.
Only suppressed or vanishing response of derived four-dimensional deficit
angles and intrinsic boundary curvature could support a pseudo-gauge
interpretation.  A separated result likewise does not identify gravitons.

No outcome derives a multi-tick instability, continuum dispersion, limiting
speed, absolute tick, `G`, Planck scale or particle mass.  Run only the new
targeted verifier and the static registry guard.

