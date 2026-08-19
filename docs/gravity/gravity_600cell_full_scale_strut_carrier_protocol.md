# Protocol: target-disclosed full scale--strut carrier audit

Date: 2026-08-19

This protocol is committed after the formula was disclosed and before any
new rational representative or accepted-background `240`-column matrix is
evaluated.  It is a falsification protocol, not a blind discovery protocol.

## Frozen provenance

| input | SHA-256 |
|---|---|
| exploratory disclosure | `e3dba59118e35cc2370beec4b081cc18fdfd1753fda74a6ca2b0a013d86bd473` |
| prior-art gate | `3fc6c3e75ad92c3c20bb420d97e26e73fdd62b69a9aac44162fa95c74c29219a` |
| universal local-lift source | `4389861a4b64d043325e0661ae9c2340f61e5c8eb50399c9fd2083a360dadbc1` |
| universal local-lift artifact | `0a569e48189c56bc081efcee33f7826fedd52afb93b6135ddb2fec385b56fbdf` |
| universal local-lift result | `646972a19450f1734ef522cb0b9693cc809b19d7895eb21823b20332a958d56d` |
| corrected strut source | `80f0a17960adee496fe7d51678ea99849280ecd3fca6254efc8acd3753aad348` |
| corrected strut artifact | `e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7` |
| two-frustum adversarial source | `b7a1f63e193aad50783929c8448ce99c18f1b50dc8e5ea27e3ed1102ec9dfa26` |
| two-frustum adversarial artifact | `0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542` |
| two-frustum consolidated result | `b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3` |
| exact face-equation source | `4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e` |
| schedule geometry source | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| accepted lapse artifact | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |

The verifier must also check the exact face-equation source against the hash
embedded in the frozen universal local-lift verifier.  It may use schedule
combinatorics and the accepted background state, but it must not load an
action value, Hessian, strong-equation matrix, sector target, continuum label
or desired nullity.

## Disclosed candidate

For a lower-to-upper diagonal `i->j`, the candidate raw response is

```text
delta d_(i->j)^2 = A sigma_i+B sigma_j+C s_i+D s_j,

A = 6 - 2 tau^2/(lambda-1)^2,
B = 2 + 2 tau^2/(lambda-1)^2,
C = -1/(lambda-1),
D =  lambda/(lambda-1).
```

The complete hypotheses and action-coordinate conversion are frozen in the
disclosure.  No coefficient may be changed after seeing a result.

## Primary method: global exact controls not used in discovery

The exploratory derivation used one symbolic cell followed by one shared
face.  The primary verifier must instead reconstruct the complete 600-cell
face-equation system and solve its universal `6 x 8` local block at three new
exact representatives:

```text
(lambda,tau)=(4,7), (5,13), (7,17).
```

These points were selected before computation and satisfy every disclosed
nondegeneracy condition.  They were not used by the old local-lift artifact.

The solver must be mechanically different from the frozen sparse affine
solver:

1. scan the complete affine rows modulo the fixed prime `1000003` only to
   select 48 candidate-independent rows;
2. solve that `48 x 48` subsystem exactly over `Q`;
3. verify the solution against **every** exact global affine row;
4. reconstruct the physical `16 x 8` displacement response;
5. differentiate all twelve oriented lower-to-upper cross diagonals directly
   and compare every coefficient with the disclosed formula.

The modular scan may select rows but cannot certify the result.  Exact
rational invertibility and zero residual on every row are required.  Record
the selected row indices, determinant, block hash, full residual count and
all 96 cross-response coefficients at each representative.

The two old representatives `(2,5)` and `(3,11)` are controls only.  They
must still reproduce the frozen physical responses, but they carry no new
evidential weight.

## Required local and corruption controls

Independently rebuild one exact one-frustum squared-length Jacobian and
verify that the endpoint ansatz yields exactly

```text
A+B=8, C+D=1
```

and no further independent equation.  This prevents the verifier from
silently promoting local underdetermination to local uniqueness.

At each new representative replace `D` by `D+1`.  The corrupted candidate
must disagree with at least one directly reconstructed cross coefficient.
Also require the frozen two-frustum artifact to retain its one-dimensional
unrestricted face stabilizer and diagonal-only constrained motion.

## Accepted-background target-blind census

Instantiate one `1560 x 240` candidate matrix for each staircase parity in
the exact row order

```text
840 internal edges = 720 diagonals + 120 poles,
720 upper-boundary edges.
```

Columns are ordered as 120 `sigma` coordinates followed by 120 logarithmic
strut-magnitude coordinates.  Require:

```text
row-support histogram       {1:120, 2:720, 4:720},
scale-column support        24 for every vertex,
strut-column support        13 for every vertex,
literal pole identity       rank 120,
upper unsigned incidence    rank 120,
complete carrier rank       240.
```

The rank proof must use the pole identity plus the connected non-bipartite
600-cell graph, not a floating threshold.  Both schedule stabilizers must
intertwine the support and coefficient roles exactly.

Require the two collective identities

```text
sum_v scale_v:
    pole=0, diagonal=L0^2/q_diag, upper edge=2/lambda;

sum_v strut_v:
    pole=1, diagonal=-rho/q_diag, upper edge=0.
```

Removing one source endpoint coefficient from the lexicographically first
diagonal must break exact equivariance and one collective identity.

For diagnosis only, report the full binary64 singular spectrum and condition
number twice, using direct SVD and an independent Gram-eigenvalue route.
Their rank decisions cannot override the exact combinatorial proof.  Large
conditioning is a coordinate warning, not a physical effect.

## Outcome hierarchy

1. `FULL_SCALE_STRUT_CONTROL_FAILED`: provenance, old controls,
   nondegeneracy, exact geometry, local underdetermination, schedule
   combinatorics or corruption control fails.
2. `FULL_SCALE_STRUT_FINITE_DISAGREEMENT`: any new exact global solution
   disagrees with the disclosed coefficients.
3. `FULL_SCALE_STRUT_CURVED_CARRIER_FAILED`: finite exact controls agree but
   the accepted-background incidence, identities, equivariance or exact-rank
   proof fails.
4. `FULL_SCALE_STRUT_FINITE_CONTROLS_CORROBORATE`: all checks pass.
5. `FULL_SCALE_STRUT_NUMERICALLY_OPEN`: all exact checks pass but either
   diagnostic spectral computation is nonfinite or the two diagnostics
   disagree beyond `1e-8` relatively.

Outcome 4 corroborates the candidate at three previously unseen rational
points and freezes the accepted-background matrix, but **does not prove the
generic formula**.  Under Rule 4, a mechanically different symbolic
adversarial verifier remains mandatory before the formula is accepted and
before any action/Hessian pullback is run.

## Interpretation firewall

This protocol tests a kinematic coordinate map.  It does not test dynamics,
canonical stationarity, gauge, propagating modes, a clock, a tick, `c`, `G`,
Planck scales or particle masses.  A rank-240 result is not a count of
physical degrees of freedom.

Only the new verifier and static registry guards may be run.
