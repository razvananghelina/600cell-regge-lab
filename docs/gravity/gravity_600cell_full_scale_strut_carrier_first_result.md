# First result: full scale--strut carrier exact controls pass, precision remains open

Date: 2026-08-20

## Provenance

| stage | commit | status |
|---|---|---|
| exploratory disclosure and prior-art gate | `bfce559` | formula disclosed before testing |
| target-disclosed protocol | `19dd489` | three unseen rational controls frozen |
| registered implementation | `6b9eaf2` | no new representative evaluated yet |
| first artifact | `90efa45` | `18/18`, numerically open |

First artifact SHA-256:

```text
6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d
```

Only the targeted verifier was run.  No full-suite run was performed.

## Exact result

**DERIVED EXACT AT THREE NEW REPRESENTATIVES.** At

```text
(lambda,tau)=(4,7), (5,13), (7,17),
```

the complete 600-cell face equations produced `51,320` exact affine
constraints on the universal `6 x 8` local block.  A modular scan selected
48 rows, their exact rational determinant was nonzero, and the resulting
solution had zero residual on every constraint.

Direct differentiation of all twelve oriented cross diagonals gave:

| `(lambda,tau)` | `(A,B,C,D)` | mismatches / 96 | corrupted `D+1` mismatches |
|---|---|---:|---:|
| `(4,7)` | `(-44/9,116/9,-1/3,4/3)` | 0 | 12 |
| `(5,13)` | `(-121/8,185/8,-1/4,5/4)` | 0 | 12 |
| `(7,17)` | `(-181/18,325/18,-1/6,7/6)` | 0 | 12 |

These are exactly the disclosed values

```text
A = 6 - 2 tau^2/(lambda-1)^2,
B = 2 + 2 tau^2/(lambda-1)^2,
C = -1/(lambda-1),
D =  lambda/(lambda-1).
```

The controls also reconfirmed that one frustum alone imposes only

```text
A+B=8, C+D=1.
```

Thus the coefficient selection is genuinely a gluing/global-closure effect,
not a local identity hidden by notation.

## Accepted curved-slab candidate

At the accepted background

```text
lambda = 0.9999968839452782140983335686805886346459,
rho    = 0.0001040396296959628201408729245721710645837,
```

both staircase parities have the preregistered exact support census:

```text
rows:          120 support-1 + 720 support-2 + 720 support-4,
scale column:  support 24 at every vertex,
strut column:  support 13 at every vertex.
```

The 120 pole rows are a literal identity in the strut columns.  The 720
upper-edge rows contain the unsigned incidence matrix of the connected,
non-bipartite 600-cell graph, whose exact modular rank is 120.  Consequently
the complete matrix has exact rank 240 without a floating threshold.

All 24 schedule stabilizer elements preserve coefficient roles exactly.  The
homogeneous scale and strut sums reproduce their analytic derivatives with
zero 100-digit error.  Removing one endpoint coefficient breaks exact
equivariance.

**STRUCTURAL.** This is a rank-240 kinematic coordinate carrier.  It is not a
count of physical degrees of freedom and has not been passed through the
action.

## Why the first outcome is still `NUMERICALLY_OPEN`

The direct binary64 SVD gives condition numbers

```text
even  1.10394725810516514e7,
odd   1.10394725841135643e7.
```

The Gram route first forms `G^T G`.  Its smallest singular values disagree
with direct SVD by `2.85%` and `2.93%`, exceeding the frozen `1e-8` criterion.
The resulting first outcome is therefore honestly

```text
FULL_SCALE_STRUT_NUMERICALLY_OPEN.
```

There is a quantitative numerical explanation:

```text
epsilon_binary64 * kappa(G)^2 ~= 0.0271,
```

which is the same scale as the observed discrepancy.  Forming normal
equations squares the condition number, so the binary64 Gram diagnostic is
not expected to resolve the smallest singular value here.  This explanation
does not retroactively change the preregistered verdict.  A separate
high-precision audit is required.

## Post-result literature check

The technical terms learned in the run were `canonical boundary data`,
`global face closure`, `ill-conditioned Gram singular values` and `normal
equations`.

- Hoehn's [*Canonical linearized Regge Calculus: counting lattice gravitons
  with Pachner moves*](https://arxiv.org/abs/1411.5672) distinguishes lapse
  and shift variables from gauge-invariant lattice gravitons in flat
  linearized Regge evolution.  It supports the interpretation firewall: a
  kinematic carrier is not yet a graviton count.
- Dittrich and Hoehn's [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974) states that initially free data
  can become fixed by constraints from subsequent moves.  This is the known
  general mechanism closest to the observed local-to-global selection.
- The [LAPACK Users' Guide](https://www.netlib.org/lapack/lug/) and its SVD
  drivers treat singular values and condition estimation directly; the
  smallest singular value should not be certified here through unscaled
  binary64 normal equations.

The post-result search still found no primary source containing the explicit
four-coefficient nonhomogeneous tetrahedral response above.  Search absence
is not a novelty proof; external novelty remains **OPEN**.

## Status ledger

| Claim | Status |
|---|---|
| One cell leaves two endpoint freedoms | **DERIVED EXACT** |
| Three unseen global rational systems select the disclosed coefficients | **DERIVED EXACT FINITE CONTROLS** |
| Generic symbolic coefficient theorem | **OPEN pending adversarial replication** |
| Curved-slab support, equivariance and rank 240 | **DERIVED EXACT for the disclosed candidate** |
| Binary64 Gram diagnosis of the weakest singular values | **OPEN / inadequate precision** |
| Gauge, constraint or physical content | **NOT EVALUATED** |
| Tick, `c`, `G`, Planck scale, graviton or particle mass | **NOT EVALUATED** |
| External novelty | **OPEN** |

## Next gates

1. Preregister and run a high-precision conditioning audit that compares two
   direct LAPACK SVD drivers with an arbitrary-precision Gram spectrum and a
   known ill-conditioned control.
2. Independently prove or refute the generic coefficient formula using a
   symbolic two-cell construction that does not call the complete global
   solver.
3. Only if both gates pass may the carrier be frozen and pulled through the
   action/Hessian or strong-equation map.

