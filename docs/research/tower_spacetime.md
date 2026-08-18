# Inflation tower as a dynamical spacetime direction

Date: 2026-07-24

## Audit provenance correction

The criteria are imported unchanged from the current holographic verifier,
but the repository contains no immutable history proving that either those
criteria or the tower weights were fixed before results were inspected.
Accordingly “frozen” below means fixed within the executable rerun, not
externally preregistered.  The negative result is reproducible conditional on
the stated construction and thresholds.

## Verdict

**DERIVED negative under the frozen finite protocol:** neither registered
exponential weight produces a qualifying four-dimensional plateau at any of
the three truncations `N=8,12,16`.  The golden choice
`c_n=phi^-n` and the McKay/PF choice `c_n=2^-n` fail in both the counting and
heat estimators.  The single-floor control retains the calibrated spatial
verdict: its counting plateau has width `0.795` decade and
`d_N=3.0688`, while it has no 4D plateau.

The uniform path control produces a 4D counting plateau and a 4D heat plateau
at `N=8`, but neither survives at `N=12` or `N=16`.  Thus the registered
partial verdict, “dimension comes from the product, not from the derived
inflation weights,” is not earned in an `N`-stable form either.  The actual
verdict is:

`NO_N_STABLE_4D_PLATEAU_FINITE_SIZE_NEGATIVE`.

This is an honest negative for this operator and these truncations, not a
proof that every dynamical tower construction fails.

## Frozen construction

The spatial operator is the verified Kähler--Dirac operator
`D3=d+d*` on all `2640` cochains of the 600-cell boundary.  With form parity
`gamma`, the product operator is

`D4 = D3 tensor I + gamma tensor D_tower`.

The tower operator is the real symmetric `N x N` Jacobi matrix

`D_tower[n,n+1] = D_tower[n+1,n] = c_n`.

The tested roster was frozen at the top of
`reproducible/verify_tower_spacetime.py` before evaluation:

| label | weights | status in this audit |
|---|---|---|
| w1 | `c_n=phi^-n` | registered golden-inflation scale candidate |
| w2 | `c_n=2^-n` | registered McKay/PF scale candidate |
| w3 | `c_n=1` | flat control only |

There is a scope qualification.  `bratteli_inflation.md` derives the
Fibonacci PF ratio `phi` and McKay PF ratio `2`, but explicitly records that
a Dirac coefficient sequence `c_n` was not derived by the AF matter audit.
Accordingly w1 and w2 are the mission-registered candidates obtained from
derived tower scale ratios; their use as hopping coefficients remains a
**STRUCTURAL** dictionary choice.  The negative result does not erase the
derived PF ratios.

All plateau constants and estimators are imported unchanged from
`verify_holographic_dimension.py`: minimum width `0.50` decade, target
tolerance `0.35`, counting RMSE at most `0.08`, counting local-dimension
standard deviation at most `0.35`, heat standard deviation at most `0.35`,
and heat log-derivative at most `1.00`.  Counting uses distinct positive
levels with multiplicities; heat retains zero modes.  Both estimators use
the `D^2` factor-two convention.

## Exact product identity

Because `D3` changes form degree,

`gamma D3 + D3 gamma = 0`.

Consequently the mixed terms cancel:

`D4^2 = D3^2 tensor I + I tensor D_tower^2`.

The verifier supplies an exact symbolic finite-dimensional certificate with
a rectangular integer coboundary block and symbolic tower hops `c0,c1`.
Both anticommutation and the complete matrix difference vanish entry by
entry.  This is an algebraic identity, not a floating-point inference.

The numerical spectrum is therefore formed only from pair sums
`lambda_3+lambda_t`, with multiplied multiplicities.  No tensor-product
matrix is diagonalized.  Multiplicity closure gives `2640 N` states for
every dataset, and `Tr(D_tower^2)=2 sum_n c_n^2` is checked independently.

## Frozen plateau results

`NONE` means no contiguous interval passes every frozen gate.

| N | weights | tower count 1 | tower heat 1 | product count 4 | product heat 4 |
|---:|---|---|---|---|---|
| 8 | w1 golden | NONE | `(1.636,0.6503,0.1515,0.9201)` | NONE | NONE |
| 12 | w1 golden | NONE | `(1.642,0.6513,0.1518,0.7649)` | NONE | NONE |
| 16 | w1 golden | NONE | `(1.622,0.6541,0.1504,0.7298)` | NONE | NONE |
| 8 | w2 McKay/PF | NONE | NONE | NONE | NONE |
| 12 | w2 McKay/PF | NONE | `(0.513,0.6514,0.0346,0.6464)` | NONE | NONE |
| 16 | w2 McKay/PF | NONE | NONE | NONE | NONE |
| 8 | w3 uniform | `(1.467,0.7909,0.0737,0.0000)` | `(0.866,1.1146,0.0355,0.9544)` | `(0.620,3.7132,0.0674,0.0585)` | `(0.775,4.3352,0.0554,0.9549)` |
| 12 | w3 uniform | `(1.732,0.7879,0.0730,0.1090)` | `(1.187,1.0972,0.0493,0.9496)` | NONE | NONE |
| 16 | w3 uniform | `(1.929,0.7951,0.0743,0.1737)` | `(1.438,1.0861,0.0568,0.9807)` | NONE | NONE |

Tuple entries are `(width, dimension, residual/std, stability)`.  The
exponential tower heat intervals near `0.65` pass the formal target-1
tolerance only at some truncations; they are plainly not 1D-like in central
value.  No tower-product dataset has a qualifying 3D interval either.  The
full counting and heat arrays are computed for every row.  They can be
printed without hidden window selection using:

```text
PYTHONPATH=/tmp/science-python-deps \
python3 reproducible/verify_tower_spacetime.py --full-curves
```

## Why exponential scale does not add one Weyl dimension here

The flat path has the familiar low-energy accumulation needed to act
approximately like an extra geometric direction, although the product
plateau is too finite and unstable here.  Exponentially decreasing Jacobi
hops instead produce a strongly inhomogeneous, hierarchical spectrum.  In
this finite unrescaled direct product, that hierarchy does not supply a
stable extra power of the spectral counting law.  Calling the coordinate
“RG scale” or “hyperbolic” does not override this measured failure.

No derived-weight 4D window appeared, so the preregistered geometry-type
comparison was not activated.  In particular there is no basis here for
fitting flat `R^4` versus hyperbolic `H^4`, and no discrete Euclidean-AdS4
claim is made.  The `H4`/hyperbolic-4 naming coincidence remains irony, not
evidence.

## Index reconciliation

The conditional index-reconciliation test was also not activated by a
derived 4D window.  Algebraically, this finite odd product does not rescue
the old vertex `Box` index `9-13+1-1=-4`: that index belongs to a different
operator hierarchy, while the Kähler--Dirac factor has two harmonic modes
and Witten index zero.  No natural equality with an index of `D4` or a
boundary index has been derived.  The old `-4` remains **OPEN /
UNEXPLAINED AS SPACETIME DIMENSION**, rather than being forced into the
product construction.

## Registered scale anchors

The frozen exponent list was `{5,6,25,35}`: `a1=5`, `b1=6`, the
electroweak exponent `25`, and the neutrino exponent `35`.

Across all tower spectra and truncations, all `134` unordered ratios of
distinct positive tower eigenvalues were tested for exact
`phi^k` agreement at log tolerance `10^-10`.  There were no hits.

For w1 hopping coefficients themselves, a separation of `k` levels has the
tautological ratio `phi^-k`.  Thus separations 5 and 6 occur whenever the
truncation is long enough; 25 and 35 do not occur in the tested roster.
Every permitted integer separation occurs by the same rule, so 5 and 6 are
not distinguished features.  Selecting them would incur the full set of
level-pair comparisons and is rejected as numerology.

## Status ledger

**Strengthened / DERIVED**

- The odd-product factorization is exact.
- Pair-sum spectra, multiplicities, tower moments, and full curves are
  reproducible without giant diagonalizations.
- The static single-floor 3D verdict survives unchanged.
- Neither exponential scale candidate produces a 4D plateau at any tested
  truncation.

**Downgraded**

- The proposed exponential-scale coordinate does not yield spectral 4D for
  this finite Jacobi product.
- Even the flat-product 4D signal is only an `N=8` finite-size occurrence,
  not an `N`-stable control.
- No H4-like/AdS4 structural interpretation is supported.

**OPEN**

- A preregistered finite-size scaling limit with a geometrically justified
  rescaling of spatial slices or tower measure.
- A tower Dirac derived from an AF spectral triple rather than a registered
  PF-ratio hopping ansatz.
- A dynamical time construction independent of the scale tower.
- Reconciliation, or explicit separation, of the old vertex index `-4`.
