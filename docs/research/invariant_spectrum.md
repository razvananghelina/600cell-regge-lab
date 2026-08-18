# Invariant Kähler--Dirac spectrum of the 600-cell

Date: 2026-07-24

## Result

The carrier spaces remain gauge.  The spectrum below is obtained without
choosing a carrier: a central class sum labels the nine `2I` isotypic
components and is scalar on each irrep.  The computation closes at dimension
2640 and has a two-dimensional kernel.

The complete row-level canonical table is emitted between
`SPECTRUM_ROWS_BEGIN/END` by
`reproducible/verify_invariant_spectrum.py`.  Each row is

`D^2 eigenvalue, multiplicity, irrep, form degree`.

There are 52 distinct positive eigenvalues and 158 nonzero
`(eigenvalue,irrep,degree)` rows.  This executable table, rather than a choice
of eigenvectors, is the canonical `U(m)`-invariant object.  Equal eigenvalues
on several rows are genuine coincidences between the listed sectors/degrees;
they are not merged because the form-degree content would be lost.

Representative exact values (and the numerical values used to identify the
full table) are:

| `D^2` eigenvalue | decimal | sectors/degrees at the bottom of the spectrum |
|---|---:|---|
| `0` | 0 | `rho0:C0(1),C3(1)` |
| `phi^-4=(7-3 sqrt(5))/2` | 0.145898033750 | `rho1:C2(4),C3(4)` |
| `phi^-2=(3-sqrt(5))/2` | 0.381966011250 | `rho2:C2(9),C3(9)` |
| `5-2 sqrt(5)` | 0.527864045000 | `rho0:C1,C2(3 each); rho2:C1,C2(3 each)` |
| `(5-sqrt(13))/2` | 0.697224362268 | `rho3:C2(16),C3(16)` |
| numerical algebraic root | 1.07457708222 | `rho4:C2(25),C3(25)` |
| `1+phi^-4` | 1.14589803375 | `rho1,rho3:C1,C2(8 each)` |

The remaining 151 rows, through the maximum
`6+6 phi = 15.7082039325`, are printed by the verifier.  Values recognized
exactly include integers, elements of `Q(sqrt(5))`,
`(5±sqrt(13))/2`, `(11±sqrt(13))/2`, and `5±sqrt(2)`.
Values not certified in such a radical form are deliberately
reported as numerical algebraic eigenvalues, not guessed exact expressions.

Numerical diagonalization is symmetric double precision (observed moment
residual below `1e-11`).  The zeta and theta values below inherit about 13
reliable decimal digits.  The requested 50-digit eigenvalue evaluation was
not used for the unresolved algebraic roots; those numbers must therefore not
be treated as 50-digit claims.  Exact sparse integer identities and exact
moments supply independent certification of the principal anchors.

## Theorem-level anchors

- **DERIVED:** `dim ker D=2`, with `rho0:C0(1)` and `rho0:C3(1)`.
- **DERIVED:** the Witten index is
  `120-720+1200-600=0`.
- **DERIVED:** every positive `D^2` eigenspace has equal total even/odd
  dimension.  Hence `D` has eigenvalues `±sqrt(lambda)`, each with half the
  total `D^2` multiplicity.
- **DERIVED:** all rows in the Frobenius--Schur-negative sectors
  `rho1,rho3,rho5,rho7` have even multiplicity.  This is the Kramers check.
- **DERIVED:** the first excited SUSY multiplet is
  `lambda=phi^-4`, in `rho1`, with four states in `C2` and four in `C3`.
  Thus `D=±phi^-2`, multiplicity four for each sign.

## Frozen-target confrontation

The registry was committed to the verifier before the baseline and before
spectral computation.  The comparison universe contains all
`binom(52,2)=1326` unordered ratios of distinct positive eigenvalues.  A
half-integer hit requires
`|log_phi(ratio)-n/2|<1e-10`.

Ten ratios pass:

`phi^2` occurs eight times,
`phi^4=(6+6phi)/(12-6phi)` once, and
`phi^6=(5+2sqrt(5))/(5-2sqrt(5))` once.

These are exact `Q(sqrt(5))` identities after recognition.  None of the nine
registered mass exponents `(0,5,3,11,16,19,26,17,11)` occurs between distinct
levels.  Exponent zero would only compare a level with itself and was excluded
by preregistration of the comparison universe.  **DERIVED negative:** the
nine-level mass ladder is not the Kähler--Dirac eigenvalue-ratio spectrum.

The numerically resolved full gap is

`gap(D^2)=0.145898033750315...`, numerically matching `phi^-4`,

whereas the registered founding value is

`lambda_1=1/(2 phi^2)=0.190983005625053...`.

If that recognition is accepted, their ratio is
`2/phi^2=3-sqrt(5)`.  The verifier supplies no exact characteristic-polynomial
or eigenvector certificate for the gap, so `gap=phi^-4` is **NUMERICALLY
IDENTIFIED**, not an exact-arithmetic theorem.  The complex therefore refines but
does not reproduce that registered gap.  The `phi^2`, `phi^4`, and `phi^6`
ratio hits are registered-family matches, but selecting one after 1326
comparisons would be invalid; their full multiplicity is reported.

## Spectral functions and action

With the two zero modes omitted from zeta:

| quantity | value |
|---|---:|
| `zeta_D2(3/4)` | 1003.86581055511 |
| `zeta_D2(1)` | 811.356982259231 |
| `zeta_D2(2)` | 852.378548683763 |
| `Tr exp(-5 D^2)` (zero modes included) | 10.7515022667804 |
| `Tr exp(-phi D^2)` (zero modes included) | 65.1892773657095 |

No registered equality is found at these special points.

The same operator reproduces the spectral-action coefficients exactly in
their finite-moment meaning:

`c0=Tr(1)=2640`,

`c1=Tr(D^2)=14880`,

`c2=(1/2)Tr(D^4)=55920`.

Division by 240 gives `(11,62,233)`.  Thus `c1` and `c2` do arise from this
Kähler--Dirac operator; they are not imported from a different graph
operator.  These are finite spectral moments, not continuum asymptotic
Seeley--DeWitt coefficients.

## Holographic lens

Applying `HOLOGRAPHIC_RG_PROTOCOL.md` to all 2638 positive modes gives:

- counting window near 3D: `d=2.9614`, RMSE `0.0643` (fails its RMSE gate);
- counting window near 4D: `d=3.9951`, RMSE `0.0158` (passes);
- heat window near 3D: mean `2.9934`, standard deviation `0.1635` (passes);
- heat window near 4D: mean `4.0165`, standard deviation `0.2220` (fails);
- most stable heat window: `d=3.4806±0.0362`.

The protocol verdict is **INCONCLUSIVE**, not a dimension flow.  The apparent
4D counting window and 3D heat window do not overlap as a stable common
plateau.  Because the spectrum is finite with only 52 distinct levels,
`log rho(lambda)` is a shell-degeneracy statistic, not an asymptotic
Cardy law.  Fits to 2D/3D/4D CFT exponents are therefore cutoff- and
bin-dependent; no CFT-like density claim survives.

## Weird lenses

### Ihara

For the 12-regular 600-cell graph (`n=120,m=720`), Bass's exact formula is

`Z_I(u)^-1=(1-u^2)^600 product_j(1-a_j u+11u^2)`,

where `a_j=12-lambda_j(Delta_0)` with graph multiplicity.  The reciprocal-root
pairing in every quadratic is the regular-graph functional equation after
completion.  Comparing these poles with all 52 `D^2` levels yields no
registered equality; the shared vertex eigenvalues are tautological because
`Delta_0` is one block of `D^2`.

### Theta and modular sweep

The two registered theta evaluations are tabulated above.  Their decimals do
not match a frozen target.  The leading multiplicities do not factor as a
standard eta quotient or a recognizable McKay--Thompson series.  Any proposed
fit after inspecting 52 shells would be an **UNREGISTERED PATTERN** with a
large, undefined modular-form search space; none is claimed.

### SUSY quantum mechanics

`D` is the supercharge, `D^2` the Hamiltonian, and form parity the fermion
number.  There are two supersymmetric ground states of opposite parity and
zero Witten index.  Every excited energy is paired.  The first excitation is
the exact quaternionic `rho1` multiplet described above.

## Synthesis

**DERIVED:** the two harmonics, SUSY pairing, Witten index, Kramers parity,
and all three finite spectral moments.  The row spectrum and first multiplet
are double-precision numerical certificates; exact radical labels require
separate algebraic certificates.

**REGISTERED MATCH:** exact golden powers occur in ten of 1326 ratios, with
the full look-elsewhere count stated.  The spectral-action triple is exactly
reproduced.

**NUMERICAL EXHAUSTIVE NEGATIVE at tolerance `1e-10`:** none of the
registered nontrivial mass exponents occurs among all 1326 ratios of the 52
numerically clustered levels.  This is exhaustive for that numerical table,
not an exact-algebraic proof for unresolved roots.  The registered gap is not
the full-complex gap; zeta/theta special values give no target
identity; holographic tests are mixed; no Cardy or modular structure is
supported.

**OPEN:** exact minimal polynomials and 50+-digit values for every unresolved
algebraic root, plus a preregistered modular-form family if that direction is
ever revisited.
