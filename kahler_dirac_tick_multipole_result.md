# Leading angular multipole of the Kähler--Dirac tick

Date: 2026-08-11

Preregistration commit: `82ca6c9`

Targeted verifier:
`reproducible/verify_kahler_dirac_tick_vertex_isotropy.py`

Targeted result: **24/24 PASS**.  The full suite was not run, by explicit
request.

## Question and provenance

After the original vertex-isotropy calculation, an LLM-generated observation
reported unusually small values denoted (A_6) at even ticks.  It supplied no
definition.  This is therefore not independent evidence.  It is useful only
as a hypothesis generator, and the audit was explicitly preregistered as
confirmatory rather than blind.

The interval in that report turns out to identify the convention: its (A_6)
is the **squared harmonic power** called (A_6^2) in the preregistered
coordinate-free normalization, not the harmonic amplitude (A_6).

## Why degree six is not fitted

For the degree-ℓ spherical harmonics, averaging the (SO(3)) character over
the icosahedral rotation group gives the trivial-representation
multiplicities

\[
(m_0,\ldots,m_{12})=(1,0,0,0,0,0,1,0,0,0,1,0,1).
\]

Thus degree six is the first nonconstant angular harmonic permitted by the
vertex stabilizer.  This is **STRUCTURAL**, not a degree selected from a scan.
The verifier independently checks the character average, the uniform
12-neighbour icosahedral shell, and a deliberately distorted-shell control.

The previously derived covariance isotropy tests only degrees one and two.
A nonzero degree-six multipole is therefore compatible with an exactly
isotropic covariance: it is the leading lattice anisotropy that the covariance
cannot see.

## Frozen definitions

For tangent directions (u_i), probabilities (p_i), radii (r_i), and
Legendre polynomial (P_\ell), the normalized power is

\[
A_\ell(q)^2=
\frac{\sum_{i,j}q_iq_jP_\ell(u_i\cdot u_j)}{(\sum_iq_i)^2}.
\]

Three variants were frozen before computation:

1. conditional angular: (q_i=p_i), after removing the origin;
2. unconditional angular: the same numerator without division by the moving
   probability;
3. solid-harmonic radial: (q_i=p_i r_i^\ell).

The third is essential: it tests the actual degree-ℓ spatial multipole
rather than discarding all radial information.

## Result

| tick | support | conditional (A_6^2) | unconditional (A_6^2) | radial (A_6^2) |
|---:|---:|---:|---:|---:|
| 1 | 12 | 0.440000 | 0.440000 | 0.440000 |
| 2 | 43 | 0.00386719 | 0.000927738 | 0.273802 |
| 3 | 114 | 0.416357 | 0.416357 | 0.0538388 |
| 4 | 255 | 0.000391270 | 0.000390365 | 0.202774 |
| 5 | 336 | 0.324272 | 0.324272 | 0.0246874 |
| 6 | 687 | 0.00958036 | 0.00357824 | 0.0332835 |
| 7 | 588 | 0.254010 | 0.254010 | 0.00799160 |
| 8 | 1109 | 0.0199053 | 0.00122998 | 0.00136736 |

The reported (4\times10^{-4}) to (2\times10^{-2}) interval is precisely
the conditional even-tick column, up to rounding.  Thus the numerical clue was
real and has been reproduced from the theory's own tick.

However, the stronger wording does not survive:

- the conditional power suppression relative to the preceding odd tick is
  about factors 114, 1064, 34 and 13 at ticks 2, 4, 6 and 8;
- hence a three-order suppression occurs near tick 4 only, not throughout the
  even subsequence;
- the even sequence is nonmonotone and is already rising from tick 4 through
  tick 8;
- the solid-harmonic radial even and odd ranges overlap strongly.  Its
  preregistered separation ratio in amplitude is 5.85, in the wrong direction
  for a universal even-tick suppression.

The angular variants separate all four even values from all four odd values,
but under the frozen amplitude criterion their ratios are 0.280 (conditional)
and 0.119 (unconditional), so both are labelled **PATTERN: weak parity
separation**, not a derived strong effect.  The radial variant is a **DERIVED
NEGATIVE** for even/odd separation.

All degrees one through five vanish at all eight ticks and in all three
weightings at the numerical precision floor (maximum amplitude
(1.71\times10^{-10})).  To make that cancellation reliable despite the
ten-decimal storage in `commons.cell600`, the verifier restores the unique
ℚ(√5) coordinate alphabet and evaluates the Legendre recurrence in
extended precision.  The original double-precision calculation was correctly
rejected by the preregistered (10^{-9}) gate after square-root amplification
of roundoff.

## Interpretation ledger

- **DERIVED / STRUCTURAL:** ℓ=6 is the first icosahedrally permitted
  anisotropy, while ℓ=1,…,5 vanish.
- **DERIVED NUMERICAL (fixed complex, ticks 1…8):** the conditional angular
  powers reproduce the reported even-tick interval exactly.
- **PATTERN:** the angular distribution has a pronounced even/odd
  cancellation, strongest at tick 4.
- **DERIVED NEGATIVE:** the cancellation is not robust to the preregistered
  radial weighting, and the claimed three-order separation is not uniform.
- **OPEN:** whether any even-tick suppression persists at longer times or
  under a geometrically selected refinement family.
- **OPEN:** any relation to Lorentz invariance, continuum rotational
  invariance, a physical clock, or suppression of observable dispersion.

## Physical verdict

There is a genuine finite-complex kinematic pattern worth retaining: the
bipartite Kähler--Dirac tick alternates between cochain parities, and the
leading allowed **angular** lattice harmonic cancels unusually well on the
even parity during the first few ticks.  But it is not yet a physical
suppression theorem.  Once radius is included, the supposed parity effect is
absent.  The honest headline is therefore:

> **A real ℓ=6 angular cancellation was found, but the fixed 600-cell does
> not yet show robust suppression of leading spatial anisotropy.**

