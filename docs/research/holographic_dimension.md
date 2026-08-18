# Calibrated holographic-dimension audit

Date: 2026-07-24

**2026-08-09 reconciliation:** see `dimension_reconciliation.md`.  The
independently reconstructed Kaehler--Dirac heat flow has global maximum
`3.295663771`, so it never reaches four.  Its target-free 2% shoulder is too
short for this note's half-decade heat-plateau gate.  The old `-4` is an
alternating nullity of a different, non-cochain `Box_p` hierarchy and its
spacetime reading is rejected.

## Verdict

**DERIVED negative under the stated fixed criteria:** the earlier rolling
`d_N=3.9951` window is not a plateau under the frozen, full-curve criterion.
The weighted 600-cell counting curve has no accepted 4D plateau.  It has one
accepted 3D interval, of width `0.795` decade, with fitted
`d_N=3.0688`, log-count residual `0.0659`, and local-dimension standard
deviation `0.2303`.

The heat curve has no accepted 3D or 4D plateau.  The 4D control also has no
accepted 4D plateau in either estimator.  Consequently the registered
holographic decision rule fails at conditions (i) and (iii).  The proper
classification is:

**ARTIFACT / INCONCLUSIVE-CALIBRATION:** the particular old 4D counting
window was a window-selection artifact, so there is no calibrated
holographic 4D-counting anomaly.  The small `T^4` control is itself too
cutoff-dominated to validate a positive plateau claim at the frozen
half-decade standard.  This does not turn the failed test into evidence for
4D.

The static complex is a triangulated `S^3`, and this audit supplies no
spectral derivation of a fourth static dimension.  The theory's registered
`d_ST=4` from the vertex spectral index `-4` must therefore be re-examined.
If four-dimensional spacetime is retained, its fourth direction must enter
dynamically—for example as RG scale, an inflation/Bratteli direction, or
physical time—not as a Weyl dimension already established by this finite
static Kähler--Dirac spectrum.  This is a precise tension, not deletion of
the vertex-index statement.

## Conventions

The spectrum is that of `D^2`.  For a first-order operator `D`, continuum
counting is `N_D(mu) ~ mu^d`.  Setting `Lambda=mu^2` gives

`N_D2(Lambda) ~ Lambda^(d/2)`.

Accordingly the two dimension estimators used here are

`d_N(Lambda) = 2 d log N(Lambda) / d log Lambda`,

`d_s(t) = -2 d log Tr exp(-t D^2) / d log t`.

The heat derivative is evaluated analytically, with the harmonic modes in
the heat-trace denominator:

`d_s(t)=2t sum(lambda exp(-t lambda))/(b+sum(exp(-t lambda)))`.

Zero modes are omitted from positive-spectrum counting and retained in the
heat trace exactly as its definition requires.  With these conventions,
“3D” has exactly the same meaning in both curves.  Omitting the factor two
would instead report the exponent of the `D^2` counting law, not spacetime
dimension.

## Stated criteria and controls; provenance limit

The repository has no version-control history or external timestamped
registry proving that these thresholds preceded spectral inspection.
`HOLOGRAPHIC_RG_PROTOCOL.md` predates this calculation, but it does not contain
the numerical thresholds below.  Thus “frozen” is a code-level immutability
claim for reruns, not a verifiable preregistration claim.  The numerical
result remains reproducible under these stated criteria, but its protection
against post-hoc tuning is **UNVERIFIED**.

These values occur at the top of
`reproducible/verify_holographic_dimension.py`, before any spectral
construction:

| item | frozen value |
|---|---:|
| minimum contiguous width | `0.50` decade |
| target tolerance | `|d-d_target| <= 0.35` |
| counting log-fit RMSE | `<= 0.08` |
| counting local-`d_N` standard deviation | `<= 0.35` |
| heat `d_s` standard deviation | `<= 0.35` |
| heat `max |d d_s/d log10(t)|` | `<= 1.00` |
| heat grid | 241 logarithmic points |
| local counting regression | centered 7 distinct levels |
| heat range | `0.01/lambda_max` through `100/lambda_min` |

All contiguous candidate intervals at least half a decade wide are tested.
An accepted result reports the widest qualifying interval and its residual
and stability statistic.  `NONE` means that every candidate failed at least
one frozen condition; it does not mean that a favorable shorter segment
could not be selected by eye.

The controls are the 5-cell, 16-cell, and triangulated 24-cell boundaries;
the periodic `4^3` Freudenthal triangulation of `T^3`; and the genuinely
four-dimensional periodic `3^4` Freudenthal triangulation of `T^4`.
Translation/Fourier blocks compute the torus spectra without constructing a
`12150 x 12150` dense matrix.  Kernel checks give `b(T^3)=8` and
`b(T^4)=16`, as required.

## Plateau results

The entries below are for multiplicity-weighted spectra.  Values in
parentheses are `(width, dimension, residual, stability)`, where residual is
log-fit RMSE for counting and standard deviation for heat.

| complex | states | positive levels | count 3 | count 4 | heat 3 | heat 4 |
|---|---:|---:|---|---|---|---|
| 600-cell boundary | 2640 | 52 | `(0.795,3.0688,0.0659,0.2303)` | NONE | NONE | NONE |
| 5-cell boundary | 30 | 1 | NONE | NONE | NONE | NONE |
| 16-cell boundary | 80 | 4 | NONE | NONE | NONE | NONE |
| triangulated 24-cell boundary | 432 | 210 | NONE | NONE | NONE | NONE |
| Freudenthal `T^3`, `n=4` | 1664 | 35 | NONE | NONE | NONE | NONE |
| Freudenthal `T^4`, `n=3` | 12150 | 173 | NONE | NONE | NONE | NONE |

The small S3 controls correctly do not manufacture a 4D plateau, but mostly
have insufficient dynamic range to pass a positive 3D plateau.  More
importantly, the small 4D control does not pass both positive 4D gates.
Thus this finite-size calibration can reject the old selected window but
cannot validate the plateau method as a positive dimension detector at
these lattice sizes.

## Full curves

No fitting window is hidden.  The verifier prints the complete tables for
every complex between:

`COUNTING_CURVE_BEGIN lambda N d_N` / `COUNTING_CURVE_END`

and

`HEAT_CURVE_BEGIN t d_s` / `HEAT_CURVE_END`.

Run:

```text
PYTHONPATH=/tmp/science-python-deps python3 reproducible/verify_holographic_dimension.py
```

This is the canonical full-curve table: 52 counting rows and 241 heat rows
for the 600-cell, with the corresponding complete rows for all five
controls.  The output also reports full-curve extrema and every interior
heat maximum.  The tables are executable rather than rounded plots, so all
candidate widths and residuals can be recomputed from the printed values.

For reference, the full 600-cell counting table begins:

| `Lambda` | `N` | local `d_N` |
|---:|---:|---:|
| 0.145898034 | 8 | 2.663746 |
| 0.381966011 | 26 | 2.731302 |
| 0.527864045 | 38 | 2.830114 |
| 0.697224362 | 70 | 2.897383 |
| 1.074577082 | 120 | 3.166295 |
| 1.145898034 | 152 | 3.232537 |
| 1.481801301 | 224 | 2.972401 |
| 1.763932023 | 272 | 2.895896 |

and ends:

| `Lambda` | `N` | local `d_N` |
|---:|---:|---:|
| 9.472135955 | 2458 | 0.334111 |
| 12 | 2508 | 0.302125 |
| 14 | 2580 | 0.282566 |
| 14.472135955 | 2598 | 0.291090 |
| 15 | 2630 | 0.390875 |
| 15.708203932 | 2638 | 0.406612 |

The finite-spectrum endpoint behavior is explicit: because harmonic modes
are retained, `d_s` tends to zero at both ends.  Its interior maxima are
`3.29566` for the 600-cell, `2.68715`, `3.11996`, and `3.05051` for the
regular S3 controls, `2.98985` for `T^3`, and `3.79843` for `T^4`.  A maximum
is not a plateau; none has the frozen half-decade width and stability.

## Degeneracy specificity

The same analysis was repeated after replacing every positive-shell
multiplicity by one.  For the 600-cell, the stripped spectrum has no
accepted 3D or 4D plateau in either estimator.  The stripped controls also
have no accepted plateaus.

This establishes a limited **DERIVED** statement: multiplicity weighting is
necessary for the sole accepted 600-cell 3D counting interval.  It does not
establish that `2I` degeneracy mimics an extra dimension, because the
registered 4D anomaly does not survive in the full weighted spectrum.
Attributing an extra dimension to multiplicities would therefore be a
rejected post-hoc **PATTERN**, not a derived mechanism.

## Status ledger

**Strengthened / DERIVED**

- Dimension conventions for `D` versus `D^2` are now consistent.
- Full weighted and degeneracy-stripped curves are reproducible.
- The old selected 4D counting fit does not satisfy the frozen plateau test.
- The 600-cell has a qualifying 3D counting interval under that same test.
- All control cochain dimensions and torus Betti kernels close.

**Downgraded**

- The former 4D counting “signal” is downgraded from INCONCLUSIVE tension to
  a rejected window claim: no calibrated holographic anomaly is present.
- The vertex-index `d_ST=4` is not a dimension derived from this static
  Kähler--Dirac spectrum.

**OPEN**

- Positive calibration on larger `T^3` and `T^4` triangulations with at least
  several decades of spectral range.
- A finite-size scaling sequence, frozen before evaluation, rather than a
  single-complex plateau test.
- A dynamical derivation of the fourth direction from RG, tower, or time
  evolution.
- Reconciliation—or explicit separation—of the vertex spectral index `-4`
  from the static cochain Weyl/heat dimensions.
