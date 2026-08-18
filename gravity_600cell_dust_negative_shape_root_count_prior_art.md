# Prior-art gate: certified root count for the negative-shape recurrence

Date: 2026-08-18

## Exact object and hypotheses

The inherited calculation selects, without a multiplier target, two
`15`-dimensional negative-stiffness shape carriers in each of the two real
parities.  It supplies four derivative audits per parity and sector.  On each
carrier the frozen centered recurrence is

```text
q2 - 2 q1 + q0 + Gamma (q2 - q0) + Omega q1 = 0,
```

and therefore the quadratic matrix polynomial to be audited is

```text
Q(z) = (I + Gamma) z^2 + (-2 I + Omega) z + (I - Gamma).
```

The carrier, the midpoint matrices and their operator-norm error balls are
inherited from the committed stiffness and centered-Jacobi calculations.  No
continuum harmonic, desired root count, speed, particle, Planck or refinement
target is an input.

The statement tested here is deliberately local.  It concerns the roots of a
single frozen two-step recurrence.  It assumes neither time-translation
invariance of later slabs nor a stationary background.  Even a resolved root
outside the unit circle would establish local frozen hyperbolicity, not a
long-time Lyapunov exponent for the nonlinear cosmology.

## Structural audit before choosing the method

The tempting reciprocal-pair shortcut is unavailable.  On the complete
selected source sectors the centered coefficients `M`, `N`, and `V` are each
Hermitian-consistent, while `N` is not anti-Hermitian:

```text
relative Hermitian defect of N       7.19e-14
relative anti-Hermitian defect of N  2.00e+00
relative *-palindromic defect        3.34e-01
```

Here the last line measures

```text
||(M+N)^* - (M-N)|| / max(||M+N||, ||M-N||).
```

Consequently the hypotheses that would make the leading and trailing
coefficients adjoints are false on the inherited operator.  A reciprocal
spectrum must not be assumed.  The restricted midpoint may display additional
numerical patterns, but they are not structural input to the test.

## Primary literature

**KNOWN.**  Gohberg and Sigal generalized the argument principle and
Rouche's theorem to analytic operator-valued functions:
<https://doi.org/10.1070/SM1971v013n04ABEH003702>.

For finite matrix polynomials, the applicable statement is especially simple.
If `P(z)` is nonsingular on a closed contour and

```text
||P(z)^(-1) DeltaP(z)|| < 1
```

there, then `det(P+DeltaP)` and `det(P)` have the same number of enclosed
zeros, with algebraic multiplicity.  Noferini, Sharify and Tisseur state this
finite-dimensional form explicitly as Theorem 2.1 and use it to localize
matrix-polynomial eigenvalues:
<https://doi.org/10.1137/14096637X>.

**KNOWN.**  Quadratic matrix polynomials, their linearizations, conditioning
and structured spectra are reviewed by Tisseur and Meerbergen:
<https://doi.org/10.1137/S0036144500381988>.  A backward-stable midpoint
linearization does not by itself propagate the inherited coefficient balls.

**KNOWN.**  Contour-integral algorithms compute eigenvalues of nonlinear
matrix functions within a contour, but ordinary numerical quadrature is not
automatically a certificate for all coefficient perturbations.  See Brennan,
Embree and Gugercin:
<https://doi.org/10.1137/20M1389303>.

No located primary source studies this particular `600`-cell Regge carrier or
supplies its root count.  External novelty remains **OPEN**.

## CONTROL, OPEN, and proposed difference

- **CONTROL:** the free recurrence `q2-2q1+q0=0` has only the double unit root,
  although its companion matrix has largest singular value `1+sqrt(2)`.
  Singular amplification alone is therefore not an instability certificate.
- **DERIVED UPSTREAM:** all `32` Gamma/Omega invariance tests are
  invariant-consistent, and every restricted forward coefficient `I+Gamma`
  is regular-resolved.
- **PATTERN UPSTREAM:** midpoint linearizations show `15` roots inside and
  `15` roots outside the unit circle in every selected cell.
- **OPEN:** whether that count survives every matrix in the inherited
  coefficient balls.
- **OPEN:** whether any resolved frozen hyperbolicity persists for independently
  solved later slabs or converges under refinement.

The new falsifiable step is a continuous-contour Rouche certificate.  For
`z=exp(i theta)`, the correlated perturbation is

```text
DeltaQ(z) = DeltaGamma (z^2 - 1) + DeltaOmega z,
```

so its pointwise operator-norm bound is

```text
delta(theta) = epsilon_Gamma |z^2-1| + epsilon_Omega.
```

The verifier will cover the complete unit circle by an adaptive deterministic
interval subdivision and lower-bound

```text
sigma_min(Q(z)) - safety * delta(theta)
```

on every interval using an analytic Lipschitz bound.  Failure to cover the
circle is an honest `OPEN`, not permission to quote the midpoint roots.  The
root count of the midpoint will be cross-checked by a generalized eigenvalue
linearization and by the winding of `det Q` before it is transferred through
Rouche's theorem.

