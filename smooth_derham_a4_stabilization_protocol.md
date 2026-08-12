# Preregistered protocol: smooth de Rham `A4` and conformal stabilization

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- NO REGISTERED `A4` RESULT**

## 1. Provenance

This is a post-recognition hostile protocol. Before writing it, hand algebra
suggested that the smooth ordinary full-de Rham `A4` Hessian factors by
`(lambda-3)` and is positive on scalar harmonics with `lambda>=8`. That
reconnaissance is not evidence. The exact trace reduction, conformal expansion
and cutoff conclusion below must be reproduced mechanically.

No finite-Regge `A4`, measured target or desired physical scale has been
inspected or fitted.

## 2. Complete scope

1. The geometry is the **smooth closed unit round three-sphere**.
2. The operator is the ordinary, ungraded de Rham Laplacian on all forms.
3. Heat coefficients use

   ```text
   Tr exp(-t D^2)
    ~(4*pi*t)^(-3/2)[A0+t*A2+t^2*A4+...].
   ```

4. The universal integrated closed-manifold Laplace-type formula is used with
   `P=-(nabla^2+E)`:

   ```text
   A4=(1/360) integral tr[
        60 R E + 180 E^2 + 30 Omega_ij Omega_ij
        +(5R^2-2|Ric|^2+2|Riem|^2) I
      ].
   ```

   Total divergences vanish because the manifold is closed. This external
   input is equation (4.28) of Vassilevich's standard heat-kernel review,
   <https://arxiv.org/abs/hep-th/0306138>.
5. Scale is removed with

   ```text
   A4_hat=Vol^(1/3) A4,
   ```

   because `A4` scales as inverse length in dimension three.
6. Conformal paths and positive Laplacian conventions are exactly those of
   protocol `b366a72`:

   ```text
   g_epsilon=(1+epsilon*f)^4 g_round,
   integral f=0,
   L f=lambda f,
   integral |grad f|^2=lambda integral f^2.
   ```

7. Scale (`l=0`) is excluded by the mean-zero condition; `l=1`, `lambda=3`
   is the conformal diffeomorphism zero sector; physical scalar tests begin at
   `l=2`, `lambda=8`.

The result will not be transferred to the singular finite-Regge geometry:
its next heat coefficients can contain different skeleton contributions.

## 3. Frozen trace reduction

For the complete exterior representation in three dimensions the verifier
must derive, not assume,

```text
rank=8,
tr E=-2R,
tr E^2=2|Ric|^2,
tr Omega_ij Omega_ij=-2|Riem|^2.
```

Substitution in the universal formula and the three-dimensional identity

```text
|Riem|^2=4|Ric|^2-R^2
```

must give

```text
A4=integral[(7/15)|Ric|^2-(1/10)R^2].
```

At unit round `S3`, where `Ric=2g` and `R=6`, the local density must be `2`
and `A4=4*pi^2`.

## 4. Frozen conformal calculation

For `u=1+epsilon*f`, use the exact three-dimensional transformations

```text
Ric_(u^4g)
 =Ric-2u^-1 Hess(u)+6u^-2 du tensor du
    +(2u^-1 L u-2u^-2|du|^2)g,

R_(u^4g)=u^-5(8 L u+6u),
dVol_(u^4g)=u^6 dVol_g.
```

The verifier must independently reduce all integrated Hessian terms using

```text
integral |Hess f|^2=lambda(lambda-2) F,
F=integral f^2,
```

and check the unnormalized second-order coefficients

```text
[epsilon^2] integral |Ric|^2 dVol
  =(24 lambda^2-104 lambda+36)F,

[epsilon^2] integral R^2 dVol
  =(64 lambda^2-288 lambda+108)F.
```

Including the `Vol^(1/3)` normalization must yield

```text
delta^2 A4_hat(f,f)
 = (48/5) V^(1/3)(lambda-3)(lambda-10/9) F.
```

Mandatory controls:

- `lambda=3` is exactly zero;
- `lambda=8` is strictly positive;
- every scalar harmonic `l>=2` is positive;
- the formula must be recovered independently for the explicit
  `f=x1^2-x2^2` using `F=V/6` and `lambda=8`.

## 5. Frozen cutoff comparison

Use first the positive Gaussian cutoff, with no fitted parameter value,

```text
chi_a(v)=exp(-a v^2), a>0,
tau=a/Lambda^2.
```

Its exact expansion weights are

```text
Lambda^3 a^(-3/2) A0
+Lambda a^(-1/2) A2
+Lambda^(-1) a^(1/2) A4+... .
```

Combining the already certified

```text
delta^2 A2_hat
 =-(32/3)V^(-1/3)(lambda-3)F
```

with the frozen `A4` prediction gives a truncated conformal threshold

```text
tau > 10 / [(9 lambda-10)V^(2/3)].
```

The hardest mode is `lambda=8`, hence the complete non-gauge scalar sector of
the `A2+A4` truncation is positive exactly when

```text
tau > tau_star=5/[31 V^(2/3)].
```

This threshold is **not** a selected cutoff and is **not** a conclusion about
the full heat trace at finite `tau`; omitted `A6,A8,...` are then not bounded.

For every fixed admissible positive cutoff with positive `A2` moment and a
uniform smooth heat expansion, the `Lambda*A2` Hessian dominates all lower
orders as `Lambda -> infinity`. The explicit positive Hopf and negative
conformal directions must therefore persist as a saddle for sufficiently
large cutoff. No numerical value of “sufficiently large” may be invented
without a remainder bound.

## 6. Preregistered decisions

- **DERIVED `A4` CONFORMAL STABILIZER:** the trace reduction and conformal
  Hessian are exact and positive for every `l>=2` scalar harmonic.
- **REFUTED `A4` STABILIZER:** the exact Hessian is nonpositive for some
  `l>=2` mode.
- **DERIVED UV SADDLE / FINITE CUTOFF OPEN:** `A4` is positive but suppressed
  by `Lambda^-2` relative to `A2`, so the full spectral action remains a
  saddle asymptotically; the truncated finite threshold depends on unselected
  `tau` and cannot prove full-action stability.
- **OPEN/FAILED:** trace conventions, conformal identities or sign controls
  do not close exactly.

Even success does not derive a physical cutoff, Newton's constant, Planck
units, Lorentzian dynamics or a graviton. Only the targeted verifier and a
static registry audit will be run; no full suite and no PDF build.
