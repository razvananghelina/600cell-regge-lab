# Smooth full-de Rham `A4`: conformal stabilization and its limit

Date: 2026-08-12

Preregistered protocol commit: `23f0c4d`

Registered verifier:
`reproducible/verify_smooth_derham_a4_stabilization.py`

Machine-readable result:
`reproducible/smooth_derham_a4_stabilization.json`

## Headline

Two statements are simultaneously true:

> **DERIVED `A4` CONFORMAL STABILIZER.** On smooth round `S3`, the ordinary
> full-de Rham coefficient `A4` has a strictly positive normalized Hessian on
> every non-gauge scalar conformal harmonic.

> **DERIVED CONDITIONAL UV SADDLE / FINITE CUTOFF OPEN.** In a standard
> spectral action `A4` is suppressed by `Lambda^-2` relative to the unstable
> positive-weight `A2`. Consequently the complete smooth action remains a
> saddle for sufficiently high cutoff, under the usual uniform heat
> expansion. An `A2+A4` truncation can reverse the conformal sign at finite
> heat time, but that time is unselected and the omitted higher orders are
> then uncontrolled.

The targeted verifier passes `19/19`. No full-suite run was performed.

## 1. Exact smooth coefficient

For a Laplace-type operator `P=-(nabla^2+E)` on a closed manifold, the
universal integrated coefficient is

```text
A4=(1/360) integral tr[
  60 R E + 180 E^2 + 30 Omega_ij Omega_ij
  +(5R^2-2|Ric|^2+2|Riem|^2)I
].
```

This is the standard formula collected in equation (4.28) of
[Vassilevich's heat-kernel review](https://arxiv.org/abs/hep-th/0306138).
Total divergences integrate to zero here because `S3` is closed.

For the complete exterior algebra in dimension three, explicit Hodge blocks
give

```text
rank=8,
tr E=-2R,
tr E^2=2|Ric|^2,
tr Omega_ij Omega_ij=-2|Riem|^2.
```

Using `|Riem|^2=4|Ric|^2-R^2` in three dimensions reduces the complete
coefficient to

```text
A4=integral[(7/15)|Ric|^2-(1/10)R^2].
```

At the unit round metric, `Ric=2g` and `R=6`, so the density is exactly `2`
and

```text
A4(round)=4*pi^2.
```

Unlike the finite matrix moments previously stripped of the misleading
Seeley--DeWitt label, this is a genuine smooth heat-asymptotic coefficient.

## 2. Exact conformal Hessian

Use

```text
g_epsilon=(1+epsilon*f)^4 g_round,
integral f=0,
L f=lambda f,
F=integral f^2,
V=Vol(S3).
```

The verifier expands the conformal Ricci and scalar-curvature formulas
independently. The unnormalized second-order terms are

```text
[epsilon^2] integral |Ric|^2 dVol
  =(24 lambda^2-104 lambda+36)F,

[epsilon^2] integral R^2 dVol
  =(64 lambda^2-288 lambda+108)F.
```

Because `A4` scales as inverse length, remove scale with

```text
A4_hat=Vol^(1/3) A4.
```

The exact Hessian factors as

```text
delta^2 A4_hat(f,f)
 = (48/5)V^(1/3)(lambda-3)(lambda-10/9)F.
```

Consequently:

- `lambda=3` is exactly zero, as required for conformal diffeomorphisms;
- every physical scalar harmonic `l>=2`, with `lambda=l(l+2)>=8`, is strictly
  positive;
- for the explicit `f=x1^2-x2^2`, `lambda=8` and `F=V/6`, giving

  ```text
  delta^2 A4_hat=(496/9)V^(4/3)>0.
  ```

Thus the same operator does contain the right higher-derivative sign to oppose
the `A2` conformal instability.

## 3. Gaussian cutoff and the exact truncated threshold

For the entire positive Gaussian family

```text
chi_a(v)=exp(-a v^2),  a>0,
tau=a/Lambda^2,
```

the heat expansion weights are

```text
Lambda^3 a^(-3/2) A0
+Lambda a^(-1/2) A2
+Lambda^(-1) a^(1/2) A4+... .
```

Combining the exact `A2` and `A4` Hessians shows that the **truncated** scalar
mode is positive precisely when

```text
tau > 10/[(9 lambda-10)V^(2/3)].
```

This decreases with `lambda`, so `l=2`, `lambda=8`, is load-bearing. The
whole non-gauge scalar sector of the truncation is positive exactly when

```text
tau > tau_star=5/[31 V^(2/3)]
                  =0.0220828513895...  on unit S3.
```

The number is a boundary between two signs, not a prediction of heat time,
cutoff or Planck scale. The family leaves `a` arbitrary, the theory leaves
`Lambda` arbitrary, and only their ratio occurs.

## 4. Why this does not yet repair the theory

At fixed cutoff shape `a`, the relative `A4/A2` weight is

```text
a/Lambda^2 -> 0  as Lambda -> infinity.
```

The already certified `A2` has both a positive Hopf direction and a negative
non-gauge conformal direction. Under a uniform smooth heat expansion, lower
orders cannot change those two signs once `Lambda` is sufficiently large.
Therefore the complete smooth spectral action is a saddle in its
asymptotic/UV regime. This conclusion does not require knowing the numerical
onset `Lambda_0`, and no such value follows without a remainder bound.

At `tau` near or above the truncated threshold, `A6,A8,...` are not bounded by
this calculation. Reporting the positive `A2+A4` truncation as a stable full
vacuum would therefore be the same peak-versus-plateau error in a new form:
one would be reading an uncontrolled finite-order crossing as the complete
observable.

The singular finite-Regge `A4` is also not supplied by the smooth formula.
Additional skeleton terms can occur, so this result cannot be attached to the
finite `(569,0,150)` Hessian by analogy.

## 5. Status ledger

| Claim | Status |
|---|---|
| Full-exterior trace identities entering smooth `A4` | **DERIVED** |
| `A4=(7/15)|Ric|^2-(1/10)R^2` in smooth 3D | **DERIVED FROM UNIVERSAL FORMULA** |
| Normalized conformal Hessian factorization | **DERIVED** |
| `A4` is positive on every scalar `l>=2` mode | **DERIVED `A4` CONFORMAL STABILIZER** |
| Gaussian `A2+A4` threshold | **DERIVED FOR THE TRUNCATION** |
| `tau_star` is selected by the theory | **REFUTED** |
| `A2+A4` positivity proves finite-cutoff full-action stability | **REFUTED / UNDERCONTROLLED** |
| Complete smooth action is a saddle for sufficiently high cutoff | **DERIVED CONDITIONAL UV SADDLE** |
| A numerical UV onset follows without a uniform remainder | **REFUTED** |
| Singular finite-Regge `A4` has the same formula | **OPEN / NOT LICENSED** |
| Higher full heat orders stabilize at a selected finite cutoff | **OPEN** |
| Newton/Planck normalization or Lorentzian gravity follows | **OPEN** |

## 6. Physical consequence

This is neither a rescue nor a dead end. It identifies the exact structure of
the remaining problem:

> The theory's own next smooth coefficient opposes the conformal instability,
> but only with a relative strength equal to an unselected heat time. In the
> genuine high-cutoff limit the destabilizing `A2` wins.

The next admissible gate is therefore the complete heat trace as a function of
`tau`, with a target-independent mechanism selecting `tau`, or an independent
derivation of a constraint/contour. Merely choosing `tau>tau_star` after this
calculation would be fitting.

The separate inventory `gravity_hamiltonian_constraint_gap_result.md` finds no
current Hamiltonian metric constraint in the repository: the local tick acts
on fixed cochains, while the Whitney constrained systems have no physical
first-class directions at the checked levels. Consequently "a constraint may
remove the conformal mode" remains **OPEN**; it is not an available rescue of
the present smooth saddle.

## 7. Reproduction history

The first execution passed the 13 geometric identities reached before SymPy
declined to solve a parameterized but already factored linear inequality. The
generic solver was replaced by exact extraction of the sole numerator root.
Two subsequent predicate failures came from asking SymPy to infer
`lambda>=8` from a symbol declared only positive and from omitting Markdown
markers in a textual guard. Their hypotheses were made explicit; no formula,
threshold or decision criterion changed. The final targeted run passed
`19/19`.

No PDF was built.
