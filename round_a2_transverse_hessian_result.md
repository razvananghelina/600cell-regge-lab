# Transverse conformal Hessian of ordinary de Rham `A2`

Date: 2026-08-12

Preregistered protocol commit: `b366a72`

Prior homogeneous input commit: `c3ec3c9`

Registered verifier:
`reproducible/verify_round_a2_transverse_hessian.py`

Machine-readable result:
`reproducible/round_a2_transverse_hessian.json`

## Headline

> **DERIVED SMOOTH SADDLE.** On the smooth unit round three-sphere, the
> normalized ordinary full-de Rham coefficient `A2` has a strictly negative
> Hessian on the first non-gauge conformal scalar mode (`l=2`), while the
> independently certified homogeneous Hopf trace-free sector has five
> positive directions. The full smooth-metric Hessian is therefore
> indefinite: round `S3` is a saddle of `A2`, not a complete local minimum.

This is a hostile negative for promoting the earlier homogeneous and
one-dimensional Regge results to a full gravitational vacuum theorem. It
does **not** invalidate either restricted theorem.

The targeted verifier passes `17/17`. No full-suite run was used as evidence.

## 1. Complete hypotheses and scope

The statement assumes:

1. the smooth unit round `S3`, rather than the finite 600-cell edge-length
   configuration space;
2. the ordinary, ungraded de Rham operator `D=d+d*` on all differential forms;
3. the already derived smooth coefficient

   ```text
   A2(g)=-(2/3) integral R_g dVol_g;
   ```

4. removal of overall scale through

   ```text
   A2_hat(g)=Vol(g)^(-1/3) A2(g);
   ```

5. smooth conformal perturbations `g_epsilon=(1+epsilon*f)^4 g_round`;
6. the positive scalar-Laplacian convention
   `integral |grad f|^2=lambda integral f^2`.

The normalization is not an inserted stabilization. It is the
scale-invariant form of the functional. Every conformal path also has an
exactly equal-volume representative obtained by a global rescaling, and the
normalized functional is unchanged by that rescaling.

## 2. Exact calculation

For `u=1+epsilon*f`, the normalized Einstein--Hilbert/Yamabe functional is

```text
Y(u)=
  [8 integral |grad u|^2 + 6 integral u^2]
  / [integral u^6]^(1/3).
```

Let

```text
V=Vol(S3),  integral f=0,
F=integral f^2,
integral |grad f|^2=lambda F.
```

Expanding the numerator and volume denominator independently gives

```text
Y(1+epsilon*f)
=6 V^(2/3)
 + 8(lambda-3)V^(-1/3)F epsilon^2
 + O(epsilon^3),

delta^2 Y(f,f)=16(lambda-3)V^(-1/3)F.
```

Multiplication by the fixed ordinary de Rham factor `-2/3` yields

```text
delta^2 A2_hat(f,f)
=-(32/3)(lambda-3)V^(-1/3)F.
```

There is no numerical fit or tolerance in this calculation.

## 3. Explicit first physical conformal mode

Use the degree-two harmonic polynomial on `S3` in `R4`

```text
f=x1^2-x2^2.
```

The verifier checks independently and exactly:

```text
Delta_R4 f=0,
degree(f)=2,
integral f=0,
integral f^2=V/6,
integral |grad_S3 f|^2=4V/3=8 integral f^2.
```

Thus `lambda=8`, and

```text
delta^2 A2_hat(f,f)
=-(160/3)V^(-1/3)F
=-(80/9)V^(2/3) < 0.
```

The `l=1`, `lambda=3` sector gives zero, as required for conformal
diffeomorphisms of the round sphere. The `l=2` direction is not gauge or
scale: the conformal scalar-curvature variation is

```text
delta R=8(lambda-3)f=40f,
```

which is nonzero and nonconstant. An infinitesimal diffeomorphism gives
`delta R=L_X(6)=0`, while a scale variation gives a constant change.

## 4. Why the Hessian is indefinite

The independent homogeneous certificate gives

```text
delta^2 A2_hat(H,H)=(8/3) Tr(H^2)>0
```

on five nonzero trace-free Hopf directions. The conformal `l=2` direction
above is strictly negative. A quadratic form taking both signs is indefinite,
irrespective of possible cross terms. Therefore the round point is a smooth
metric saddle.

This does not contradict the global homogeneous theorem: a function can have
a unique minimum on a restricted submanifold while decreasing in a transverse
direction. It also does not contradict the continuous affine round--Regge
path certificate, which controls one different path.

## 5. Framing attack

The result is not evidence for a new physical instability unique to the
600-cell. The conformal sign is the familiar structural obstruction carried
by an Einstein--Hilbert-type Euclidean functional. What is new inside this
repository is that the obstruction has now been applied to the exact
ordinary-de Rham coefficient and placed against the previously positive Hopf
and Regge results under a preregistered decision boundary.

Nor does this prove that the fixed finite 600-cell Regge Hessian is
indefinite. A smooth conformal mode is not automatically an admissible vector
in a particular finite edge-length space. That discrete transverse Hessian,
and its behavior under refinement, remain **OPEN**.

The negative also cannot be repaired by changing the sign after the fact.
The positive spectral cutoff gives `A2` a positive asymptotic weight, and the
ordinary de Rham curvature multiplier is already fixed at `-2/3`. A formal
`+1` multiplier reverses the conformal sign, but it is a different operator
functional rather than a result of this theory.

## 6. Status ledger

| Claim | Status |
|---|---|
| Exact conformal expansion of normalized Einstein functional | **DERIVED** |
| Ordinary de Rham conformal Hessian formula | **DERIVED** |
| `l=1` conformal zero sector | **DERIVED** |
| Explicit `l=2` mode has `lambda=8` and `F=V/6` | **DERIVED** |
| `l=2` is outside diffeomorphism and scale directions | **DERIVED** |
| Ordinary de Rham `A2` is negative on `l=2` | **DERIVED** |
| Homogeneous Hopf sector is positive | **DERIVED PRIOR INPUT** |
| Full smooth-metric Hessian at round is indefinite | **DERIVED SMOOTH SADDLE** |
| Round is a local minimum of `A2` over all smooth metrics | **REFUTED** |
| Round is the unique minimum in the left-invariant fixed-volume class | **DERIVED; UNAFFECTED** |
| Round is the unique minimum on the certified affine Regge path | **DERIVED; UNAFFECTED** |
| Fixed finite 600-cell transverse Regge Hessian is indefinite | **DERIVED LATER AT THE EQUILATERAL REGGE POINT** |
| Higher spectral terms stabilize the conformal direction | **OPEN** |
| Complete finite-cutoff action has a stable round vacuum | **OPEN** |
| Lorentzian dynamics, `G`, Planck mass or Planck time follows | **OPEN** |

## 7. Physical consequence

The strongest honest statement is now:

> The theory's ordinary geometric operator supplies a real smoothing signal
> along the certified homogeneous and round--Regge directions, but its leading
> curvature coefficient alone is not a stable gravitational action on the
> full smooth metric space.

Consequently higher spectral orders, a derived constraint/contour, or a
different dynamical principle are not optional refinements; at least one is
load-bearing if the programme is to obtain a stable gravitational vacuum.
Selecting one of them after seeing this sign would be fitting unless it is
fixed by a new target-independent protocol.

The later preregistered audit `finite_regge_a2_hessian_result.md` finds the
same obstruction intrinsically at the equilateral finite-Regge point: after
scale is removed, the full 719-dimensional edge Hessian has inertia
`(569,0,150)`, and its frozen nine-dimensional discrete conformal carrier is
negative. This is structural sign agreement, not yet a refinement theorem.

## 8. Reproduction note

The first execution passed `16/17`; its sole failure was a textual scope guard
because the searched phrase crossed a Markdown line break. Only that string
predicate was corrected; no formula, sign, input or decision boundary
changed. The targeted rerun passed `17/17`.

No PDF was built.
