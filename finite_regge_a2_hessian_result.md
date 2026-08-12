# Full finite-Regge Hessian of ordinary de Rham `A2`

Date: 2026-08-12

Preregistered protocol commit: `70afdbc`

Registered verifier: `reproducible/verify_finite_regge_a2_hessian.py`

Machine-readable result: `reproducible/finite_regge_a2_hessian.json`

## Headline

> **DERIVED COMPUTATIONAL FINITE SADDLE.** In the complete local
> configuration space of all 720 edge lengths, the equilateral 600-cell is a
> stationary saddle of the exact equal-volume ordinary full-de Rham conical
> `A2`. After removing only global scale, the 719-dimensional Hessian has
>
> ```text
> inertia (positive, zero, negative) = (569, 0, 150).
> ```

The full nine-dimensional vertex-quadratic discrete conformal carrier is
strictly negative. Thus the conformal instability found on smooth round `S3`
survives in a canonical finite probe, although no refinement theorem yet
identifies the two Hessians.

The targeted verifier passes `18/18`. No full-suite run was performed.

## 1. Framing correction

The smooth round metric is not an element of the space of Euclidean Regge
edge lengths on the fixed triangulation. The finite point tested here is the
canonical equilateral, singular Regge metric. Therefore this is not literally
the “finite Hessian at round”. It answers the well-posed neighboring question:
whether the equilateral 600-cell is stable under arbitrary intrinsic Regge
edge deformations.

The continuous round--Regge path uses piecewise-smooth interior metrics and
leaves the pure Euclidean-Regge edge space. Its monotone smoothing direction
is consequently compatible with exact stationarity inside the pure Regge
space.

## 2. Complete functional and configuration space

Use squared edge lengths `x_e=l_e^2`. Every assignment in a sufficiently
small neighborhood of `x_e=1` is admissible when all 600 local tetrahedral
Gram matrices remain positive definite. Shared faces automatically have the
same intrinsic metric because their three edge lengths are shared.

The exact coefficient is

```text
V(x)=sum_t Vol_t(x),
beta_e(x)=sum_(t contains e) theta_(t,e)(x),
C(beta)=16*pi^2/(3*beta)+8*beta/3-8*pi,

Ahat(x)=V(x)^(-1/3) sum_e sqrt(x_e) C(beta_e(x)).
```

The positive omitted factor `(2*pi^2)^(1/3)` normalizes the volume to that of
unit round `S3` and does not change signs. With it restored, the implementation
reproduces the independently committed endpoint

```text
A2_Regge=-78.871998592694... .
```

There is no face term for an intrinsic piecewise-flat Regge metric: each
shared Euclidean face is totally geodesic from both incident tetrahedra. Open
edges carry the conical `A2`, while vertices first enter the next heat order.

## 3. Stationarity is structural

At the equilateral point the `H4` symmetry is transitive on all 720 edges.
Consequently every component of the gradient is the same constant `c`.
Meanwhile

```text
Ahat(s*x)=Ahat(x)
```

because lengths scale as `sqrt(s)` and volume as `s^(3/2)`. Differentiating
in the all-ones direction gives `720c=0`, hence the gradient vanishes exactly.

The numerical automatic-differentiation audit found

```text
max |gradient component| < 4.8e-16,
RMS |Hessian * scale|    < 9.1e-15.
```

Only this exact scale mode was removed. No prospective vertex-displacement
or continuum diffeomorphism direction was discarded.

At a stationary point, changing from squared lengths to lengths gives

```text
H_l=4 H_x
```

at `l_e=1`. It is a nonsingular congruence, so Sylvester inertia is unchanged.
The verifier checks both coordinate descriptions.

## 4. Full Hessian result

Every local tetrahedral volume and dihedral angle was differentiated to
second order by forward automatic differentiation before global assembly.
The calculation includes the second derivatives of angles and volumes,
`C''(beta) d beta tensor d beta`, and all mixed derivatives from
`F V^(-1/3)`.

On the scale quotient:

```text
minimum eigenvalue             = -0.0691635459709...
maximum eigenvalue             = +4.16247741411...
smallest absolute eigenvalue   =  0.0000695284957411...
inertia                        = (569,0,150).
```

The seven negative numerical symmetry clusters are

```text
-0.0691635459709... x  9
-0.0675136709728... x 16
-0.0246592846802... x 25
-0.00380138765069.. x 36
-0.00190991618123.. x 24
-0.00139801604707.. x 16
-0.000968159727765. x 24
                         ---
                         150
```

The remaining 569 eigenvalues are positive. There are no extra numerical
null modes after scale is removed.

The sign count is stable when local analytic constants are evaluated at 50,
80 and 120 decimal digits. Symmetric diagonalization and an independent
pivoted `LDL^T` congruence both give `(569,0,150)`. The smallest claimed
nonzero magnitude is almost 700 times the preregistered `10^-7` acceptance
threshold.

## 5. Discrete conformal sector

Before calculation the protocol froze every trace-free quadratic form `Q` on
the ambient quaternion coordinates and the induced edge variation

```text
f_Q(v)=v^T Q v,
h_Q(ij)=f_Q(v_i)+f_Q(v_j).
```

These variations span nine dimensions and are orthogonal to scale. The
Hessian compression to this whole carrier is

```text
Q_conf^T H Q_conf
  = -0.0207468181369... I_9,
```

so every nonzero preregistered discrete conformal direction is negative.
The explicit analogue of the smooth mode `f=v_1^2-v_2^2` has the same value.

An important hostile correction was made during review: this compression is
**not** a full Hessian eigenspace. Although the carrier is symmetry-invariant,
the 720-edge representation contains enough multiplicity for an equivariant
Hessian to mix it with equivalent copies. The measured off-carrier norm is
`0.0986205...`. Thus the displayed scalar is a negative restricted quadratic
form, not an eigenvalue of the complete spectrum. Its only required
consequence remains exact: the negative index is at least nine.

## 6. Independent controls

The verifier additionally checks:

- the complete f-vector `(120,720,1200,600)` and five tetrahedra per edge;
- the regular angle by both projection and Cayley--Menger inverse formulas;
- covariance of all local volume/angle first and second derivatives under all
  24 tetrahedral vertex relabellings;
- the full order-20 edge stabilizer and all 62 relative-position orbits;
- direct nonlinear five-point second derivatives on 61 edge-orbit contrasts
  and all nine conformal directions.

The 70 nonlinear controls reproduce the automatic Hessian to maximum relative
error `2.23e-7`. They independently display both signs:

```text
smallest positive normalized orbit contrast > +0.4265,
largest conformal second derivative          < -0.0207.
```

Thus the saddle verdict does not rely only on diagonalizing one assembled
matrix.

## 7. Status ledger

| Claim | Status |
|---|---|
| Exact finite conical functional on all 720 edges | **DERIVED FROM PRIOR ANALYTIC INPUT** |
| Equilateral point is stationary | **DERIVED STRUCTURALLY** |
| Global scale is an exact Hessian null mode | **DERIVED** |
| No other quotient null modes occur | **DERIVED COMPUTATIONAL** |
| Quotient inertia is `(569,0,150)` | **DERIVED COMPUTATIONAL** |
| Equilateral 600-cell is a local minimum of `A2` | **REFUTED** |
| Equilateral 600-cell is a finite-Regge saddle | **DERIVED COMPUTATIONAL FINITE SADDLE** |
| Entire frozen nine-dimensional conformal carrier is negative | **DERIVED COMPUTATIONAL** |
| That carrier is a full Hessian eigenblock | **REFUTED** |
| Smooth and finite conformal signs agree | **STRUCTURAL AGREEMENT** |
| Finite Hessian converges to the smooth Hessian under refinement | **OPEN** |
| Positive `A2` alone supplies a stable finite or smooth vacuum | **REFUTED** |
| Higher spectral terms stabilize all negative modes | **OPEN** |
| Complete finite-cutoff action is stable | **OPEN** |
| Lorentzian gravity, `G`, Planck scale or a graviton follows | **OPEN** |

## 8. Physical consequence

The combined smooth and finite result closes one tempting route:

> With the standard positive spectral-action sign, the leading ordinary
> de Rham curvature coefficient is not a stable gravitational vacuum selector.
> It is a saddle both on the full smooth metric space and at the canonical
> finite Regge point.

This does not kill spectral dynamics. It makes a higher-order term,
constraint/contour or different derived action load-bearing. No such term may
now be chosen because it repairs these 150 directions; it must be frozen by a
new target-independent protocol and then tested against the complete negative
subspace.

## 9. Reproduction history

The first targeted execution found `(569,0,150)` and passed the 12 checks then
present. Hostile review noticed that the conformal compression value was not
in the full spectrum. Work stopped before documentation. The discrepancy was
resolved as equivariant copy mixing, and explicit symmetry-invariance,
off-carrier leakage, coordinate-change and nonlinear two-sign controls were
added. None changed the Hessian or its inertia. The final targeted verifier
passes `18/18`.

No PDF was built.
