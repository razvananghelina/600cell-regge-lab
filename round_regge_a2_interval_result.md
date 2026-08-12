# Continuous round--Regge path certificate

Date: 2026-08-12

Preregistered protocol commit: `8b322b1`

Verifier: `reproducible/verify_round_regge_a2_interval.py`

Machine-readable result: `reproducible/round_regge_a2_interval.json`

## Headline

Under the complete hypotheses below, the finite-grid pattern has been
promoted to a computer-assisted continuous theorem:

> **DERIVED PATH SELECTION.** On the frozen affine metric path from the
> fixed regular 600-cell Regge metric (`u=0`) to the unit-round metric on
> `S^3` (`u=1`), the equal-volume ordinary complete-exterior de Rham heat
> coefficient satisfies
>
> `A2_eq'(u)<0` for every real `0<=u<1`.
>
> Consequently `A2_eq(u)>A2_eq(1)` for every `0<=u<1`; the round endpoint
> is the unique minimizer on this path.

This is one Riemannian path and one heat coefficient. It is not a global
metric theorem and does not derive time, Einstein dynamics, `G`, or a Planck
scale.

## 1. Complete hypotheses

The statement assumes all of the following:

1. the carrier is the radially identified boundary of the regular,
   unit-circumradius 600-cell;
2. on each tetrahedron
   `g_u=(1-u)I+u(I/r^2-yy^T/r^4)`;
3. the operator is the ordinary complete-exterior Hodge--de Rham Laplacian;
4. faces carry the previously derived transmittal domain, continuous in the
   form and in `(d+d*)`;
5. edges carry the exact five-sector tangent-cone Hodge coefficient;
6. the assembled coefficient is normalized to volume `2*pi^2` and no local
   relative weight is fitted;
7. the standard stratified heat formula and the first-variation formula for
   the distributional scalar-curvature functional are applicable to this
   continuous piecewise-smooth path.

The local operator/domain hypotheses are those audited before the interval
protocol. The seventh item is an analytic theorem input, not a numerical
finding and not a hidden claim that the path is smooth for `u<1`.

## 2. Exact coefficient and endpoint lemma

Put

```text
a^2=(7+3*sqrt(5))/16,
r^2=a^2+|y|^2,
q=1-u+u/r^2,
p=1-u+u*a^2/r^4.
```

The frozen coefficient is

```text
A2_raw=B+F+E,

B=-(2/3) integral R dV,
F=-(4/3) integral_faces (L+ + L-) dA,
E=integral_edges [-(4/3)delta+4*delta^2/(3*beta)] dl,

delta=2*pi-beta,
A2_eq=(2*pi^2/V)^(1/3) A2_raw.
```

Thus the bulk, face and linear-deficit terms are exactly `-2/3` times the
distributional total scalar curvature: the codimension-one mean-curvature
jump and codimension-two deficit each have the required relative factor two.
The remaining heat-specific cone term is exactly quadratic in `delta`.

At `u=1`, all jumps and deficits vanish and the metric is the unit-round
metric. For a covariant variation `h` in dimension three,

```text
delta integral R dV
 = integral (R*g/2-Ric):h dV
 = integral tr(h) dV
 = 2 delta V,
```

because `R=6` and `Ric=2g`. Since `(integral R dV)/V=6`, the first variation
of `(integral R dV)*V^(-1/3)` is exactly zero. Internal face and linear corner
variations are the boundary/corner completion of this formula, while
`4*delta^2/(3*beta)` has zero first variation at `delta=0`. Therefore

```text
A2_eq'(1)=0
```

exactly. This is the endpoint fact used on the final interval.

A hostile independent differentiation of the explicit radial formulas gave

```text
B'(1)+F'(1)+E'(1)
  = -4.1764083149295052060343290057878568555...,
(4/3)V'(1)
  = +4.1764083149295052060343290057878568555...,
residual
  = 5.4e-79.
```

That numerical cancellation is a cross-check of every component, not the
proof of the exact endpoint lemma.

## 3. Continuous certificate

The protocol covers `[0,1]` by the 20 closed cells
`I_j=[j/20,(j+1)/20]`. At each midpoint it computes the formal Taylor series
of the normalized coefficient through degree 18 using Arb balls.

Every spatial integral uses composite two-node Gauss--Legendre quadrature.
The order-16 enclosure includes a cellwise fourth-derivative remainder; the
same test is repeated at the preregistered orders 20 and 24. The omitted
Taylor tail is bounded by Cauchy's estimate on a complex disk of radius
`1/10`, with cell half-width `1/40` and the preregistered majorant `M=1000`.

The complex audit obtained the stronger bounds

```text
min |q|                         > 0.9832815
min |p|                         > 0.8424390
min Re(p)                       > 0.8424390
min |face denominator|          > 0.8599457
min |V(z)|                      > 15.3191493
max |delta|                     < 0.1520674
min |beta|                      > 6.1319899
max |A2_eq(z)|                  < 213.448
```

so the frozen Cauchy majorant has a factor greater than four of slack. The
square-root and inverse-cosine arguments stay on their analytic branches.

At every one of the three spatial orders:

- `A2_eq'(I_j)` is strictly negative for `j=0,...,18`;
- `A2_eq''(I_19)` is strictly positive.

The weakest required bounds were:

| Spatial order | upper bound for `A2_eq'(I_18)` | lower bound for `A2_eq''(I_19)` |
|---:|---:|---:|
| 16 | `-0.00354591866452637...` | `+0.247887278451632...` |
| 20 | `-0.0114408331638871...` | `+0.264914880440868...` |
| 24 | `-0.0131907048499377...` | `+0.268668405156668...` |

On the last cell, strict convexity and `A2_eq'(1)=0` imply
`A2_eq'(u)<0` for `0.95<=u<1`. Together with the other 19 cells this proves
the headline inequality.

The exact radial reduction also encloses all six values previously computed
by an independent three-dimensional Duffy implementation at `u=1/2`.

Targeted verifier result:

```text
RESULT: 37/37 checks passed
DERIVED PATH SELECTION
elapsed: 3348.6 s with four workers, plus a sub-second exact algebra post-check
```

The integral run initially exited `34/37` because three exact rational
identities were compared with `Arb ==` rather than exact rational/symbolic
predicates. All 34 analytic/numerical checks had passed. Replacing only those
three predicates by `Fraction`/SymPy and running the dedicated algebra audit
gave `4/4`; no interval, Taylor coefficient, bound or decision criterion was
changed. The stored `37/37` is this transparent aggregate, not a claim that
the pre-fix process returned zero.

No full repository suite was run, following the user's explicit instruction.

## 4. Framing attack and literature boundary

The standard spectral-action positivity axiom was subsequently audited in
`round_regge_spectral_action_sign_result.md`. It fixes the asymptotic `A2`
weight to be positive, so that sign preserves the preference. The cutoff
moment, scale, higher coefficients and finite-action remainder remain
unselected and can still prevent complete-action minimization.

The theorem also controls only one affine direction in metric space. A curve
can be minimizing while a transverse direction is unstable. The next
geometric gate is therefore the complete canonical transverse Hessian,
modulo scale and diffeomorphism directions.

Primary literature already uses the 600-cell as a Regge approximation to a
closed FLRW spatial slice and studies its evolution, for example
[Barrett et al.](https://arxiv.org/abs/gr-qc/9411008),
[De Felice--Fabri](https://arxiv.org/abs/gr-qc/0106077), and
[Tsuda--Fujiwara](https://arxiv.org/abs/2011.04120). The general transmittal
heat input is also established independently
([Gilkey--Kirsten--Vassilevich](https://arxiv.org/abs/hep-th/0101105)). A
targeted search did not find the particular complete-de-Rham,
transmittal-plus-cone, equal-volume monotonicity theorem proved here. That is
not proof of bibliographic novelty; a specialist literature review is still
required before any novelty claim.

## 5. Status ledger

| Claim | Status |
|---|---|
| Exact radial/Duffy changes of variables | **DERIVED** |
| Uniform complex analyticity on all Taylor disks | **DERIVED** |
| Validated spatial quadrature remainders at orders 16/20/24 | **DERIVED** |
| Continuous negativity of `A2_eq'` on `[0,0.95]` | **DERIVED** |
| Exact stationarity at the round endpoint under the stated stratified-variation hypothesis | **DERIVED** |
| Continuous negativity of `A2_eq'` on `[0,1)` | **DERIVED** |
| Unique round minimum on this one path | **DERIVED PATH SELECTION** |
| Unique minimum among all admissible 600-cell metrics | **OPEN** |
| Positive asymptotic `A2` sign in the standard spectral action | **DERIVED LATER** |
| Dominance in the complete finite-cutoff spectral action | **OPEN** |
| Lorentzian evolution, causal speed, universal gravity, `G`, or Planck units | **OPEN** |
| Bibliographic novelty | **OPEN** |

## 6. Ordered continuation

1. derive, rather than assume, the sign and weight of this `A2` term in the
   full admissible spectral action;
2. enumerate the full coefficient-free canonical metric-deformation space
   and compute the constrained Hessian at round;
3. if that Hessian is positive modulo gauge/scale, test a local transfer law
   between selected spatial slices for finite propagation speed;
4. seek an independently selected dimensional scale before discussing
   Newton, Planck mass, or Planck time.

Skipping step 4 would be dimensionally empty: the current construction is
scale covariant and cannot select an absolute mass or time by itself.
