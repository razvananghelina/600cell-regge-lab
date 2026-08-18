# The ordinary de Rham coefficient on the full round--Regge path

Date: 2026-08-12

Protocol commit: `033c935` (written before any interior value was evaluated)

Verifier: `reproducible/verify_round_regge_interior_a2.py`

## Headline

The complete preregistered numerical audit found no interior counterexample:

> **PATTERN:** on every point `u=j/200`, at all five frozen quadrature
> orders, the equal-volume ordinary complete-de Rham coefficient decreases
> strictly from the fixed Regge metric toward the unit-round metric.

This note records the preregistered finite-grid stage. It was subsequently
superseded by the continuous Arb certificate in
`round_regge_a2_interval_result.md`, which proves strict monotonicity on every
real `0<=u<1` under the same stated operator/domain hypotheses. Global path
selection is therefore **DERIVED on this one path**, while selection among
all admissible metrics remains OPEN.

The result is nevertheless stronger than the previous two-endpoint
comparison. It uses the theory's de Rham operator and a coefficient fixed by
the metric. There is no fitted weight between the volume, face and edge terms.

## 1. Complete hypotheses

The result is conditional on all of the following, with none left implicit:

1. the carrier is the radially identified boundary of the regular
   unit-circumradius 600-cell;
2. on every facet,
   `g_u=(1-u)I+u(I/r^2-yy^T/r^4)`, exactly as preregistered;
3. the operator is the ordinary, not graded, complete-exterior Hodge--de Rham
   Laplacian;
4. open faces use the de Rham transmittal domain selected by continuity of
   both the form and `(d+d*)`;
5. open edges use the exact tangent-cone Hodge coefficient with the derived
   five-sector `C5` link;
6. the assembled coefficient is rescaled only to volume `2*pi^2`;
7. minimizing this single `A2` coefficient is treated only as a conditional
   selector, not as an already derived physical action.

The heat convention is

```text
Tr exp(-t Delta_u)
  ~ (4*pi*t)^(-3/2) [A0(u)+t*A2(u)+...].
```

## 2. Independent local derivation

Put `s^2=|y|^2`, `r^2=a^2+s^2` and

```text
a^2=(7+3*sqrt(5))/16,
q=1-u+u/r^2,
p=1-u+u*a^2/r^4.
```

The metric has radial and transverse eigenvalues `p,q,q`, hence
`sqrt(det g)=q*sqrt(p)`. Reducing it to the warped product

```text
g=p ds^2+q s^2 dOmega^2
```

gives the exact scalar curvature

```text
R(s,u)=8*u*a^2/(q*p*r^6)-2*u/(q^2*r^4).
```

The verifier symbolically recovers this expression from the warped-product
curvature formula. It gives `R=0` at `u=0` and `R=6` at `u=1`.

For `g_ij=q delta_ij+B y_i y_j`, `B=-u/r^4`, the only connection term that
survives tangential projection onto the representative face is

```text
Gamma_radial = 2*u*(1-u)/(r^6*q*p).
```

This proves that the face curvature has the mandatory factor `u(1-u)`.
It is not optional and is nonzero in the interior.

The sign convention was checked against the primary transmittal theorem,
which defines

```text
L_ab^+/- = (nabla^+/-_{e_a} e_b, nu^+/-),
```

with inward normals. That is exactly the convention implemented here. For
the complete exterior algebra in dimension three, direct tracing of the
de Rham transmittal endomorphism gives `Tr(U)=4(L_aa^++L_aa^-)`, while
`Tr(I)=8`. The face coefficient is therefore

```text
A2_face=-(4/3) integral_face (L_aa^+ + L_aa^-) dA.
```

Source: [Gilkey--Kirsten--Vassilevich, Theorem 2.3 and Lemma
6.1](https://arxiv.org/abs/hep-th/0101105).

At a representative edge the two inward face covectors determine the sector
angle intrinsically:

```text
cos(theta_u)=-(n1^T g_u^-1 n2)
             /sqrt((n1^T g_u^-1 n1)(n2^T g_u^-1 n2)),
beta_u=5*theta_u.
```

The exact full-de Rham tangent-cone density is

```text
16*pi^2/(3*beta_u)+8*beta_u/3-8*pi.
```

At this heat order a codimension-two density must be dimensionless;
extrinsic curvature has inverse-length dimension and first enters one order
later. The `A2` edge density therefore depends on the tangent-cone angle,
not on a fitted seam coefficient. This uses the standard locality of the
stratified heat expansion; a dedicated theorem for this exact mixed
piecewise-smooth family has not been reproduced from first principles here.
That analytic input is stated rather than hidden.

## 3. Complete assembled coefficient

With the derived multiplicities of 600 tetrahedra, 1200 faces and 720 edges,
the raw result is

```text
A2_raw=B+F+E,

B=-(2/3)*600 integral_tetra R dV,
F=-(4/3)*1200 integral_face (L_aa^+ + L_aa^-) dA,
E=720 integral_edge
  [16*pi^2/(3*beta)+8*beta/3-8*pi] dl.
```

There is no freedom to tune the relative coefficients. Equal-volume
comparison uses

```text
A2_eq=(2*pi^2/V(u))^(1/3) A2_raw.
```

The face term is positive on this path. Omitting it would be a different,
incorrect operator-domain calculation.

## 4. Numerical result and hostile controls

Tensor Gauss--Legendre quadrature after the frozen Duffy maps was run at
orders `16,24,32,40,48`; the edge rule uses twice each order. The last three
orders agree over the full 201-point grid within

```text
volume      3.91e-14
bulk        1.42e-13
face        8.66e-15
edge        6.32e-13
raw A2      6.25e-13
normalized  6.25e-13
```

Representative order-48 values are:

| `u` | volume | bulk `B` | face `F` | edge `E` | equal-volume `A2` |
|---:|---:|---:|---:|---:|---:|
| 0.000 | 16.6925267711 | 0 | 0 | -74.5853661323 | -78.8719985927 |
| 0.100 | 16.9893449514 | -8.1486905483 | 0.4863393158 | -67.3727618291 | -78.8827865058 |
| 0.250 | 17.4378947560 | -20.2608108764 | 1.0079574615 | -56.4534225500 | -78.9000199696 |
| 0.500 | 18.1942438235 | -40.1620264568 | 1.3332757773 | -37.9836200041 | -78.9277542590 |
| 0.750 | 18.9614029535 | -59.7217512512 | 0.9928683783 | -19.1688730940 | -78.9486493941 |
| 0.900 | 19.4268190093 | -71.3007286570 | 0.4747412264 | -7.7107235722 | -78.9554401126 |
| 0.990 | 19.7078941334 | -78.1934643411 | 0.0521084403 | -0.7736898693 | -78.9568207529 |
| 0.995 | 19.7235493692 | -78.5752116038 | 0.0261827383 | -0.3869179556 | -78.9568315878 |
| 1.000 | 19.7392088022 | -78.9568352087 | 0 | 0 | -78.9568352087 |

Every order separately gives strict grid monotonicity. The weakest adjacent
decrease is `-3.621e-6`; the smallest interior margin above round is
`+3.621e-6`. Both are more than six orders of magnitude above the measured
quadrature spread.

Near round, the sampled values show cancellation of the linear contributions
between bulk, face, edge and volume normalization. On the grid the quotient

```text
[A2_eq(u)-A2_eq(1)]/(1-u)^2
```

stays positive in `[0.0848366,0.144837]`. This is useful evidence that the
endpoint approach is quadratic rather than a missed linear crossing, but it
is still a sampled statement, not an interval proof.

Both exact endpoint decompositions pass at all five orders:

```text
u=0: B=F=0; E gives the exact Regge cone result,
u=1: V=2*pi^2; B=-8*pi^2; F=E=0.
```

## 5. Framing attack

The calculation does **not** derive gravity.

- **DERIVED:** the metric path, local bulk formula, transmittal coefficient,
  conic edge coefficient, endpoint values and quadrature convergence.
- **PATTERN:** all preregistered interior samples prefer round and do so
  monotonically.
- **DERIVED LATER:** the continuum inequality on the entire path is proved by
  the preregistered Arb certificate in `round_regge_a2_interval_result.md`.
- **OPEN:** selection among all `H4`-invariant metrics. One affine path is
  not the whole metric space.
- **OPEN:** why a physical spectral action should minimize this coefficient
  with this sign and dominate over its other heat orders.
- **OPEN:** Lorentzian time, a kinetic term, diffeomorphism constraints,
  universal stress-energy coupling, Newton's constant and a Planck scale.

The honest physical statement is therefore limited but useful:

> the same operator-derived curvature coefficient that removes homogeneous
> Hopf anisotropy also behaves as a smoothing energy all along the tested
> Regge-to-round interpolation. It supplies a coherent Euclidean shape
> preference, not yet gravitational evolution.

## 6. Status ledger

| Claim | Status |
|---|---|
| Radial scalar-curvature formula | **DERIVED** |
| Face term has exact `u(1-u)` factor | **DERIVED** |
| GKV sign and full-exterior trace give `-4/3 L_sum` | **DERIVED** |
| Both exact endpoints and all quadrature controls pass | **DERIVED** |
| Last three frozen orders agree below `6.3e-13` | **DERIVED NUMERICAL** |
| All 200 interior grid values lie above round | **PATTERN** |
| Frozen grid is strictly monotone | **PATTERN** |
| Round globally minimizes `A2` on every real point of this path | **DERIVED later; see `round_regge_a2_interval_result.md`** |
| Round is selected among every admissible metric | **OPEN** |
| `A2` is the complete physical gravitational action | **OPEN** |
| Gravity, time, `G` or Planck units have been derived | **OPEN** |
