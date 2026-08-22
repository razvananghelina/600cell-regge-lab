# Frozen method correction: integral Taylor remainders

Date: 2026-08-22.

Parent protocol commit: `559ca7a`.

Status: frozen before the primary verifier exists or any invariant-region
inequality is evaluated.

## Unchanged content

The half-strip

```text
0<m<=2/5,
x>=125
```

and every exact, physical, global-root, control and outcome gate in the parent
protocol remain unchanged.  No threshold or sign target is altered.

## Reason for the correction

The parent protocol requested degree-12 Taylor models with high-order
derivative remainders.  Before implementation, rewriting in

```text
u=t^2
```

showed that the removable quotients admit lower-order *exact integral Taylor
remainders*.  These avoid unnecessary dependency growth in a fourteenth
derivative and are mechanically stronger.  This is a pre-run method
correction, not a response to an observed interval failure.

## Corrected primary enclosure

Let `M(u)` and `P(u)` denote the regular functions in the parent protocol and
set

```text
N(u)=P(u)-P(0)+2*pi*u*M(u),
C(u)=[P(u)-P(0)]/u+4*pi*M(u).
```

Use the exact identities

```text
W(u)=[P(u)-P(0)]/u
    =integral_0^1 P'(s*u) ds,

Bbar(u)=N(u)/u^2
       =integral_0^1 (1-s)*N''(s*u) ds,

W'(u)=integral_0^1 s*P''(s*u) ds,

-C'(u)=-W'(u)-4*pi*M'(u).
```

Consequently outward-rounded Arb bounds on `P'`, `P''`, `N''`, `M` and `M'`
over the one rational interval

```text
0<=u<=4/390625
```

rigorously enclose `W`, `Bbar` and `-C'`; the weights integrate to `1`,
`1/2` and `1/2`, respectively.  Formal even coefficients through degree 12
in `t` remain required as a non-load-bearing algebraic control, but the exact
integral identities replace the high-order remainder as the load-bearing
certificate.

With

```text
U=z*M,
B=u*Bbar,
C=W+4*pi*M,
r=2/U-1,
```

the decisive normalized gap must be derived exactly as

```text
G/(m^2*z^2)
 =4*(1-U)*Bbar/M^2
  -z^2*(1-r^(-2))*Cbar,
```

where `Cbar` is the mean value of `-C'(u)` between
`u_plus=(m_plus*z)^2` and `u=(m*z)^2`.  Thus the same Arb interval for
`-C'` encloses `Cbar`.  This formula includes both axes continuously and is
the load-bearing no-grid proof of the same-`x` bracket.

The remaining inequalities may be certified directly from these one-variable
bounds.  Deterministic bisection at 192 digits remains a permitted fallback
only if a full-interval Arb enclosure is inconclusive.  Maximum depth 28 and
the frozen splitting order remain unchanged.

## Falsification rule

If any full-interval or subdivided Arb bound includes the wrong sign and
cannot be resolved within the frozen depth, the outcome is `OPEN`.  If a
strict wrong sign is certified, the candidate half-strip is refuted.  No
smaller domain or alternative coefficient construction may replace it in this
mission.
