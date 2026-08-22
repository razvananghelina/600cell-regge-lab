# Result: the internal kernel is the lapse-constraint tangent

Date: 2026-08-22.

Status: **DERIVED EXACT/COMPUTATIONAL; BOUNDED NEGATIVE FOR FIXED-INPUT
EVOLUTION.**

## Complete hypotheses

Use the first positive-height homogeneous 600-cell dust slab from `v=3/2`,
zero cosmological constant, conserved dust, positive height and endpoint
scales, either certified staircase parity, and the exact rank-240
scale-plus-strut carrier.  Hold the old spatial boundary and incoming
canonical momentum fixed.

The replicated internal-rank result gives

```text
rank R_p=239,
ker R_p is one-dimensional and homogeneous.
```

The present result identifies that line and asks whether it survives the
fixed incoming momentum equation.  It does not study a forced response to
varying incoming data.

## Reproducible verdict

The registered exact reconciliation reports

```text
INTERNAL_KERNEL_IS_LAPSE_CONSTRAINT_TANGENT_FIXED_INPUT_REMOVES_IT
14/14 PASS
```

Implementation commit: `8387926`.

Artifact commit: `bf0e694`.

Artifact SHA-256:

```text
81ec0379247023451e82ab42f5beb026ee2d1b083aa5e2553e42b894554266f6
```

Only this verifier and static registry checks ran.  The full suite did not
run.

## Exact coordinate identity

For normalized fixed incoming data `(m,pi)`, the homogeneous one-slab
equations are

```text
C(h,q)=8*pi[mu(q)-m]+4*pi*h*q*mu(q),
P(h,q)=p(q)-pi-2*pi*h*mu(q).
```

`C=0` is the internal lapse constraint; `P=0` fixes the incoming canonical
momentum.  The already certified state identity

```text
4*pi*mu'(q)+q*p'(q)=0
```

gives

```text
det partial(C,P)/partial(h,q)=8*pi^2*h*mu(q)^2>0.
```

The present carrier coordinates satisfy

```text
sigma=delta lambda,
c=delta log rho,
lambda=1+h*q,
rho=h^2.
```

Therefore

```text
delta h=h*c/2,
delta q=sigma/h-q*c/2,

partial(h,q)/partial(sigma,c)
  =[[0,h/2],[1/h,-q/2]],

det=-1/2.
```

The combined canonical determinant is consequently

```text
det partial(C,P)/partial(sigma,c)
  =-4*pi^2*h*mu(q)^2 !=0.
```

The exact `C`-tangent ratio is

```text
c/sigma=-2*C_q/[h*(h*C_h-q*C_q)].
```

At the frozen state it is

```text
c/sigma=0.4589898592210244392413761746604434209605...
```

and the fixed-momentum response for `sigma=1` is

```text
dP=2.1296728633240540701... !=0.
```

Symbolically, on the unnormalised tangent generator used by the verifier,

```text
dP=4*pi^2*h*mu(q)^2.
```

## Independent projector closure

The analytic 240-vector

```text
(1 repeated 120 times, (c/sigma) repeated 120 times)
```

was compared only after its exact derivation with four stored projectors:

| Construction | Even distance | Odd distance |
|---|---:|---:|
| 180-digit orbit-sector primary | `9.82e-15` | `4.84e-15` |
| complete-real-space adversarial | `2.29e-11` | `2.26e-11` |

The agreement gate was `9.33e-10`.

Two hostile conventions were strongly rejected:

```text
sigma=delta log lambda   projector distance 0.4850...
c=delta rho              projector distance 0.3997...
```

Omitting the factor `1/2` in `delta h=h*c/2` changes the exact coordinate
determinant from `-1/2` to `-1` and also fails its symbolic gate.

## Logical conclusion

Because the complete internal kernel is exactly the homogeneous `dC=0`
line and `dP` is nonzero on it,

```text
ker R_p intersect ker dP = {0}
```

for both schedule parities.

**DERIVED BOUNDED NEGATIVE:** the sole internal survivor is not a free tick,
gauge-protected physical mode or fixed-input evolution direction.  It is the
tangent to the lapse constraint, and the fixed incoming canonical momentum
removes it.

## Framing correction: zero fixed-input kernel is not zero dynamics

This result does **not** say that the one-slab dynamics has no response.  A
locally regular canonical relation should give an isolated output when its
incoming canonical data are held fixed.  The positive determinant is the
same local-regularity property that permits implicit evolution.

The physically meaningful linear object is therefore a forced derivative

```text
(delta incoming geometry, delta incoming momentum)
    -> (delta outgoing geometry, delta outgoing momentum),
```

obtained from the action's implicit canonical/Jacobi equations.  The
homogeneous finite-height map already supplies its minisuperspace analogue.
The nonhomogeneous version has not been derived.

Thus this result closes a false interpretation rather than the entire Regge
programme.  It also explains why searching the fixed-input kernel for
gravitons was the wrong observable.

No physical tick, limiting speed, graviton, wave equation, continuum limit,
`G`, Planck scale, particle mass or Standard-Model sector is derived.

## Exact next gate

Formulate, before computing a spectrum, the geometry-selected incoming
canonical perturbation space and the forced linearized one-slab equations.
Then test existence, uniqueness and schedule independence of the implicit
map from incoming to outgoing data.  Do not enlarge the 240-dimensional
configuration carrier by fitted directions; the cotangent/input variables
must come from the action's boundary Legendre transform.
