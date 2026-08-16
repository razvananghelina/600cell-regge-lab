# Exact regular-lapse identity of the 600-cell dust slab

Date: 2026-08-16

Prior-art gate: `898d123`

Frozen protocol: `cf492b9`

Implementation registered before evaluation: `55d953c`

Preserved implementation-control failures: `6d5330c`, `4e7acaa`

Exactness-only corrections: `70d4477`, `edb2604`

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_regular_lapse_identity.py`

Artifact:
`reproducible/gravity_600cell_dust_regular_lapse_identity.json`

Artifact SHA-256:
`5079428fade247f730ebc07e5e2eae388b48045cd5201e84afb3186bfc248a51`

## Verdict

**DERIVED:** the targeted verifier passes `13/13` and returns

```text
REGULAR_LAPSE_IDENTITY_PROVED
```

For both derived order-24 schedule parities, retain the complete corrected
Lorentzian Regge plus dust action and set

```text
q_old = q_new = L^2,
x_diagonal = L^2-rho,
x_pole physical square = -rho,
tau=sqrt(rho)>0.
```

On the entire connected interval `0<rho<=rho0`, exactly

```text
S_grav =  720 epsilon3 L tau,
S_dust = -720 epsilon3 L tau,
S_total = 0,

g_internal = 0                         on all 35 orbit coordinates,
p_pre,e  = -epsilon3 L tau/4           on every old boundary edge,
p_post,e = +epsilon3 L tau/4           on every new boundary edge,
```

where `epsilon3=2*pi-5*acos(1/3)` and
`M=(90/pi)*epsilon3*L` is the already derived dust normalization.

This is an exact identity, not an interpolation through the 41 previously
known continuation points.

## Proof audit

### 1. Product geometry and Lorentzian branch

The regular squared distances are realized explicitly by a regular spatial
tetrahedron at times `0` and `tau` in a Minkowski product:

```text
same-layer spatial distance^2   = L^2,
different-vertex cross distance^2 = L^2-tau^2,
same-vertex pole distance^2       = -tau^2.
```

The 2,400 staircase 4-simplices reduce to twenty ordered and eleven marked
angle types.  After normalizing `L=1` and writing `u=rho/L^2`, every simplex
has signed volume square

```text
-u/1152.
```

The ten exact leading-minor forms have one sign change and no zero on
`0<u<1/2`, so every simplex has inertia `(3,1)`.  Facet, hinge and angle
arguments are also nonzero there.

The physical interval lies inside this domain by elementary exact bounds:
`rho0=(51/5000)^2<1`; `pi>3` and `sqrt(2)>1` give `zeta>1/2`; `pi<10/3`
gives `R0>4`; hence `L=zeta R0>2`, so `u0<1/4<1/2`.

### 2. Complete hinge census

Both parities have exactly the same census:

```text
2,400 boundary SSS triangles,
2,400 internal DDS triangles,
1,440 internal vertical DPS triangles.
```

For all fifteen complete angle-incidence patterns, exact products plus the
branch anchored at `u -> 0+` give

```text
boundary SSS curvature = 0,
internal DDS deficit   = 0,
vertical DPS deficit   = epsilon3.
```

Each vertical triangle has five incident 4-simplices, each with angle
argument

```text
1/3 - (2 sqrt(2)/3) i,
```

so its angle is `-acos(1/3)`.  The eleven radical argument forms remain
positive-real, negative-imaginary or in the fixed fourth quadrant throughout
the proof domain; none crosses the logarithm cut.

### 3. Action

A vertical triangle has Lorentzian area

```text
A_DPS = i L tau/2.
```

There are 1,440 of them, so the corrected action convention gives

```text
-i * 1440 * A_DPS * epsilon3 = 720 epsilon3 L tau.
```

The five order-24 pole orbits encode all 120 dust worldlines.  Their term is

```text
-8 pi M tau = -720 epsilon3 L tau,
```

and cancels gravity exactly.

### 4. Unrestricted internal equations

The action restriction alone would not prove stationarity.  The exact
edge-incidence and area derivatives close that gap.

Every staircase diagonal meets two curved vertical triangles, but on the
product geometry

```text
partial A_DPS / partial log(diagonal square) = 0.
```

Every pole meets twelve vertical triangles.  The per-pole logarithmic
derivatives are

```text
gravity = +3 epsilon3 L tau,
dust    = -3 epsilon3 L tau.
```

Thus all 720 diagonal-edge equations and all 120 pole-edge equations vanish;
expansion into the 35 stored orbit coordinates is exact in both schedules.

### 5. Boundary momentum

Every one of the 1,440 old/new boundary spatial edges meets exactly one
vertical triangle.  Its logarithmic derivative is

```text
g_boundary = +epsilon3 L tau/4.
```

With the canonical sign convention `p_pre=-g_old`, while
`p_post=+g_new`, the registered momentum formula follows.  Uniformity makes
it invariant under the independently derived old-to-final orbit permutation.

## Independent numerical control

The original, non-symbolic 100-digit geometric evaluator was run at the
preregistered ratios

```text
rho/rho0 in {1, 3/4, 1/2, 1/4, 1/16, 1/256}
```

for both parities.  The largest errors were

```text
internal-gradient absolute error  6.354e-95
pre-momentum relative error        1.121e-90
post-momentum relative error       5.775e-94
action relative error              1.832e-92
imaginary contamination            9.253e-97.
```

These evaluations independently control the symbolic convention and signs;
they are not the proof of the identity.

## Exact consequence for the failed continuation

At the published datum,

```text
P p_post(rho0) = -p_pre(rho0).
```

Therefore the committed homotopy target is exactly

```text
p(lambda)=(1-2 lambda) p_pre(rho0).
```

The positive-lapse regular solution is consequently

```text
tau(lambda)=tau0(1-2 lambda),
rho(lambda)=rho0(1-2 lambda)^2,     0<=lambda<1/2.
```

At `lambda=1/2`, `tau=rho=0`; the 4-simplices degenerate and this is not an
accepted Lorentzian frame.  For `lambda>1/2`, the target momentum has the
wrong sign for every `tau>0` on this orientation branch.  Thus the earlier
failure near `lambda=21/64` was a Newton/Armijo artifact, not the physical
endpoint.

**DERIVED NEGATIVE, scoped:** the same-orientation regular branch cannot
reach the desired `lambda=1` forward momentum and cannot produce a
nondegenerate next spatial frame.  A nonregular branch that bifurcates away
from this exact family is not excluded, but it would require a new selection
rule after loss of uniqueness.  The opposite temporal orientation is also a
separate branch and remains **OPEN**.

The homotopy parameter `lambda` is not physical time.  The identity proves a
lapse family, not emergent time.

## Relation to the original self-consistent theory

**DERIVED:** the normalization

```text
M=(90/pi) epsilon3 L
```

is an exact balance of intrinsic 600-cell curvature against dust on a static
product slab.  It selects a mass-to-radius relation but leaves `tau`
arbitrary.  The original self-referential construction can therefore remain
a consistency condition for one spatial slice; it does not close the
dynamical map or select a tick.

For a nonstatic slab with `L_old != L_new`, conserved `M` cannot cancel the
intrinsic term independently of the scale change.  Extrinsic-curvature terms
must enter.  This is the correct next route to a discrete Friedmann relation.
It does not follow that `tau` becomes physical: an exact reparameterization
symmetry could leave lapse gauge, while a finite Regge pseudo-constraint
could fix it as a discretization effect.  Those alternatives must be tested.

## Post-result primary-source audit

The general setting is known.  De Felice and Fabri evolve dust-filled
600-cell cosmologies and analyze a causality-breaking endpoint
(`<https://arxiv.org/abs/gr-qc/0009093>`,
`<https://arxiv.org/abs/gr-qc/0106077>`).  Liu and Williams report null-strut
endpoints in closed Regge FLRW models (`<https://arxiv.org/abs/1501.07614>`).
Canonical simplicial momenta and pseudo-constraints are treated by Dittrich,
Hoehn, Bahr and collaborators (`<https://arxiv.org/abs/1108.1974>`,
`<https://arxiv.org/abs/0912.1817>`,
`<https://arxiv.org/abs/0905.1670>`).

The search did not locate this complete finite-carrier factorization, the
eleven exact branch-controlled angle types or the per-edge momentum identity
for the present schedule pair.  External novelty remains **OPEN** pending a
dedicated literature review; a targeted search cannot prove novelty.

## Next gate

The next mission should preregister the first nonstatic homogeneous slab with

```text
M fixed,
L_old fixed at the published turning datum,
L_new != L_old,
the complete internal equations and canonical pre-momentum retained.
```

It must decide whether the equations select only the invariant velocity
`(L_new-L_old)/tau`, select the lapse itself through a pseudo-constraint, or
admit no connected Lorentzian solution.  Only after that distinction is made
can the result be called a discrete Friedmann tick or rejected as a lapse
artifact.

