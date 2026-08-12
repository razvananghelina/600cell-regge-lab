# The exact all-form cone coefficient narrowly prefers the round endpoint

Date: 2026-08-12  
Protocol commit: `c3d2500`  
Registered verifier: `reproducible/verify_regge_de_rham_cone_selector.py`  
Machine-readable result: `reproducible/regge_de_rham_cone_selector.json`

## Headline

The analytic gate is passable, and the exact result goes the opposite way
from the naive smooth-Regge estimate:

> **DERIVED CONDITIONAL ENDPOINT SELECTION.**  For the ordinary Hodge--de
> Rham heat coefficient on the complete exterior algebra, after normalizing
> the two metrics to equal three-volume, the unit round `S3` endpoint has a
> strictly lower `A2` than the fixed piecewise-flat 600-cell endpoint.

Numerically,

```text
A2(round)                    = -78.9568352087148689507...
A2(fixed Regge, equal volume)= -78.8719985926940468466...
Regge - round                =  +0.0848366160208221041...
```

The relative separation is only `0.1074468%`, but its sign is stable at 40,
80 and 140 decimal digits.  It is not a floating-point crossing.

This result is evidence for spectral smoothing of the fixed endpoint.  It is
not yet emergent gravity: it compares one asymptotic coefficient at two
metrics, not a complete action on every metric, and supplies neither time nor
an absolute scale.

## 1. The operator and domain are the same analytic construction

The 600-cell boundary was rebuilt from its 120 quaternion vertices.  Its
complete f-vector is

```text
(vertices, edges, triangles, tetrahedra)=(120,720,1200,600).
```

Every open edge has normal link `C5`, so its neighborhood is the product of
the edge with a two-dimensional cone.  Every vertex link is the icosahedral
triangulation of `S2`, with f-vector `(12,30,20)` and vanishing middle
cohomology.  Thus the local link condition in Cheeger's Hodge theory is met:
no extra ideal-boundary datum is required.  The generalized Dirichlet and
Neumann Hodge extensions coincide, matching the closed-Hilbert-complex
extension selected by the existing fixed-Regge Whitney limit.

The protocol used “Friedrichs” as shorthand.  The precise certified statement
for forms is the one above: Cheeger's generalized Dirichlet and Neumann Hodge
extensions coincide.  Calling every degree separately “the Friedrichs
extension” is unnecessary and potentially ambiguous, so that stronger
terminology is not used in the result.

Cheeger's piecewise-flat expansion places the coefficient of
`t^(-3/2+2/2)=t^(-1/2)` on the `(3-2)=1` skeleton.  Therefore open edges
supply this `A2`; the vertex links first enter the constant term.  In
dimension three the higher-dimensional logarithmic complication does not
contaminate this order.  This closes the earlier concern that intersecting
edge strata might add an unknown term to the same coefficient.

The relevant primary result is [Cheeger, *Spectral Geometry of Singular
Riemannian Spaces*, Sections 7.1 and
7.5](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/cheeger.pdf).

## 2. Exact full-exterior conical coefficient

Let `beta` be the cone angle and `L` the length of a codimension-two
stratum.  For the scalar Hodge Laplacian the singular contribution is

```text
S(beta,L)=beta/6 * ((2*pi/beta)^2-1) * L.
```

For one-forms in ambient dimension three, Fursaev and Miele obtain

```text
V(beta,L)=3*S(beta,L)+2*(beta-2*pi)*L.
```

These are coefficients for the Hodge--de Rham operator, not for an unrelated
rough vector Laplacian.  Our cone satisfies `beta<2*pi`, the range treated
directly in the mode calculation rather than only by analytic continuation.
See [Fursaev and Miele, *Cones, Spins and Heat Kernels*, equations (i6a) and
(i10)](https://arxiv.org/abs/hep-th/9605153).

There is also an independent domain-sensitive recovery.  On a closed conic
two-sphere, Hodge decomposition for the Cheeger extension gives the exact
identity

```text
K_1(t)=2*K_0(t)+b1-2*b0=2*K_0(t)-2.
```

Conic Gauss--Bonnet converts the constant `-2` into the local one-form term
`2*S+2*(beta-2*pi)` at each tip.  Cheeger's Kunneth formula on
`C_beta x R` then adds the tangent scalar component and gives precisely
`V=3*S+2*(beta-2*pi)`.  Thus the coefficient is not being imported from a
different self-adjoint domain by analogy.

Hodge duality gives equal coefficients in degrees `0/3` and `1/2`.  The
ordinary full-exterior contribution is consequently

```text
C_full(beta,L)
  = 2*S + 2*V
  = [16*pi^2/(3*beta) + 8*beta/3 - 8*pi] * L.
```

Writing the deficit as `delta=2*pi-beta` exposes the crucial correction:

```text
C_full
  = -(4/3)*delta*L + 4*delta^2*L/(3*beta).
```

The first term is exactly what the smooth distributional-curvature rule

```text
A2=-(2/3) integral R,
integral R=2*sum(delta_e L_e)
```

would give.  The second term is positive and genuinely conical.  It cannot be
discarded at the finite 600-cell deficit.

## 3. Frozen 600-cell comparison

The exact metric data are

```text
edge length = 1/phi,
beta        = 5*arccos(1/3),
delta       = 2*pi-beta,
V_Regge     = 600*(1/phi)^3/(6*sqrt(2)),
V_round     = 2*pi^2.
```

All 720 coordinate edges and all 600 facet volumes were recomputed from the
quaternion coordinates.  The Regge coefficient scales as length, so its
equal-volume value is multiplied by

```text
c=(V_round/V_Regge)^(1/3)=1.0574728352587844587...
```

The exact cone formula then gives the headline values.  A particularly useful
hostile control is what happens if the quadratic cone term is omitted:

```text
linearized smooth-Regge A2 = -80.5523086694673133363...
round A2                   = -78.9568352087148689507...
exact conical Regge A2     = -78.8719985926940468466...
```

The linear approximation says that Regge wins by `1.59547...`; the exact
conical correction adds `1.68031...` and reverses the ordering.  This is a
direct demonstration that applying the smooth `integral R` formula to the
singular endpoint would have produced the wrong physical narrative.

## 4. What this establishes about emergent gravity

There are now two compatible results for the same ordinary all-form
coefficient:

1. among every smooth left-invariant metric of fixed volume, the round metric
   is the unique global minimum;
2. against the fixed 600-cell Regge endpoint at the same volume, it also has
   the lower value.

This is a **STRUCTURAL POSITIVE** for an emergent-gravity route: the theory's
geometric operator supplies a curvature-sensitive term that removes both
homogeneous smooth anisotropy and the particular polyhedral endpoint.  The
term behaves like a smoothing energy rather than a fitted endpoint score.

The stronger statements remain **OPEN**:

- the interior of the full round--Regge metric family has not been evaluated;
- `A2` alone is not a complete spectral action and its physical sign still
  has to be supplied by that action;
- no cutoff function or dimensionless heat parameter is selected;
- no absolute length, Newton constant or Planck scale is derived;
- no Lorentzian evolution, Hamiltonian constraint, diffeomorphism gauge
  redundancy, massless graviton or universal coupling to stress energy has
  been obtained.

Accordingly the correct headline is not “gravity derived”.  It is:

> the first common spectral curvature functional now agrees on which of the
> tested smooth and singular shapes is preferred, but the dynamics and scale
> that would turn that preference into gravity are still missing.

## 5. Status ledger

| Claim | Status |
|---|---|
| 600-cell open-edge cone angle is `5*acos(1/3)<2*pi` | **DERIVED** |
| Every edge link is `C5`; every vertex link is icosahedral `S2` | **DERIVED** |
| Vertex strata cannot alter the `t^(-1/2)` coefficient | **DERIVED, using Cheeger's heat theorem** |
| Full ordinary cone term is `2*S+2*V` | **DERIVED, using the primary Hodge coefficients** |
| Exact-minus-linear term is `4*delta^2*L/(3*beta)>0` | **DERIVED** |
| Equal-volume round endpoint has lower `A2` than fixed Regge | **DERIVED CONDITIONAL ENDPOINT SELECTION** |
| The smooth Regge curvature sum is accurate enough for this ordering | **REFUTED** |
| Round is the minimum over every interior `g_u` | **OPEN** |
| This coefficient alone defines a physical gravitational action | **OPEN** |
| Gravity, `G`, Planck scale or Lorentzian time are derived | **OPEN** |
