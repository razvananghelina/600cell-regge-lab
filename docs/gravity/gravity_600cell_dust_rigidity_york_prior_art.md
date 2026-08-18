# Prior-art gate: rigidity image, self-stress carrier and a possible discrete York split

Date: 2026-08-18

## Exact object, carrier and hypotheses

Use only the literal regular 600-cell embedding

```text
x_v in S^3 subset R^4,   v=1,...,120,
```

its 720 edges, and the fixed logarithmic squared-edge carrier `Q=R^720`.
For an ambient displacement `z=(z_v)` define the rigidity tangent, up to one
common nonzero edge factor, by

```text
(R z)_uv = (x_u-x_v) dot (z_u-z_v).
```

Because every edge has the same length, this has the same image as the
variation of `log ell_uv^2`.

Split each ambient displacement uniquely into radial and tangential parts,

```text
z_v = sigma_v x_v + y_v,   x_v dot y_v = 0.
```

The radial part maps to the already certified conformal image:

```text
R(sigma_v x_v) = constant * C sigma.
```

Let `D` denote `R` restricted to the `360`-dimensional tangent-vertex
carrier.  Then

```text
im R = im C + im D.
```

The Euclidean edge-orthogonal complement

```text
S = ker R^T
```

is the equilibrium self-stress space of the embedded framework.  It is a
geometry-selected carrier, but the word "stress" here is framework
terminology; it is not matter stress-energy.

The proposed project-specific question is whether the committed centered
Regge recurrence dynamically decouples this carrier.  That is a separate
question from its dimension.

## What is already known before calculation

[Whiteley's infinitesimal-rigidity theorem for convex 4-polytopial
frameworks](https://doi.org/10.1090/S0002-9947-1984-0752486-6) applies because
the regular 600-cell is strictly convex and all its two-faces are triangles.
Consequently its Euclidean bar-and-joint framework is infinitesimally rigid.
The kernel of `R` consists only of four translations and six rotations, so

```text
rank R = 4*120 - 10 = 470,
dim ker R^T = 720 - 470 = 250.
```

Thus the numbers `470` and `250` are **KNOWN CONTROL**, not discoveries and
not evidence for gravity.

The same theorem fixes the radial/tangential ranks.  A rigid Euclidean motion
which is tangent to the unit sphere at all 120 spanning vertices cannot
contain a translation, so

```text
ker D = the six rotations,
rank D = 360 - 6 = 354.
```

Since `rank(C,D)=rank R=470` and `rank C=120`,

```text
dim(im C intersect im D) = 120 + 354 - 470 = 4.
```

These four directions are supplied explicitly by projecting the four ambient
translations tangentially to `S^3`; their lost radial parts are conformal.
Together with six rotations they are the ten reducibility directions familiar
from conformal Killing geometry of the round three-sphere.

Calling the remaining dimension

```text
720 - 470 = 250
```

"two polarizations per vertex plus ten" would be only a **PATTERN**.  In
rigidity theory it is simply left nullity/self-stress.  No physical quotient
has been derived from this count.

## Regge prior art and the load-bearing warning

- [Bahr--Dittrich](https://arxiv.org/abs/0905.1670) show that curved Regge
  solutions generally have no exact discrete diffeomorphism gauge symmetry;
  constraints become pseudo-constraints.
- [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817) derive exact local first-
  class vertex-displacement constraints only for linearized Regge calculus
  on a flat background and show that nonlinear corrections break them.
- [Hoehn](https://arxiv.org/abs/1411.5672) identifies lattice gravitons only
  after using those flat-background constraints and gauge-invariant curvature
  data.
- [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974) give the general
  action-generated pre/post-constraint formalism for simplicial evolution.

The present dust `S^3` background is curved.  Moreover, the repository has
already certified an invertible full boundary Legendre cross block in every
sector.  Therefore the fixed one-slab map has no action-derived primary
pre/post constraint whose kernel is `im D`.  Subtracting the 354 tangent
directions by hand and announcing 250 physical modes would be illegitimate.

## KNOWN / CONTROL / OPEN

- **KNOWN:** `rank R=470`, `rank D=354`, intersection dimension `4`, and
  self-stress dimension `250` follow from convex-framework rigidity.
- **CONTROL:** `C` is the unique literal conformal map and has exact rank
  `120`.
- **CONTROL:** the centered recurrence and all of its sector matrices were
  committed before this carrier comparison.
- **CONTROL / NEGATIVE:** the full discrete Legendre map is regular, so this
  coarse curved model has no exact canonical gauge quotient selected by a
  null cross block.
- **OPEN:** the exact sector multiplicities of `im R` and `ker R^T` under the
  frozen binary-tetrahedral action.
- **OPEN:** whether `ker R^T` is decoupled by `M,N,V`, or equivalently closed
  under the normalized recurrence operators `Gamma,Omega`.
- **OPEN:** whether any different action-derived pseudo-constraint carrier
  approaches `im D` under refinement.
- **OPEN:** any physical tensor/graviton interpretation.

## Framing attack and licensed test

The `250` count is forced by a theorem that knows nothing about the dust
action or evolution.  Reproducing it numerically would have zero evidential
weight.  The only new finite question worth computing here is dynamical
compatibility.

Construct `R` from the literal coordinates and freeze the unique self-stress
carrier `S=ker R^T` before loading the centered operators.  Then, in every
minimal sector and both schedules, report without fitting:

1. the known rigidity ranks as controls;
2. both cross blocks between `im R` and `S` for `M,N,V`;
3. leakage of `S` under `Gamma=M^-1N` and `Omega=M^-1V`;
4. schedule robustness and complete error envelopes.

Exact or error-consistent zero cross blocks would select a closed
stress/tensor-like recurrence sector.  Resolved leakage would be a clean
**STRUCTURAL NEGATIVE** for this particular York route on the fixed coarse
carrier.  Neither outcome alone supplies Hamiltonian/diffeomorphism
constraints, two physical polarizations, or a continuum graviton.
