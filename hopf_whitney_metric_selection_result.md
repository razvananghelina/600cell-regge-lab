# Symmetry and Whitney refinement do not select round versus Regge metric

Date: 2026-08-12  
Protocol commit: `8565d19`  
Registered verifier: `reproducible/verify_hopf_whitney_metric_selection.py`  
Machine-readable result: `reproducible/hopf_whitney_metric_selection.json`

## Headline

**DERIVED METRIC-SELECTION NO-GO.**  The positive Hopf curvature response was
derived on the unit round metric `g_0`, whereas the exact Whitney tower was
certified for the pushed-forward piecewise-flat Regge metric `g_R`.  These
metrics are genuinely different, but both satisfy every frozen kinematic
selection criterion:

- full 600-cell symmetry;
- construction from the existing radial identification, without labels;
- positivity and face compatibility;
- the same Whitney subcomplex and exact refinement maps;
- uniform, level-independent metric equivalence.

In fact all metrics in the injective continuum

```text
g_u=(1-u)g_R+u g_0,   0<=u<=1,
```

pass those criteria.  Symmetry plus refinement therefore does not select the
metric on which an emergent gravitational action should live.

This does not kill emergent gravity.  It kills the inference that the round
Hopf heat Hessian has already transferred to the exact fixed-Regge Whitney
theory.

## 1. The two endpoints are not the same metric

Let `R(x)=x/|x|` identify the polyhedral boundary with `S3`.  On a flat facet,

```text
R*g_R = g_flat,
R*g_0 = I/r^2 - x x^T/r^4.
```

The exact squared support distance of a unit-circumradius 600-cell facet is

```text
a^2=(7+3 sqrt(5))/16.
```

Relative to the flat facet metric, the three eigenvalues of `R*g_0` are

```text
a^2/r^4, 1/r^2, 1/r^2.
```

At the midpoint of the facet-centroid-to-vertex segment, `r^2>a^2`, so the
parallel and transverse eigenvalues are exactly unequal.  Thus `g_0` and
`g_R` are not even pointwise proportional there.  This is not merely a
different global normalization.

Their local curvature also differs:

```text
R(g_0)=6,
R(g_R)=0 in every open tetrahedral facet.
```

The global Regge curvature is concentrated on the lower-dimensional
skeleton.  Consequently a smooth round heat coefficient cannot be assigned
to `g_R` by changing coordinates.

## 2. Both endpoints have the required symmetry

For every orthogonal transformation `Q`, radial normalization obeys

```text
R(Qx)=Q R(x).
```

The verifier establishes the stronger `O(4)` covariance identity

```text
M(Qx)=Q M(x) Q^T,
M(x)=I/r^2-xx^T/r^4,
```

using exact plane rotations and reflection.  Hence the round pullback is
invariant under the full orthogonal 600-cell symmetry group.  The flat
facet metric is invariant under the same group, so its pushforward is also
invariant.

On a shared triangular face, the flat tangential Gram matrix is the same
from either regular tetrahedron.  The round pullback uses only the ambient
point and face-tangent vectors and contains no parent-facet datum.  Both
metrics therefore satisfy the Regge tangential compatibility condition, as
does every convex combination.

This sharpens an earlier phrase.  Averaging the six Hopf projectors selects
an isotropic tensor within that local anisotropy ansatz.  It does **not**
distinguish the globally round metric from an `H4`-invariant polyhedral Regge
metric.  Both are locally isotropic in their own natural frames and globally
symmetry-equivariant.

## 3. There is a whole admissible metric continuum

The exact round/flat comparison gives the uniform bounds

```text
(7+3 sqrt(5))/16 <= g_0/g_R <= 28-12 sqrt(5).
```

Since the interval contains one, every `g_u` is positive and obeys the same
bounds.  At the explicit interior point `g_0-g_R` is nonzero, so

```text
g_u=g_v  implies  u=v.
```

The family is therefore genuinely continuous and not a list of coordinate
descriptions of one metric.  No value of `u` was compared with an action or
physical target.

Calling this entire family canonical would be misleading.  The family is
canonically available; no individual interior `u` is selected.  Its role is
to prove that the frozen criteria are insufficient.

## 4. Whitney induction cannot break the tie

The verifier independently subdivides a tetrahedron across an edge and uses
a generic positive diagonal metric

```text
G=diag(mx,my,mz).
```

It constructs the exact Whitney Gram matrices, including the metric-induced
inner products on every exterior degree.  With the metric-independent
inclusions obtained from barycentric minors, it proves symbolically

```text
P_p^T M_f,p(G) P_p=M_c,p(G),   p=0,1,2,3.
```

Thus exact Galerkin refinement isometry is not a selector for the metric.  It
holds because the inherited coarse Whitney form is literally the same
differential form integrated with the same metric on the same domain.

The middle-degree mass matrices for `diag(1,1,1)` and `diag(2,1,1)` are not
proportional, so the family is not invisible to the operator.  Positivity is
checked on a generic anisotropic control.  Refinement compatibility and
operator equality are distinct statements.

## 5. Norm equivalence is not curvature transfer

Uniform equivalence was sufficient for the earlier FEEC compactness and
spectral-convergence theorem.  It is not sufficient to identify adjoints,
spectra or heat coefficients.  The verifier includes an exact finite
control where two uniformly equivalent metrics give different minimal right
inverses of the same differential.

Therefore the valid earlier theorem remains:

> the exact Whitney spectra converge to the fixed-Regge continuum.

It does not become:

> the fixed-Regge spectra or curvature Hessian equal those of the round
> continuum.

The latter statement is refuted by the endpoint curvature difference before
any spectral calculation.

## 6. Attack on the framing

There is a clean conditional fork.

1. If exact flatness on every original tetrahedron and the existing affine
   mass matrices are declared fundamental, then `g_R` is selected by that
   extra premise.  The round Hopf heat Hessian does not apply, and a separate
   Regge/singular-metric response must be derived.
2. If the ambient unit sphere and restored local isotropy license `g_0`, then
   a canonical round Whitney construction exists, but the kinematic axioms
   do not explain why nature chooses it instead of `g_R` or `g_u`.

The only honest way to break this fork is an independently selected dynamic
functional or an additional metric axiom frozen before inspecting the
desired gravitational response.

## 7. Status ledger

| Claim | Status |
|---|---|
| `g_R` and `g_0` are exactly distinct | **DERIVED** |
| Both are full-symmetry equivariant | **DERIVED** |
| Both are face-compatible and uniformly equivalent | **DERIVED** |
| `g_u`, `0<=u<=1`, is a nontrivial admissible continuum | **DERIVED** |
| Exact Whitney refinement isometry selects one metric | **REFUTED** |
| Norm equivalence transfers adjoints/heat Hessians | **REFUTED** |
| Round Hopf heat Hessian is already a fixed-Regge result | **REFUTED** |
| Current exact Whitney theory uses `g_R` | **DERIVED SCOPE** |
| A round-metric Whitney theory is canonically constructible | **STRUCTURAL POSITIVE** |
| An existing axiom selects round over Regge | **OPEN / absent** |
| An independently derived action could select a metric | **OPEN** |
| Emergent gravity is impossible | **NOT CLAIMED** |

## 8. Next strict gate

The next noncircular question is dynamical:

> Does one already licensed ordinary spectral functional select an endpoint
> or a unique `u` throughout the refinement tower, with the same answer over
> its full allowed scale range?

A single chosen heat time or one finite level would be a fit.  The admissible
family and complete heat-scale interval must be registered first.  If no
functional is already licensed, inventing one after seeing this ambiguity
would merely rename the missing metric axiom.

Even a positive Riemannian selector would still be static.  Lorentzian time,
gauge constraints and universal source coupling remain later gates.

## Reproduction

Run only the targeted verifier:

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_hopf_whitney_metric_selection.py
```

Result: `18/18` checks passed.  The full suite was not run under the user's
current instruction.
