# Result: local dust weights are canonical inside the P1 ansatz

Date: 2026-08-19

## Headline

On both canonical projected carriers, continuous nodal `P1` dust has unique
vertex-only affine-exact local weights

```text
w_v = sum_(tetrahedra t incident on v) Vol(t)/4.
```

After normalization to the already blind-selected total dust mass,

```text
m_v = M*w_v/sum_u w_u,
```

all local masses are positive, exactly conservative, cubically covariant under
spatial scaling and equivariant under all 241 certified spatial actions.

This is **KNOWN / DERIVED CONDITIONAL, COMPUTATIONALLY VERIFIED AND
ADVERSARIALLY CORROBORATED**.  The condition is load-bearing: the continuous
nodal `P1` density ansatz is **STRUCTURAL**, not derived from the 600-cell.

## Provenance

| stage | commit |
|---|---|
| prior-art and framing gate | `ba7de6c` |
| primary protocol | `9156c0f` |
| registered primary verifier | `694e97c` |
| frozen primary artifact | `de0c220` |
| adversarial protocol | `0257417` |
| registered adversarial verifier | `ffacaf5` |
| frozen adversarial artifact | `f0e3e29` |

Primary artifact:

```text
reproducible/gravity_600cell_projected_rank_edgewise_local_dust.json
SHA-256 53463e5271301ae41eb26564875d26991ddea8024a9e09ae3c302d428ad39779.
```

Adversarial artifact:

```text
reproducible/gravity_600cell_projected_rank_edgewise_local_dust_adversarial.json
SHA-256 7e3466504eddb98ca464d2b8b3b9b64e8b256f00caa5f7af5e0324a85c41aa15.
```

The targeted verifiers pass `11/11` and `11/11`.  No full suite was run.

## 1. Exact conditional theorem

On the reference tetrahedron, the four barycentric basis functions integrate
exactly to

```text
(1/24,1/24,1/24,1/24) = Vol(t)*(1/4,1/4,1/4,1/4).
```

A linear quadrature supported only at the four vertices is exact on every
affine density if and only if it is exact on these four basis functions.
Therefore its four coefficients are uniquely `Vol(t)/4`.

The verifier reconstructs this with exact SymPy rational integration and an
exact four-equation solve.  This is a standard `P1` finite-element fact; no
external novelty is claimed.

## 2. Finite canonical carriers

The primary assembly gives:

| carrier | vertices | weight range | H4 orbits |
|---|---:|---:|---:|
| `P(sd K_600)` | 2,640 | `0.00398915` -- `0.0398915` | 4 |
| `P(Esd_2(sd K_600))` | 19,680 | `0.000503286` -- `0.00515803` | 10 |

The four base orbits have sizes

```text
120, 600, 720, 1200,
```

the four face-rank populations of the barycentric chamber complex.  The fine
carrier has orbit sizes

```text
120, 600, 720, 1200, 1440,
2400, 2400, 3600, 3600, 3600.
```

These counts were printed target-free.  They are geometry diagnostics, not
particle multiplets.

For the selected unit-volume-radius backgrounds:

```text
P(sd K_600):
  sum w = 19.147932918312847
  sum m = 2.365802630146490
  m_v in [0.000492876,0.00492876]

P(Esd_2(sd K_600)):
  sum w = 19.583480465413963
  sum m = 2.358675470499193
  m_v in [0.0000606168,0.000621244].
```

Maximum weight residuals over the 241 left/right/conjugation actions are
`1.74e-12` and `2.29e-13`.  Cubic scaling residuals are at most `1.30e-15`.
For homogeneous proper time, the assembled action collapses to

```text
-8*pi*sum_v m_v*tau = -8*pi*M*tau
```

within `1.88e-16` relative.

## 3. Independent consistent-matrix audit

The adversarial path uses NetworkX maximal cliques, Cayley--Menger volumes and
the consistent local `P1` mass matrix

```text
M_t = V_t/20 * (I + all_ones).
```

Its symbolic row sums are exactly `(V_t/4,V_t/4,V_t/4,V_t/4)`.  Global sparse
row sums reproduce the primary pointwise weights to

```text
1.53e-15  on P(sd K_600),
1.73e-15  on P(Esd_2(sd K_600)).
```

Both global matrices are symmetric at roundoff and all local eigenvalues are
positive.  The independently reconstructed Gram-weight byte digests equal the
frozen primary digests exactly.

## 4. Why symmetry alone does not select the weights

The negative control assigns the globally uniform value

```text
sum_t Vol(t) / number_of_vertices
```

to every vertex.  It is positive, exactly conserves total volume and respects
every spatial symmetry.  Nevertheless it fails the exact integral of a global
`P1` hat function by as much as

```text
81.8% on the base carrier,
97.7% on the fine carrier.
```

Thus positivity, conservation and `H4` equivariance do not suffice.  Exact
affine `P1` quadrature is the substantive selection hypothesis.

## 5. Scientific status

- **KNOWN / DERIVED CONDITIONAL:** `Vol(t)/4` is the unique vertex-only linear
  quadrature exact on affine tetrahedral densities.
- **DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED:** the assembled finite
  weights are positive, conservative and spatially equivariant on both
  canonical carriers.
- **STRUCTURAL:** representing dust as a continuous comoving nodal `P1`
  density.
- **OPEN:** why this approximation class, rather than discontinuous,
  higher-order, circumcentric, optimized-dual or independently specified
  particle masses, is physically selected.
- **OPEN:** density perturbations and their initial data.
- **OPEN:** a canonical Lorentzian carrier supporting independent vertex
  lapses on the refined mesh.

The post-result literature check confirms rather than elevates the result.
Dittrich, Gielen and Schander use fixed-mass particles with action `-m` times
proper time in simplicial cosmology (arXiv:`2109.00875`, DOI
`10.1088/1361-6382/ac42ad`).  Barycentric tetrahedral mass lumping and its
alternatives are standard; Jacobson explicitly compares them in
arXiv:`2406.08647`, DOI `10.1111/cgf.15133`.

No located source selects this `P1` ansatz as fundamental 600-cell matter.
Search absence is not a novelty proof.

## 6. Consequence for the gravity programme

One ambiguity is now bounded: if we use continuous nodal `P1` dust, local
masses no longer carry adjustable coefficients.

The next blocker is temporal, not material.  A local lapse equation requires a
four-dimensional Lorentzian slab in which same-vertex strut lengths may vary
independently.  The homogeneous cellular frustum used for the acceleration
calibration has only one common `rho`; inserting one `rho_v` into that formula
would be mathematically invalid.

Before forming any refined Hessian, the next mission must determine whether
the rank/chromatic structure of the canonical spatial carrier selects a
globally compatible staircase triangulation of `K x I`, or whether temporal
scheduling freedom returns.  If no choice-free slab exists, the local
gravity route remains blocked even though the dust weights are clean.
