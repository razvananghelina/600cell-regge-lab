# The dual constraint hierarchy stays bounded on three edgewise levels

Date: 2026-08-12

Preregistration commit: `36a56b7`

Inherited dual-resolution result commit: `799966f`

Targeted verifier:
`reproducible/verify_whitney_dual_resolution_edgewise.py`

Targeted result: **9/9 PASS**.  The verifier is registered exactly once.  No
spectrum, metric or phenomenological target was used.  The full suite was not
run, by explicit user request.

## Headline

The locality failure seen under repeated barycentric subdivision does not
occur on the three frozen levels of the actual rank-edgewise control tower.

For

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}\partial\Delta^4),
 \qquad k=1,2,4,
\]

the complete dual-cell constraint hierarchy has the identical maximum
incidence data at every level:

\[
 \boxed{a_0=24,\qquad a_1=6,\qquad r_3=14.}
\]

Here `a0` is the maximum tetrahedron multiplicity of a vertex, `a1` the
maximum tetrahedron multiplicity/link-cycle length of an edge, and `r3` the
maximum number of dual faces in a vertex dual-volume boundary.

> **PATTERN TOWARD BOUNDED DUAL LOCALITY:** every preregistered incidence
> maximum is nonincreasing from `k=2` to `k=4`; in fact all three are exactly
> constant from `k=1` through `k=4`.

This is a strong finite structural positive, not yet an all-`k` theorem.

## Why the earlier negative does not transfer

The previous audit compared the 600-cell directly with its first
barycentric subdivision.  There the maximum vertex multiplicity grew
`20 -> 120`, edge-cycle length `5 -> 10`, and dual-volume boundary
`12 -> 62`.

That transition is not the chosen refinement law.  The chosen tower performs
one barycentric subdivision to fix a canonical rank order and thereafter
uses direct edgewise refinement.  The new audit tests precisely this second
mechanism.  Treating the barycentric negative as a verdict on the edgewise
tower would therefore have been a category error.

## Exact link gates

The dual resolution from commit `799966f` transfers only if the local links
have the required topology.  The verifier exhausts every link at all three
levels and proves combinatorially:

- every triangle has exactly two parent tetrahedra;
- every edge link is one connected cycle;
- every vertex link has connected one-skeleton, every link edge has two link
  faces, and Euler characteristic exactly two;
- every tetrahedron-occurrence graph for every simplex is connected;
- neighbour-constraint occurrence degrees remain exactly `(3,2,1,0)`.

Thus every triangle has a dual interval, every edge a dual disk, and every
vertex a dual 3-ball.  Under the already-proved signed dual-cell theorem, the
complete relation hierarchy is exact on each frozen control.  **DERIVED ON
THE CONTROLS.**

No floating rank or support threshold is involved.

## Complete finite data

| `k` | f-vector | max tetrahedra/vertex | max tetrahedra/edge | max edges/vertex |
|---:|---:|---:|---:|---:|
| 1 | `(30,150,240,120)` | 24 | 6 | 14 |
| 2 | `(180,1140,1920,960)` | 24 | 6 | 14 |
| 4 | `(1320,9000,15360,7680)` | 24 | 6 | 14 |

The full vertex multiplicity histograms are:

| `k` | histogram: tetrahedra per vertex |
|---:|---|
| 1 | `12:20, 24:10` |
| 2 | `12:20, 16:30, 24:130` |
| 4 | `12:20, 16:90, 24:1210` |

The edge multiplicity/link-cycle histograms are:

| `k` | histogram: tetrahedra per edge |
|---:|---|
| 1 | `4:90, 6:60` |
| 2 | `4:540, 6:600` |
| 4 | `4:3960, 6:5040` |

The vertex degrees are:

| `k` | histogram: incident edges per vertex |
|---:|---|
| 1 | `8:20, 14:10` |
| 2 | `8:20, 10:30, 14:130` |
| 4 | `8:20, 10:90, 14:1210` |

More strongly, only three vertex-link f-vectors occur through `k=4`:

\[
 (8,18,12),\qquad(10,24,16),\qquad(14,36,24),
\]

and only edge cycles of length four or six occur.  This is exact finite
enumeration, not a sampled neighbourhood census.

## Relation to the edgewise construction

The original Edelsbrunner--Grayson construction subdivides every
`d`-simplex into `k^d` equal-volume simplices, uses only finitely many shape
classes in fixed dimension, agrees on shared faces, makes repeated
subdivision equivalent to increasing `k`, and has translation-equivalent
interior vertex neighbourhoods.  These properties explain why a uniform
local bound is plausible for a fixed finite coarse complex:

- Herbert Edelsbrunner and Daniel R. Grayson,
  [*Edgewise Subdivision of a Simplex*](https://doi.org/10.1007/s004540010063),
  Discrete & Computational Geometry 24 (2000), 707--719.

But the paper's general properties do not by themselves prove the sharp
global bounds `(24,6,14)` after all coarse-face gluings.  A self-contained
all-`k` link-type classification for this particular rank-ordered carrier is
still required before changing the label from **PATTERN** to **DERIVED
UNIFORM**.

## What this advances physically

The result rescues one necessary locality condition:

> the complete canonical constraint-redundancy hierarchy can be represented
> with a bounded-looking incidence stencil on the first three levels of the
> selected refinement family.

It does not supply the missing dynamics.  The copy constraints are still
second class on the original phase space.  The hierarchy contains no
auxiliary symplectic bracket, positive metric, BRST Hamiltonian or
gauge-invariant physical embedding.  Exact physical dressing in the minimal
conversion still contains the global Gram inverse.

Therefore this is a kinematic locality positive, not a derivation of a causal
unitary tick.

## Status ledger

- **DERIVED:** every frozen edge and vertex link passes the exact topology
  gates.
- **DERIVED:** the complete signed dual resolution transfers to all three
  controls.
- **DERIVED:** maxima `(24,6,14)` are identical at `k=1,2,4`.
- **DERIVED:** only three vertex-link types and two edge-link lengths occur on
  those controls.
- **PATTERN TOWARD BOUNDED DUAL LOCALITY:** no maximum grows from `k=2` to
  `k=4`.
- **OPEN:** a proof and preferably sharp bounds for every integer `k`.
- **STRUCTURAL POSITIVE:** choice-free bounded-stencil kinematics on the
  tested tower.
- **STRUCTURAL NEGATIVE:** no physical first-class gauge theory or local
  Hamiltonian follows from this incidence result.
- **NOT CLAIMED:** spectrum, time, causality, inertia, mass or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_dual_resolution_edgewise.py
```

Expected result: `9/9`.

## Subsequent all-resolution upgrade

The open all-`k` classification in this note is closed by
`whitney_dual_resolution_all_k_result.md`.  Using the finite partition
classification of edgewise links together with the exact barycentric
coface-gluing multiplicities, it proves the sharp identities

\[
 a_0(q)=24,\qquad a_1(q)=6,\qquad r_3(q)=14
\]

for every positive integer `q`.  The finite **PATTERN** label above is
therefore superseded by **DERIVED UNIFORM** under the hypotheses stated in
that subsequent note.
