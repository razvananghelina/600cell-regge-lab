# A canonical shape-regular refinement exists

Date: 2026-08-11

Preregistration commit: `58fa9fc`

Targeted verifier:
`reproducible/verify_whitney_rank_edgewise_refinement.py`

Targeted result: **9/9 PASS**.  The verifier is registered exactly once.  The
full suite was not run, by explicit user request.

## Headline

The apparent choice between canonicity and mesh quality was too narrow.  A
third construction passes both requirements:

\[
 K_n=\operatorname{Esd}_{2^n}(\operatorname{sd}K).
\]

Apply barycentric subdivision once, order every resulting chamber by the
dimensions `(0,1,2,3)` of the faces in its flag, and thereafter use direct
edgewise subdivisions at resolutions `1,2,4,8,...`.

> **DERIVED POSITIVE:** on the piecewise-Euclidean 600-cell carrier, this is a
> canonical, conforming, nested, uniformly shape-regular refinement tower.

This result repairs the numerical continuum carrier.  It is not yet a
physical continuum limit and selects no dimensional constant.

## Why the order is geometric rather than fitted

A chamber of the barycentric subdivision is a complete flag

\[
 \sigma_0\subset\sigma_1\subset\sigma_2\subset\sigma_3,
 \qquad \dim\sigma_i=i.
\]

Its vertices are the barycentres of those four faces.  Their order is fixed by
face dimension.  Every simplicial automorphism preserves inclusion and face
dimension, so it preserves this order.  There is no global vertex numbering,
edge selection, octahedral diagonal, or numerical optimization in the rule.

The exact verifier constructs all 24 barycentric chambers of one tetrahedron,
applies the rule, and checks the resulting complex under all 24 elements of
`S4`:

| edgewise resolution | top tetrahedra in one subdivided tetrahedron | failed `S4` actions |
|---:|---:|---:|
| 2 | 192 | 0 |
| 3 | 648 | 0 |

The count is `24 k^3`, as required.  The local argument transfers to every
600-cell symmetry because every 600-cell automorphism restricts to a
simplicial map on its tetrahedral cells.

## Why it is uniformly shape regular

For a regular reference tetrahedron, the rank-ordered barycentric chamber is
an exact orthoscheme.  Its consecutive shape vectors have Gram matrix

\[
 Q^*Q=\operatorname{diag}\left(2,\frac23,\frac13\right).
\]

The independently implemented color-scheme construction verifies exactly
that every edgewise child has consecutive shape vectors

\[
 \frac1k(q_{\pi(1)},q_{\pi(2)},q_{\pi(3)})
\]

for a permutation `pi`.  Therefore scale changes with `1/k`, but shape can
take only finitely many values.  Reversal gives a congruent tetrahedron, so
there are at most `3!/2 = 3` congruence classes for every `k`.

This is an all-level proof, not an extrapolation.  The finite enumeration is a
control:

| `k` | top cells | vertices | strict interior overlaps | normalized shape classes |
|---:|---:|---:|---:|---:|
| 1 | 1 | 4 | 0 | 1 |
| 2 | 8 | 10 | 0 | 3 |
| 3 | 27 | 20 | 0 | 3 |
| 4 | 64 | 35 | 0 | 3 |

The minimum mean-ratio quality over the three classes is approximately
`0.524454`; it does not decay with level.  The maximum reported affine
condition is approximately `7.10073`.  These numbers diagnose the fixed
classes; no threshold was selected from them.

This agrees with the primary theorem of Edelsbrunner and Grayson, which proves
face compatibility, composition in scale, and at most `d!/2` congruence
classes independent of resolution: [Edgewise Subdivision of a
Simplex](https://doi.org/10.1007/s004540010063).

## Conformity and nesting

The face order induced by ranks is the same from either tetrahedron adjacent
to a shared face.  Exact two-tetrahedron controls give identical boundary
triangulations:

| `k` | boundary triangles from left/right |
|---:|---:|
| 2 | 24 / 24 |
| 3 | 54 / 54 |

For nesting, every tetrahedron of `Esd_4` lies in exactly one tetrahedron of
`Esd_2`.  Each of the eight coarse tetrahedra contains exactly eight fine
ones, with exact rational containment and volume checks.  In general the
edgewise composition theorem gives

\[
 \operatorname{Esd}_{\ell}
 \operatorname{Esd}_{k}
 =\operatorname{Esd}_{k\ell}.
\]

Thus the powers-of-two sequence is a genuine nested tower.  On the full
600-cell its number of top tetrahedra is

\[
 600\cdot24\cdot(2^n)^3=14400\cdot8^n.
\]

## Hostile controls: why the simpler routes fail

### Direct edgewise subdivision

Without the one-time barycentric rank, ordered `Esd_2` has exactly three
distinct tetrahedral subdivisions.  The rotational stabilizer `A4` acts on
each with orbit size three; no variant is fixed.

> **DERIVED NEGATIVE:** direct edgewise refinement of an unranked regular
> tetrahedron is shape regular but not selected by the tetrahedral geometry.

This reproduces, rather than hides, the ordering freedom described in the
primary paper.

### Symmetrizing the three midpoint choices

Adding the tetrahedron centroid and coning it to the eight faces of the
central octahedron gives a fully `S4`-invariant 12-split.  It contains a
repeatable central child with affine transform

\[
 T=\frac14
 \begin{pmatrix}
 1&-1&-1\\
 -1&1&-1\\
 -1&-1&1
 \end{pmatrix},
 \qquad
 \sigma(T)=\left\{\frac12,\frac12,\frac14\right\}.
\]

Hence the condition number of `T^n` is exactly `2^n`.

> **DERIVED NEGATIVE:** the obvious fully symmetric 12-split is canonical but
> not uniformly shape regular.

### Bisection

All six edges of a regular tetrahedron are tied, and the tetrahedral rotation
group is transitive on them.  Longest-edge or newest-vertex bisection can
control shapes, but on this carrier its initial edge or vertex is extra data.
This is a structural exclusion of an unlabelled canonical rule, not a
criticism of bisection algorithms.  Shape-control results for such labelled
algorithms are standard; for example, Liu and Joe prove a level-independent
quality bound for their bisection scheme: [Quality Local Refinement of
Tetrahedral Meshes Based on Bisection](https://doi.org/10.1137/0916074).

## Attack on the framing

The earlier framing suggested that symmetry and shape regularity might be
intrinsically incompatible.  That was false.  It considered rules applied
identically to every unlabelled child and missed a functorial source of labels:
the dimensions in a face flag.

The construction does add one derived layer, `sd K`.  This is defensible here
because it is functorial, contains no choice, and the project already uses the
barycentric carrier.  If a future physical axiom forbids changing from `K` to
`sd K`, this route would be excluded by that new axiom; no such axiom is
currently derived.

The proof is for the affine, piecewise-Euclidean realization used by the
Whitney operators.  It does not by itself prove an identical constant after a
nonlinear radial projection onto the round `S^3`.  Such a projection is a
separate metric choice and must be tested separately if introduced.

## What this advances

- **DERIVED:** the refinement rule is selected by face ranks.
- **DERIVED:** exact `S4` equivariance and shared-face conformity.
- **DERIVED:** an all-level finite-shape bound, hence uniform shape regularity.
- **DERIVED:** a nested powers-of-two tower.
- **DERIVED NEGATIVE:** direct edgewise subdivision has three unselected
  variants.
- **DERIVED NEGATIVE:** the canonical symmetric 12-split degenerates.
- **STRUCTURAL:** this removes mesh-shape degeneration as a contaminant in the
  next Whitney scaling experiment.

## What remains open

- **OPEN:** assemble the exact Whitney metric and trace stiffness on this new
  tower and test spectral convergence.
- **OPEN:** determine whether degreewise stiffness ratios approach a common
  factor once shape degeneration is removed.
- **OPEN:** whether any continuum limit is physically selected rather than
  merely mathematically available.
- **OPEN:** an absolute length, tick, mass, inertia, Lorentzian signature,
  causal speed, Planck scale, or quantum dynamics.

In particular, this is useful mathematical infrastructure, not yet a new law
of physics.  It gives the theory a clean place to ask its next question.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_rank_edgewise_refinement.py
```

Expected result: `9/9`.

