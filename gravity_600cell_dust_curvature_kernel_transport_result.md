# Result: the curvature-kernel line is not transported to the next slab

Date: 2026-08-17

## Provenance

```text
prior-art gate                    ed51853
target-disclosed protocol         8fe4176
registered implementation         a1e0c94
passing artifact                  9f4d7f1
```

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_curvature_kernel_transport.py`.

Artifact:

```text
reproducible/gravity_600cell_dust_curvature_kernel_transport.json
SHA-256 d638977e8c44fb04278892b31169ae011f6d6eda97b20effa41673f48876ebfe
```

Only this targeted verifier and its direct 43-control geometry import were
run.  The full suite was not run.  Two complete runs were byte-identical:

```text
12/12 PASS
CURVATURE_KERNEL_TRANSPORT_REFUTED
```

## Background and seam control

The calculation uses the two previously accepted fixed-mass slabs, not a
fitted continuation.  The complete first-post/second-pre momentum seam has
maximum residual

```text
6.0468e-44
```

in each schedule, versus inherited bound `3.6514e-21`.  The second slab uses
the same conserved total dust mass; it does not recompute mass from its later
scale.

All base and displaced simplices retain the certified Lorentzian branch.  The
65-dimensional pre-Legendre matrices exclude zero in Flint ball arithmetic.

## A fresh second-slab kernel exists

Both schedules independently reproduce

```text
rank F_1 = 59, nullity F_1 = 1,
rank F_2 = 59, nullity F_2 = 1.
```

The restrictions to the fixed homogeneous plane have singular values

```text
slab 1: (1.458e-34, 0.19202874541896045),
slab 2: (1.353e-34, 0.38405877995082597).
```

Thus each full response has a unique homogeneous kernel.  Their normalized
coordinates in the preregistered `(q,p)` basis are

```text
K_1 = (-0.0034313802072921571, 0.9999941127976069258),
K_2 = (-0.0068627812467832903, 0.9999764508394979196).
```

The two independent schedule constructions identify `K_1`, `K_2` and
`T_1K_1` separately.  Their cross-schedule line distances are respectively
`3.90e-45`, `2.15e-45` and `3.10e-33`, all far inside their frozen errors.

**DERIVED COMPUTATIONAL:** the instantaneous one-dimensional curvature-blind
direction persists as a *count* on the second slab, but its direction changes
with the background.

The roughly doubled `q/p` slope is a **POST-RESULT PATTERN** only.  It was not
preregistered and is not used in the verdict.

## Decisive transport test

The target was disclosed before constructing `F_2`:

```text
T_1 K_1 ?= K_2.
```

For both schedules the projective separation is

```text
sin angle(T_1 K_1, K_2) = 0.0034314093660378607,
epsilon                  = 1.0042e-16.
```

The separation is approximately `3.42e13` uncertainty units.  As an
independent test, the freshly reconstructed second-slab curvature response is

```text
||F_2 T_1 K_1|| / ||T_1 K_1|| = 0.0013178628946323380,
epsilon                         = 7.205e-13,
```

which is nonzero by approximately `1.83e9` uncertainty units.

**DERIVED COMPUTATIONAL NEGATIVE:** the first slab's unique
curvature-preserving phase direction is not mapped into the second slab's
unique curvature-preserving direction.  This remains false after allowing the
kernel to depend on the evolving background.

## Consequence for the `119+1` interpretation

Two increasingly charitable hypotheses have now been falsified:

```text
fixed line:       T_1 K_1 != K_1,
moving line:      T_1 K_1 != K_2.
```

Therefore the one instantaneous curvature-kernel line cannot supply the
missing dynamically preserved `+1` in the observed `119+1=120` count.  That
numerical count remains an algebraic pattern, not a canonical decomposition
of the tangent dynamics.

This closes the specific proposal that the curvature kernel itself is the
homogeneous gauge/lapse/time mode.  The negative is stronger than mere failure
to hit an eigenvector: even the correctly time-dependent target fails.

## Interpretation firewall

- **DERIVED:** each of the first two slabs has one unique homogeneous
  curvature-blind tangent direction.
- **DERIVED NEGATIVE:** neither a fixed nor a transported bundle of those
  directions is preserved by the canonical evolution.
- **STRUCTURAL:** this is strong evidence that the line is an instantaneous
  degeneracy of the chosen internal-curvature diagnostic, not an exact
  canonical gauge generator.
- **OPEN:** `K_n` is a tangent-vector kernel, whereas a Hamiltonian constraint
  is a function/covector.  This calculation does not enumerate every possible
  pseudo-constraint.
- **OPEN:** physical scalar/vector/tensor classification of the remaining
  perturbations, matter perturbations, refinement and continuum propagation.

The curved-background literature is consistent with this outcome.  Exact
vertex-displacement constraints are controlled in flat linearized Regge
calculus, while curvature generically produces background-dependent
pseudo-constraints: <https://arxiv.org/abs/1411.5672>,
<https://arxiv.org/abs/0905.1670>, and
<https://arxiv.org/abs/0912.1817>.  No post-result search located the present
600-cell response kernel or its transport.  External novelty is **OPEN**.

## Correct next physical mission

Stop trying to promote the curvature kernel into a gauge/time direction.  The
background is a non-autonomous multi-tick trajectory, so individual one-step
eigenvalues are not physical modes either.  The next clean object is the
second slab's complete canonical tangent `T_2` and the two-step cocycle

```text
T_2 T_1
```

on the full 1,440-dimensional boundary phase carrier, blockwise in the seven
binary-tetrahedral sectors.  Its blind, schedule-controlled census would be
the first anisotropic propagation result across two dynamically solved slabs.
Only after that should curvature-changing sectors be compared with continuous
`S^3` perturbations.
