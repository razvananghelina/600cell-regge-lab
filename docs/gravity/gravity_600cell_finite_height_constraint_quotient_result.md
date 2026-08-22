# Bounded no-go: no exact local constraint quotient on the accepted two-slab map

Date: 2026-08-22

## Verdict

**DERIVED EXACT/COMPUTATIONAL, BOUNDED NO-GO UNDER THE COMPLETE SCOPE
BELOW.**

On the frozen branch-B history at incoming representative `v=3/2`, the first
and second finite-height Regge-plus-conserved-dust slabs define regular local
canonical maps on the complete `1440`-dimensional boundary phase carrier in
both staircase parities.  Their local pre- and post-images are open and have
constraint codimension zero.  Therefore the current action does not select a
nontrivial exact local boundary gauge quotient through pre/post-Legendre
degeneracy on these two slabs.

The targeted verifier returned

```text
FINITE_HEIGHT_TWO_STEP_EXACT_CONSTRAINT_QUOTIENT_BOUNDED_NO_GO
14/14 PASS
```

Artifact SHA-256:

```text
bb3735fd475c08068f8240ed1f85a0b205a529228cb46f2d6d8a788114b9e49f
```

Verifier SHA-256:

```text
1dd2a14559c11b65b012528f772ff12823188a7396aaf3f1a5ceb474d334d78d
```

Only the targeted verifier was run.  The full suite was not run.

## Frozen provenance

| stage | commit |
|---|---|
| prior-art and repository gate | `cbbf7b6` |
| target-free theorem protocol | `286c34b` |
| registered verifier before first execution | `de719a2` |
| first result artifact | `2bbc4d7` |

No tangent eigenvalue, desired physical mode count, continuum dispersion
target, limiting speed, `G`, Planck unit or particle datum was inspected.

## Complete hypotheses

The no-go requires all of the following:

1. the fixed regular 600-cell and its complete one-slab carrier of `720` old
   boundary, `840` internal and `720` new boundary variables;
2. the frozen zero-`Lambda` Lorentzian Regge action with conserved homogeneous
   dust;
3. the accepted branch-B history at incoming representative `v=3/2`;
4. only the first and second accepted finite-height slabs;
5. both frozen even/odd staircase schedules;
6. logarithmic signed-squared-edge variables and action-derived pre/post
   canonical momenta;
7. the accepted numerical error envelopes and 10/100 classifier;
8. smooth irreducible local constraints selected by degeneracy of the same
   discrete Legendre relation;
9. no deletion of a nonzero direction by a numerical threshold.

It does not cover a third or later slab, an infinite history, refinement, an
enlarged dust carrier, a different/perfect action, or a separately derived
continuous momentum-map reduction.

## Route A: regularity and open images

For fixed old canonical data `(O,p_minus)`, define

```text
F(O,p_minus;X,N) = (S_X, -S_O-p_minus).
```

The derivative with respect to the `1560` variables `(X,N)` is

```text
J = [[ S_XX,  S_XN],
     [-S_OX, -S_ON]].
```

The two independent dense calculations already certify every required `J` as
regular.  The new verifier rechecked every frozen rank margin:

```text
first slab, sigma_min/gate       1119.69 ... 1119.69
second audit, first physical     4478.67 ... 4478.67
second normalized               73369.68 ... 73369.68
second physical                 86495.30 ... 86495.31
```

Thus the implicit-function theorem gives a unique local `(X,N)` for every
incoming datum in an open subset of the complete old phase space.  The exact
action-generated evolution is canonical; its accepted tangent and all four
two-step products pass the symplectic controls.  A symplectic square tangent
is invertible, so the outgoing image is also open in the complete phase space.

Consequently both local constraint codimensions are zero.

## Route B: independent constraint-covector contradiction

Assume that a smooth irreducible post-constraint `C` vanishes on the image of
the local map `T`.  Differentiation gives

```text
dC * DT = 0.
```

Because `DT` is invertible, `dC=0`, contradicting irreducibility.  Applying
the same argument to `T^-1` excludes a local pre-constraint.  This proof does
not enumerate or fit candidate covectors.

The verifier implemented this argument with exact rational linear algebra,
separately from the artifact rank bookkeeping.  Both routes returned

```text
pre-constraint codimension   0
post-constraint codimension  0
```

for each first/second parity realization and all four parity-ordered
products.

## Adversarial controls

The exact regular quadratic control has an invertible combined Legendre
matrix, a determinant-one symplectic tangent and zero constraint codimension
by both routes.

The deliberately singular control has rank one and produces one independent
pre-constraint and one independent post-constraint.  The no-go classifier
refused to fire on it.  Thus the test space could have produced a positive
constraint signal.

An exact mixed coefficient `10^-100` has rank two in rational arithmetic.  A
hostile `10^-12` cutoff misclassifies it as rank one.  The verifier explicitly
rejected the threshold quotient.  This is the finite-dimensional version of
why resolved pseudo-constraint directions cannot be deleted merely because
they are weak.

Finally, a regular first move followed by a deliberately singular future
move pulls the later constraint covector back from `(1,1)` to `(1,2)`.  This
control demonstrates that regularity of the present two moves cannot be
promoted to an infinite-history theorem.

## Why pseudo-constraints do not rescue the quotient

Exact pre/post constraints are relations on one phase-space slice and arise
from a singular Legendre image.  Curved-background pseudo-constraints depend
also on neighbouring-step data and can fix lapse-like variables.  They are
nonzero dynamics, not generators of exact gauge orbits.

Therefore a `pseudo-constraint quotient` is not a neutral operation on the
fixed finite theory.  Choosing a small-singular-value cutoff and deleting the
corresponding directions would manufacture a different theory.  The old full
anisotropic control makes this concrete: its `120` weak candidates are
resolved nonzero while the complete rank remains `1560/1560`.

## Reconciliation with the existing refined H4 result

The repository already contains a genuine singular control on a different
construction.  On the stationary barycentric product

```text
K0=P(sd K_600), tau0=0.0102,
```

the homogeneous `H4` invariant sector has one exact product-duration internal
null.  It couples to a nonzero boundary compatibility row and yields a
schedule-independent constrained response on an eleven-dimensional
hyperplane.

This does not contradict the present no-go.  It uses a different stationary
seed, a projected barycentric carrier, selected rank masses and only twelve
homogeneous boundary coordinates.  It is one refined construction, not a
nested scaling family connecting the current finite-height full boundary map
to a continuum constraint kernel.

It is nevertheless an important repository control: exact singularity is
not forbidden by the code or formalism.  It also shows why the next route must
compare matched actions and backgrounds across resolutions rather than rerun
the same coarse matrix.

## Primary-literature reconciliation after the result

The result matches the established discrete canonical distinction:

- Dittrich and Hoehn derive regular versus singular discrete Legendre maps,
  pre/post-constraint surfaces and region-dependent propagated constraints in
  [arXiv:1303.4294](https://arxiv.org/abs/1303.4294).
- Bahr and Dittrich show that curved Regge backgrounds generically replace
  exact gauge constraints with pseudo-constraints in
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670).
- Bahr and Dittrich discuss restoration through perfect actions in
  [arXiv:0909.5688](https://arxiv.org/abs/0909.5688) and construct improved or
  perfect Regge actions in special settings in
  [arXiv:0907.4323](https://arxiv.org/abs/0907.4323).
- Bahr, Dittrich and He develop coarse graining with gauge symmetry in the
  linearized setting in
  [arXiv:1011.3667](https://arxiv.org/abs/1011.3667).

The general mechanism is therefore **KNOWN**, not a discovery.  The
project-specific application to the accepted 600-cell two-slab artifacts is
reproducible, but its external novelty remains **OPEN**.

## Status ledger

- **DERIVED COMPUTATIONAL:** all accepted first/second complete
  pre-Legendre systems are regular by strict frozen margins.
- **DERIVED EXACT/STRUCTURAL:** regular local canonical evolution has no
  positive-codimension pre/post-constraint image.
- **DERIVED BOUNDED NO-GO:** the present two-slab action does not select a
  nontrivial exact local boundary constraint quotient.
- **REFUTED AS FITTING:** deleting weak nonzero pseudo-constraint directions
  by a numerical threshold.
- **OPEN:** a constraint propagated backward by a later singular slab.
- **OPEN:** exact gauge restoration under a coherent refinement family.
- **OPEN:** an enlarged matter action or a different/perfect action.
- **OPEN:** physical graviton modes, a wave operator, `c`, `G`, Planck units
  and particle physics.

## Physical meaning and next gate

This is a negative result for the coarse model, not progress toward local GR.
At this resolution the action evolves the full boundary carrier, including
directions which continuum GR would treat as gauge.  The accepted symplectic
map is mathematically real but cannot yet be interpreted as a physical
graviton phase space.

The next legitimate gate is **constraint restoration under a coherent
refinement or improved-action family**.  Before any new Hessian is evaluated,
that mission must freeze:

1. nested carriers and maps between their boundary and vertex-displacement
   spaces;
2. the same physical background, matter normalization and action prescription
   across levels;
3. the geometric candidate weak carrier before opening singular values;
4. a target-free scaling classifier for its coupling and singular values;
5. a negative control which stays regular and a positive gauge control which
   becomes exactly singular.

Only convergence toward an exact kernel, not smallness at one level, can
reopen the physical quotient route.
