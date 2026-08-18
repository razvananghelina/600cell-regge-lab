# Prior-art gate: nonlinear canonical continuation of the 600-cell dust slab

Date: 2026-08-16

Status: **completed before any nonlinear evaluation away from the published
canonical datum**.

Upstream rank result: `715b6ad` (ledger update `961803f`).

This is a targeted primary-source map, not proof of external novelty.

## 1. Exact new object and complete hypotheses

For each already derived even/odd order-24 schedule carrier, retain the same
complete Lorentzian Regge plus De Felice--Fabri dust action, angle branch and
95 logarithmic coordinates as in the canonical-rank census.

Hold the old spatial geometry fixed at the published regular boundary.  For
the 65 unknowns

```text
y = (log x_internal[35], log q_new[30])
```

solve

```text
F(y;p_target) = (g_internal[35], p_pre(y)-p_target[30]) = 0,
p_pre = -g_old.
```

There are two distinguished targets, neither fitted:

```text
p_target = p_pre(published)    reproduction control,
p_target = P p_post(published) candidate forward junction datum.
```

The gluing verifier derived the orbit map `P` and found it to be the identity
in its independent quotient order, but the implementation must reconstruct
and apply the stored map rather than silently assume identity.

The new object is the connected nonlinear solution branch from the
reproduction datum to the forward datum.  An unrelated root reached from a
chosen seed is not canonically selected and does not count.

## 2. KNOWN structure

Dittrich and Hoehn's discrete canonical formalism uses Hamilton's principal
function to generate pre/post momenta and canonical evolution.  It permits
constraints and variable phase-space dimension and ensures agreement with
the covariant action equations:

- <https://arxiv.org/abs/1108.1974>.

Their covariant-to-canonical analysis shows that nonlinear curvature can
break exact Regge gauge symmetries and replace constraints by
background-dependent pseudo-constraints:

- <https://arxiv.org/abs/0912.1817>.

Implicit Sorkin evolution, including an early 600-cell Friedmann application,
is known:

- <https://arxiv.org/abs/gr-qc/9411008>.

Dust-filled 600-cell evolution and its causal endpoint are also prior art:

- <https://arxiv.org/abs/gr-qc/0009093>;
- <https://arxiv.org/abs/gr-qc/0106077>.

Therefore neither implicit Regge evolution, a 600-cell cosmology, nor the use
of pre/post momenta is new here.

## 3. Upstream controls and the new obstruction

**DERIVED UPSTREAM.**  The complete pre-Legendre Jacobian has calibrated rank
65/65 in both parity carriers.  Hence a unique local branch exists near the
published canonical datum in this reduced space.

**DERIVED UPSTREAM.**  The repeated-slab junction residual is nonzero and is
entirely homogeneous in its 30 pre-momentum equations.

**DERIVED DESIGN DIAGNOSTIC, before nonlinear solve.**  Applying the committed
linear inverse to the complete forward residual gives, in both parities,

```text
mean Delta log(staircase square) = +5.4091038624e-5,
mean Delta log(pole square)      = -4.0000000000,
mean Delta log(q_new)            = approximately 0,
```

with negligible within-sector spread.  This is the tangent to the already
known collective lapse family: decreasing the common pole magnitude while
adjusting staircase diagonals so `slant_square + pole_square = l0^2`.

Thus the seam is a covector in equation space, not a configuration
displacement.  The inverse canonical Hessian rotates its first linear response
almost entirely into lapse, refuting the framing that the seam direction by
itself selects a spatial-scale Newton step.  The size `-4` also lies outside a
trustworthy one-step linear neighborhood.

## 4. Symmetric reduction: use and limitation

The published datum, the two target momenta and the action are invariant
under the full regular-boundary symmetry.  Since the full pre-Legendre map is
locally unique at the base, its locally connected solution must remain in the
fixed subspace.  In that subspace the 65 unknowns reduce to

```text
one common staircase square,
one common pole square,
one common final-boundary square,
```

and the 65 equations reduce to one representative equation of each type.

This three-variable reduction is **STRUCTURAL and exact only while local
uniqueness persists**.  It is a predictor and a symmetry control, not a proof
that no nonsymmetric branches exist after a bifurcation.  Every reduced root
must be substituted into all 65 equations and, at accepted endpoints,
corrected/verified by the complete 65-variable system.

## 5. Framing attack

The affine interpolation

```text
p(lambda) = p_pre + lambda (P p_post-p_pre),  0 <= lambda <= 1,
```

is a numerical homotopy, not physical intermediate time.  Only `lambda=0`
and `lambda=1` have the canonical meanings above.

Full rank at `lambda=0` proves only local solvability.  It does not prove that
the connected branch reaches `lambda=1`.  In particular the linear direction
suggests a possible lapse-collapse obstruction near the point where the
interpolated momentum changes sign.  A failure at that point cannot be
repaired by choosing an unrelated initial seed and calling the resulting root
the evolution.

Conversely, a symmetric endpoint would be a discrete-Friedmann tick, not yet
a test of general gravitational degrees of freedom.  The latter requires
perturbed initial data after the background tick exists.

## 6. OPEN difference and next gate

No located source supplies the connected canonical branch for this exact
order-24 carrier, two parity schedules, complete Lorentzian branch and dust
action.  External novelty is **OPEN**.

The next preregistration must freeze:

1. the reproduction residual and recovery tolerance;
2. a deterministic continuation and backtracking rule with no seed search;
3. Lorentzian branch and positivity gates along the branch;
4. full 65-equation substitution gates;
5. outcomes distinguishing a reached forward target, lapse collapse,
   bifurcation/numerical openness and control failure;
6. scale/lapse/shape decomposition of any reached endpoint.

No calculation in this gate establishes a next frame, expansion, a physical
clock, inflation, `c` or Planck scales.
