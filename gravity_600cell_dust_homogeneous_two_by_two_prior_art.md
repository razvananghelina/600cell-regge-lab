# Prior-art gate: homogeneous two-dimensional phase reduction

Date: 2026-08-17

## Exact object, operator, carrier and hypotheses

Use the already accepted non-static fixed-mass 600-cell Regge--dust slab and
its two independently constructed order-24 staircase schedules (`even` and
`odd`).  On the trivial binary-tetrahedral sector, the action-generated
canonical tangent is

```text
T : (delta q_old[30], delta p_pre[30])
      -> (delta q_new[30], delta p_post[30]).
```

The variables are logarithmic squared boundary lengths and their conjugate
covectors.  The proposed homogeneous carrier is the *fixed*, geometry-defined
real plane

```text
U = span{u_q, u_p},
u_q = (1_30, 0_30)/sqrt(30),
u_p = (0_30, 1_30)/sqrt(30).
```

No basis direction may be learned from the tangent or from its spectrum.  The
internal-curvature response is the already defined map

```text
F : C^60 -> C^160,
```

obtained after solving the linearized pre-Legendre equations for the 35
internal and 30 final quotient variables and then differentiating all 160
internal triangle-deficit orbit values.  Its calibrated kernel has already
been shown to be one-dimensional and contained in `U`.

The new questions are exactly:

1. does `T U` lie in `U`, rather than merely having one nearly invariant
   vector there;
2. what is the calibrated `2 x 2` compression `U* T U`;
3. what is the one-dimensional kernel of the calibrated `160 x 2` response
   `F U`;
4. is that line invariant under the full `T`, and if so which eigenvalue of
   the `2 x 2` compression it carries;
5. are the answers identical under the literal, target-independent
   even-to-odd orbit identification already derived?

The calculation is local and linearized at one curved homogeneous slab.  It
does not assume a continuum limit, a physical proper-time normalization, a
gauge interpretation or nonlinear integrability.

## Primary prior art

- Dittrich and Hoehn derive action-generated pre/post canonical maps and the
  linearized tent-move dynamics of Regge data.  On flat backgrounds their
  Hessian null vectors generate lapse/shift; beyond that regime the
  constraints become background-dependent pseudo-constraints:
  <https://arxiv.org/abs/0912.1817>.
- De Felice and Fabri evolve a homogeneous dust-filled 600-cell with the
  Sorkin algorithm and compare the evolution with closed Robertson--Walker
  dynamics: <https://arxiv.org/abs/gr-qc/0009093>.
- Their generalized 600-cell evolution allows more free variables and studies
  the resulting causality-breaking endpoint:
  <https://arxiv.org/abs/gr-qc/0106077>.
- Barrett et al. give the earlier implicit Sorkin evolution of a homogeneous
  dust-filled 600-cell: <https://arxiv.org/abs/gr-qc/9411008>.

These establish homogeneous 600-cell Regge evolution and canonical
linearized Regge dynamics as prior art.  No located primary source computes
the present fixed-mass, 30-orbit boundary tangent, its exact constant
position/momentum plane, or its intersection with the 160-orbit internal
curvature response.  A search is not a novelty proof; external novelty is
**OPEN**.

## KNOWN / CONTROL / OPEN

- **KNOWN:** a regular discrete action is a generating function for a
  canonical pre/post evolution map.
- **KNOWN:** on flat Regge backgrounds, Hessian null directions can represent
  vertex-displacement gauge; curvature generically replaces exact constraints
  by pseudo-constraints.
- **CONTROL:** both schedule tangents are symplectic inside their calibrated
  derivative/ball errors.
- **CONTROL:** the full curvature response has rank `1439/1440`, and its
  unique null line lies in the trivial sector.
- **CONTROL:** in each schedule the trivial-sector null line is contained in
  `U` but is separated from both coordinate axes.
- **CONTROL:** the two schedule kernel lines agree after the literal
  orbit-set permutation.
- **OPEN:** whether the whole plane `U` is invariant under `T`.
- **OPEN:** whether the curvature-kernel line alone is invariant under `T`.
- **OPEN:** whether a resulting multiplier is exactly `-1`, merely close to
  it, or separated from it.
- **OPEN:** whether the line is a gauge/lapse direction, a linearized
  constraint direction, or a physical homogeneous mode.
- **OPEN:** nonlinear continuation, proper time and refinement stability.

## Framing attack

Homogeneity of the boundary state does **not** prove invariance of `U`.  The
staircase triangulation is certified only under the order-24 action used for
the quotient, not under the full 600-cell symmetry acting transitively on all
boundary edges.  Consequently, `T` may send a uniform perturbation into
zero-sum shape components.  A `2 x 2` compression is an actual restriction
only after this leakage is independently resolved as zero.

Even if the kernel line is an eigenline with multiplier `-1`, that fact alone
does not make it time, lapse or gauge.  Gauge requires the relevant canonical
constraints; time requires a clock/proper-time interpretation; physical
evolution requires nonlinear integrability or a multi-slab propagation test.

The new computation therefore tests a finite linear-algebra statement and
nothing broader.  Only its targeted verifier will be run; the full suite is
excluded by instruction.
