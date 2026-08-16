# Prior-art gate: stationary-root handoff to the second canonical dust tick

Date: 2026-08-16

## Exact object and hypotheses

The carrier is the fixed order-24 staircase slab between two homothetic
600-cell boundaries.  The Lorentzian Regge action, five dust-pole orbits,
angle branch, boundary permutation and canonical momentum signs are those of
the already certified first tick.  The total dust mass is held fixed at

```text
M=(90/pi)*(2*pi-5*acos(1/3))*L0,
```

where `L0` is the original published boundary scale.  It is not recomputed
from either later boundary scale.

At the accepted first output `a1=log(L1/L0)` and lapse `r1=log(rho1/rho0)`,
the target-independent stationary-root enumeration has already been committed
at `caaf1f1`.  The next object is the full comparison of both committed
stationary pre-momenta with the canonically mapped post-momentum of the first
tick.  No stationary root may be changed during that comparison.

## Framing correction: what mass conservation does and does not imply

The exact identity

```text
S_grav=+720 epsilon3 L tau,
S_dust=-720 epsilon3 L tau
```

holds on the static product family with one common scale `L=L0` and with the
above normalization of `M`.  It does not assert

```text
M=(90/pi)*epsilon3*L_boundary
```

separately on every evolving boundary.  Reimposing that relation after every
scale change would change the physical mass from slab to slab and is excluded
by the present hypotheses.

With fixed `M` and unequal lower/upper scales, the lapse/pole equation is not
an identity.  This is the first branch of the dichotomy proposed in the
current audit.  It is already the mechanism used by the first canonical tick.
It does not yet derive an absolute clock: the original `tau0=0.0102` remains
externally supplied, while the pole equation and momentum seam determine
relative scale and lapse changes.

## Primary prior art

- De Felice and Fabri evolve a dust-filled regular 600-cell with the Sorkin
  algorithm and report the later causality-limiting endpoint.  This establishes
  that fixed-matter 600-cell evolution is **KNOWN**, but it does not provide
  the present target-independent stationary-root firewall:
  <https://arxiv.org/abs/gr-qc/0009093>.
- Their generalized calculation allows more variables and discusses how
  matter is introduced in the Regge model:
  <https://arxiv.org/abs/gr-qc/0106077>.
- Dittrich and Höhn derive canonical discrete evolution from Hamilton's
  principal function, with pre/post constraints and data that may be fixed by
  later consistency conditions.  This supplies the **KNOWN** canonical
  framework for the momentum seam:
  <https://arxiv.org/abs/1108.1974>.
- Their covariant-to-canonical analysis shows that curvature can break exact
  discrete gauge symmetry and turn constraints into background-dependent
  pseudo-constraints:
  <https://arxiv.org/abs/0912.1817>.

## KNOWN / CONTROL / OPEN

- **KNOWN:** action-generated pre/post momenta and consistency across a
  discrete hypersurface are the canonical Regge evolution condition.
- **KNOWN:** dust-filled 600-cell evolution and finite-step causal obstructions
  have been studied.
- **CONTROL:** the first accepted tick, the gluing permutation and its target
  momentum are committed independently of the present comparison.
- **CONTROL:** the complete multiset of two stationary roots and their
  pre-momenta was committed at `caaf1f1` before the target is parsed.
- **DERIVED:** the implementation keeps the original total mass constant; it
  does not enforce a changing `M proportional to L` at later boundaries.
- **OPEN:** whether either stationary root equals the second canonical target.
- **OPEN:** if neither equals it, whether a preregistered local two-variable
  correction on the structurally contracting branch produces a complete
  second tick.
- **OPEN:** absolute lapse/time selection, continuum convergence and external
  novelty.

A literature search cannot prove novelty.  Only the targeted verifier for the
current comparison will be run; the full repository suite will not be run.
