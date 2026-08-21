# Prior-art gate: boundary cotangent of the curvature-matched refined seed

Date: 2026-08-21

Status: written after acceptance of the local curvature-mass identity and
before evaluating any refined boundary derivative.

## Exact object and complete hypotheses

Use the now on-shell static product over

```text
K0=P(sd K_600)
```

with exact projected rank geometry, supplied `tau0=0.0102`, all 24 staircase
schedules and conserved vertex masses

```text
m_v=K_v/(8*pi).
```

The previous theorem supplies six refined old-boundary and six refined
new-boundary squared-edge orbits.  Let `G^-_i` and `G^+_i` be the total action
derivatives with respect to their logarithms, in orbit order

```text
(01,02,03,12,13,23).
```

With the repository's canonical signs,

```text
P^-_i=-G^-_i,   P^+_i=+G^+_i.
```

For common logarithmic length scale `s`, every squared edge logarithm has
derivative two.  The already proved cotangent pullback is therefore

```text
p^-_s=2 sum_i P^-_i,
p^+_s=2 sum_i P^+_i.                              (1)
```

The narrow questions are:

1. do all 24 schedules produce the same complete six-component boundary
   covectors, not merely the same scalar pullback;
2. does (1) satisfy the static product identity

   ```text
   p^-_s=-tau0*K_fine/2=-4*pi*tau0*M_fine,
   p^+_s=+tau0*K_fine/2=+4*pi*tau0*M_fine;         (2)
   ```

3. how does this compare with the accepted unrefined regular 600-cell
   momentum at the same unit volume radius?

No Hessian, spectrum, unequal-boundary solve, tick or physical constant is in
scope.

## The comparison must distinguish two meanings of "same state"

At unit volume radius the regular coarse 600-cell has

```text
zeta^3=pi^2*sqrt(2)/50,
epsilon3=2*pi-5*acos(1/3),
K_coarse=720*zeta*epsilon3,
M_coarse=K_coarse/(8*pi).
```

The projected refined carrier has its separately certified `K_fine` and
`M_fine=K_fine/(8*pi)`.  These are not exactly equal at finite resolution:
the already committed spatial actions give approximately

```text
K_fine/K_coarse = M_fine/M_coarse = 0.9841903774... .
```

Therefore fixed radius, fixed mass and exact raw momentum equality cannot all
be assumed before a perfect/improved action relates the two regulators.  The
protocol must report both:

- **fixed-radius raw comparison:** `p_fine/p_coarse`;
- **mass-normalized comparison:** `(p_fine/M_fine)/(p_coarse/M_coarse)`.

Calling the expected finite-curvature mismatch a canonical-transport failure
would be too strong.  Conversely, using mass normalization to hide the raw
mismatch would also be misleading.

## KNOWN from primary literature

- A discrete action or Hamilton principal function generates pre/post
  canonical data by its boundary derivatives: Marsden and West, *Discrete
  mechanics and variational integrators*, DOI
  `10.1017/S096249290100006X`, and Dittrich and Hoehn, *Canonical simplicial
  gravity*, [arXiv:1108.1974](https://arxiv.org/abs/1108.1974).
- Dynamically faithful relations between different Regge resolutions require
  coarse graining of the action; equality of bare actions is not automatic:
  Bahr and Dittrich, *Improved and Perfect Actions in Discrete Gravity*,
  [arXiv:0907.4323](https://arxiv.org/abs/0907.4323), and Bahr, Dittrich and
  He, [arXiv:1011.3667](https://arxiv.org/abs/1011.3667).
- Local and global Regge variations can differ at finite resolution, and
  refinement changes the accuracy of closed FLRW models: Liu and Williams,
  [arXiv:1501.07614](https://arxiv.org/abs/1501.07614).

No located source computes these six projected-600-cell boundary momentum
orbits.  Search absence is not proof; external novelty is **OPEN**.

## Repository controls

- **DERIVED EXACT:** for the regular coarse static slab,
  `p^-_s=-tau0*K_coarse/2` and `p^+_s=+tau0*K_coarse/2`.
- **DERIVED EXACT / STRUCTURAL:** the refined cotangent pullback (1) has rank
  one and a five-dimensional inverse fiber.
- **DERIVED COMPUTATIONAL:** the curvature-matched refined seed is internally
  on shell for all schedules.
- **CONTROL:** direct analytic boundary derivatives must agree with centered
  finite differences of the complete action on a schedule and its reverse.
- **CONTROL:** exchanging old and new boundaries must exchange and reverse the
  canonical covectors according to the frozen sign convention.

## Framing attack and acceptance boundary

Equation (2) alone would determine only the scalar pullback and could hold for
many points in the five-dimensional lift fiber.  The evidential content is
the complete six-component vector supplied by the action and its equality
across all 24 schedules.

The route advances only if:

1. the seed is on shell under the curvature masses;
2. every schedule selects the same six-component old/new covectors;
3. their pullbacks satisfy (2);
4. the normalized coarse/fine relation is exact within the certified
   numerical envelope.

A positive result selects one action-generated refined covector for the
curvature-matched branch.  It does not establish a perfect action because the
fixed-radius raw mismatch remains.  Schedule dependence closes the claimed
canonical lift even if the scalar pullbacks agree.

## Next admissible calculation

Preregister one all-schedule direct-action verifier.  Print the complete
six-component covectors before comparing their pullback with either coarse
target.  Include a corrupted boundary-gradient control and an artificial
uniform six-vector having the correct pullback but the wrong components, so
the verifier cannot cite a scalar match as evidence for vector selection.
