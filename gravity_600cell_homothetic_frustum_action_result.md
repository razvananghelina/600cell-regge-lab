# Homogeneous 600-cell Regge action is exactly subdivision independent

Date: 2026-08-17

Prior-art gate: `f133de0`

Frozen protocol: `cd7a2ea`

Registered implementation: `b6ae341`

Lorentzian branch correction after the preserved first failure: `16389c3`

Independent closed-form control: `e55ba6e`

Targeted verifier:
`reproducible/verify_gravity_600cell_homothetic_frustum_action.py`

Artifact:
`reproducible/gravity_600cell_homothetic_frustum_action.json`

Artifact SHA-256:
`c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d`

## Verdict

**DERIVED:** the targeted verifier passes `16/16` and returns

```text
HOMOTHETIC_FRUSTUM_ACTION_INVARIANT
```

On the complete regular 600-cell slab, the direct cellular action of 600
flat Lorentzian tetrahedral frusta is exactly the restriction of both
committed global staircase actions to

```text
q_old      = L_minus^2,
q_new      = L_plus^2,
q_diagonal = L_minus L_plus-rho,
q_strut    = -rho,
L_minus,L_plus,rho>0.
```

The equality includes the Regge boundary term and the same conserved-dust
term.  Every artificial-diagonal equation vanishes on this complete
homothetic family.  Therefore the reduced old/new momenta and lapse
derivative are also independent of the even/odd staircase choice.

This removes a subdivision objection to the four already accepted
homogeneous ticks.  It does not derive new dynamics by itself and does not
settle the anisotropic theory.

## 1. Direct cellular action

Write

```text
Delta = L_plus-L_minus,
h     = sqrt(rho+Delta^2/4),
c_t   = (Delta^2+2 rho)/(2(Delta^2+3 rho)),
b     = Delta/sqrt(8(Delta^2+3 rho)),
epsilon_t = 2 pi-5 acos(c_t).
```

The exact gravitational action on the selected causally regular branch is

```text
S_grav = 360 (L_minus+L_plus) h epsilon_t
       + 600 sqrt(3) (L_minus^2-L_plus^2) asinh(b).
```

The first term is the contribution of the 720 timelike lateral trapezoids.
The second is the oriented lower/upper spatial-boundary contribution from
the 1,200 triangular hinges on each boundary.  With conserved total dust
mass `M`,

```text
S_total = S_grav-8 pi M sqrt(rho).
```

At `L_minus=L_plus=L`, the boost term vanishes and this reduces exactly to

```text
S_grav = 720 L sqrt(rho) [2 pi-5 acos(1/3)],
```

recovering the already proved regular-lapse identity.

## 2. Independent angle reconstruction

No cellular angle was copied from the staircase evaluator.  The verifier
reconstructed the intrinsic Lorentzian metric of the common frustum and used
the inverse metric on its six facet conormals.

The exact squared conormal relations reduce to

```text
cos(theta_lateral)^2
  = [(Delta^2+2 rho)/(2(Delta^2+3 rho))]^2,

cos(theta_base)^2
  = -Delta^2/[8(Delta^2+3 rho)].
```

The branch is fixed from the static anchors and the repository's lower-side
Lorentzian logarithm convention.  Thus the lower and upper base cosines are
opposite imaginary values, while the lateral angle is real.  This gives the
real `asinh` boundary term above.

### Preserved first implementation failure

The first run used the principal square root on the negative real axis.  It
passed all carrier and area gates and matched the static action to about
`1e-99`, but reversed the nonstatic boundary boost.  It consequently failed
four of fifteen checks.

That was not evidence for subdivision dependence.  The negative outcome was
not admissible under the protocol because the independent boundary-angle
convention gate itself had failed.  The correction in `16389c3` uses

```text
sqrt(-x-i0) = -i sqrt(x),
```

consistent with the already fixed
`log(-x-i0)=log(x)-i pi` branch.  No test point, coefficient or expected
action value was altered.  The corrected direct action then agreed at every
pre-registered point.

## 3. Complete carrier census

Both global parities independently give

```text
spatial boundary triangles       2,400
cell-facet subdivision triangles 2,400
trapezoid subdivision triangles  1,440
total                            6,240.
```

The cellular supports are sharp:

- every spatial triangle belongs to the two tetrahedral frusta adjacent to
  its 600-cell face;
- every facet-subdivision triangle belongs to the two frusta sharing the
  corresponding triangular-prism facet;
- every vertical triangle belongs to the five frusta meeting at the
  corresponding 600-cell edge.

The 1,440 vertical triangles pair exactly into the 720 lateral trapezoids,
and every one of the 720 staircase diagonals occurs in precisely the two
triangles of one such trapezoid.  There are no unclassified or mixed
supports.

Since the earlier geometric theorem proves that each local staircase is a
nonfolded triangulation of the same flat frustum, these support facts meet the
hypotheses of the flat-cell Regge subdivision theorem:

- subdivision-only hinges have zero deficit;
- cellular wedge angles add across their simplex pieces;
- signed hinge areas add.

This is the exact all-domain proof.  The finite evaluations below only
control branch and implementation conventions.

## 4. Why the artificial diagonals carry no equation

Let `d^2` be the squared diagonal of one lateral trapezoid.  Holding its two
spatial edges and two struts fixed, the two Lorentzian triangle-area
derivatives at the geometric value

```text
d^2=L_minus L_plus-rho
```

are exactly

```text
partial A_lower/partial d^2
 = +i (L_plus-L_minus)/(4 sqrt(4 rho+Delta^2)),

partial A_upper/partial d^2
 = -i (L_plus-L_minus)/(4 sqrt(4 rho+Delta^2)).
```

They cancel.  Every other triangle touching such a diagonal is a
subdivision-only hinge with zero curvature.  Hence all 720 physical diagonal
derivatives, or all 30 stored diagonal-orbit gradients, vanish identically on
the homothetic family.

This is the missing condition that makes the canonical statement valid:
the chain rule from the cellular action to a staircase action contains no
hidden diagonal term.  Therefore the cellular derivatives really equal the
sums of the 720 old-edge, 720 new-edge and 120 strut derivatives.

## 5. Arbitrary-precision controls

The independently constructed cellular-normal action, the closed real
formula, and both staircase actions were evaluated at the six frozen points

```text
(1,1,1/16),       (1,3/4,1/16),
(1,5/4,1/16),     (1,1/2,1/8),
(1,3/2,1/8),      (2,3,1/2).
```

All 2,400 simplices per parity remained Lorentzian with inertia `(3,1)`.
The maximum errors were

```text
cellular versus staircase action       5.0783e-99 relative
closed form versus cellular normals    1.2559e-99 relative
collective derivative                  2.0690e-70 relative
individual artificial-diagonal gradient 4.2878e-100 absolute
imaginary contamination                2.0421e-97 absolute.
```

The derivatives used both preregistered centered steps and Richardson
extrapolation; they are controls, not a fitted reconstruction of the action.

## 6. Physical status

### What advances

**DERIVED:** schedule ambiguity is gauge/subdivision data throughout the
regular homothetic sector, at the level of metric, complete action and
collective canonical data.  The old even/odd agreement is now explained by
one cellular Regge action rather than by a numerical coincidence.

**DERIVED:** the relevant homogeneous dynamics can now be studied with the
three physical variables `(L_minus,L_plus,rho)` and the exact action above;
the 30 artificial diagonal orbits need not be treated as physical degrees of
freedom on this sector.

### What does not advance

The equality is largely an instance of a known flat-block subdivision
principle.  It is a consistency theorem, not new gravitational dynamics.

It does not prove:

- equality of the unrestricted `30+35+30` variable actions;
- flat-frustum embeddability after anisotropic edge perturbations;
- stability or a lattice-graviton spectrum;
- a selected nonzero lapse or emergent time;
- a causal limiting speed;
- Planck time, Planck mass or particle masses.

Those remain **OPEN**.

## 7. Post-result prior-art audit

The technical terms learned in the calculation strengthen rather than weaken
the prior-art classification.

- Tsuda and Fujiwara, [*Oscillating 4-Polytopal Universe in Regge
  Calculus*](https://doi.org/10.1093/ptep/ptab074), explicitly use the
  tetrahedral 4-frustum, its polygon/trapezoid hinges and the 600-cell counts;
  they also state flat-block triangulation invariance of the Regge action.
- Dittrich and Steinhaus, [*Path integral measure and triangulation
  independence in discrete gravity*](https://arxiv.org/abs/1110.6866), and
  Bahr, Dittrich and Steinhaus, [*Discretization independence implies
  non-locality in 4D discrete quantum
  gravity*](https://arxiv.org/abs/1404.5288), delimit the result: generic 4D
  triangulation independence is not supplied by this flat-sector identity.
- Jercher and Steinhaus, [*Cosmology in Lorentzian Regge
  calculus*](https://arxiv.org/abs/2312.11639), make the Lorentzian frustum
  branch and oriented boundary angles explicit in a modern cosmological
  setting.
- Dittrich and Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974), establish why equality of the
  complete one-step action, including its boundary term, controls the
  canonical momenta.

Thus the cellular action and its subdivision invariance are **KNOWN**.  The
exact reconciliation of this repository's two chromatic 600-cell schedules,
including the artificial-diagonal cancellation and its inherited complex
branch, was not located explicitly.  External novelty of that narrow audit
remains **OPEN**.

## 8. Next falsification test

Use the closed cellular action, not either 2,400-simplex representation, to
derive the exact homogeneous equations

```text
partial S_total/partial rho     = 0,
p_pre  = -partial S_total/partial L_minus^2,
p_post = +partial S_total/partial L_plus^2.
```

Then compare their solution branches with the four committed ticks and with
the closed Friedmann dust equation.  The decisive questions are whether the
lapse remains gauge, whether conserved dust selects it as a
pseudo-constraint, and whether the accepted numerical branch is unique.

Only after that should the programme perturb individual boundary edges.  In
that anisotropic test, loss of flat-frustum embeddability or schedule
dependence would be physical rather than a coordinate artifact.
