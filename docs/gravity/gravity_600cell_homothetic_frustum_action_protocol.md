# Protocol: cellular versus staircase action on the homothetic 600-cell slab

Date: 2026-08-17

Prior-art commit: `f133de0`

Status: frozen before implementation and before evaluating any new cellular
frustum action.

## 1. Registered question

On the complete regular 600-cell slab and under all hypotheses stated in
`gravity_600cell_homothetic_frustum_action_prior_art.md`, compare:

1. the direct cellular Regge action of 600 flat tetrahedral frusta;
2. the existing even global staircase action;
3. the existing odd global staircase action.

The comparison includes the additive spatial-boundary term and the same 120
conserved-dust worldlines.  It is restricted to the regular homothetic family

```text
q_old       = L_minus^2,
q_new       = L_plus^2,
q_diagonal  = L_minus L_plus-rho,
q_strut     = -rho,
L_minus,L_plus,rho>0.
```

No accepted tick, root or fitted parameter is used to choose the control
points or the outcome.

## 2. Independent cellular construction

Use the already certified intrinsic coordinates

```text
0 <= z <= 1,
y_i >= 0,
sum_i y_i = R_minus+(R_plus-R_minus)z,
R_+-R_- = phi(L_plus-L_minus),
```

and independently reconstruct their exact Lorentzian metric.  The six
facets are the lower tetrahedron, upper tetrahedron and four lateral
three-frusta.

For each codimension-two cellular hinge, compute the cell dihedral angle from
the inverse intrinsic metric and the two exact facet covectors.  Fix the
branch by continuity from `L_plus=L_minus`, where:

- the lateral-lateral internal angle is `acos(1/3)`;
- each base-lateral internal angle is `pi/2`;
- the one-slab spatial boundary exterior curvature is zero.

The direct action uses:

```text
720  lateral trapezoidal hinges,
1200 lower spatial triangular hinges,
1200 upper spatial triangular hinges,
600  flat frusta,
120  dust struts.
```

No angle or curvature may be copied from the staircase evaluator.  That
evaluator may be called only after the cellular quantities have been built.

## 3. Exact structural gates

The verifier must decide each of the following without a floating tolerance.

### 3.1 Carrier and hinge classification

1. Reconstruct both committed global schedule carriers and retain their
   existing certificate counts.
2. Recover exactly 600 cellular frusta and the 600-cell f-vector
   `(120,720,1200,600)`.
3. Classify every staircase triangle by its minimal cellular support:
   spatial hinge, lateral-trapezoid subdivision, or subdivision-only hinge.
4. Require no unclassified or multiply classified triangle.
5. Pair the 1,440 vertical triangles into 720 lateral trapezoids in each
   parity, two triangles per trapezoid with the same geometric diagonal.

### 3.2 Flat subdivision terms

1. Every subdivision-only hinge must have zero total Regge curvature as a
   geometric consequence of the common flat-frustum realization.
2. The simplex angles incident on a spatial boundary triangle must sum to
   the independently computed cellular base-lateral angle contribution.
3. The simplex angles incident on either triangle of a lateral trapezoid must
   give the independently computed cellular lateral-lateral angle.
4. The two complex Lorentzian triangle areas must add exactly to the signed
   trapezoid area.
5. Spatial triangle areas must equal `sqrt(3)L^2/4` on the corresponding
   boundary.

The geometric subdivision theorem may be used only after these hypotheses
have been checked for the actual carriers; it may not replace the carrier
audit.

### 3.3 Artificial-diagonal equation

For a trapezoid with fixed lower edge, upper edge and two equal struts,
differentiate the sum of its two triangle areas with respect to the squared
diagonal while holding those four sides fixed.  At

```text
d^2=L_minus L_plus-rho
```

the derivative must vanish exactly.  Together with zero curvature on every
subdivision-only hinge, this must imply that every one of the 30 stored
diagonal-orbit gradients in both parities vanishes on the complete
homothetic family.

This gate is essential.  Without it, equality of restricted actions cannot
be promoted to equality of canonical boundary momenta.

## 4. Functional action and derivative comparison

After the exact gates pass, form the cellular action

```text
S_cell(L_minus,L_plus,rho,M)
```

from the three cellular hinge contributions plus the dust term.  Prove as a
symbolic geometric identity, rather than by interpolation, that each
staircase sum reduces to `S_cell` on the stated connected branch.

Then require the chain-rule identities:

1. the cellular derivative with respect to common old boundary squared
   length equals the sum of the 720 per-edge old staircase derivatives;
2. the analogous identity holds for the 720 new edges;
3. the cellular `rho` derivative equals the sum over the 120 strut
   derivatives after including the zero diagonal contribution;
4. even and odd parities agree in all three derivatives.

The pre/post sign convention remains

```text
p_pre  = -partial S/partial q_old,
p_post = +partial S/partial q_new.
```

Claims are only for the homogeneous collective momenta.  No unrestricted
anisotropic Hessian identity is tested.

## 5. Frozen arbitrary-precision controls

The symbolic proof is independently controlled at 100 decimal digits at the
following six rational points, listed as
`(L_minus,L_plus,rho)`:

```text
(1,   1,   1/16),
(1,   3/4, 1/16),
(1,   5/4, 1/16),
(1,   1/2, 1/8),
(1,   3/2, 1/8),
(2,   3,   1/2).
```

At every point:

- all 2,400 four-simplices in each parity must have inertia `(3,1)`;
- direct, even and odd gravitational actions must agree relatively below
  `1e-70`;
- the complete actions including the already fixed dust mass must agree
  relatively below `1e-70`;
- centered log-coordinate derivatives, evaluated independently at steps
  `1e-18` and `3e-18` with Richardson agreement, must match below `1e-50`;
- imaginary contamination of the final real action and momenta must be below
  `1e-70`.

These finite evaluations are convention controls only.  They cannot prove
the functional identity.

## 6. Registered outcomes

### `HOMOTHETIC_FRUSTUM_ACTION_INVARIANT`

Report **DERIVED** only if every exact carrier, geometry, area, curvature,
diagonal-gradient and symbolic action gate passes, and all arbitrary-
precision controls agree.  Then the four previously accepted homogeneous
ticks are protected against the even/odd staircase-subdivision objection.

### `HOMOTHETIC_FRUSTUM_ACTION_SUBDIVISION_DEPENDENT`

Report **DERIVED NEGATIVE** only if the independently built cellular action
or derivative disagrees on the same certified Lorentzian branch, after the
boundary convention and angle branch controls pass.  Identify the first
nonzero term and whether it is bulk, boundary or dust.

### `HOMOTHETIC_FRUSTUM_ACTION_OPEN`

Use **OPEN** if a branch, boundary convention, carrier classification or
exact structural implication cannot be certified.  Numerical agreement
alone is `PATTERN`, never `DERIVED`.

## 7. Scope ledger fixed before evaluation

- Equality of homogeneous actions and collective momenta: tested.
- Equality of all local schedules that form a conforming global
  triangulation: follows only if the checked flat-subdivision hypotheses are
  schedule-independent.
- Existence of arbitrary independent global schedule choices: not tested.
- Unrestricted 30+35+30-variable action, anisotropic perturbations and
  graviton modes: **OPEN**.
- Selection of a nonzero lapse or physical clock: **OPEN**.
- Emergent `c`, Planck time, Planck mass and particle masses: **OPEN**.
- External novelty: **OPEN**.
