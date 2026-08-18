# Whitney copy equality is second-class rigidity, not gauge freedom

Date: 2026-08-11

Preregistration commit: `5efcfb6`

Targeted verifier:
`reproducible/verify_whitney_constraint_dirac_bergmann.py`

Targeted result: **10/10 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Result

For the exact constrained first-order action

\[
S=\int dt\left[
\frac{i}{2}(u^*M_{\rm loc}\dot u-\dot u^*M_{\rm loc}u)
-u^*A_{\rm loc}u-\lambda^*Cu-u^*C^*\lambda
\right],
\]

the complex copy-equality constraints are

\[
\phi=Cu=0.
\]

Their Dirac bracket Gram matrix is

\[
G=CM_{\rm loc}^{-1}C^*.
\]

Because (M_{\rm loc}>0), exactly

\[
y^*Gy=(C^*y)^*M_{\rm loc}^{-1}(C^*y),
\]

and hence

\[
\ker G=\ker C^*,
\qquad
\operatorname{rank}G=\operatorname{rank}C.
\]

Therefore every independent copy-equality constraint is second class.  There
are no physical first-class directions.

> **DERIVED NEGATIVE:** the cycle redundancy of the local Whitney constraints
> is not a Maxwell-like Gauss gauge symmetry.  It is freedom only in the
> redundant multiplier coordinates.

## Complete counts

| level | constraint rows | independent complex constraints | real second-class constraints | physical first-class | multiplier-only redundancy |
|---|---:|---:|---:|---:|---:|
| base 600-cell | 8,400 | 6,360 | 12,720 | 0 | 2,040 |
| first barycentric | 201,600 | 153,120 | 306,240 | 0 | 48,480 |

The real physical dimensions remaining after the second-class reduction are

\[
2\times9000-12720=5280=2\times2640

\]

and

\[
2\times216000-306240=125760=2\times62880,

\]

exactly twice the conforming complex dimensions.  Thus the constraint count
neither removes nor adds hidden physical modes.

## Constraint preservation

The equations are

\[
iM_{\rm loc}\dot u=A_{\rm loc}u+C^*\lambda,
\qquad Cu=0.
\]

Differentiating (Cu=0) once gives

\[
G\lambda=-CM_{\rm loc}^{-1}A_{\rm loc}u.
\]

On independent rows, (G) is positive definite and fixes the multiplier
uniquely.  The descriptor is therefore index two.  Redundant multiplier
shifts in ​(ker C^*) do not act on (u).

Eliminating the multiplier yields

\[
i\dot u=P_D M_{\rm loc}^{-1}A_{\rm loc}u,

\]

where

\[
P_D=I-M_{\rm loc}^{-1}C^*G^+C.

\]

The verifier checks on an independent small triangulated 3-sphere that

\[
CP_D=0,
\qquad
P_DJ=J,
\qquad
P_D^2=P_D,
\qquad
P_D^*M_{\mathrm{loc}}=M_{\mathrm{loc}}P_D
\]

as exact rational matrix identities.

It also checks exactly

\[
P_D M_{\mathrm{loc}}^{-1}A_{\mathrm{loc}}J
=J(J^TM_{\mathrm{loc}}J)^{-1}(J^TA_{\mathrm{loc}}J).
\]

Thus Dirac reduction reproduces the same assembled Whitney evolution; it does
not produce a new local generator.

## Independent control

On the boundary of a 4-simplex:

- the independent constraint matrix has shape `45 x 75` and exact rank 45;
- (G) has exact rank 45;
- its smallest numerical eigenvalue is 1.90036;
- all four projector identities pass exactly;
- the reduced vector field equals the assembled one exactly;
- the exact projector contains 1,000 nonzero entries between distinct
  tetrahedron blocks.

The cross-block count is only a small-mesh control, not a proof of the full
support.  Full-complex nonlocality is independently supported by the exact
Whitney inverse and refinement certificates.

## Causal interpretation

The local KKT pencil remains a valid and useful spectral result.  But exact
time preservation requires solving the positive constraint Gram system for
(lambda) at every instant.  After reduction, the increasingly complex
assembled inverse reappears.

Second-class status alone would not prove acausality—some local relativistic
systems have second-class constraints.  Here the additional load-bearing
facts are:

1. the multiplier is determined by the (G) solve rather than by a local
   gauge choice;
2. the exact reduced vector field equals the previously certified assembled
   Whitney generator;
3. the inverse-polynomial lower bounds grow from `(8,21,26)` to at least
   `(116,511,511)` after one refinement.

The honest verdict is therefore:

> **The constrained formulation saves uniformly local spectral geometry, but
> exact copy rigidity does not yet give a bounded-depth causal tick.**

## Remaining routes

Two nontrivial continuations remain:

1. give copy differences finite stiffness and let them propagate as genuine
   local modes; conformity then becomes approximate and introduces a scale;
2. add new variables and convert the constraints into a genuinely first-class
   system, with a derived gauge symmetry and physical embedding.

The first is outside exact Whitney assembly and was already proved incapable
of exact finite-κ invariance.  The second is mathematically possible in broad
constraint-conversion formalisms, but an arbitrary extension would be pure
fitting.  Geometry must select it.

## Status ledger

- **DERIVED:** all independent copy-equality constraints are second class.
- **DERIVED:** zero physical first-class directions at both levels.
- **DERIVED:** multiplier redundancies `(2040,48480)` act only on ​(lambda).
- **DERIVED:** descriptor differentiation index two.
- **DERIVED:** exact Dirac projector and reduced Whitney vector field.
- **DERIVED NEGATIVE:** no Gauss-law gauge rescue in the canonical system.
- **STRUCTURAL:** exact conformity as a rigid/instantaneous limit.
- **OPEN:** a geometry-selected first-class extension.
- **OPEN:** causal finite-stiffness approximate conformity and its continuum
  scaling.
- **NOT CLAIMED:** Lorentzian time, inertia, mass or (c).

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_constraint_dirac_bergmann.py
```

Expected result: `10/10`.
