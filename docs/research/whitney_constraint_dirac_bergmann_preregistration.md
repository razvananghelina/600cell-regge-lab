# Preregistration: Dirac--Bergmann class of the local Whitney constraints

Date: 2026-08-11

## Exact system and complete hypotheses

Let (M_{\rm loc}>0) be the block-local Whitney metric,
(A_{\rm loc}=A_{\rm loc}^*) the local weak Kähler--Dirac matrix, and (C)
the canonical degree-at-most-three neighbour constraint from protocol
`a819a52`.  Consider the complex first-order constrained action

\[
S=\int dt\left[
\frac{i}{2}(u^*M_{\rm loc}\dot u-\dot u^*M_{\rm loc}u)
-u^*A_{\rm loc}u
-\lambda^*Cu-u^*C^*\lambda
\right].
\]

Its equations are

\[
iM_{\rm loc}\dot u=A_{\rm loc}u+C^*\lambda,
\qquad
Cu=0.
\]

This action is selected by the already-derived local pencil.  No multiplier
kinetic term, stiffness, gauge coupling or new scale is added.

## Frozen Dirac--Bergmann test

The unconstrained Poisson bracket is

\[
\{u_a,u_b^*\}=-i(M_{\rm loc}^{-1})_{ab}.
\]

For the complex equality constraints ​(phi=Cu), define

\[
G=CM_{\rm loc}^{-1}C^*.
\]

The real constraint-bracket matrix has off-diagonal blocks ​(pm iG).
Classify only independent constraint rows; redundant rows are retained in the
canonical local presentation but quotiented for the Dirac class.

Because (M_{\rm loc}) is positive definite,

\[
y^*Gy=(C^*y)^*M_{\rm loc}^{-1}(C^*y),

\]

so

\[
\ker G=\ker C^*,
\qquad
\operatorname{rank}G=\operatorname{rank}C.
\]

The equality is a theorem under the stated hypotheses, not a numerical rank
guess.  Combine it with the exact graph-incidence ranks already recomputed at
both levels:

- base: `rows(C)=8400`, `rank(C)=6360`;
- first barycentric: `rows(C)=201600`, `rank(C)=153120`.

Record:

- complex independent constraint rank (r);
- real second-class count (2r);
- multiplier-only redundancy ​`rows(C)-r`;
- physical first-class count, which is zero iff the quotient bracket is
  nonsingular.

## Preservation and reduced vector field

Constraint preservation gives

\[
G\lambda=-CM_{\rm loc}^{-1}A_{\rm loc}u.

\]

On the independent row quotient, (G) is positive definite and fixes
(lambda) uniquely.  With (G^+) denoting the inverse on
(operatorname{ran}C), the reduced evolution is

\[
i\dot u=
\left[I-M_{\rm loc}^{-1}C^*G^+C\right]
M_{\rm loc}^{-1}A_{\rm loc}u.

\]

Freeze the following identities:

1. the bracket projector
   ​(P=I-M^{-1}C^*G^+C) satisfies (CP=0), (PJ=J), (P^2=P), and is
   (M)-self-adjoint;
2. on (u=Jx), the reduced vector field equals

   \[
   J(J^TM_{\rm loc}J)^{-1}(J^TA_{\rm loc}J)x;
   \]

3. the semi-explicit descriptor has differentiation index two: the algebraic
   constraint must be differentiated once to determine the multiplier.

## Independent small calibration

On the boundary of a 4-simplex, use a full-row-rank constraint basis so no
pseudoinverse or multiplier gauge remains.  Compute (G), its eigenvalues and
the projector directly.  Verify all four projector identities and equality
of the reduced and assembled Whitney vector fields.  Record exact/nonzero
support across distinct tetrahedron blocks; do not infer the full-complex
support solely from one small mesh.

## Decision boundary

- **DERIVED SECOND-CLASS CONSTRAINTS:**
  `rank(G)=rank(C)` at both levels and the complete small-control identities
  pass.  Then there is no physical first-class gauge freedom in copy equality;
  row-cycle redundancy acts only on multipliers.
- **FIRST-CLASS REFUTATION OF THE FRAMING:** (G) has a kernel beyond
  ​(ker C^*).  That would provide genuine gauge directions and invalidate
  the second-class claim.

## Mandatory scope

A second-class result does not by itself prove acausality; local relativistic
theories can contain second-class constraints.  The relevant additional fact
is that multiplier preservation requires solving (G), and the already
certified reduced Whitney projector contains the assembled inverse whose
depth grows under refinement.

Therefore a pass licenses only:

- **DERIVED NEGATIVE:** the current constraint redundancy is not a Maxwell-like
  first-class gauge rescue;
- **STRUCTURAL:** exact copy equality behaves as a rigid/instantaneous limit;
- **OPEN:** a different first-class extension with new auxiliary variables;
- **OPEN:** finite-stiffness local waves that approximate, rather than exactly
  impose, Whitney conformity.

Only the new targeted verifier will be run.
