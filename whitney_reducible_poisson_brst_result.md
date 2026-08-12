# The dual hierarchy gives local reducible Poisson--BRST kinematics, not dynamics

Date: 2026-08-12

Protocol commit: `6a0b208`

Inherited all-resolution locality commits: `8d0c557`, `7ba2b7b`

Registered verifier:
`reproducible/verify_whitney_reducible_poisson_brst.py`

Targeted result: **14/14 PASS**.  The full suite was not run, by explicit
user request.

## Verdict

The complete signed dual-cell resolution does canonically produce a
reducible, first-class, nilpotent and uniformly local Poisson--BRST
*kinematic complex* for the Whitney copy constraints.

> **DERIVED LOCAL REDUCIBLE POISSON--BRST KINEMATICS.**

This is a genuine advance over the earlier minimal conversion: it retains
all canonical neighbour rows and uses the dual two- and three-cell relations
instead of choosing an independent row basis or a spanning tree.

It is not yet a physical gauge theory:

- the inverse symplectic form on the auxiliary leaf is not local in the
  canonical Moore--Penrose representation;
- the obvious fully local nondegenerate realisation leaves one complete
  extra conforming spectator sector;
- the BRST kinematics permits inequivalent coefficient-free local
  Hamiltonians and therefore selects none;
- exact recovery of the already established Whitney Hamiltonian retains the
  assembled inverse metric.

The resulting label is

> **DERIVED KINEMATIC ADVANCE / PHYSICAL DYNAMICS GATE STILL CLOSED for the
> stated linear constructions.**

No claim is made against every nonlinear or differently embedded gauge
theory.

## 1. Fixed exact resolution

Degree by degree, let

\[
 0\longrightarrow Z_3\mathrel{\mathop{\longrightarrow}^{R_3}}Z_2
 \mathrel{\mathop{\longrightarrow}^{R_2}}Z_1
 \mathrel{\mathop{\longrightarrow}^{C^*}}Z_0
 \mathrel{\mathop{\longrightarrow}^{J^*}}W\longrightarrow0
\]

be the canonical occurrence/dual-cell resolution.  Here `J:W->Z0` copies a
global cochain value into every tetrahedral occurrence and `C:Z0->Z1`
records every neighbour difference.  Exactness gives

\[
 J^*C^*=0,qquad C^*R_2=0,qquad R_2R_3=0.
\]

On the exact rational boundary-of-a-four-simplex control, the dimensions and
ranks are

| degree | `(Z0,Z1,Z2,Z3,W)` | `(rank C,rank R2,rank R3)` |
|---:|---:|---:|
| 0 | `(20,30,20,5,5)` | `(15,15,5)` |
| 1 | `(30,30,10,0,10)` | `(20,10,0)` |
| 2 | `(20,10,0,0,10)` | `(10,0,0)` |
| 3 | `(5,0,0,0,5)` | `(0,0,0)` |

Every composition and rank identity is verified exactly over the rationals.

## 2. Canonical degenerate Poisson conversion

Let `M=M_loc>0` be the block-diagonal exact Whitney metric and define

\[
 G=CM^{-1}C^*.
\]

Retain every canonical constraint row and introduce `eta in Z1` with

\[
 \{\eta,\eta^*\}=+iG,
 \qquad
 \Phi=Cu+\eta.
\]

The original bracket contributes `-iG`, so

\[
 \{\Phi,\Phi^*\}=-iG+iG=0.
\]

This cancellation fixes the coefficient.  If the auxiliary bracket is
`i alpha G`, first-classness requires

\[
 (\alpha-1)G=0.
\]

Since `G` is nonzero, exactly

\[
 \alpha=1.
\]

Thus the conversion has no adjustable scale.

Because the row set is redundant, `G` is singular.  Positivity of `M` and
exactness of the dual resolution give

\[
 \ker G=\ker C^*=\operatorname{im}R_2,
\]

and, by transposition,

\[
 \operatorname{im}G=\operatorname{im}C=\ker R_2^*.
\]

The zero-Casimir symplectic leaf is therefore

\[
 L_0=\operatorname{im}G.
\]

It is not arbitrarily imposed: `Phi=0` forces
`eta=-Cu in im C=im G`.  Other affine symplectic leaves contain no point of
the constraint surface.

## 3. Reducible nilpotent differential

With ghosts `c1`, ghosts-for-ghosts `c2` and second-stage ghosts `c3`, define

\[
 \begin{aligned}
 su&=-iM^{-1}C^*c_1,\\
 s\eta&=+iGc_1,\\
 sc_1&=R_2c_2,\\
 sc_2&=R_3c_3,\\
 sc_3&=0.
 \end{aligned}
\]

Then

\[
 s\Phi=(-iG+iG)c_1=0,
\]

while

\[
 s^2u=-iM^{-1}C^*R_2c_2=0,
\]

\[
 s^2\eta=iGR_2c_2=0,
 \qquad s^2c_1=R_2R_3c_3=0.
\]

All identities pass exactly in all four cochain degrees.  This is the
minimal coordinate/ghost differential.  A full quantum BFV measure,
antighost sector and gauge-fixing fermion have not been constructed.

## 4. Correct physical quotient

The gauge transformation of the occurrence field is generated along

\[
 \operatorname{im}(M^{-1}C^*).
\]

Define the local assembled covector

\[
 y=J^*Mu.
\]

It is invariant because

\[
 \delta y=-iJ^*C^*\epsilon=0.
\]

Moreover,

\[
 \ker(J^*M)=\operatorname{im}(M^{-1}C^*),
\]

so `y` separates the gauge orbits; it is not merely one invariant among
many.  On the constraint surface `eta=-Cu`, the count is

\[
 \begin{array}{c|c}
 \text{space}&\text{complex dimension}\\ \hline
 u&n\\
 \eta\text{ leaf}&r\\
 \Phi=0&-r\\
 \text{gauge quotient}&-r
 \end{array}
 \quad\Longrightarrow\quad n-r=\dim W.
\]

For the exact control, `(n,r,dim W)=(75,45,30)` and the reduced dimension is
exactly `30`.

On a conforming representative `u=Jw`,

\[
 y=J^*MJw=M_Ww,qquad M_W=J^*MJ.
\]

The reduced Poisson bracket is correspondingly

\[
 \{y,y^*\}=-iM_W,
\]

which is precisely the covector form of the assembled Whitney bracket.

## 5. Uniform locality of the differential

The new all-resolution theorem is load bearing.  For every

\[
 K_q=\operatorname{Esd}_q(\operatorname{sd}\partial\Delta^4),
 \qquad q\geq1,
\]

the exact support bounds are independent of `q`:

| map | degree 0 | degree 1 | degree 2 | degree 3 |
|---|---:|---:|---:|---:|
| max row support of `G` | 21 | <=21 | 7 | 0 |
| max row support of `M^-1 C*` | 12 | 12 | 4 | 0 |
| max column support of `M^-1 C*` | 8 | 12 | 8 | 0 |
| max row support of `R2` | 2 | 1 | 0 | 0 |
| max column support of `R2` | 6 | 6 | 0 | 0 |
| max row support of `R3` | 1 | 0 | 0 | 0 |
| max column support of `R3` | 14 | 0 | 0 | 0 |

The `G` bounds follow because a constraint has two parent tetrahedra and the
local Whitney inverse mixes only the fixed number of forms inside either
tetrahedron.  The `R2/R3` column bounds are exactly the all-resolution edge
and vertex incidence bounds proved previously.

Thus every nonzero appearing in `s` has a uniformly bounded local stencil.
This is the precise positive result.

## 6. Where nonlocality reappears

A local Poisson tensor is not the same thing as a local symplectic form.  On
the leaf `L0`, the symplectic form is the inverse of `G` restricted to
`im G`.  Its canonical counting-metric representative is the
Moore--Penrose inverse `G+`.

On the exact control, all four Moore--Penrose identities hold, but `G+` has
300 nonzero entries joining constraints with disjoint tetrahedron endpoint
sets.  Therefore the local Poisson description has not supplied a local
canonical action in the auxiliary variables.

This support count is representation-sensitive off the leaf: adding terms
that vanish on `im G` can alter an ambient generalized-inverse matrix.  The
result therefore closes the canonical Moore--Penrose realisation, not every
possible local symplectic embedding.  A universal no-go would require
proving that no local generalized inverse or extra-field resolution induces
the same leaf form.

## 7. The obvious local symplectic realisation doubles the physics

There is a simple way to avoid the degenerate bracket: add another complete
occurrence field `v in Z0` with the opposite nondegenerate local bracket and
set

\[
 \Phi=C(u+v).
\]

The constraints commute and retain the same reducibility.  But the physical
dimension is

\[
 2n-2r=2(n-r)=2\dim W.
\]

On the control this is `60`, not `30`.  In variables `u+v` and `u-v`, one
conforming sector is the desired one and another complete conforming sector
survives as a spectator.  Removing it requires additional constraints or
gauge transformations not selected by the present geometry.

This is a **DERIVED NEGATIVE for the obvious occurrence-double symplectic
realisation**, not a theorem against every auxiliary carrier.

## 8. Hamiltonian obstruction and nonselection

The exact assembled Whitney Hamiltonian is

\[
 H_W(w)=w^*A_Ww,
 \qquad A_W=J^*A_{\rm loc}J.
\]

Since `w=M_W^-1 y`, its unique expression in the invariant coordinate is

\[
 H_W(y)=y^*M_W^{-1}A_WM_W^{-1}y.
\]

With the reduced bracket `-iM_W`, the flow is exactly

\[
 i\dot y=A_WM_W^{-1}y,
\]

equivalent to the assembled Whitney equation.  The verifier proves both the
energy transformation and the flow identity exactly.

The boundary-of-a-four-simplex is too small to expose remote entries of this
particular inverse: its exact remote count is zero.  That control therefore
cannot support a nonlocality claim.  The independent full-600-cell
certificate is the relevant evidence: the assembled inverse reaches graph
distances `(5,10,15)` in degrees zero, one and two, and the inverse depth
grows under refinement.  We keep these logically separate.

More basically, BRST invariance does not select the Whitney Hamiltonian.
Both

\[
 y^*y
 \quad\text{and}\quad
 y^*A_W^2y
\]

are coefficient-free, positive, gauge-invariant and local (the second has
bounded radius two), yet they are inequivalent: their exact ranks on the
control are `30` and `28`.  Hence the kinematic complex alone cannot decide
which dynamics is physical.

## 9. Attack on the framing

The phrase “a local BRST complex solves the nonlocality problem” would be
false.  Three distinct notions are involved:

1. **local differential:** achieved here;
2. **local nondegenerate symplectic action with the correct physical
   quotient:** still open;
3. **selected local Hamiltonian reproducing the desired physical theory:**
   still open, and not implied by BRST nilpotency.

The present result establishes only the first, plus the correct abstract
Poisson quotient.  It gives a choice-free mathematical arena in which the
remaining physical question can be asked without spanning trees.  It does
not derive a gauge force or a clock.

## Status ledger

- **DERIVED:** canonical degenerate auxiliary Poisson tensor `+iG`.
- **DERIVED:** first-class cancellation uniquely fixes its normalization.
- **DERIVED:** `ker G=im R2` and `im G=ker R2*`.
- **DERIVED:** the zero-Casimir leaf is selected by nonempty constraint
  surface.
- **DERIVED:** complete two-stage reducible nilpotent differential.
- **DERIVED:** correct physical quotient and assembled Poisson bracket.
- **DERIVED UNIFORM:** every map in the differential has bounded all-`q`
  incidence support.
- **DERIVED CONTROL:** canonical `G+` has 300 remote entries.
- **DERIVED SCOPED NEGATIVE:** the obvious local nondegenerate occurrence
  double leaves one full spectator sector.
- **DERIVED:** exact Whitney dynamics in `y` contains `M_W^-1`.
- **DERIVED NONSELECTION:** at least two inequivalent positive local
  coefficient-free Hamiltonian forms pass BRST invariance.
- **OPEN:** a canonical local symplectic realisation with no spectators.
- **OPEN:** a physical principle selecting one Hamiltonian.
- **NOT CONSTRUCTED:** full quantum BFV measure and gauge fixing.
- **NOT CLAIMED:** time, causality, inertia, mass, `c`, `hbar`, Newton's `G`
  or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_reducible_poisson_brst.py
```

Expected result: `14/14`.
