# Protocol: reducible Poisson--BRST conversion of Whitney copy constraints

Date: 2026-08-12

This protocol is committed before constructing or testing the candidate
conversion below.  No spectrum, particle target, coupling, mass, `a1=5`,
clock rate or Planck quantity may be consulted.

The purpose is to separate two questions that must not be conflated:

1. does the canonical dual resolution produce a local nilpotent reducible
   first-class/BRST *kinematic complex*?
2. does that complex also select a local symplectic physical theory and a
   local Hamiltonian with the original Whitney dynamics?

A positive answer to the first question is not evidence for the second.

## Fixed input and complete hypotheses

Work degree by degree and then take the direct sum over cochain degrees.
Let

\[
 0\longrightarrow Z_3\mathrel{\mathop{\longrightarrow}^{R_3}}Z_2
 \mathrel{\mathop{\longrightarrow}^{R_2}}Z_1
 \mathrel{\mathop{\longrightarrow}^{C^*}}Z_0
 \mathrel{\mathop{\longrightarrow}^{A}}W\longrightarrow0
\]

be the exact signed dual-cell resolution already derived for the canonical
neighbour copy constraints.  Thus

\[
 AC^*=0,\qquad C^*R_2=0,\qquad R_2R_3=0,
\]

and the sequence is exact at every term.  Set `J=A^*`, so `CJ=0` and
`im J=ker C`.

Let `M=M_loc>0` be the block-diagonal exact tetrahedral Whitney metric and
let `A_loc=A_loc^*` be the fixed local weak Kähler--Dirac form.  The original
complex field has bracket

\[
 \{u,u^*\}=-iM^{-1}.
\]

The all-resolution theorem in commits `8d0c557` and `7ba2b7b` is inherited:
the nonzero degrees of `C^*`, `R_2` and `R_3` are uniformly bounded on

\[
 K_q=\operatorname{Esd}_q(\operatorname{sd}\partial\Delta^4),
 \qquad q\geq1.
\]

No independent constraint basis, spanning tree or fitted Schur coefficient
is allowed.

## Frozen Poisson conversion

Retain every canonical constraint row and define

\[
 G=CM^{-1}C^*\quad\text{on }Z_1.
\]

Introduce a row-space auxiliary coordinate `eta in Z_1` with the constant
Poisson tensor

\[
 \{\eta,\eta^*\}=+iG,
 \qquad \{u,\eta^*\}=0,
\]

and define

\[
 \Phi=Cu+\eta.
\]

Because `C` is redundant, `G` is singular.  The candidate physical
auxiliary arena is the zero-Casimir symplectic leaf

\[
 L_0=\operatorname{im}G.
\]

The audit must prove, not assume,

\[
 \ker G=\ker C^*=\operatorname{im}R_2,
 \qquad
 \operatorname{im}G=\operatorname{im}C=\ker R_2^*.
\]

On this leaf the constraints must be reducible identities:

\[
 R_2^*\Phi=0.
\]

## Frozen minimal BRST differential

Introduce ghosts `c_1 in Z_1`, ghosts-for-ghosts `c_2 in Z_2`, and, when
present, `c_3 in Z_3`, with successive ghost numbers one, two and three.
Freeze the linear differential

\[
 \begin{aligned}
 su&=-iM^{-1}C^*c_1,\\
 s\eta&=+iGc_1,\\
 sc_1&=R_2c_2,\\
 sc_2&=R_3c_3,\\
 sc_3&=0.
 \end{aligned}
\]

Signs may be corrected only if the simultaneous bracket convention requires
it; a correction must precede numerical evaluation.  The required exact
identities are

\[
 s\Phi=0,qquad s^2u=s^2\eta=s^2c_1=s^2c_2=0.
\]

This is a certificate for the minimal coordinate/ghost differential, not by
itself a full quantum BFV measure or gauge-fixed path integral.

## Physical quotient gate

On `Phi=0` and `eta in L_0`, gauge transformations are

\[
 \delta u=-iM^{-1}C^*\epsilon,
 \qquad
 \delta\eta=+iG\epsilon.
\]

The proposed complete invariant coordinate is the assembled covector

\[
 y=J^*Mu.
\]

The audit must prove:

1. `y` is gauge invariant;
2. `ker(J^*M)=im(M^{-1}C^*)`, so `y` separates gauge orbits;
3. the quotient has complex dimension `dim W=n-rank(C)`;
4. on the conforming representative `u=Jw`,
   `y=M_W w` with `M_W=J^*MJ`.

Recovering `w` from `y` requires `M_W^{-1}`.  This fact must be reported and
may not be hidden by calling `y` a local copy of `w`.

## Locality and scale gates

1. Every nonzero of `G` must join constraints incident to a common
   tetrahedron.  Its row support must therefore have an all-`q` bound derived
   from the local 15-dimensional tetrahedral block and the already proved
   incidence bounds.
2. `s` must use only `M^{-1}` inside one tetrahedron and the uniformly local
   maps `C^*`, `R_2`, `R_3`.
3. Replacing the auxiliary bracket by `i alpha G` must force `alpha=1` from
   first-class cancellation.  No free conversion scale is allowed.
4. The inverse of `G` on `L_0` is the symplectic form of that leaf.  Its
   support must be audited separately from the local Poisson tensor.

## Symplectic-realisation negative control

Test the obvious fully local nondegenerate realisation: introduce a second
occurrence field `v in Z_0` with opposite bracket

\[
 \{v,v^*\}=+iM^{-1},
 \qquad \Phi=C(u+v).
\]

The constraints commute and are reducible, but the reduced complex dimension
is predicted to be

\[
 2n-2r=2(n-r)=2\dim W,
\]

so a complete extra conforming spectator sector remains.  Confirm or refute
this.  If confirmed, the occurrence-field realisation does not provide the
desired physical theory unless further independently selected constraints or
gauge structure are added.

## Hamiltonian gate

Let

\[
 M_W=J^*MJ,qquad A_W=J^*A_{\rm loc}J.
\]

The original Whitney physical Hamiltonian in the conforming coordinate is

\[
 H_W(w)=w^*A_Ww.
\]

In the local invariant covector `y=M_Ww`, exact equivalence requires

\[
 H_W(y)=y^*M_W^{-1}A_WM_W^{-1}y,
\]

and the strong flow contains `M_W^{-1}`.  Verify the transformation exactly
on an independent rational control.

Also exhibit at least two inequivalent coefficient-free local Hermitian
quadratic forms in `y` allowed by the BRST kinematics, if they exist.  Their
existence proves that the incidence/BRST complex alone does not select a
physical Hamiltonian.  Do not compare their spectra with a desired target.

## Frozen controls

1. Use the boundary of a 4-simplex for exact rational Poisson, quotient,
   nilpotency and Hamiltonian identities.
2. Use every cochain degree, including the nontrivial two-stage reducibility
   in degree zero.
3. Use the all-`q` incidence theorem for the infinite locality statement; do
   not infer it from a finite support fit.

## Decision boundaries

**Kinematic acceptance:** the Poisson leaf, first-class cancellation,
complete reducibility, nilpotent differential, correct quotient and uniform
locality all pass without a free coefficient.  Label this **DERIVED LOCAL
REDUCIBLE POISSON--BRST KINEMATICS**.

**Physical acceptance:** in addition, a canonical nondegenerate symplectic
realisation with no extra physical modes and a selected local Hamiltonian
reproducing the Whitney physical dynamics must be obtained without `G^+`,
`M_W^{-1}`, a basis choice or a free scale.

**Relocation verdict:** kinematic acceptance passes but every exact physical
realisation/dynamics tested either contains a global inverse, adds spectator
modes, or leaves multiple Hamiltonians unselected.  Label this **DERIVED
KINEMATIC ADVANCE / PHYSICAL GATE STILL CLOSED**, with the scope of the tested
linear constructions stated completely.

**Refutation:** any Poisson, reducibility, nilpotency, quotient-dimension or
relative-locality identity fails.

No successful outcome here derives time, causality, inertia, mass, `c`,
`hbar`, Newton's `G` or Planck units.
