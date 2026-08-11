# Preregistration: can conservation select the broken-FEEC projection?

Date: 2026-08-11

This protocol is committed before evaluating any new moment-preservation
residual.  It uses no low eigenvalue, continuum target, phenomenological
number or preferred candidate spectrum.

## Question with complete hypotheses

On the fixed closed controls

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}\partial\Delta^4),
 \qquad k=1,2,4,
\]

let `V_p` be the duplicated elementwise Whitney `p`-cochains, let `W_p` be
the conforming cochains embedded by the signed occurrence injection `J_p`,
and equip `V_p` with the exact block-local Whitney mass `M_p`.

The two previously preregistered strict-occurrence-local projections are

\[
 P^C_p=J_pL^C_p
 \quad\hbox{and}\quad
 P^D_p=J_pL^D_p,
\]

where `C` is equal counting and `D` is diagonal-Whitney weighting.  Both
satisfy `LJ=I`, have image `W_p`, form exact broken-FEEC complexes, and give
the same Betti vector `(1,0,0,1)`.  They are distinct in degrees one and two
after nontrivial refinement.

The frozen question is whether either of two independently intelligible
conservation requirements uniquely selects one of these local projections.

## Conservation notion I: every conforming moment

For every broken field `v` and every conforming field `w`, require

\[
 \langle Pv,w\rangle_M=\langle v,w\rangle_M.
\]

In matrices this is

\[
 P^TMJ=MJ.
\]

This is deliberately stronger than preservation of topology.  Before the
numerical audit, the following finite-dimensional lemma fixes its logical
content.

**Lemma.**  If `P^2=P` and `im(P)=im(J)`, preservation of every conforming
moment is equivalent to `P` being the `M`-orthogonal projector onto
`im(J)`.  Hence it uniquely fixes

\[
 P^A=J(J^TMJ)^{-1}J^TM.
\]

**Proof frozen before computation.**  The moment identity says
`<Pv-v,w>_M=0` for every `w` in `im(J)`.  Since `Pv` lies in `im(J)`, it is
the unique orthogonal projection of `v` onto that subspace.  Conversely the
orthogonal projection has this defining property.  Positive definiteness of
`M` and injectivity of `J` make `J^TMJ` invertible.

The verifier will evaluate the defect

\[
 \Delta^X_p=(L^X_p)^T(J_p^TM_pJ_p)-M_pJ_p,
 \qquad X\in\{C,D\},
\]

without forming a dense local-by-local projector.  It will record the
maximum absolute defect and its nonzero count at `k=1,2,4` for degree zero,
which already suffices to falsify an all-degree claim.  It will also rebuild
over exact rationals the previously observed `k=1,p=0` off-occurrence
coefficient of `L^A`, rather than merely importing its value.

Decision:

- exact nonzero off-occurrence support of `L^A` gives a **DERIVED LOCALITY
  CONFLICT** for full-moment selection;
- a local candidate passes full-moment conservation only if its defect is
  below `1e-11` on every frozen control, but this numerical gate alone will
  not be promoted to an all-level theorem;
- the algebraic uniqueness lemma is **DERIVED**, independent of those
  controls.

## Conservation notion II: only harmonic/topological moments

Let `H_p` be the conforming harmonic space for the assembled Whitney
complex.  Require the same identity only for `w` in `J_pH_p`:

\[
 \langle Pv,Jh\rangle_M=\langle v,Jh\rangle_M
 \quad\hbox{for all }h\in H_p.
\]

The exact topology of every frozen control is `S^3`, so

\[
 (\dim H_0,\dim H_1,\dim H_2,\dim H_3)=(1,0,0,1).
\]

The test is therefore completely enumerated rather than sampled:

1. in degree zero, use the global constant cochain and evaluate
   `(P^TM-M)J1` for both candidates at `k=1,2,4`;
2. degrees one and two contain no harmonic test vector;
3. in degree three, verify that every top simplex has one local occurrence,
   so both projections are the identity and preserve the entire space, not
   just its one harmonic line.

No numerical eigensolver is needed to choose a harmonic vector.

Decision:

- if exactly one candidate passes all nonvacuous harmonic tests, report a
  **STRUCTURAL SELECTION BY THE PROPOSED AXIOM**, not a derived physical law;
- if both pass, report a **DERIVED NEGATIVE FOR TOPOLOGICAL-MOMENT
  UNIQUENESS ON THE FROZEN TOWER**;
- if neither passes, report that harmonic-moment conservation does not
  validate either local candidate;
- the `1e-11` residual gate is frozen, and exact combinatorial identities
  will be used wherever available.

## Hodge/Poincare framing audit

The current primal cochain construction does not itself contain a
primal-to-dual Hilbert space, a dual metric, or a chain-level Hodge star.
Moreover its primal dimensions are not paired degreewise.  A Hodge-star
selection argument must therefore state new dual-carrier and metric data.

The result note will distinguish:

- homological Poincare duality, already expressed by the Betti vector and
  unable by itself to distinguish projections with the same cohomology;
- a metric Hodge star, which may impose stronger conditions but is **OPEN**
  until its dual carrier and metric are derived rather than selected after
  seeing a spectrum.

This is not a theorem that every possible primal-dual construction fails.
It is a scope statement about the data currently present.

## Attack on the framing

The phrase "conservation law" can overstate both tests.  Full conforming
moment preservation is the variational characterization of an orthogonal
projection; it is not yet derived from time translation or a physical
Noether symmetry.  Harmonic-moment preservation protects global topological
charges, but says little about positive-energy dynamics.  A successful
matrix test therefore cannot by itself supply time, causality, inertia,
mass, or Planck units.

The strongest admissible outcome is a precise selection trilemma:

1. full metric information selects uniquely but may be nonlocal;
2. topological information may remain local but may not select uniquely;
3. stronger local selection needs an additional independently derived axiom.

## Frozen outputs and scope

The registered verifier will write
`reproducible/whitney_projection_selection.json` and the result note will
record every residual, exact witness and label.

Excluded:

- no candidate spectrum;
- no fit or target comparison;
- no continuum-limit theorem from three finite controls;
- no newly chosen Hodge star or dual mesh;
- no full-suite run, by explicit user request.
