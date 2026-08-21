# Prior-art and framing gate: canonical selection after two-slab branching

Date: 2026-08-21.

Input result: `docs/gravity/gravity_600cell_finite_height_composition_result.md`.

Status: written after the two physical continuations were known, but before
evaluating any proposed selector on either continuation.  This is therefore a
post-result falsification audit, not blind preregistration relative to the
existence or numerical location of the roots.

## 1. Question with complete hypotheses

On the fixed positive-Lorentzian homogeneous tetrahedral-frustum action, with
zero cosmological constant, conserved global dust, the committed canonical
pre/post momentum convention and the accepted first state `v=3/2`, two
isolated positive-height second slabs solve the same incoming canonical data.

Ask whether any selector already implied by the action and its geometric
domain admits exactly one of these solutions:

1. future time orientation;
2. causal regularity of the tetrahedral frustum;
3. local regularity of the discrete Legendre relation;
4. membership in the connected real branch on which the cellular action was
   derived.

No numerical cutoff, least-change norm, preferred coordinate, chosen lapse,
minimum-action rule or comparison with a desired cosmology is admissible.

## 2. Primary prior art

- Marsden and West, [*Discrete mechanics and variational
  integrators*](https://doi.org/10.1017/S096249290100006X), define the
  discrete pre/post Legendre transforms and obtain the evolution equation as
  momentum matching.  Their regularity condition supplies local inversion;
  a global inverse is stronger data.  West's thesis states the local
  isomorphism condition explicitly in Theorem 2.10:
  <https://thesis.caltech.edu/2492/>.
- Dittrich and Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974), formulate additive Regge actions
  as generating functions for canonical evolution.  They allow data to be
  fixed by later consistency conditions, but do not provide a general rule
  selecting one of several isolated global solutions of the same discrete
  Legendre relation.
- Dittrich and Hoehn, [*From covariant to canonical formulations of discrete
  gravity*](https://arxiv.org/abs/0912.1817), show that broken discrete gauge
  symmetry produces lapse-dependent pseudo-constraints.  Local nondegeneracy
  and global uniqueness are different questions in that framework.
- Jercher and Steinhaus, [*Cosmology in Lorentzian Regge
  calculus*](https://arxiv.org/abs/2312.11639), distinguish a Lorentzian
  top-dimensional frustum from causally regular lower-dimensional cells.  In
  their cubical carrier, causal regularity constrains central height relative
  to spatial edge change.  That inequality is carrier-specific and may not be
  imported into the present tetrahedral carrier.

Thus a causal selector is a legitimate question, while treating a nonzero
Jacobian as global uniqueness is not.

## 3. Geometry-specific correction

The present variable

```text
rho=h^2
```

is the positive square of the proper length of each timelike same-vertex
strut, not the squared central-coordinate height used in the cubical formula.
The already-certified tetrahedral embedding has

```text
R=phi*L,
T^2=rho+(R_plus-R_minus)^2,
strut^2=-rho.
```

For `L_plus=L_minus+h*q`, the physical coordinate speed of a vertex is

```text
beta(q)^2 = phi^2*q^2/(1+phi^2*q^2) < 1
```

for every finite real `q`.  Consequently, copying the cubical-looking bound
`|q|<constant` would double-count a causal restriction already built into
the proper-strut parametrisation.

The direct cellular action and all its angles were previously certified on
the connected domain

```text
L_minus>0, L_plus>0, rho>0.
```

Any stronger causal exclusion must therefore be newly derived from the
tetrahedral light-cone geometry; it cannot be inferred from the magnitude of
`q` alone.

## 4. Exact regularity certificate to test

For arbitrary normalized incoming data `(m,pi)`, the second-slab equations
are

```text
C=8*pi[mu(q)-m]+4*pi*h*q*mu(q),
P=p(q)-pi-2*pi*h*mu(q).
```

The already-derived state functions obey

```text
4*pi*mu'(q)+q*p'(q)=0.
```

Before inserting either root, expand the two-variable determinant and test
the exact identity

```text
det partial(C,P)/partial(h,q) = 8*pi^2*h*mu(q)^2.
```

If true, every positive-height root is locally regular.  The existence of
two such roots then proves that the discrete Legendre relation is locally
regular on both sheets but not globally injective at the frozen state.

## 5. Frozen outcome hierarchy

### `ACTION_DERIVES_CANONICAL_BRANCH_SELECTOR`

Use only if a condition already implied by the action or exact tetrahedral
causal geometry admits exactly one physical branch and rejects the other,
without a numerical cutoff or fitted convention.  Because the candidate was
formulated after the roots were known, initially label a positive hit
**PATTERN** pending a mechanically independent derivation.

### `STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH`

Use **DERIVED NEGATIVE, selector-scoped** if both roots have future
orientation, lie in the same certified causal/action domain and have positive
exact regularity determinant.  This excludes only the four selectors listed
in Section 1.  It does not prove that no additional matter, global
extendibility condition, path-integral contour or nonhomogeneous equation can
select a branch.

### `SELECTOR_AUDIT_OPEN`

Use **OPEN** if the causal domain, determinant identity or provenance cannot
be certified without adding a convention.

Only a targeted verifier may be run.  The full suite is excluded.
