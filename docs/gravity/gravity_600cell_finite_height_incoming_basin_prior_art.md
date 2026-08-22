# Prior-art gate: incoming-state branch and invariant-basin classification

Date: 2026-08-22.

Status: completed before preregistering or evaluating the global incoming-state
classification.

## Exact proposed object

Use the fixed homogeneous tetrahedral-frustum 600-cell action at zero
cosmological constant with conserved global dust and the committed canonical
pre/post momentum convention.

The already certified one-slab physical incoming domain is

```text
I=(v_A,v_star) union (v_star,v_C),
```

where `v_A`, `v_star` and `v_C` are intrinsic roots from the exact one-slab
classification, not fitted numerical cutoffs.

For every `v` in `I`, let the unique first physical slab define its exact
outgoing normalized canonical state

```text
F(v)=(m1(v),pi1(v)).
```

From that state, retain every physical root of the complete second-slab
equations.  Repeat canonical momentum matching on every surviving branch.

The proposed object is the complete branch diagram over `I`, including all
values at which:

- the number of all-real or physical successors changes;
- a branch is born, merges, loses positive height or reaches zero endpoint
  scale;
- a surviving branch enters the already certified invariant half-strip

  ```text
  D={(m,x): 0<m<=2/5, x>=125},
  x=m*q;
  ```

- a branch remains unresolved instead of either dying or entering `D`.

Once a branch enters `D`, the accepted invariant theorem supplies exactly one
physical successor at every later finite step.  The requested basin is thus
the subset of original incoming `v` whose branch tree contains a path into
`D`.  This is not a scan of selected representatives and not a fit to the
known `v=3/2` history.

## Primary sources checked

1. Dittrich and Höhn, *Canonical simplicial gravity*, Class. Quantum Grav. 29
   (2012) 115009, [arXiv:1108.1974](https://arxiv.org/abs/1108.1974).
   Their canonical framework uses a discrete action as a generating function
   and treats pre-/post-constraints, changing phase-space dimension and data
   that can be fixed a posteriori by later moves.  Branching and later
   admissibility are therefore **KNOWN structural possibilities** in
   simplicial canonical dynamics.  They do not classify the present
   homogeneous dust map or its incoming-state basin.

2. De Felice and Fabri, *The Friedmann universe of dust by Regge Calculus:
   study of its ending point*,
   [arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093), and
   *Singularities of the closed RW metric in Regge Calculus: a generalized
   evolution of the 600-cell*,
   [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077).
   They evolve dust 600-cell geometries with a Sorkin scheme and find a causal
   stopping point.  This is a load-bearing warning that forward existence
   cannot be inferred from finite samples.  Their evolution variables and
   equations are not the exact frustum canonical relation tested here.

3. Jercher and Steinhaus, *Cosmology in Lorentzian Regge calculus: causality
   violations, massless scalar field and discrete dynamics*,
   [arXiv:2312.11639](https://arxiv.org/abs/2312.11639).
   Their spatially flat cuboidal-frustum model has matter-dependent existence
   inequalities and monotone expanding/contracting branches.  It supplies a
   close structural control, but uses different cells, scalar matter and a
   different canonical map.

4. Jercher, Simão and Steinhaus, *(2+1) Lorentzian quantum cosmology from
   spin-foams: chances and obstacles for semi-classicality*,
   [arXiv:2411.08109](https://arxiv.org/abs/2411.08109).
   They report discontinuous dependence of a single-frustum classical/quantum
   observable on scalar mass and emphasize causal-sector dependence.  This is
   another control for parameter-dependent branch changes, not a
   classification of the present classical 600-cell dust relation.

## Boundary after the search

### KNOWN

- Discrete action principles generate canonical relations rather than
  automatically single-valued maps.
- Later pre-constraints can restrict data that were free at an earlier move.
- Symmetry-reduced Lorentzian Regge cosmologies can have matter-dependent
  existence regions, branch changes and causal endpoints.
- Published dust 600-cell evolutions can terminate; indefinite existence is
  not generic by assumption.

### CONTROL

- The exact one-slab domain `I` and its unique first physical update.
- The representative two-slab nonuniqueness at `v=3/2` and `v=3`, and the
  no-successor control at `v=20`.
- The complete branch-A death / branch-B survival result at `v=3/2`.
- The rigorous invariant half-strip `D` and its unique-successor theorem.

These controls may falsify a classifier after it is constructed; they may not
be used to choose its thresholds or branch labels.

### OPEN

- The complete number of physical second and later branches as a function of
  every `v` in `I`.
- All intrinsic bifurcation values and their multiplicities.
- The exact incoming basin of branches that enter `D`.
- Whether some incoming interval dies, branches forever outside `D`, or
  requires arbitrarily many pre-entry steps.
- External novelty of this exact branch diagram.

### Proposed difference

The proposed result is a complete-domain, interval-certified branch and basin
classification for one exact normalized dust canonical relation.  It is
stronger than the existing three-point representative census and narrower
than a theorem for general Regge gravity.

No checked primary source states the present equations, intrinsic thresholds
or invariant-basin classification.  Search failure is not evidence of
novelty; external novelty remains **OPEN** pending a dedicated expert review.

## Post-result search

After the 36-signature, 50-cell skeleton was known, the search was repeated
using the more specific terms `multivalued canonical relation`, `branch
birth/merger`, `a-posteriori constraint`, `causal stopping`, `600-cell dust`
and `linearized lattice gravitons`.

The same four sources above remained the closest matches. De Felice and Fabri's
generalized 600-cell evolution, [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077),
reinforces that adding variables does not remove causal stopping by
assumption. Dittrich and Höhn's later linearized analysis, *Canonical
linearized Regge Calculus: counting lattice gravitons with Pachner moves*,
[arXiv:1411.5672](https://arxiv.org/abs/1411.5672), supplies the closest
primary framework for the later nonhomogeneous physical gate, not for the
present homogeneous branch atlas.

No checked source contains the normalized scalar relation `E(m,pi,q)=0`, the
36-signature diagram or the 50 candidate cells. This is still a search result,
not proof of external novelty; novelty remains **OPEN**.

## Framing warning

A successful basin theorem would remove the arbitrary representative-state
scope of the current complete history.  It would not make infinite
extendibility a local law and would still contain no nonhomogeneous degrees of
freedom.  A failed or unresolved global classification is a first-class
result and must not be replaced by a denser plot or another chosen `v`.
