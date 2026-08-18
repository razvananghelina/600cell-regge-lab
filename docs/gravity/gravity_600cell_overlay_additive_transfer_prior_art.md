# Prior-art gate: additive transfer on the universal staircase overlay

Date: 2026-08-17

Status: written before constructing or ranking the transfer matrix.

## 1. Exact object and complete hypothesis

Use only the certified local universal overlay of `Delta^3 x I` from commit
`94c42ad`, containing 148 open full-dimensional chambers, and the 96 coarse
four-simplices indexed by

```text
(vertex order o, staircase position k),
o in S4, k in {0,1,2,3}.
```

The hypothesis tested here is deliberately narrow:

> A scalar quantity is represented by one additive weight `x_C` on each fine
> top-dimensional chamber.  Its value on a coarse staircase simplex is the
> sum of the weights of the fine chambers contained in that simplex.

Define the exact incidence/aggregation matrix

```text
R[(o,k),C] = 1  if chamber C lies in simplex (o,k),
             0  otherwise.
```

Then `y=R x` is the vector of all 96 coarse totals.  The question is whether
compatibility with all 24 staircases, positivity and the complete local
`S4 x C2` symmetry determine `x` uniquely.

This is a test of piecewise-constant **additive top-cell transfer**.  It is not
a model of the full Regge action: deficit angles live on hinges and couple
neighbouring simplices.  A negative result here cannot prove that every
nonlocal or hinge-based dynamics is noncanonical.

## 2. Analytic nonuniqueness bounds before computation

There are 148 fine variables and at most 96 coarse equations, so

```text
nullity(R) >= 148-96 = 52.
```

For each vertex order, its four coarse simplices partition the 148 chambers.
The sum of those four rows is therefore the same all-ones row for all 24
orders.  Comparing orders with one fixed order gives 23 independent row
relations, hence

```text
rank(R) <= 96-23 = 73,
nullity(R) >= 148-73 = 75.
```

Under `S4 x C2`, the fine chambers form 14 certified orbits.  The 96 coarse
simplex labels form only two orbits: the outer positions `{0,3}` and the inner
positions `{1,2}`, with time reflection exchanging `k` and `3-k`.  Therefore
the invariant restriction maps a 14-dimensional fine space to a space of
dimension at most two:

```text
rank(R_invariant) <= 2,
nullity(R_invariant) >= 12.
```

Thus uniqueness under the stated hypotheses is already **DERIVED NEGATIVE**.
The computation will audit the premises, determine the exact ranks and emit
an explicit positive invariant pair with identical 96 coarse totals.

## 3. Primary prior art

- A common refinement (or supermesh) supports conservative transfer by
  integrating over intersections of source and target elements.  Jiao and
  Heath formulate accurate conservative transfer in
  [*Common-refinement-based data transfer between non-matching meshes in
  multiphysics simulations*](https://doi.org/10.1002/nme.1147), International
  Journal for Numerical Methods in Engineering 61 (2004), 2402--2427.
- Menon and Schmidt extend the supermesh construction to cell-centred
  polyhedral finite-volume variables in
  [*Conservative interpolation on unstructured polyhedral
  meshes*](https://doi.org/10.1016/j.cma.2011.04.025), Computer Methods in
  Applied Mechanics and Engineering 200 (2011), 2797--2804.  Their transfer
  also needs an accuracy principle, not conservation alone.
- Bahr and Dittrich construct improved/perfect actions through dynamical
  coarse graining in [*Improved and Perfect Actions in Discrete
  Gravity*](https://arxiv.org/abs/0907.4323).  A common carrier by itself is
  not that construction.
- Dittrich and Steinhaus analyze the additional path-integral measure needed
  for triangulation independence in
  [*Path integral measure and triangulation independence in discrete
  gravity*](https://arxiv.org/abs/1110.6866).
- Dittrich, Kaminski and Steinhaus prove a four-dimensional obstruction to a
  local triangulation-independent measure in linearized Regge calculus in
  [*Discretization independence implies non-locality in 4D discrete quantum
  gravity*](https://arxiv.org/abs/1404.5288).  This makes it especially
  important not to infer a local action merely from a common refinement.

No cited source claims that inclusion, conservation, positivity and finite
symmetry uniquely select fine weights.  Standard transfer methods add a norm,
shape functions, cell volumes or a variational principle.  Those are extra
data, not consequences of the incidence matrix.

## 4. KNOWN / CONTROL / OPEN

### KNOWN

- The 148-chamber overlay and all 24 unique staircase assignments are
  certified in the source artifact.
- Common-refinement incidence gives a canonical aggregation map `R` for
  piecewise-constant additive quantities.
- The dimension and partition identities above rule out an injective `R`.

### CONTROL

- Reconstruct all 96 rows from the frozen sign words rather than trusting
  printed assignment counts.
- Require each order's four rows to be a disjoint partition of 148 columns,
  with row counts `(19,55,55,19)`.
- Compute rank over the rationals and cross-check it modulo several primes.
- Reconstruct the 14 fine and two coarse `S4 x C2` orbits and rank the
  invariant restriction exactly.
- Construct two distinct, strictly positive, invariant rational weight
  vectors with identical 96 coarse totals.

### OPEN

- The exact two ranks and nullities before the registered computation.
- Whether a metric volume, Galerkin norm, path-integral measure or dynamical
  variational principle selects one point in the affine family.
- A transfer law for hinge curvature and the nonlinear Regge--dust action.
- Locality or unavoidable nonlocality of a perfect action on the global
  600-cell cylinder.

## 5. Framing consequence

The combinatorial aggregation map may be unique while its inverse is highly
nonunique.  Conflating those statements would turn the certified common
carrier into false evidence for selected dynamics.  The mission must report
both separately.

## 6. Post-result structural identification

The exact rank exposed a simpler Boolean-lattice factorization.  For a chamber
`C`, write

```text
b_A(C)=1 if h_A(C)>0, else 0,
b_empty=1, b_full=0.
```

A staircase top cell is an edge `L subset U`, `|U\L|=1`, on a maximal chain
of the Boolean lattice `B4`.  Monotonicity of subset sums gives the exact
identity

```text
R[(o,k),C] = b_L(C)-b_U(C).
```

Across all 24 orders, the 96 labelled rows reduce to the 32 Hasse edges of
`B4`: eight boundary edges occur six times and 24 internal edges occur twice.
Consequently every row lies in the 15-dimensional span of the constant
function and the 14 internal `b_A`.  The computed rational rank is exactly 15,
so this elementary upper bound is saturated.

This factorization is **DERIVED** and explains why the full 96-row data carry
far less information than their raw count suggests.  The post-result search
located standard work on positive threshold functions and Boolean lattices,
but not this particular transfer-matrix rank statement.  External novelty of
the observation remains **OPEN**.
