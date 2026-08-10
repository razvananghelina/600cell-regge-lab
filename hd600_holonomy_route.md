# HD-600: canonical spin holonomy on the 600-cell spatial slice

Date: 2026-08-10  
Verifier: `reproducible/verify_hd600_holonomy.py`

## Decision

**STRUCTURAL ADVANCE, NOT A PHYSICS GATE.**  If the embedded boundary of the
600-cell is interpreted as a sampling of the unit round three-sphere, its
short geodesic edges select a Levi--Civita spin transport.  This transport has
nonzero curvature, is equivariant under the orientation-preserving 600-cell
symmetry group, generates the full spin-fibre algebra
`H tensor_R C = M2(C)`, and is consistent with the first canonical geodesic
refinement.

Two stronger readings are false:

- **DERIVED NEGATIVE:** the simpler group-difference links are pure gauge.
- **DERIVED NEGATIVE:** the curved Levi--Civita links are not a flat
  local-system twist of the certified cochain differential.  A covariant
  Kähler--Dirac operator would be a new operator, not the old `d+d*` with an
  unnoticed flat coefficient system.

This does not yet realise the holonomy--diffeomorphism proposal.  It gives one
distinguished background connection, whereas that proposal is built on a
configuration space of connections and a Dirac-type operator on that
configuration space.  It also does not select a Standard-Model internal
algebra.

The first representation gate is already negative: on the doubled spin fibre
`S+S`, ordinary spin geometry selects the diagonal action `U -> diag(U,U)`.
The resulting algebra is still only `M2(C)`, has a four-complex-dimensional
commutant, and cannot distinguish states which differ only by their copy
label.  Enlarging it requires a new representation choice rather than more
600-cell geometry.

The projective connection-space gate is mixed.  Normalised Haar measure is
canonical and cylindrically consistent, but edge subdivision leaves the scale
of every genuinely new connection mode free.  Consequently, projective
compatibility alone does not select a configuration-space Dirac operator.
Round Whitney `L2` geometry fixes that scale and nests exactly, but its
constant edge Gram matrix fails local gauge invariance.  No positive diagonal
product-link replacement can retain tetrahedral symmetry and exact
cylindrical pullback.  The remaining load-bearing object is therefore a
connection-dependent covariant metric.  A positive local candidate is now
constructed below; its curved refinement compatibility remains open.

## Complete hypotheses

Every canonicity claim below assumes all of the following.

1. The 120 vertices are the usual unit-quaternion embedding of the regular
   600-cell.
2. The carrier is the unit round `S^3`, with the metric inherited from the
   ambient Euclidean four-space.
3. An edge is traversed on its unique short round geodesic.  Its length is
   `theta=arccos(phi/2)=pi/5`, so there is no antipodal ambiguity.
4. Spin transport uses the unique lift continuous from the identity along
   that short geodesic (the principal quaternion square root).
5. Spinors are written in a global left-invariant orthonormal frame.

Hypothesis 2 is **STRUCTURAL**, not derived from the abstract face poset alone.
The embedded regular polytope supplies it naturally, but the current theory
does not prove that this metric is dynamically selected.  Subject to these
hypotheses the connection is rigid; there are no fitted edge parameters.

The left-frame choice is not physical freedom.  The right-invariant frame
gives a vertexwise gauge-equivalent connection.  More generally, under
`q -> a q b^-1`, the left-frame links change by the constant spin-fibre gauge
`b`.

## Two candidate transports

For unit quaternions `q_i,q_j`, column-vector transport from `i` to `j` is
composed from right to left.

### Pure-gauge control

The inverse-convention group-difference link is

```text
U_flat(i,j) = q_j^-1 q_i.
```

The product around every triangular face is exactly the identity.  Therefore
this construction contains no curvature and no local geometric signal.  Its
failure is a useful control: group multiplication alone is insufficient.

### Levi--Civita spin transport

In the global left frame, the round Levi--Civita transport is

```text
U_LC(i,j) = sqrt(q_i^-1 q_j)^-1,
```

where `sqrt` is the principal unit-quaternion square root.  The right-frame
formula is gauge related by

```text
U_R(i,j) = q_j U_LC(i,j) q_i^-1.
```

Both reversal by inversion and this frame relation are checked directly on
all 720 edges.

## Exact face holonomy

Each of the 1200 faces is a round equilateral spherical triangle of side
`pi/5`.  If `alpha` is its interior angle, the spherical cosine rule gives

```text
cos(alpha) = 1/sqrt(5),
E = 3 alpha - pi,
```

where `E` is the spherical excess.  Tangent holonomy rotates by `E`; its spin
lift has half-angle `E/2`.  Every face holonomy consequently obeys

```text
Re(H)^2   = (25 + 11 sqrt(5))/50,
|Im(H)|^2 = (25 - 11 sqrt(5))/50.
```

Numerically, `Re(H)=0.995959313990` and
`|Im(H)|=0.089805595319`.  All faces have this conjugacy class and none has
identity holonomy.

**DERIVED:** two face holonomies fail to commute.  The real algebra generated
by them has quaternion rank four, and its complexification is `M2(C)`.

**Scope warning:** this is the algebra acting on a two-component spin fibre.
It is not a selected finite internal algebra, and it is not evidence for
`C + H + M3(C)`.  Calling it an internal sector would conflate spatial spin
transport with particle gauge structure.

## First refinement

Radially normalised barycentres of all cells give

```text
120 + 720 + 1200 + 600 = 2640
```

distinct points.  Their cover-relation graph has
`2*720 + 3*1200 + 4*600 = 7440` edges.  This is the Hasse carrier used by the
oriented cochain construction; it must not be confused with the full
barycentric one-skeleton, which also joins non-consecutive comparable cells.

Two refinement checks are **DERIVED**:

1. On every coarse edge, transport factors through its normalised geodesic
   midpoint with maximum residual `4.8e-12`.
2. Splitting every face into its six congruent geodesic barycentric triangles
   gives 7200 small faces.  Every small holonomy has spin angle `E/12`, hence
   carries exactly one sixth of the coarse curvature.  The measured values are
   `Re(H_small)=0.999887685165` and
   `|Im(H_small)|=0.014987230528`.

This is one-level cylindrical consistency for the fixed Levi--Civita
background.  It is not convergence of a quantum configuration-space measure.

## Relation to the existing Kähler--Dirac operator

The certified operator in `verify_kahler_dirac.py` is the untwisted cochain
operator `D=d+d*`, with `d^2=0`.  A flat local system preserves this complex
property.  The Levi--Civita spin links do not: on each triangular two-cell the
covariant square contains `H_f-I`, which is nonzero.

This is not a defect of differential geometry: a covariant exterior
derivative with curvature normally satisfies `d_A^2=F_A`.  It is instead a
boundary on what has been established here:

- **DERIVED NEGATIVE:** the new links cannot be inserted while claiming that
  the already-certified cochain complex is unchanged.
- **OPEN:** construct and certify a genuinely covariant discrete Dirac or
  Kähler--Dirac operator, including its Hilbert-space weights, adjoint and
  refinement maps.
- **DERIVED:** because every face holonomy is a nonidentity `SU(2)` element,
  it has no eigenvalue one; the round connection has no parallel spinor on
  this carrier.

## Projective connection-space gate

### Hypotheses

The following negative result assumes only:

1. one compact-group element per unoriented graph edge;
2. the standard vertex gauge action;
3. edge subdivision projected by multiplication
   `m(g1,g2)=g1*g2`;
4. an `Ad`-invariant metric on each Lie-algebra component;
5. symmetry under exchanging the two equal subedges;
6. cylindrical isometry for modes pulled back from the coarse edge.

It does not assume a particle target or choose a spectral scale.

The connected 600-cell graph has cycle rank

```text
E-V+1 = 720-120+1 = 601.
```

For `SU(2)`, the edge connection space has dimension 2160, the based quotient
has dimension 1803, and the full generic quotient has dimension 1800.  For
`U(2)` the corresponding numbers are 2880, 2404 and 2401.  The last number is
not 2400 because the constant central `U(1)` is a continuous generic
stabiliser.  Thus passing from the selected spatial `SU(2)` spin connection
to the paper's `U(2)` configuration space is a real enlargement.

### Measure: positive result

**DERIVED:** normalised Haar measure is projectively consistent.  The
push-forward of `dg1 dg2` by multiplication is `dg`, by invariance and
uniqueness of normalised Haar measure.  The verifier includes an exhaustive
finite-group control: every element of a group of order 120 has exactly 120
preimages under multiplication.

### Metric and Dirac: underdetermination theorem

Linearise multiplication at the identity:

```text
dm(X,Y) = X+Y.
```

On one Lie-algebra component, the most general exchange-symmetric invariant
quadratic form is

```text
M = [[a,b],[b,a]].
```

Let the coarse metric coefficient be `c`.  Requiring the horizontal lift
`Z -> (Z/2,Z/2)` to be isometric imposes only

```text
a+b = 2c.
```

Writing `delta=a-b>0` leaves

```text
a = c+delta/2,
b = c-delta/2,
||(X,-X)||^2 = 2 delta |X|^2.
```

The new vertical mode `(X,-X)` is killed by `dm`, so no coarse observable can
fix `delta`.  Three explicit choices `delta/c=1,2,3` are positive,
exchange-symmetric and exactly identical on all inherited coarse modes.  The
product metric `b=0` is only the middle choice; demanding it is an additional
locality hypothesis.

For the local maximal-torus model, the Hodge--Dirac square has Fourier
eigenvalue `k^T M^-1 k`.  Therefore

```text
lambda_coarse  for k=(1,1)  = 1/c,
lambda_new     for k=(1,-1) = 2/delta.
```

All inherited spectra intertwine, while the new spectrum varies.  This proves
the precise statement:

**DERIVED NEGATIVE:** subdivision, symmetry and cylindrical compatibility do
not select the configuration-space Dirac operator.  They leave at least one
positive scale for each new refinement sector.

This agrees with the earlier projective spectral-triple construction, which
states explicitly that its Dirac operators form a large class labelled by a
sequence of real parameters `{a_n}` and that a scaling behaviour must be
chosen ([arXiv:0802.1783](https://arxiv.org/abs/0802.1783)).  The 2025 paper
does not derive away this freedom; it assumes a gauge-covariant metric on the
configuration space and notes remaining variation in the Dirac operator.

The next subsection closes the cylindrical-consistency part of that question:
the round `L2` metric on radially pulled-back Whitney one-forms does fix the
vertical sector.  What remains open is whether the theory selects that
finite-element Hilbert space.  Merely choosing the product metric would still
be fitting by convenience.

### Whitney `L2` candidate and its limitation

There is a mathematically natural candidate once a piecewise-flat metric is
chosen: integrate the inner products of Whitney one-forms exactly.  On a
reference tetrahedron and its full barycentric subdivision, the verifier
constructs

```text
coarse:  4 vertices, 6 edges, 1 tetrahedron
fine:   15 vertices, 50 edges, 24 tetrahedra.
```

Let `P` record the integral of every coarse Whitney basis form over every fine
edge, and let `M_coarse`, `M_fine` be the exact rational mass matrices.
Direct assembly gives

```text
P^T M_fine P = M_coarse
```

exactly over `Q`.  By affine naturality the same local identity holds on each
nondegenerate piecewise-flat tetrahedron.  Thus:

- **DERIVED:** Whitney `L2` geometry supplies a cylindrically consistent
  metric, including the new vertical modes.
- **DERIVED NEGATIVE:** the counting inner product used by the current
  unweighted Kähler--Dirac operator is not cylindrically consistent.  Its
  pullback `P^T P` has two distinct eigenvalues `7/9` and `125/72`, so no
  overall rescaling repairs it.

The exact flat calculation uses **affine** barycentres inside a reference
tetrahedron.  Gate 0 used radially normalised barycentres and round geodesics
on `S^3`.  These physical point sets are distinct.  For a regular 600-cell
tetrahedron the affine barycentre has norm

```text
r = sqrt(7+3 sqrt(5))/4,
```

so radial normalisation moves it out of the affine facet hyperplane by

```text
1-r = 0.074385206589...
```

**DERIVED:** the affine-Whitney and round-geodesic refinements are distinct.
That distinction, however, is **not** an obstruction.  Define a smooth round
simplex map on the reference tetrahedron by

```text
F(lambda) = (sum_i lambda_i q_i) / |sum_i lambda_i q_i|.
```

Its edges are the short great-circle arcs.  Restricting `F` to the 24 affine
barycentric sub-tetrahedra maps their barycentres to exactly the radially
normalised points used at Gate 0.  The verifier proves pointwise, on every
reference subcell, that every coarse Whitney one-form is the `P`-weighted sum
of the fine Whitney forms.  This identity is independent of the metric.
Pulling all forms through the common map `F` and integrating with the round
metric therefore preserves

```text
P^T M_fine(round) P = M_coarse(round).
```

**DERIVED ADVANCE:** round Levi--Civita refinement and cylindrically nested
Whitney spaces are compatible when both levels are defined as pullbacks from
one smooth radial simplex map.  The earlier apparent affine/round bifurcation
was a framing error, now corrected.

This construction matches the modern intrinsic-manifold FEEC framework:
finite-element forms are pulled from Euclidean reference simplices through a
smooth triangulation, and commuting projections provide a stable discrete de
Rham complex ([Licht, arXiv:2310.14276](https://arxiv.org/abs/2310.14276),
updated 2026).

**STRUCTURAL/OPEN:** selecting lowest-order Whitney forms still changes the
theory's current unweighted cochain Hilbert space.  The current axioms do not
yet establish why this finite-element family, rather than a higher-order or
different discrete Hodge star, is fundamental.  Lowest order is the minimal
choice with exactly one degree of freedom per edge and no added nodes, but
minimality is not yet a proved physical selection principle.

### Local-gauge obstruction

The Whitney result above concerns the metric on Lie-algebra-valued one-forms.
It does not automatically give a gauge-invariant metric on the finite product
of link holonomies.  State the additional hypotheses explicitly:

1. a graph connection is represented by group elements on oriented edges;
2. tangent vectors are left-trivialised at each edge;
3. the Whitney Gram matrix is constant over the connection space;
4. vertex gauge transformations act independently.

In left trivialisation, an edge tangent based at source `s` transforms by
`Ad(g_s)`.  A constant cross term between edges with different sources is
therefore not invariant under independent `g_s`.

The exact tetrahedral Whitney mass contains the witness

```text
M[(1,2),(2,3)] = -1/120.
```

Rotating the Lie-algebra vector at source 2 from `Z` to `-Z`, while leaving
source 1 fixed, changes the quadratic norm by `1/30`.

Could one retain gauge invariance by using a positive diagonal product metric
on the 50 fine links?  Tetrahedral symmetry leaves six weights, indexed by the
dimensions of the nested faces joined by a fine edge.  Solving

```text
P^T diag(w) P = I_6
```

exactly forces

```text
w_(1,3) = -9 w_(1,4)/16 - w_(2,3)/4
          -9 w_(2,4)/32 - w_(3,4)/16.
```

It is negative whenever the other free weights are positive.  Averaging any
hypothetical positive diagonal solution over the tetrahedral group would give
a positive orbit-constant solution, so this excludes nonsymmetric positive
diagonal solutions as well.

**DERIVED NO-GO, with complete hypotheses:** a connection-independent
positive product-link metric cannot be simultaneously locally gauge
invariant, tetrahedrally symmetric and exactly cylindrical with the Whitney
cochain injection.

This does not kill the configuration-space route.  It proves that the metric
must be **connection dependent**, using parallel transport to compare tangent
vectors attached to different sources.  Such a fibrewise gauge-covariant
metric is precisely an input assumed in the 2025 HD construction.

There is a local tetrahedrally symmetric construction.  For each tetrahedron
and each of its four possible base vertices `r`:

1. transport every edge tangent from its source to `r` along the unique direct
   tetrahedral edge;
2. evaluate the Whitney quadratic form on the transported vectors;
3. average the four results.

Each summand is positive because the exact Whitney mass has spectrum

```text
1/6 (x1), 1/24 (x4), 1/60 (x1)
```

and adjoint transport is orthogonal.  The average is therefore positive and
does not privilege a vertex.  Under a local gauge transformation all vectors
transported to `r` rotate by the same `Ad(g_r)`, proving gauge invariance.
Forty deterministic numerical convention checks give a maximum residual
`1.78e-15`; the algebraic argument is the evidence, not the sampling.

At the identity connection this metric reduces exactly to the Whitney mass.
Gauge covariance extends the result to every pure-gauge/flat configuration,
so exact cylindrical consistency is **DERIVED on the flat sector**.

**OPEN:** on a curved configuration, coarse holonomy is a noncommutative
product of fine holonomies and the tangent injection acquires adjoint factors.
It has not been proved that the basepoint-averaged metric is cylindrically
consistent there.  Passing this test, rather than mere gauge invariance, is
the next acceptance boundary.

## Why this route was tested

Aastrup and Grimstrup's 2025 construction starts with a configuration space of
gauge connections on a spatial three-manifold and an algebra generated by
parallel transports along flows.  In a semiclassical limit it produces a
spatial Dirac operator and an almost-commutative structure; importantly, its
finite factor depends on the representation and on the state/localisation
point ([arXiv:2504.03391](https://arxiv.org/abs/2504.03391)).

That makes it relevant to the repository's independently established
three-dimensional spatial geometry, but it does not supply the missing data
for free.  In particular, our fixed Levi--Civita connection is only one point
of the required configuration space, and `M2(C)` above is only the spin-fibre
algebra.

## First HD representation gate

The 2025 construction represents its holonomies on `S+S`.  Its own
Stone--Weierstrass discussion notes that the diagonal `U(2)` representation on
the two copies never separates all pure states, while a four-dimensional
irreducible representation may do so.  The latter is explicitly a
representation choice, not a conclusion forced by the underlying
three-manifold.

The 600-cell calculation makes the obstruction concrete:

```text
natural action       U -> diag(U,U)
generated algebra    {diag(A,A): A in M2(C)}
complex dimension    4              (M4(C) would have 16)
commutant dimension  4
```

For example, a fixed spin-up state in copy one and the same spin-up state in
copy two have identical expectation values on every face holonomy.  Hence the
natural doubled algebra does not separate these two pure states.

**DERIVED NEGATIVE:** the natural geometric representation fails the paper's
state-separation gate.  This does not refute the HD framework, because a
different representation can be imposed.  It does refute the claim that the
600-cell Levi--Civita geometry has already selected the larger finite factor.

There is a second missing input.  Levi--Civita spin transport is `SU(2)`-valued
(`det U=1`); the paper specialises to the configuration space of `U(2)`
connections.  A `U(1)` connection can be added, but none has been selected by
the construction tested here.  The repository's Hopf-bundle data do not by
themselves prove that this `U(1)` should be combined with the spatial spin
connection in the HD configuration space.

## Status ledger

| Claim | Status | Result |
|---|---|---|
| Unit-quaternion 600-cell and round short edges | DERIVED | f-vector `(120,720,1200,600)`, edge length `pi/5` |
| Round metric is selected by the full theory | OPEN | inherited from the embedding for this test |
| Group-difference transport has curvature | DERIVED NEGATIVE | pure gauge on all faces |
| Levi--Civita spin transport is canonical under the stated hypotheses | STRUCTURAL | unique principal lift, no fitted link coefficients |
| Nonzero, face-uniform curvature | DERIVED | exact spherical-excess holonomy on 1200/1200 faces |
| Spin-fibre algebra | DERIVED | `H` over `R`, `M2(C)` after complexification |
| Natural diagonal action on `S+S` | DERIVED NEGATIVE | algebra dimension 4, not 16; pure states not separated |
| Alternative four-dimensional irreducible action | OPEN | allowed by the HD paper but not selected here |
| `U(1)` factor needed for `U(2)` connection space | OPEN | absent from Levi--Civita `SU(2)` transport |
| Standard-Model internal algebra selected | OPEN | no comparison made; spin and internal roles are distinct |
| First geodesic refinement preserves transport and curvature | DERIVED | 720 edges and 7200 small faces pass |
| Existing `d+d*` accepts this as a flat twist | DERIVED NEGATIVE | `d_A^2` contains nonzero face curvature |
| Finite graph connection space | DERIVED | `SU(2)^720` and generic gauge quotient dimensions certified |
| Configuration-space Dirac operator | OPEN | a compatible family exists, but no member is selected |
| Normalised Haar measure under edge subdivision | DERIVED | exactly projectively consistent |
| Horizontal configuration metric | DERIVED | fixed by the coarse metric |
| New vertical refinement scale | DERIVED FREE | arbitrary positive `delta` |
| Dirac selected by projective compatibility alone | DERIVED NEGATIVE | new eigenvalue `2/delta` varies |
| Round Whitney `L2` geometry fixes vertical scales | DERIVED CONDITIONED | exact if the Whitney Hilbert space is adopted |
| Affine Whitney `L2` mass under barycentric refinement | DERIVED | exact local pullback identity |
| Current unweighted cochain metric under refinement | DERIVED NEGATIVE | eigenvalues `7/9`, `125/72` |
| Round physical points equal flat affine points | DERIVED NEGATIVE | centre displacement `0.074385...` |
| Round radial Whitney spaces nest under refinement | DERIVED | pointwise restriction identity on all 24 subcells |
| Round `L2` mass is cylindrical | DERIVED | follows from pointwise nesting and partitioned integration |
| Theory selects lowest-order Whitney Hilbert space | OPEN/STRUCTURAL | minimal, but not forced by current axioms |
| Constant Whitney link metric is locally gauge invariant | DERIVED NEGATIVE | cross-source coefficient `-1/120` changes the norm |
| Positive diagonal link metric can replace it | DERIVED NO-GO | exact six-orbit system forces a negative weight |
| Basepoint-averaged covariant Whitney metric | DERIVED LOCAL | positive and gauge invariant; no preferred tetra vertex |
| Covariant metric on flat/pure-gauge refinement sector | DERIVED | reduces gauge-covariantly to exact Whitney nesting |
| Covariant metric on curved refinement sector | OPEN | noncommutative tangent pullback not certified |
| Representation-independent finite factor | OPEN | explicitly not guaranteed by the motivating paper |
| Continuum/refinement state or rigging map | OPEN | one background and one refinement level are insufficient |
| Time/fourth spacetime direction | OPEN | this route starts from a spatial `S^3`; none is generated here |

## Next falsifiable gate

The finite graph connection space, Haar measure and a cylindrically nested
round Whitney metric are now available.  The next calculation must define,
without target fitting:

1. why the connection group is `SU(2)` or `U(2)`;
2. curved-sector cylindrical consistency of the positive basepoint-averaged
   covariant Whitney metric;
3. the resulting global Clifford module and self-adjoint Dirac-type operator;
4. refinement embeddings for spinors/states, not only for one-forms;
5. a semiclassical state and a representation of the HD algebra selected
   before examining its finite factor.

Acceptance requires these objects to be fixed by geometry and functoriality
before inspecting an almost-commutative target.  If the operator, measure or
finite factor can only be obtained by selecting representation-dependent
coefficients after the fact, the route has reproduced the same fitting
freedom that closed earlier matter routes.
