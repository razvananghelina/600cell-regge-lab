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

## Status ledger

| Claim | Status | Result |
|---|---|---|
| Unit-quaternion 600-cell and round short edges | DERIVED | f-vector `(120,720,1200,600)`, edge length `pi/5` |
| Round metric is selected by the full theory | OPEN | inherited from the embedding for this test |
| Group-difference transport has curvature | DERIVED NEGATIVE | pure gauge on all faces |
| Levi--Civita spin transport is canonical under the stated hypotheses | STRUCTURAL | unique principal lift, no fitted link coefficients |
| Nonzero, face-uniform curvature | DERIVED | exact spherical-excess holonomy on 1200/1200 faces |
| Spin-fibre algebra | DERIVED | `H` over `R`, `M2(C)` after complexification |
| Standard-Model internal algebra selected | OPEN | no comparison made; spin and internal roles are distinct |
| First geodesic refinement preserves transport and curvature | DERIVED | 720 edges and 7200 small faces pass |
| Existing `d+d*` accepts this as a flat twist | DERIVED NEGATIVE | `d_A^2` contains nonzero face curvature |
| Configuration space of connections and its Dirac operator | OPEN | not constructed |
| Representation-independent finite factor | OPEN | explicitly not guaranteed by the motivating paper |
| Continuum/refinement state or rigging map | OPEN | one background and one refinement level are insufficient |
| Time/fourth spacetime direction | OPEN | this route starts from a spatial `S^3`; none is generated here |

## Next falsifiable gate

The next calculation must not promote one background connection into a theory
of connections.  It should define, without target fitting:

1. the finite graph connection space (including the choice `SU(2)` versus
   `U(2)` and the vertex gauge quotient);
2. a measure and Hilbert space;
3. a Dirac-type operator on that configuration space;
4. refinement embeddings under which the operator and chosen states are
   cylindrically consistent.

Acceptance requires these objects to be fixed by geometry and functoriality
before inspecting an almost-commutative target.  If the operator, measure or
finite factor can only be obtained by selecting representation-dependent
coefficients after the fact, the route has reproduced the same fitting
freedom that closed earlier matter routes.
