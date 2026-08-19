# Prior-art gate: fixed-frame gluing of two cellular frusta across one face

Date: 2026-08-19

## Exact object and complete hypotheses

Take two nondegenerate spacelike tetrahedra sharing one triangular face.
Develop their union into one Minkowski frame, so the lower vertices are

```text
p0,p1,p2  shared face,
p3        left apex,
p4        right apex.
```

The apices lie on opposite sides of the shared face.  Use the homothetic
upper placement

```text
q_i = lambda p_i + tau n,
```

with timelike unit normal `n` and `tau!=0`.  For each of the two tetrahedral
frusta separately, hold fixed its six upper intrinsic lengths and its four
corresponding strut lengths.  The accepted local theorem gives a
six-dimensional relative-Poincare kernel `K_L` and `K_R`.

This mission fixes the background face development and imposes physical
vertex gluing:

```text
delta q_i^L = delta q_i^R,   i=0,1,2.
```

It asks:

> Does any relative local motion survive on the shared face, or do all
> compatible pairs reduce to one common six-dimensional motion?

This is a fixed-frame metric-gluing question.  It does not introduce an
independent face holonomy, connection or frame-transition variable.

## Disclosed algebraic warning

In the full ten-dimensional Poincare algebra, a nondegenerate spacelike
triangle has a one-dimensional pointwise stabilizer: the Lorentz boost in
the two-plane normal to the triangle, with the translation that places its
fixed plane on the triangle.  Therefore a correct implementation must find

```text
dim{Poincare Killing fields vanishing on the face} = 1
```

as a positive control.

However, the strut-preserving homothetic kernel is smaller.  For
`lambda!=1`, its elements obey the tetrahedron-independent relation

```text
b(A) = tau/(lambda-1) A n
```

in a common homothety-centred development.  At `lambda=1`, it contains only
`A n=0` and `<b,n>=0`.  These additional conditions may have trivial
intersection with the pointwise face stabilizer.  The disclosed prediction
is therefore:

```text
full Poincare face stabilizer                 dimension 1,
strut-preserving relative face stabilizer     dimension 0,
glued pair space                              dimension 6.
```

No two-frustum gluing matrix has been evaluated while writing this gate.

## What the primary literature establishes

### KNOWN

- First-order and connection phase spaces contain variables beyond length
  Regge data.  Gluing/metricity constraints are required to recover a Regge
  geometry: B. Dittrich and J. P. Ryan,
  [*Phase space descriptions for simplicial 4d
  geometries*](https://arxiv.org/abs/0807.2806).
- Twisted geometries attach independent normals and extrinsic-angle data to
  shared faces; without shape matching, adjacent polyhedra do not glue to a
  common Regge metric: L. Freidel and S. Speziale,
  [*Twisted geometries: A geometric parametrisation of SU(2) phase
  space*](https://arxiv.org/abs/1001.2748).
- In Lorentzian cellular decompositions, secondary simplicity constraints
  can impose shape matching and recover discrete extrinsic geometry from a
  holonomy.  This requires dynamical constraints, not just equality of
  lengths: F. Anza and S. Speziale,
  [*A note on the secondary simplicity constraints in loop quantum
  gravity*](https://arxiv.org/abs/1409.0836).
- Poincare 2-group/higher-gauge approaches retain Lorentz and translational
  geometric data as distinct parts of a larger connection framework:
  S. K. Asante et al.,
  [*Quantum geometry from higher gauge
  theory*](https://arxiv.org/abs/1908.05970), and
  A. Mikovic and M. Vojinovic,
  [*Poincare 2-group and quantum gravity*](https://arxiv.org/abs/1110.4694).

These works make the scope boundary sharp: failure of a relative mode in
the present fixed-frame length kernel does not prove that connection
formulations have no face holonomy.  It proves only that such a holonomy is
not already hidden in these six cellular flexes.

### CONTROL

- the accepted one-frustum six-dimensional Poincare stratification;
- the one-dimensional pointwise stabilizer of a spacelike triangle in the
  full Poincare algebra;
- exact agreement between parameter-space gluing and direct equality of all
  three shared vertex displacements.

### OPEN

- the dimension of the glued pair kernel after both frusta's strut
  constraints;
- whether every compatible pair is diagonal in a common developed frame;
- whether varying the frame transition as a new variable restores a
  connection/holonomy mode;
- closure around an edge, torsion-free conditions, shape matching and the
  global 600-cell carrier;
- action, symplectic structure, dynamics and continuum physics.

No primary source found in the 2026-08-19 pre-computation search states this
exact two-homothetic-frustum kernel.  External novelty remains **OPEN**.

## Acceptance and kill boundary

- **Candidate hidden face mode:** a nonzero relative direction survives all
  length and shared-vertex equations and is uniquely the full-Poincare
  pointwise face stabilizer.
- **Local connection interpretation killed:** the full-Poincare control has
  dimension one but its intersection with the two strut-preserving kernels
  is zero, so the glued space is only the six-dimensional diagonal.
- **Underdetermination worse:** more than one relative mode survives.
- **Control failure:** the full-Poincare stabilizer, local kernels or direct
  shared-displacement equations are not reconstructed exactly.

A diagonal-only result would not kill first-order gravity.  It would close
only the claim that the already-found six flexes themselves supply an
independent face connection without adding new variables.
