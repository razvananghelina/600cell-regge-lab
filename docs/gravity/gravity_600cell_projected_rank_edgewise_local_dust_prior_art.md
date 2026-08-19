# Prior-art gate: local conserved dust on the canonical projected carrier

Date: 2026-08-19

## Exact object and hypotheses

Let `K` be either certified spatial carrier

```text
P(sd K_600),
P(Esd_2(sd K_600)).
```

Represent a comoving dust density on one spatial slice by continuous
piecewise-linear (`P1`) nodal data.  On every chordal Euclidean tetrahedron
`t`, require a local quadrature with the following complete hypotheses:

1. it uses only the four nodal density values of `t`;
2. it is linear in those values;
3. it is exact for every affine density on `t`;
4. global nodal masses are assembled additively from incident tetrahedra;
5. masses are then held fixed along the corresponding comoving worldlines.

For a total dust mass `M`, define

```text
w_v = sum_(t incident on v) Vol(t)/4,
m_v = M*w_v/sum_u(w_u).
```

The candidate local matter action for independently assigned proper strut
lengths `tau_v` is

```text
S_dust = -8*pi*sum_v m_v*tau_v.
```

This mission asks only whether the weights are uniquely selected within the
stated `P1` ansatz, positive, exactly conservative and equivariant on the two
finite carriers.  It does not yet construct a Lorentzian slab with
independent `tau_v`.

## The conditional uniqueness argument

Let `lambda_i` be the four barycentric coordinate functions on one
tetrahedron.  Every affine density is

```text
rho(x)=sum_i rho_i lambda_i(x).
```

Symmetry of the reference simplex and `sum_i lambda_i=1` give

```text
integral_t lambda_i dV = Vol(t)/4.
```

Therefore exactness on the four basis functions `lambda_i` forces every
vertex-only linear quadrature coefficient to be `Vol(t)/4`.  Conversely those
weights integrate every affine density exactly.  There is no free coefficient
inside this hypothesis class.

This is a standard finite-element fact, not a new theorem.

## Primary prior art

Dittrich, Gielen and Schander formulate simplicial Lorentzian cosmology with
dust particles; conserved particle number leaves fixed masses and each
particle contributes minus its mass times proper time:

- *Lorentzian quantum cosmology goes simplicial*, arXiv:`2109.00875`, DOI
  `10.1088/1361-6382/ac42ad`, dust discussion in the shell model.

Equal division of tetrahedron volume among its four vertices is the standard
barycentric dual/lumped `P1` mass construction.  A modern comparison of
barycentric, circumcentric and optimized tetrahedral dual volumes is:

- A. Jacobson, *Optimized Dual-Volumes for Tetrahedral Meshes*,
  arXiv:`2406.08647`, DOI `10.1111/cgf.15133`, especially Section 3.1 and
  equation (15).

That comparison is also a warning against claiming unconditional uniqueness:
other dual-volume constructions exist, and their properties differ on
irregular or non-Delaunay tetrahedra.

Mass lumping and barycentric dual grids are long-standing finite-element and
finite-volume techniques.  Their use here is a discretization choice, not a
new physical principle.

## KNOWN / CONTROL / OPEN

### KNOWN

- A vertex-only linear quadrature exact on affine tetrahedral fields has the
  unique weights `Vol(t)/4`.
- Barycentric mass lumping is positive on every nondegenerate tetrahedron and
  exactly integrates constant and affine `P1` densities.
- A conserved dust particle contributes `-m` times its proper worldline
  length, up to the repository's fixed `8*pi` normalization.

### CONTROL

- The two spatial f-vectors and chordal volumes must reproduce the frozen
  carrier artifact.
- Every local tetrahedral contribution and every assembled vertex weight must
  be strictly positive.
- `sum_v w_v` must equal the total chordal volume to roundoff.
- Under all certified `H4` spatial automorphisms, assembled weights must map
  exactly within the frozen numerical tolerance.
- Scaling all coordinates by `s` must scale all weights by `s^3`.
- On the homogeneous strut `tau_v=tau`, the local action must collapse exactly
  to the previous global action `-8*pi*M*tau`.

### OPEN

- Whether `P1` nodal dust is selected by fundamental physics rather than used
  as a standard consistent discretization.
- Inhomogeneous dust-density degrees of freedom and their initial data.
- A canonical Lorentzian product triangulation for independent vertex lapses
  on the refined spatial carrier.
- The local Regge constraints and the constraint-reduced gravitational
  Hessian.

## Framing attack

The formula `Vol(t)/4` is not uniquely forced among all imaginable local
mass assignments.  It is forced only after choosing continuous nodal `P1`
density and exact affine quadrature.  Higher-order elements, discontinuous
density, circumcentric duals, optimized duals or freely specified point-particle
masses evade the theorem.

Accordingly:

- the quadrature theorem is **KNOWN / DERIVED CONDITIONAL**;
- using the `P1` dust ansatz in this theory is **STRUCTURAL**;
- positivity, conservation, symmetry orbits and collapse to the global dust
  action on these carriers can be **DERIVED COMPUTATIONAL**;
- no local lapse equation follows from the spatial weights alone.

The mission advances the gravity route only if it supplies a coefficient-free
matter discretization under a stated approximation class.  It does not by
itself authorize the inhomogeneous Hessian.
