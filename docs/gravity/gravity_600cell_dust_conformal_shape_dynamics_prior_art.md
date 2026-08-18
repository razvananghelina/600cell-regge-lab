# Prior-art gate: dynamic closure of the conformal/shape supermetric split

Date: 2026-08-18

## Exact object, carrier and hypotheses

Use the fixed regular 600-cell boundary carrier with `120` vertices, `720`
spatial edges and logarithmic squared-edge variations

```text
E = C^720.
```

The literal vertex-conformal map is

```text
C : C^120 -> E,
(C sigma)_uv = sigma_u + sigma_v.
```

Its image `K = im C` has exact dimension `120`.  The already committed
centered two-slab Jacobi equation is

```text
M (q_2 - 2 q_1 + q_0)
+ N (q_2 - q_0)
+ V q_1 = 0,
```

with `M` invertible and normalized operators

```text
Gamma = M^-1 N,
Omega = M^-1 V.
```

For each of the two frozen schedule parities, seven minimal
binary-tetrahedral sectors and four frozen derivative variants, define the
Hermitian kinetic form

```text
H = (M + M*)/2
```

and the action-relative shape complement

```text
S_H = ker(C* H).
```

The preceding certified result gives a nondegenerate restriction of `H` to
`K`, with minority sign on all `120` conformal directions and the opposite
sign on the `600`-dimensional complement.  Hence, sector by sector,

```text
E = K direct-sum S_H
```

is a uniquely defined `H`-orthogonal direct sum.  The present question is:

> Do both normalized action-generated operators preserve both factors of
> this direct sum, in every frozen audit?

Equivalently, do `Gamma` and `Omega` define separate conformal and shape
recurrences, or do they mix them?

There is one geometry/action-selected split per frozen matrix audit.  The
schedule and derivative variants are robustness checks, not candidates from
which a favorable split may be chosen.  No rotation of either factor after
looking at `Gamma` or `Omega` is allowed.

## What is already known

### KNOWN

1. Vertex conformal variations of piecewise-flat metrics and their curvature
   response are established discrete differential geometry.  See
   Glickenstein,
   [*Discrete conformal variations and scalar curvature on piecewise flat two
   and three dimensional manifolds*](https://arxiv.org/abs/0906.1560), and
   Champion--Glickenstein--Young,
   [*Regge's Einstein-Hilbert Functional on the Double
   Tetrahedron*](https://arxiv.org/abs/1007.0048).
2. The continuum conformal--traceless decomposition is not, by itself, a
   physical mode decomposition.  Its evolution is tied to lapse, shift,
   constraints and gauge conditions.  Brown gives an action/evolution
   treatment in
   [*Conformal invariance and the conformal-traceless decomposition of the
   gravitational field*](https://arxiv.org/abs/gr-qc/0501092).
3. The Lund--Regge supermetric can have triangulation- and point-dependent
   signature and degeneracy.  A continuum `1:5` analogy does not force a
   discrete dynamical split.  See Hartle--Miller--Williams,
   <https://arxiv.org/abs/gr-qc/9609028>.
4. On curved Regge backgrounds exact vertex-displacement gauge symmetry is
   generically broken and constraints become background-dependent
   pseudo-constraints.  See Bahr--Dittrich,
   <https://arxiv.org/abs/0905.1670>.

### CONTROL

- `C` is literal, rank `120`, exactly equivariant under the frozen order-24
  action and paired with each schedule's exact edge ordering.
- In every minimal sector of irrep dimension `d`, `K` has dimension `5d`.
- In all `56` audits the restriction `C* H C` is positive and nondegenerate,
  while the `H`-orthogonal complement carries the opposite sign.
- All `M` blocks are invertible and the stored `Gamma,Omega` identities are
  certified in ball arithmetic.
- The complete fixed-carrier recurrence has already been committed before the
  present subspace comparison.

### OPEN

- Whether `K` is invariant under `Gamma` and `Omega`.
- Whether `S_H` is invariant under `Gamma` and `Omega`.
- Whether any observed mixing survives a continuum refinement.
- Any scalar/vector/tensor or constraint interpretation of the factors.
- External novelty of this exact 600-cell comparison.

## Distinction from the previous scale/shape test

The earlier `60 x 60` quotient calculation isolated only the single uniform
scale line among `30` boundary-edge orbits and found its two-dimensional phase
plane invariant.  That statement concerns

```text
1 global scale + 29 zero-sum quotient shapes.
```

It does not test the complete conformal carrier.  Under the order-24 quotient,
the `120` vertex conformal variables reduce to five independent vertex-orbit
variables, so four inhomogeneous conformal directions were included in the
old object called "shape".  On the full edge carrier the present split is

```text
120 conformal + 600 H-orthogonal shape.
```

The proposed calculation is therefore not a repetition of the old global
scale result.

## Framing attack

The word "shape" must not be read as "transverse-traceless graviton".  The
`600`-dimensional factor can contain longitudinal, scalar, vector,
constraint-violating and discretization modes.  No Hamiltonian or momentum
constraint has reduced it to two polarizations.

Likewise, complete conformal/shape mixing would not refute general relativity.
On an FLRW background, scalar metric and matter/constraint variables can mix;
only a fully reduced scalar/vector/tensor analysis supplies physical sectors.
A negative result would instead close the shortcut

```text
kinematic 1:5 signature  =>  dynamically independent 5-component sector.
```

A positive result would also be limited: it would select a finite invariant
subsystem, not prove that this subsystem is gauge-reduced or converges to the
continuum tensor equations.

Finally, `S_H` is not geometry-only: it uses the action-derived kinetic form.
That dependence is explicit and is the reason the split is relevant to the
recurrence.  It is canonical only relative to the declared pair `(C,H)` and
the literal adjacent-slice identification; it is not asserted to be a
universal York decomposition.

## Licensed target-free test

Before inspecting any new cross residual:

1. reconstruct `K` and `S_H` independently in every sector and audit;
2. certify their dimensions, direct-sum conditioning and complete numerical
   error envelopes;
3. test `Gamma K subset K`, `Omega K subset K`,
   `Gamma S_H subset S_H`, and `Omega S_H subset S_H`;
4. classify every residual as zero-consistent, resolved nonzero or open using
   frozen error bands;
5. require every residual to be zero-consistent for dynamic decoupling.

One resolved nonzero residual refutes the universal closure claim.  Open
residuals cannot be counted as hits.  No continuum spectrum, desired mode
count, wave speed or refinement target belongs in this gate.
