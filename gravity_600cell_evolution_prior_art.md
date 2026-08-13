# Prior-art gate: 600-cell Regge evolution

Date: 2026-08-13

Status: **primary-literature map completed before the next control protocol**.

Scope: a targeted, non-exhaustive search for the same carrier, tent/staircase
evolution, Regge action, boundary canonical data and closed-universe source
terms.  This note does not prove bibliographic novelty.

## 1. Exact object being compared

The repository currently uses:

- the regular 600-cell boundary as a triangulated spatial `S3`;
- a phased Sorkin/tent advancement of all 120 vertices;
- a staircase triangulation of the product slab with 2400 four-simplices;
- ordinary Lorentzian length Regge curvature action;
- 35 independent internal squared-length orbits under the order-24
  pointwise schedule stabilizer;
- 30 independent final-boundary edge orbits and their post-momenta.

The literature question is not whether Regge calculus or the 600-cell is old.
It is which items in this exact list are old and which precise restriction is
not yet located.

## 2. Primary-source timeline

### 1973: regular closed Friedmann carrier

Collins and Williams used regular four-polytopes, including the 600-cell as a
spatial Cauchy surface, in a Regge model of Friedmann dynamics.  Their blocks
were symmetry-fixed prisms rather than the fully simplicial staircase used
here.  This is the origin of the closed-polytopal cosmology, not a result of
this repository.

### 1994/1997: simplicial tent evolution and staircase slab

Barrett, Galassi, Miller, Sorkin, Tuckey and Williams gave the parallelizable
implicit tent-move scheme and explicitly described ordering the upper vertices
to subdivide every tetrahedral prism into four 4-simplices.  Thus both local
tent advancement and the staircase product construction are **KNOWN**:

[A Parallelizable Implicit Evolution Scheme for Regge Calculus](https://arxiv.org/abs/gr-qc/9411008)

The paper also applied the scheme to a dust-filled 600-cell Friedmann model.
It imposed homogeneous edge lengths and used only a minimal subset of the
equations, explicitly noting that correctness required the ignored full
equations to admit the same solution.

That paper asserted four nonadjacent classes of 30 vertices.  This assertion
is false for the 600-cell edge graph: the repository's exact result
`alpha=24` and the four-dimensional kissing bound both exclude a 30-set.

### 2000/2001: correction to five phases and generalized dust evolution

De Felice and Fabri corrected the schedule to five classes of 24, printed one
class and generated the other four, and advanced them in five stages:

[The Friedmann universe of dust by Regge Calculus: study of its ending point](https://arxiv.org/abs/gr-qc/0009093)

[Singularities of the closed RW metric in Regge Calculus: a generalized evolution of the 600-cell](https://arxiv.org/abs/gr-qc/0106077)

They also recognized three inequivalent old-edge types between each pair of
color classes.  For each newly advanced class they then:

1. fixed the pole length as lapse;
2. imposed three shift/equality conditions on newly created slant edges;
3. equated all new connections according only to the class of their other
   endpoint;
4. solved four retained Regge equations for four lengths;
5. repeated this for all five classes.

The resulting five `4 x 4` systems break the original full symmetry but still
use a stronger ansatz than the order-24 orbit decomposition.  In particular,
the three stabilizer orbits within a phase pair are independent in the
repository's 35/65-variable action but are partly equated in the published
evolution.

The paper publishes an initial time-symmetric sandwich control:

```text
M approximately 10.202
tau_0 = 0.0102
l_0^2 = 7.69379990138304
d^2   = 7.69369586138301
l_1^2 = 7.69379990138297
```

with `d^2 = l_0^2-tau_0^2` to about fourteen significant figures.  It obtains
the mass from

```text
epsilon_3 = 2*pi - 5*acos(1/3),
M = (90/pi)*epsilon_3*l_0.
```

These values are a legitimate published control for the next verifier.  They
are not targets selected by the repository.

### 2011 onward: canonical boundary data

The statement that a discrete action is Hamilton's principal function and
generates pre/post momenta is standard canonical simplicial gravity:

[Canonical simplicial gravity](https://arxiv.org/abs/1108.1974)

Therefore the boundary-Legendre interpretation is **KNOWN STRUCTURE**.  The
repository's contribution is at most its explicit realization and rank audit
on this fixed 600-cell slab.

Closed `Lambda`-FLRW Regge models, including global-versus-local variation and
moment-of-time-symmetry constraints, were studied by Liu and Williams:

[Regge calculus models of the closed vacuum Lambda-FLRW universe](https://arxiv.org/abs/1501.07614)

Recent Lorentzian shell models use 600-cell boundaries and compare their
Hamilton-Jacobi functions with the continuum:

[Lorentzian quantum cosmology goes simplicial](https://arxiv.org/abs/2109.00875)

Thus dust, a positive cosmological constant, boundary momenta and 600-cell
Hamilton-Jacobi data all have published precedents.

## 3. KNOWN / CONTROL / OPEN boundary

### KNOWN

- 600-cell spatial slices in closed Regge cosmology;
- Sorkin/tent evolution;
- staircase subdivision of the tetrahedral product prisms;
- five nonadjacent classes of 24 and five-stage evolution;
- the existence of three inequivalent edge types between color classes;
- dust and positive-cosmological-constant closed-universe models;
- action-generated pre/post canonical momenta.

### CONTROL

The next calculation must insert the published dust mass and time-symmetric
lengths into the repository's exact full action.  It must report separately:

- the published homogeneous/equality-restricted residuals;
- all 35 internal orbit residuals before those equalities;
- the five pole equations including the dust derivative;
- the 30 diagonal equations unaffected directly by dust;
- even/odd schedule parity and all causal/branch checks.

Agreement only after averaging equations is weaker than agreement orbit by
orbit and must be labelled accordingly.

### OPEN

- whether the published time-symmetric sandwich solves every one of the 35
  independent orbit equations of the complete staircase action;
- whether the order-24 35/65-variable full-action audit has appeared in other
  literature;
- whether relaxing the published shift/equality ansatz yields a stationary
  dust, positive-`Lambda`, or pure-vacuum tick;
- whether either ordered phase parity has physical significance;
- whether the theory selects a dust mass, cosmological coefficient, lapse or
  physical duration rather than accepting it as input.

## 4. Decision

Do not start another unconstrained pure-vacuum root search yet.  First
preregister and run the published time-symmetric dust sandwich as an external
control on the same complete action.  This distinguishes three possibilities:

1. all full equations pass: the implementation reproduces a known physical
   tick and the generalized search has a trustworthy baseline;
2. only the reduced/averaged equations pass: the old solution relies on its
   extra ansatz and is not a full-action solution in the present sense;
3. even the published restricted equations fail: the action normalization,
   Lorentzian branch or slab conventions are not aligned, and new root results
   must not be physically interpreted until that mismatch is resolved.
