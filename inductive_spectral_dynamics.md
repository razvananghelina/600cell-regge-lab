# Inductive spectral dynamics on barycentric refinement paths

Date: 2026-08-10

## Decision

A nontrivial inductive spectral dynamics can be constructed canonically on
the **top-cell refinement path space** of the 600-cell.  It is not a fitted
level-weight toy:

- the 600 root states are the actual tetrahedra of the 600-cell;
- the root coupling is their connected 4-regular dual graph;
- every top simplex has exactly 24 barycentric children;
- normalized conditional expectations fix the refinement embeddings;
- the hierarchical Laplacian is a Markov generator and intertwines those
  embeddings exactly;
- inverse three-volume scaling fixes `b=24^(1/3)`, and the square-root Dirac
  has zeta abscissa three.

It nevertheless fails a stronger physical gate suggested by asking one
operator to explain inertia, mass and the limiting propagation speed.  Among
the 24 barycentric children of one tetrahedron, only 36 unordered pairs share
a face.  The level-one hierarchical generator `I-P` directly couples all 276
pairs, including the 240 non-neighbours, and gives every nonconstant child
mode the same eigenvalue.  It therefore forgets spatial wavelength and does
not supply an incidence-local propagation cone.

The honest verdict is consequently narrower: this is a **STRUCTURAL ADVANCE
AS AN ULTRAMETRIC DIFFUSION**, but a **DERIVED NEGATIVE AS LOCAL SPACETIME
DYNAMICS**.  Its uniform state is tracial and its time is Euclidean diffusion
time, not Lorentzian time.  The construction must be moved from top-cell
paths to the full glued barycentric cochain complex before the modular/KMS
question can carry physical weight.

The natural noncommutative extension is the Cuntz algebra `O_24`.  Its
Perron--Frobenius KMS state supplies a nontrivial modular flow, but the compact
Dirac has unbounded ordinary commutators with the Cuntz generators.  The
modularly twisted commutators vanish algebraically.  Combining the KMS GNS
representation, compact resolvent and the real/Lorentzian gates in one
spectral object remains **OPEN**.

Exact verifiers:

- `reproducible/verify_inductive_spectral_dynamics.py`;
- `reproducible/verify_inductive_relativistic_gate.py`;
- `reproducible/verify_local_refinement_dynamics_gate.py`.

## 1. Complete hypotheses

Every positive statement below assumes all of the following.

1. The spatial carrier is the boundary triangulation of the regular
   600-cell, with 600 tetrahedra.
2. Refinement is complete barycentric subdivision.
3. Only nested **top-dimensional cells** are retained in the path coordinate;
   shared faces and lower-dimensional incidence are not retained below the
   root level.
4. A top cell has the uniform probability measure on its 24 barycentric
   children.  This is invariant under every permutation of the four parent
   vertices and chooses no flag.
5. Fine-to-coarse Hilbert maps are conditional averaging; coarse-to-fine maps
   repeat a parent amplitude on all children and are isometric for the
   normalized measure.
6. Dirac scale is inverse isotropic linear scale inferred from the exact
   affine volume ratio `1/24`.  This uses the already-derived spatial
   dimension three.
7. The root tetrahedron radius/energy unit is fixed separately.  The
   construction selects relative scales, not a conversion to seconds or GeV.

Hypothesis 3 is the main loss of geometry.  Hypothesis 6 is canonical as a
volume law but is not an independent derivation of dimension.

## 2. The finite inductive system

Let `X_0` be the 600 tetrahedra.  A level-`n` state is

`X_n = X_0 times {1,...,24}^n`,

so `|X_n|=600*24^n`.  The branch labels are the 24 flags of a tetrahedron;
the construction uses only their uniform sum, hence is invariant under their
permutation.

Define

`A_n=C(X_n)`, `H_n=L2(X_n,mu_n)`,

with uniform probability in every child fibre.  The algebra embedding repeats
a function on all children.  The Hilbert embedding `I_n:H_n->H_(n+1)` does the
same and is exactly isometric.

Let `E_k^(n)` be conditional expectation at level `n` onto functions which
depend only on the root and their first `k` refinement choices.  Then

`E_j E_k=E_min(j,k)`

and

`W_k=E_(k+1)-E_k`

are mutually orthogonal wavelet projections of rank

`rank(W_k)=600*23*24^k`.

The exact census is

`600 + sum_(k=0)^(n-1) 600*23*24^k = 600*24^n`.

No refinement modes are omitted.

## 3. The hierarchical Markov Laplacian

Reconstructing the 600 tetrahedra from the vertex graph gives 1200 triangular
faces, each incident on exactly two cells.  Their dual graph is connected,
4-regular and has 1200 edges.  Let `L_0` be its graph Laplacian.

Put

`b=24^(1/3)`, `lambda_k=b^(2(k+1))`.

At finite level define

`L_n = lifted(L_0) + sum_(k=0)^(n-1) lambda_k W_k`,

where `lifted(L_0)` acts on the root conditional averages.  Equivalently, the
vertical part is a positive sum of generators `I-E_k`, because the sequence
`lambda_k` is increasing.  Therefore:

- `L_n` is symmetric and positive;
- its off-diagonal entries are nonpositive and its row sums vanish;
- its kernel is one-dimensional;
- `exp(-t L_n)` is a positivity-preserving, constant-preserving Markov
  semigroup;
- `L_(n+1) I_n=I_n L_n` exactly.

Functional calculus then gives the positive Dirac

`D_n=sqrt(L_n)`

and

`D_(n+1) I_n=I_n D_n`.

Doubling `H_n` and placing `D_n` off diagonal supplies an even grading if
needed.  No real structure or Lorentzian reconstruction is claimed here.

This is nontrivial dynamics: the heat operator mixes cylinder amplitudes,
first inside refinement fibres and at long scale through the root dual graph.
It is not merely a phase attached to every level.

It is not local in the actual child incidence geometry.  A barycentric child
tetrahedron is a complete flag of the four parent vertices.  Two children
share an internal face exactly when their flags differ by an adjacent
transposition.  This gives a 3-regular dual graph with 36 edges.  By contrast,
`I-P` has the off-diagonal entry `-1/24` for every distinct pair.  Its support
is the complete graph on 24 children.

For the associated wave equation the exact propagator on this fibre is

`cos(t sqrt(I-P)) = P + cos(t)(I-P)`.

Every non-neighbour therefore receives the direct amplitude
`(1-cos(t))/24 = t^2/48 + O(t^4)`.  This does not prove that no effective
metric whatsoever could be put on the path space; it proves the relevant
negative: the operator is not selected by, and is not local for, the actual
face-incidence metric of the refined tetrahedron.

## 4. The combined inertia--mass--causality gate

The physical suggestion can be stated without metaphor.  Assume a positive
energy branch `E(p)`, a low-energy momentum variable supplied by spatial
geometry, and ordinary local Lorentz kinematics.  One and the same spectral
relation must yield:

1. rest mass: `E(0)=m c^2`;
2. inertial response: `E''(0)=1/m`;
3. limiting signal/group speed: `sup_p |E'(p)|=c`.

The calibration is the mass shell

`E(p)^2=m^2 c^4+c^2 p^2`.

It satisfies all three exactly.  More importantly, the repository already
contains the algebraic mechanism which could produce it rather than impose
it.  If `D_sp` is the spatial Kahler--Dirac operator and `gamma` is form
parity, then `{D_sp,gamma}=0`.  For an internal self-adjoint mass operator
`M` commuting with the spatial factor,

`H = c D_sp tensor 1 + gamma tensor M`

obeys

`H^2 = c^2 D_sp^2 tensor 1 + 1 tensor M^2`.

Thus a spatial singular value `p` and an internal eigenvalue `mu` give
`E^2=c^2 p^2+mu^2`, with `m=|mu|/c^2`.  The same operator then encodes rest
mass, inertial curvature and limiting speed.  **STRUCTURAL:** this identifies
a coherent mathematical bridge.  **OPEN:** the geometry has not selected
`M`, the physical normalization of `c`, or a Lorentzian time evolution.

The hierarchical path generator fails before those open parameters are
reached: its 23-dimensional vertical band is flat and has no spatial
momentum/wavelength label.  Reading mass or a light cone from it would be an
invention.

### 4.1 The first local replacement and its exact obstruction

There is a canonical local alternative on top cells: use the dual graph in
which two fine tetrahedra couple only when they share a face.  In one parent
tetrahedron every barycentric flag-child has three internal neighbours and
one neighbour across a parent face.  Exactly six of the 24 children meet each
of the four parent faces.

Let `I` repeat a coarse amplitude isometrically with factor `1/sqrt(24)`, and
let `L_f`, `L_c` be the unweighted fine and coarse dual-graph Laplacians.  The
flag census gives the exact compression identity

`I* L_f I = (1/4) L_c`.

Consequently the Galerkin requirement `I* (a L_f) I=L_c` selects `a=4`
uniquely.  This is a genuine geometric scale selection, not a fitted level
coefficient.  **DERIVED ADVANCE.**

The stronger operator identity nevertheless fails.  If the parent value is
zero and the four neighbouring parent values are `(1,0,0,0)`, then

`||4 L_f I-I L_c||^2=3`.

The residual has zero fibre average and lies entirely in the new directional
child modes.  In general its four channels are `sum(q)-4q_j`; it vanishes
only when all four neighbouring values agree.  No scalar normalization can
remove it.  **DERIVED NEGATIVE:** the incidence-local top-cell Laplacians do
not form an exactly intertwining inductive operator family.

This refines the continuation.  Exact operator intertwining is too strong
for the local discretization; exact **quadratic-form compression** plus
controlled spectral convergence is the viable criterion.  The repository's
independent Whitney calculation already supplies an exact isometric nested
one-form space, `P* M_f P=M_c`, while also showing that the current counting
inner product is not cylindrical.  The surviving candidate is therefore the full
cochain/Whitney (FEEC) Kahler--Dirac family with its metric-derived mass
matrices, not the top-cell graph and not the ultrametric path generator.

## 5. Dimension and compactness

Every barycentric top simplex has affine volume `1/24` of its parent.  Under
inverse isotropic three-dimensional linear scaling,

`b=(24)^(1/3)`.

The non-base Dirac spectrum has eigenvalue `b^(k+1)` with multiplicity
`600*23*24^k`.  Hence

`zeta_vertical(s) = 600*23*b^(-s)/(1-24*b^(-s))`

and its abscissa is

`log(24)/log(b)=3`.

The finite 600-dimensional root sector does not change the abscissa.
The resolvent is compact and the heat trace is finite for every positive
time.  **DERIVED conditional on Hypothesis 6.**

This does not recover four dimensions.  It confirms that the intrinsic
top-cell refinement dynamics is spatial and three-dimensional.

## 6. Perron--Frobenius/KMS completion

The infinite refinement paths form the full 24-shift.  Its noncommutative
path completion is generated by Cuntz isometries

`S_i* S_j=delta_ij`, `sum_i S_i S_i*=1`.

Uniform cylinder weights are

`mu([w])=24^(-|w|)`.

They are exactly projectively consistent.  For gauge dynamics

`alpha_t(S_i)=exp(i epsilon t) S_i`,

the generator KMS equation is

`phi(S_i S_i*)=exp(-beta epsilon) phi(S_i* S_i)`,

so

`beta epsilon=log(24)`.

This is the standard Perron--Frobenius mechanism for graph/Cuntz--Krieger KMS
states; uniqueness under the usual irreducibility hypotheses is established
for generalized gauge actions by
[Exel](https://arxiv.org/abs/math/0110183).  Spectral triples and diffusion on
weighted path trees are established by
[Pearson--Bellissard](https://arxiv.org/abs/0802.1336), and the stationary
Bratteli/higher-rank-graph version relates the zeta measure to the canonical
path measure in
[Farsi et al.](https://arxiv.org/abs/1803.09304).

Using the geometric refinement energy

`epsilon=log(b)=log(24)/3`

gives critical inverse temperature `beta=3`.  This is an exact compatibility
between branching, spatial scale and KMS balance.  It is not a temperature in
kelvin and does not determine a physical time unit.  Without the geometric
choice of `epsilon`, KMS fixes only the product `beta epsilon`.

## 7. Ordinary versus modular spectral triple

On a word of length `n`, the exponential Dirac has size `b^n`.  A creation
operator increases the length by one, so

`[D,S_i]|_(level n) = (b-1)b^n S_i`,

which is unbounded.  Thus the Cuntz generators are not in the ordinary
Lipschitz algebra of this compact-resolvent Dirac.  **DERIVED NEGATIVE.**

There is nevertheless an exact modular relation.  On the algebraic analytic
core define

`sigma(S_i)=b S_i`, `sigma(S_i*)=b^(-1) S_i*`.

This is the imaginary modular/gauge scaling.  Then

`D S_i-sigma(S_i)D=0`,

`D S_i*-sigma(S_i*)D=0`

away from the standard rooted boundary term.  The twist preserves the
algebraic Cuntz relations, although it is not star-preserving for real
`b != 1`.  **DERIVED CONDITIONAL:** the algebraic twisted-commutator gate
passes.

What has not been constructed is one Hilbert representation which
simultaneously carries:

- the critical `O_24` KMS state and its modular operator;
- the compact-resolvent dimension-three Dirac above;
- the required twisted boundedness on a dense star algebra;
- a real structure, reflection positivity or Lorentzian causal data.

The Fock representation naturally carrying the depth Dirac is Toeplitz--
Cuntz and its critical KMS state is not a normal finite Gibbs state there.
The critical KMS GNS representation is a different, type-III/modular arena.
Conflating the two would fake the missing link.

## 8. Acceptance and kill boundaries

### What is accepted

- **DERIVED:** an exact, connected, inductively compatible Markov dynamics on
  all top-cell refinement paths;
- **STRUCTURAL:** its volume-selected ultrametric Dirac and dimension three;
- **DERIVED:** projectively consistent Perron--Frobenius cylinder state;
- **DERIVED:** nontrivial KMS modular frequency on the `O_24` completion;
- **DERIVED CONDITIONAL:** exact modularly twisted commutators.

### What is not accepted

- a dynamics on the full glued barycentric `S^3`;
- an ordinary spectral triple containing the Cuntz generators;
- a single KMS-GNS compact-resolvent modular spectral triple;
- Lorentzian time, a fourth dimension or physical units;
- an incidence-local spatial propagation cone;
- a spatial dispersion relation capable of defining inertia;
- any finite internal algebra or matter sector.

### Next falsifiable gate

Construct metric-weighted local Kahler--Dirac operators `D_n=d_n+d_n*` on the
**full glued barycentric complexes**, together with refinement maps, and test
simultaneously:

1. exact or controlled inductive compatibility;
2. incidence-local support at every level;
3. an unbounded refinement momentum spectrum;
4. a wave/Dirac evolution with one scale-independent limiting speed;
5. the graded product identity with a geometrically selected internal `M`.

Only after this gate passes is it meaningful to construct its KMS/modular
completion.  Exact intertwining has already failed for the top-cell local
Laplacian, so the remaining acceptance standard is metric-derived quadratic
forms plus a proved, controlled spectral limit.  If that cannot be obtained
without fitted edge weights, the proposed refinement spacetime route is
closed.  Passing it would still leave the finite/matter selection contract
from `missing_link_audit.md` separate.

## Status ledger

- **DERIVED:** root dual graph `(600 vertices, 1200 edges, degree 4)`.
- **DERIVED:** exact inductive Hilbert/algebra/refinement system.
- **DERIVED:** connected hierarchical Markov Laplacians and heat semigroups.
- **DERIVED NEGATIVE:** at one refinement step the vertical generator has
  240 direct couplings between non-face-adjacent child tetrahedra.
- **DERIVED NEGATIVE:** its single flat 23-mode band carries no spatial
  wavelength/momentum dispersion.
- **DERIVED ADVANCE:** face-local fine energy compresses exactly to one
  quarter of the coarse dual-graph energy, uniquely selecting scale `4`.
- **DERIVED NEGATIVE:** at that unique scale the local operators do not
  intertwine; directional refinement modes leak from the inherited sector.
- **STRUCTURAL:** the graded product `c D_sp + gamma M` is the precise bridge
  capable of unifying mass, inertia and limiting speed.
- **STRUCTURAL:** `b=24^(1/3)` from inverse affine three-volume scaling.
- **DERIVED CONDITIONAL:** compact resolvent and zeta dimension three.
- **DERIVED:** uniform `24`-adic cylinder measure and KMS balance.
- **DERIVED NEGATIVE:** KMS balance alone fixes only `beta*epsilon`.
- **DERIVED NEGATIVE:** ordinary Cuntz commutators are unbounded.
- **DERIVED CONDITIONAL:** analytic modular twist kills the algebraic
  commutators.
- **OPEN:** a refinement-compatible incidence-local Kahler--Dirac/wave
  system, followed by a single modular/KMS spectral triple and Lorentzian
  interpretation.
- **NOT CLAIMED:** four-dimensional spacetime or particle physics.
