# Prior-art gate: identify the 120+120 compatible data space

Date: 2026-08-19

## Exact object and complete hypotheses

On the frozen nonstatic complete variable-face flat-frustum 600-cell system,
let `K_E` be the modular 120-dimensional compatible upper-edge-only data
space established in commit `9b97775`.  Let `U` be the 720-by-120 unsigned
vertex-edge incidence map of the 600-cell graph.  For an edge `{u,v}`, the
derived squared-length variation from vertex scale variables is

\[
\delta q_{uv}=8\lambda(\sigma_u+\sigma_v).
\]

The nonzero factor `8 lambda` does not change the image.  The target-disclosed
question is whether `im(U)=K_E` over each frozen finite field.  The strut
candidate is the entire 120-dimensional strut coordinate space, already
implied modularly by `rank([F S])=rank(F)`.

This mission identifies boundary-data subspaces only.  It does not reuse the
refuted local cell-flex lift and does not construct an action, Hessian,
symplectic form, time step, propagation speed, or physical unit.

## KNOWN

Vertex-factor edge scalings are standard in discrete conformal geometry.  In
piecewise-flat two- and three-manifolds, Glickenstein studies conformal
variations of edge lengths and their effect on Regge scalar curvature; see
[arXiv:0906.1560](https://arxiv.org/abs/0906.1560).  Therefore the formula
`sigma_u + sigma_v` is not itself new.

Flat linearized Regge calculus has vertex-displacement gauge freedom.  In the
canonical Pachner-move analysis, a 1-4 move introduces four lapse-and-shift
variables; see Hoehn
[arXiv:1411.5672](https://arxiv.org/abs/1411.5672), especially the 1-4 move.
Thus an unrestricted modular strut sector is consistent with known
flat-background kinematics and is not evidence for a physical clock.

## CONTROL

The 600-cell graph is connected and contains triangular faces.  If
`U sigma=0`, then `sigma_u=-sigma_v` on every edge.  Propagation around a
triangle gives `sigma_u=-sigma_u`; over the odd frozen primes this forces
`sigma_u=0`, and connectedness forces every entry to vanish.  Hence `U` has
rank 120 over the rationals and both frozen fields.

For matrix blocks `[F E S]`, image inclusion of the candidate is tested
without selecting a cell-flex lift:

\[
E\,\operatorname{im}U\subseteq\operatorname{im}F
\quad\Longleftrightarrow\quad
\operatorname{rank}[F\;EU]=\operatorname{rank}F=3600.
\]

If inclusion holds, equal dimensions `dim im(U)=dim K_E=120` then prove
equality.  Neither equality of dimensions nor the existence of a differently
chosen local lift is sufficient by itself.

As a negative construction control, change exactly one incidence row from
`sigma_u+sigma_v` to `sigma_u` while keeping every other row fixed.  The
corrupted map must have rank 120 and its joint image with `U` must have rank
121.  If both distinct 120-dimensional images were compatible, the frozen
`K_E` dimension would be at least 121, contradicting the census.

## OPEN and proposed difference

The literature search did not identify the equality of this precise
600-cell flat-frustum compatibility space with the unsigned incidence image.
Search absence is not proof; external novelty remains **OPEN**.

The exact questions left open are:

1. modular equality `im(U)=K_E` under every frozen construction and
   convention;
2. equality over the rationals;
3. construction of the unique global cell-flex lift through full-column-rank
   `F`;
4. whether the action turns the arbitrary kinematic strut directions into
   gauge, constraints, or selected lapse data.

The first question is the next registered test.  The other three are not
licensed by a modular positive.
