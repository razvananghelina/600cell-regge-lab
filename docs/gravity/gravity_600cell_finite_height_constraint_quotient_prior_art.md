# Prior-art gate: exact constraint quotient of the finite-height boundary map

Date: 2026-08-22

Status: written before implementing or executing the targeted theorem verifier.

## Exact object and complete hypotheses

For each accepted finite-height slab, write the frozen Regge-plus-conserved-dust
action as

```text
S(O,X,N),
```

where `O` contains all `720` old-boundary logarithmic signed-squared edge
variables, `X` contains all `840` internal variables, and `N` contains all
`720` new-boundary variables.  For fixed old canonical data `(O,p_minus)`, the
implicit update equations are

```text
S_X(O,X,N) = 0,
-S_O(O,X,N) = p_minus.
```

Their derivative with respect to `(X,N)` is the complete `1560 x 1560`
pre-Legendre matrix

```text
J = [[ S_XX,  S_XN],
     [-S_OX, -S_ON]].
```

The proposed quotient is required to be selected by exact pre/post constraints
or an exact continuous gauge symmetry of this same action.  It must be defined
without deleting nonzero singular directions by a numerical threshold.  The
question is local and first-order at the frozen branch-B history with incoming
representative `v=3/2`, for the first and second accepted slabs, both staircase
parities, and the unreduced `1440`-dimensional boundary phase carrier.

The following are not included in the hypotheses:

- a third or later slab;
- a refinement family or continuum limit;
- independently fluctuating dust variables;
- a different or perfect action;
- a newly derived continuous momentum-map symmetry;
- a continuum Hamiltonian or diffeomorphism constraint imposed by hand.

## Repository gate

The authoritative map and the wider repository were searched under
`constraint quotient`, `pre-constraint`, `post-constraint`, `pseudo-constraint`,
`Legendre rank`, `Lagrangian two-form`, `gauge kernel`, `presymplectic`,
`reduced phase space`, and `vertex displacement`.

### ACCEPTED INPUT

The complete dense first-slab calculation independently certifies, in both
parities,

```text
classification                 REGULAR
normalized sigma_min(J)        4.0880740e-4
frozen gap gate                3.6510838e-7
complete boundary tangent      CANONICAL
```

The complete dense second-slab replication independently certifies regular
normalized and physical pre-Legendre matrices in both parities.  Its smallest
normalized singular values are approximately `3.34e-4` to `4.09e-4`, against
frozen gates from `4.56e-9` to `9.13e-8`.  All four parity-ordered two-step
products are canonical under the accepted uncertainty model.

The first accepted artifact has SHA-256

```text
ee9491b2ae5fdf3f2a9d0d78c0e837c8c2692797d87ccd8e1757efeadd8060e7
```

and the second accepted artifact has SHA-256

```text
1355f8cf339d18c1cf2855ecb1228e97e868d73f7a1ef739e4c11ce9521fcd4b
```

The internal one-dimensional lapse-constraint tangent is not a boundary gauge
direction: the accepted exact reconciliation proves that fixing incoming
momentum removes it.

### REUSABLE CONTROL

On the older stationary order-24 carrier, one collective lapse null and four
weak relative-pole directions were found.  The full `2280`-variable audit then
resolved the apparently weak `120`-dimensional carrier as nonzero and certified

```text
rank(J_full) = 1560/1560,
error-consistent nullity = 0.
```

That result is a control for the distinction between an exact kernel and weak
pseudo-constraint candidates.  It is not an input proving the finite-height
claim.

The Euclidean embedded-framework self-stress carrier is also already refuted
as a dynamically closed York proxy.  It cannot be reused as an action-derived
quotient.

## Primary prior art

Dittrich and Hoehn define regular variational discrete systems by nonsingularity
of the mixed Lagrangian two-form.  In that case the discrete Legendre
transforms are locally invertible, canonical evolution is unique and
symplectic, and constraints do not arise.  In the singular case, left and
right null vectors of the Lagrangian two-form make the Legendre transforms
non-surjective; their images are the pre- and post-constraint surfaces.  See
Sections 3.2 and 3.3 of
[*Constraint analysis for variational discrete systems*](https://arxiv.org/abs/1303.4294),
J. Math. Phys. 54, 093505 (2013), DOI `10.1063/1.4818895`.

Hoehn's quadratic classification is explicitly based on the null vectors of
the Lagrangian two-form:
[*Classification of constraints and degrees of freedom for quadratic discrete
actions*](https://arxiv.org/abs/1407.6641), J. Math. Phys. 55, 113506 (2014),
DOI `10.1063/1.4900926`.

Bahr and Dittrich distinguish exact constraints from pseudo-constraints in
curved Regge calculus.  Curvature generically breaks the exact discrete gauge
symmetry; the resulting relations depend weakly on next-step data and can fix
lapse and shift rather than constrain a single time slice.  They also identify
continuum/refinement limits and alternative actions as possible routes by
which exact symmetry may return:
[*\(Broken\) Gauge Symmetries and Constraints in Regge
Calculus*](https://arxiv.org/abs/0905.1670), Class. Quant. Grav. 26, 225011
(2009), DOI `10.1088/0264-9381/26/22/225011`.

The action-generated canonical framework for simplicial evolution is given in
Dittrich and Hoehn,
[*Canonical simplicial gravity*](https://arxiv.org/abs/1108.1974), Class.
Quant. Grav. 29, 115009 (2012), DOI `10.1088/0264-9381/29/11/115009`.

These papers establish the general mechanism.  They do not publish the present
600-cell matrices.  External novelty of the project-specific numerical
regularity result remains **OPEN**.

## Exact implication to be tested

Define

```text
F(O,p_minus;X,N) = (S_X, -S_O-p_minus).
```

If `J = partial_(X,N) F` is invertible, the implicit-function theorem gives a
unique local solution `(X,N)=Y(O,p_minus)`.  The outgoing momentum
`p_plus=S_N(O,Y,N)` then defines a local boundary map

```text
T : (O,p_minus) -> (N,p_plus).
```

If the derivative of `T` is symplectic, it is invertible and maps an open
neighbourhood of the complete old phase space to an open neighbourhood of the
complete new phase space.  Its domain and image therefore cannot be
positive-codimension pre/post-constraint surfaces.  Equivalently, there is no
nontrivial local action-derived pre/post-constraint quotient at that slab.

This implication does not say that every exact symmetry of every regular
Hamiltonian system is absent.  A separately derived continuous group action
and momentum level could permit a different symplectic reduction.  No such
action has been derived here; the finite staircase symmetry is not one.

## KNOWN / CONTROL / OPEN

- **KNOWN:** regular discrete Legendre evolution has no nontrivial
  pre/post-constraint surface generated by Legendre degeneracy.
- **KNOWN:** a symplectic square tangent is invertible.
- **CONTROL:** a deliberately singular mixed generating function must produce
  a positive-codimension constraint image.
- **CONTROL:** an exactly nonzero but arbitrarily small mixed coefficient must
  remain regular; thresholding it to zero must be rejected.
- **OPEN BEFORE THE TEST:** whether every accepted first/second slab and all
  parity products meet the frozen regularity and canonicality premises.
- **OPEN:** constraints propagated backward by a later singular move.
- **OPEN:** an exact continuous symmetry from an enlarged matter carrier.
- **OPEN:** restoration of vertex-displacement gauge symmetry under a declared
  refinement family or a different/perfect action.
- **OPEN:** a physical graviton quotient, wave equation, limiting speed, `G`
  or Planck scale.

## Framing attack

The current route name, `P-CONSTRAINT-QUOTIENT`, combines two objects that are
not interchangeable.  Exact constraints can select a presymplectic reduction.
Pseudo-constraints are nonzero dynamical consistency relations and do not
generate exact gauge orbits.  Quotienting their weak directions by a numerical
cutoff changes the finite theory.

Regularity for two slabs does not prove regularity for an infinite history.
In variational discrete systems, a later singular move can impose a secondary
condition which propagates backward.  The strongest licensed conclusion is
therefore a bounded local no-go for the accepted two-slab region, not a theorem
that Regge calculus or the 600-cell can never possess constraints.

If the premises survive the targeted verifier, the honest next gate is not a
mode spectrum on a fitted quotient.  It is a refinement/alternative-action
test asking whether a geometrically identified weak carrier approaches an
exact Lagrangian-two-form kernel with a preregistered scaling law.
