# Protocol: conditional canonical P1 dust weights

Date: 2026-08-19

Prior-art commit: `ba7de6c`

## Frozen inputs

```text
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f

reproducible/verify_gravity_600cell_projected_rank_edgewise_carrier.py
  50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23

reproducible/gravity_600cell_projected_rank_edgewise_carrier.json
  b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84

reproducible/gravity_600cell_projected_rank_edgewise_acceleration_blind.json
  2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2
```

No continuum acceleration, particle mass or spectral target is loaded.

## Frozen carriers and weights

Independently reconstruct

```text
K_0=P(sd K_600),
K_1=P(Esd_2(sd K_600))
```

using the direct rank-selected split.  For every tetrahedron `t` with chordal
Euclidean volume `V_t`, add `V_t/4` to each of its four vertices.  Assemble

```text
w_v=sum_(t incident on v) V_t/4.
```

Load only the already blind-selected total mass `M_j` at each carrier level
and set

```text
m_v=M_j*w_v/sum_u w_u.
```

Record the complete weight vector in the artifact by SHA-256 digest and record
its minimum, maximum, sum and symmetry-orbit census.  The artifact need not
store tens of thousands of redundant entries.

## Exact conditional uniqueness control

On the standard tetrahedron

```text
x>=0, y>=0, z>=0, x+y+z<=1,
```

integrate the four barycentric basis functions

```text
1-x-y-z, x, y, z
```

with exact SymPy rational arithmetic.  Require every integral to equal
`Vol/4=1/24`.

Then solve the four exactness equations for a general vertex-only linear
quadrature.  Require the unique solution `(1/24,1/24,1/24,1/24)`.  This proves
uniqueness only inside the stated `P1` ansatz.

## Frozen numerical gates

For each carrier require:

1. the exact frozen f-vector;
2. total chordal volume agreement with the carrier artifact within `5e-11`;
3. every tetrahedral volume `>0`;
4. every local and assembled vertex weight `>0`;
5. `abs(sum_v w_v - sum_t V_t)<5e-12`;
6. `abs(sum_v m_v-M_j)<5e-13`;
7. all 120 left and 120 right binary-icosahedral actions plus quaternion
   conjugation map the weights with maximum residual `<2e-10`;
8. under coordinate scales `s=(0.5,2.0)`, the maximum relative deviation from
   `w_v(s)=s^3 w_v(1)` is `<2e-12`;
9. for `tau=(0.01,0.1,1.0)`, the local dust action agrees with
   `-8*pi*M_j*tau` to relative error `<2e-14`.

The symmetry-orbit count and weight multiplicity distribution are printed
target-free.  No expected count is frozen because it has not yet been
computed.

## Output and hierarchy

Write

```text
reproducible/gravity_600cell_projected_rank_edgewise_local_dust.json.
```

Exactly one outcome is allowed:

1. `P1_LOCAL_DUST_INTERNAL_FAILURE` if any exact, topology, positivity,
   conservation, symmetry, scaling or action-collapse gate fails;
2. `P1_LOCAL_DUST_WEIGHTS_DERIVED_CONDITIONALLY` otherwise.

## Acceptance boundary

A positive outcome establishes only:

> Given a comoving continuous nodal `P1` density and vertex-only exact affine
> quadrature, the finite canonical carriers have coefficient-free, positive,
> conserved and `H4`-equivariant local dust masses.

It does not establish that Nature selects `P1` dust, does not construct
independent vertex lapses and does not authorize a physical Hessian by itself.

Before accepting the finite computation, a separate audit must assemble the
consistent `P1` tetrahedral mass matrix and recover the same weights by row
sums.  It must include the equal-global-vertex-mass rule as a negative control:
that rule conserves total mass and symmetry but must fail exact affine
element quadrature on the irregular carrier.
