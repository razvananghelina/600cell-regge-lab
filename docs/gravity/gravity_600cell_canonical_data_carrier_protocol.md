# Protocol: exact classification of the 240 canonical-data directions

Date: 2026-08-19

This target-disclosed protocol is committed after the target-blind modular
nullity 240 and before constructing any non-homogeneous rational kernel
vector.

## Frozen provenance

| input | SHA-256 |
|---|---|
| classification prior-art gate | `fd7158f80af48fadc88c121c6001e258d2bccab480d8624b709f0ee142d145af` |
| admissibility protocol | `8db29cb9af699da660b969988eeb76c5e605e67c5ec65716795ada2e34674185` |
| admissibility verifier | `4d3595fbf418fc0876dba5a1129bdbcbd49d43a68ef9e6fd5fba2f0cb6e6873e` |
| target-blind admissibility artifact | `fa45c80739ca0dda4f82c9da98a4b22f4d8a18c182a40696a2a22d1d26ec89a1` |

The target-blind verifier must reproduce `10/10`, fixed ranks 3600,
augmented ranks 4200 and modular nullities 240 at both representatives and
both primes before its exposed matrices may be used.

## Candidate carrier derived locally

At each of

```text
(lambda,tau)=(2,5),(3,11),
```

derive the same `10 x 16` local natural-length Jacobian `J`, its pivot right
inverse `P` and its rational six-dimensional kernel `K` as in the frozen
audit.

For arbitrary values `(sigma_i,nu_i)` on the four local vertices form the
`16 x 8` displacement map

```text
U_i(sigma_i,nu_i)=sigma_i p_i+nu_i n.
```

Compute

```text
Y_local=J U,
Z_local=K^left (U-P Y_local),
```

where `K^left` is obtained from an exact nonsingular six-row minor.  Require
entrywise

```text
U=P Y_local+K Z_local.
```

The familiar coefficients must be outputs of `Y_local`, not inputs:

```text
upper edge (i,j): 8 lambda on sigma_i and sigma_j, zero on nu,
strut i:          6(lambda-1) on sigma_i, -2 tau on nu_i.
```

## Global 240-column map

Assemble a rational map

```text
R : Q^(120 sigma + 120 nu) -> Q^(3600 cell flex + 840 data).
```

For each tetrahedron, place `Z_local` in its six cell-flex rows using its four
global vertex labels.  For every global edge and vertex, place the data
coefficients derived from `Y_local`.

Do not form a fitted basis of the modular kernel.  `R` is defined solely by
the vertex displacement formula and the already frozen local decomposition.

For every sparse rational row `c` of the complete augmented matrix require

```text
c R=0
```

entrywise.  Report the number of tested rows and maximum support; no floating
tolerance participates.

## Rank and exhaustion certificate

Let `B_+` be the unsigned `720 x 120` vertex-edge incidence matrix,

```text
(B_+ sigma)_(uv)=sigma_u+sigma_v.
```

Require rank 120 modulo both frozen odd primes.  Since it has only 120
columns, either nonzero minor certifies rational rank 120.  The data part of
`R` is block-equivalent to

```text
[8 lambda B_+      0       ]
[6(lambda-1) I   -2 tau I  ],
```

so require its exact rational rank to be 240.

The frozen augmented modular rank 4200 proves

```text
dim_Q ker(C_aug)<=4440-4200=240.
```

Exact inclusion of a rank-240 rational map proves the reverse inequality.
Only their conjunction may assign

```text
ker_Q(C_aug)=im_Q(R), dim=240.
```

## Coordinate and falsification controls

1. Repeat the complete inclusion with the frozen alternate right inverse
   `P_alt=P+K T`; `Z_local` may change, but the data map and exhaustion may
   not.
2. Replace the endpoint sum on every upper edge by the oriented difference
   `sigma_u-sigma_v`, while retaining the legitimate local flex map.  At
   least one exact face row must become nonzero.
3. Delete the `nu_v` contribution from the lexicographically first strut
   data row while retaining its local displacement.  At least one exact face
   row must become nonzero.
4. The two constant columns must reproduce the already frozen homogeneous
   scale and lapse controls exactly.

The two wrong carriers are falsification controls only.  Their ranks are not
physical alternatives.

## Outcome hierarchy

1. `CANONICAL_DATA_CARRIER_CONTROL_FAILED` if provenance, reproduced frozen
   ranks, local derivation, decomposition, incidence-rank, alternate-graph,
   negative-control or constant-column checks fail.
2. `CANONICAL_DATA_EXACTLY_VERTEX_SCALE_PLUS_LAPSE` if all controls pass,
   `C_aug R=0` exactly, `rank_Q(data(R))=240`, and the frozen modular upper
   bound is 240 at both representatives.
3. `CANONICAL_DATA_VERTEX_CARRIER_PROPER_SUBSPACE` only if the rational
   rank-240 inclusion passes but an explicit additional rational kernel
   vector outside `im R` is supplied.  A modular count cannot trigger this
   branch.
4. `CANONICAL_DATA_VERTEX_CARRIER_REFUTED` if correctly reconstructed rows
   give a nonzero residual for `R` while all construction controls pass.
5. `CANONICAL_DATA_CARRIER_OPEN` otherwise.

## Interpretation firewall

A positive exact classification gives two vertex scalars, not a tensor
field.  It establishes neither that either scalar is dynamical nor that
`nu` is a physical clock.  The action must still be varied on this carrier.

Because the tangent data contain no independent upper-edge shape sector, a
positive result closes the linear gravitational-wave and effective-speed
route for schedule-free **flat tetrahedral frusta with only natural length
data**.  It does not close simplicial Regge calculus, first-order connection
extensions, non-flat cells or a separately selected temporal schedule.

Only the mission-specific verifier and static registry guards may be run.
The full suite remains excluded.
