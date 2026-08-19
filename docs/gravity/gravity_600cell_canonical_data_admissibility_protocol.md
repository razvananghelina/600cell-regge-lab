# Protocol: canonical-data tangent admissibility on the regular 600-cell

Date: 2026-08-19

This protocol is committed before constructing any natural-data forcing
column or evaluating any augmented rank.  The two homothetic tangent
directions are disclosed controls; no other kernel dimension is assumed.

## Frozen provenance

| input | SHA-256 |
|---|---|
| prior-art/framing gate | `8acf8e29a809f065033d49dd03c1e858f1ed7a219d1be6acfb27851f78f1ce56` |
| complete variable-face protocol | `ed79c6a15ade377ae09854b3cad3028eb1c0f43cf8e85789d46993fe25ed1b49` |
| complete variable-face verifier | `ec44be8e4d82634e30944739d10d3f80fbb9f6fee0883ec1f612690c38d90ab6` |
| complete variable-face artifact | `61cebd1cd67fcdc56de088855b1fc7b805d0f70f9f9b3029d4a61209d7a53944` |
| local variable-face result | `2db55cb87ec1c01d537cdbc11010bc9ea740762c598108e4c2de0f3acca72cc8` |
| local variable-face verifier | `69a5d7479a5df427cead76f82db31fe62a9190c28c967f699c846881634fb0f6` |
| local variable-face artifact | `001212016553d006862e68edc4f780f37ca1476110b6e0aed3e987f52a43b5e3` |

## Exact carrier and representatives

Rebuild the golden-field 600-cell by clique incidence and require exactly

```text
f=(120,720,1200,600),
two tetrahedra per triangular face,
five tetrahedra per spatial edge,
connected four-regular dual graph.
```

Use

```text
eta=diag(1,1,1,-1),
n=(0,0,0,1),
p0=( 1, 1, 1,0), p1=( 1,-1,-1,0),
p2=(-1, 1,-1,0), p3=(-1,-1, 1,0)
```

and the two predeclared nonstatic representatives

```text
(lambda,tau)=(2,5),(3,11).
```

The lower slice is fixed.  Global data columns are, in lexicographic order,

```text
720 upper squared-edge variations,
120 corresponding-strut squared-length variations.
```

## Local natural-length Jacobian

Stack the sixteen upper-vertex displacement coordinates as
`(delta q_0,...,delta q_3)`.  Derive, rather than insert, the `10 x 16`
Jacobian `J` whose rows are

```text
delta l_ij^2 = 2 <q_i-q_j,delta q_i-delta q_j>,  i<j,
delta s_i^2  = 2 <q_i-p_i,delta q_i>.
```

Require exactly:

1. `rank_Q(J)=10`;
2. a deterministic pivot right inverse `P` with `J P=I_10`;
3. a rational kernel basis `K` of shape `16 x 6`;
4. equality of `im K` with the already derived six-dimensional constrained
   relative-Poincare displacement image;
5. an alternate right inverse

   ```text
   P_alt=P+K T,
   T_ab=(a+1)(b+1),
   ```

   which must leave all final ranks unchanged.

The right inverse is only a coordinate graph.  Any verdict depending on its
choice is a control failure.

## Face equations with forced data

For every directed shared face, independently reconstruct the same full
Lorentzian affine transition used by the complete variable-face audit.  In a
common ordering of its three global vertices, form

```text
F_z = [E_source K, -L_f E_target K],
F_y =  E_source P D_source - L_f E_target P D_target.
```

Here `D_T` maps the six local upper edges and four local struts of tetrahedron
`T` to the 840 global data columns, and `L_f` is repeated on the three vertex
displacements.

Let `s_f` be the displacement column generated on the upper face by the
unique Poincare line fixing the lower triangle pointwise.  Eliminate its
coefficient exactly with a rational row basis `Q_f` of

```text
ker(s_f^T).
```

The exact forced face equations are

```text
Q_f (F_z z + F_y y)=0.
```

Select a row basis from the combined local matrix only by exact rational
row reduction.  Record the complete local-rank census; do not assume that
the data columns leave the old rank five unchanged.  The fixed-data part
must nevertheless reproduce rank five and the old six-to-seven local
theorem on every face.

## Complete augmented matrix

Assemble

```text
C_aug : Q^(3600+840) -> Q^m,
```

where the first `3600=6*600` columns are cell-flex coordinates and the last
840 are the global data.  Let `C_z` denote the first 3600 columns.

Compute exact finite-field ranks at

```text
p1=1000003,
p2=1000033.
```

For each representative require the old fixed-data rank

```text
rank(C_z mod p)=3600
```

at both primes.  Record, without target comparison,

```text
r_p=rank(C_aug mod p),
nu_p=4440-r_p.
```

Because `C_z` has full rational column rank, the projection of a rational
augmented kernel vector to the 840 data columns is injective.  A modular
nullity alone remains only an upper bound on the rational admissible
dimension unless matched by explicit rational witnesses.

## Disclosed rational positive controls

For each representative construct the two infinitesimal homothetic motions
directly from upper displacements:

```text
scale: delta q_i=p_i,
lapse: delta q_i=n.
```

Their global data are derived with the Jacobian, not hard-coded.  Decompose
each local displacement exactly as

```text
delta q=P y_local+K z_T
```

and require both complete rational vectors to satisfy `C_aug(z,y)=0`.
Their data projections must have rank two.  Repeat the construction with
`P_alt`; the complete vectors may change in `z`, but their data columns and
null residuals may not.

These two controls prove only

```text
dim_Q admissible data >=2.
```

They do not assume that the dimension equals two.

## Falsification and convention controls

1. Reversing all face orientations must preserve both fixed and augmented
   modular ranks at `(lambda,tau)=(2,5)`.
2. The odd canonical relabelling `(0 1)` must preserve those ranks.
3. Replacing `eta` by `-eta` must preserve those ranks.
4. The deterministic alternate right inverse `P_alt` must preserve all
   augmented ranks at both representatives.
5. On the lexicographically first face, remap one target shared upper-edge
   datum to a private column while leaving the source incidence unchanged.
   At least one of the two exact homothetic controls must then fail that
   corrupted face equation.  This is an incidence-wiring negative control,
   not a rank target.

## Outcome hierarchy

1. `CANONICAL_DATA_ADMISSIBILITY_CONTROL_FAILED` if provenance, incidence,
   local Jacobian, Poincare-kernel, transition, fixed-rank, homothetic,
   right-inverse or convention controls fail.
2. `CANONICAL_DATA_ONLY_HOMOTHETIC` if controls pass and both primes give
   `nu_p=2` at both representatives.  The two rational independent controls
   and modular upper bound then prove exact rational admissible dimension
   two.
3. `CANONICAL_DATA_INTERMEDIATE_MODULAR_OPEN` if controls pass, the two
   primes agree at each representative and `2<nu_p<840`, but no complete
   rational basis of that size is certified.
4. `CANONICAL_DATA_FULL_MODULAR_OPEN` if controls pass and `nu_p=840` at both
   primes, but 840 rational independent witnesses are absent.
5. `CANONICAL_DATA_PRIME_DISAGREEMENT_OPEN` if controls pass but the prime
   nullities disagree.
6. `CANONICAL_DATA_ADMISSIBILITY_OPEN` otherwise.

No outcome called “full rational admission” exists without 840 explicit
rational data directions.  Agreement of two finite fields cannot supply
that missing upper-rank certificate.

## Interpretation firewall

`CANONICAL_DATA_ONLY_HOMOTHETIC` would close the schedule-free flat-cell
anisotropic route at the tangent level: the cellular object would support
only scale and common-lapse changes, not an inhomogeneous canonical phase
space.

An intermediate modular result remains **OPEN** until its rational carrier
is built and classified.  Full admission, if later proved, merely authorizes
an action/Hessian calculation.  None of the outcomes by itself derives a
physical tick, gravitational waves, `c`, `G` or Planck units.

Only the mission verifier and static registry guards may be run.  The full
suite is excluded by the user's standing instruction.
