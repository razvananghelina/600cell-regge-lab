# Preregistered protocol: asymmetric Lorentzian tent existence and selection gate

Date: 2026-08-12

Status at registration:
**PROTOCOL ONLY -- TARGET-FOUND ASYMMETRIC WITNESS DISCLOSED**

## 1. Why this gate is different

The symmetric Lorentzian tent theorem proves that equal final cone lengths
give an internal angle below `2*pi/5` for every pole length and final radius.
It does not prove the same bound when the 12 final cone lengths differ.

This protocol separates two questions:

1. **existence:** can any strictly admissible asymmetric final spatial star
   satisfy the zero-volume Regge pole equation?
2. **selection:** does the bare 600-cell geometry select that asymmetry and
   its pole, rather than merely permitting a fitted solution?

A positive answer to existence refutes an extension of the symmetric no-go.
It is not by itself a physical tick.

## 2. Complete Lorentzian hypotheses

Use the same canonical local carrier

```text
T_v=[v,v']*L_v,
```

with `L_v` the icosahedral link. Freeze:

1. every old cone edge and link edge has spacelike squared length `a^2`;
2. the pole has squared interval `-rho a^2`, `rho>0`;
3. the 12 final cone edges have squared lengths `q_u a^2`, `q_u>0`;
4. every final tetrahedron `[v',u_i,u_j,u_k]` must be strictly spacelike;
5. every four-simplex must have signature `(-,+,+,+)`;
6. use the ordinary real Lorentzian Regge action on timelike internal
   triangles, with zero volume/cosmological coefficient and no matter or
   higher-curvature term;
7. vary only the internal pole while holding all `q_u` fixed.

For an internal hinge `h_u=[v,v',u]`, define

```text
epsilon_u=2*pi-sum_(five incident simplices) theta,
w_u=d(A_u/a^2)/d rho>0.
```

After removing the common positive factor `d rho/d tau`, the pole equation is

```text
E(rho,q)=sum_u epsilon_u w_u=0.
```

This is one constraint on boundary data, not 12 separate zero-deficit
conditions.

## 3. Frozen general simplex derivation

For one link triangle, put `v=0`, `y=v'-v`, and use its three old edge vectors
`u_i`. In units `a=1`, set

```text
y.y=-rho,
u_i.u_i=1,
u_i.u_j=1/2,
c_i=y.u_i=(1-rho-q_i)/2.
```

The verifier must derive

```text
det G = -Q/4,
Q = c1^2+c2^2+c3^2
    +(c1-c2)^2+(c1-c3)^2+(c2-c3)^2+2rho > 0.
```

Because the link block is positive definite, every such Gram matrix has
inertia `(3 positive,1 negative)` for arbitrary real `c_i` and `rho>0`.

At the hinge spanned by `(y,u_1)`, project `u_2,u_3` into its two-dimensional
spacelike orthogonal complement. Their projected Gram matrix `P` must be
derived as

```text
P11=[4c1^2-4c1c2+4c2^2+3rho]/[4(c1^2+rho)],
P22=[4c1^2-4c1c3+4c3^2+3rho]/[4(c1^2+rho)],
P12=[2c1^2-2c1c2-2c1c3+4c2c3+rho]/[4(c1^2+rho)],

theta=acos[P12/sqrt(P11 P22)].
```

The internal timelike area and positive weight are

```text
A_u/a^2=(1/4)sqrt[4rho+(1-q_u-rho)^2],
w_u=(1+q_u+rho)/[4 sqrt(4rho+(1-q_u-rho)^2)].
```

At deterministic asymmetric points, explicit coordinates and facet normals
in `R^(1,3)` must independently reproduce the projected angle.

## 4. Disclosed target-found shell witness

Choose one vertex `u_0` of the icosahedral link. Its graph-distance shells
have sizes

```text
1,5,5,1.
```

The order-ten stabilizer of `u_0` is transitive on each shell, so a
stabilizer-invariant final metric has four values. The following values were
found exploratorily after looking for a stationary pole and are therefore
frozen as a **fitted witness**, not a prediction:

```text
rho = 1/4,
(q_0,q_1,q_2,q_3) = (x, 3/2, 4/5, 3/2),
x in [11/25,9/20] = [0.44,0.45].
```

The 20 link triangles split into four shell types, five of each:

```text
(0,1,1), (1,1,2), (1,2,2), (2,2,3).
```

Using Arb interval arithmetic, the verifier must establish on the complete
closed `x` interval:

1. all final squared edge lengths are positive;
2. every principal minor of all four final-tetrahedron Gram types is strictly
   positive;
3. every angle radicand/denominator stays on its real nonsingular branch;
4. every four-simplex remains Lorentzian and every internal hinge timelike;
5. every `w_u` is strictly positive.

At the exact rational endpoints, certify

```text
E(11/25)<0,
E(9/20)>0.
```

Continuity then proves at least one root. The verifier should additionally
certify `dE/dx>0` throughout the bracket, proving the root there is unique,
and report a high-precision enclosure/value without treating its decimals as
an exact constant.

The preliminary root is

```text
x*=0.44333089835748125745... .
```

At that root, report all four deficit and weight orbits. The expected
mechanism is cancellation between positive and negative deficits, not local
flatness at every hinge.

## 5. Mandatory canonicity attack

The full order-120 icosahedral automorphism group acts transitively on the 12
link vertices. Therefore a bare-geometry invariant scalar assignment has
only one length orbit and is exactly the symmetric family already killed.

The four-shell witness becomes invariant only after choosing one of 12 link
directions. Neither that direction nor the fitted values `1/4,3/2,4/5` are
selected by the bare local geometry. The verifier must record:

- full automorphism-invariant length dimension: `1`;
- selected-vertex-stabilizer length dimension: `4`;
- stationary equations in this ansatz: `1`;
- consequently, existence of a root cannot select the four boundary lengths
  and pole from geometry alone.

Do not claim a four-dimensional solution manifold solely from variable minus
equation counting. A local implicit-function statement requires a certified
nonzero derivative; only the root bracket derivative is required here.

## 6. Decision boundary

- **DERIVED ASYMMETRIC VACUUM EXISTENCE:** all causal/admissibility interval
  checks pass and the pole equation has the certified root.
- **REFUTED UNIVERSAL ANGLE NO-GO:** one admissible asymmetric simplex has
  angle above `2*pi/5`, so the symmetric bound cannot be promoted.
- **STRUCTURAL/FITTED NON-SELECTION:** the witness was target-found, needs an
  extra link direction, and leaves continuous boundary freedom.
- **REFUTED WITNESS:** the endpoint signs, causal domain or angle convention
  fails.
- **OPEN:** a theory-selected asymmetric state, pre/post Legendre map,
  constraint algebra, causal global schedule and continuum limit.

The result cannot be described as a derived physical clock even if the root
exists. It only proves that asymmetry can evade the local symmetric vacuum
obstruction.

Only the targeted verifier and static registry guard may run. No full suite
and no PDF build.
