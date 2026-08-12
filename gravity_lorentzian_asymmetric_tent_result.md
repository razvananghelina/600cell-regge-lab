# An asymmetric Lorentzian vacuum tent exists but is not selected

Date: 2026-08-12

Preregistered protocol commit: `4a63f66`

Registered verifier:
`reproducible/verify_gravity_lorentzian_asymmetric_tent.py`

Machine-readable result:
`reproducible/gravity_lorentzian_asymmetric_tent.json`

## Headline

> **DERIVED ASYMMETRIC VACUUM EXISTENCE.** On the canonical local carrier
> `[v,v']*L_v`, with ordinary real Lorentzian Regge calculus, zero
> volume/cosmological coefficient and the preregistered asymmetric spacelike
> boundary data, the internal pole equation has one certified real root in
> the frozen interval `x in [11/25,9/20]`.

> **STRUCTURAL/FITTED NON-SELECTION.** This is an existence theorem, not a
> derived clock. The witness was found by searching for a zero, requires a
> selected link direction, and fixes three shell lengths by hand. Bare
> icosahedral symmetry selects only the symmetric one-orbit family already
> ruled out in vacuum.

The targeted verifier passes `21/21`. No full suite was run.

The result breaks a genuine obstruction: the angle bound proved for the
symmetric Lorentzian tent does not extend to arbitrary spacelike final
boundary data. It does not yet provide a law that evolves one spatial frame
to the next.

## 1. Complete scope

The carrier is

```text
T_v=[v,v']*L_v,
```

where `L_v` is the combinatorial icosahedral link of a 600-cell vertex. The
calculation assumes:

1. old cone and link edges have squared length `a^2`;
2. the pole has squared interval `-rho a^2`, with `rho>0`;
3. the twelve final cone edges have squared lengths `q_u a^2`, with `q_u>0`;
4. every final tetrahedron is strictly spacelike;
5. every four-simplex has signature `(-,+,+,+)`;
6. timelike hinges use their real area magnitude and the ordinary rotation
   angle in the spacelike normal plane;
7. the action is the zero-volume Lorentzian Regge action, without matter or
   higher-curvature terms;
8. only the internal pole is varied; all boundary `q_u` are held fixed.

These hypotheses must accompany the existence claim. In particular this is
not stationarity under arbitrary boundary variations and not a solution of a
derived Hamiltonian constraint system.

## 2. Exact asymmetric simplex

For one link triangle, in units `a=1`, write

```text
y.y=-rho,
u_i.u_i=1,
u_i.u_j=1/2,
c_i=y.u_i=(1-rho-q_i)/2.
```

The signed four-simplex Gram determinant is

```text
det G=-Q/4,

Q=c1^2+c2^2+c3^2
 +(c1-c2)^2+(c1-c3)^2+(c2-c3)^2+2rho.
```

Thus `Q>0` for `rho>0`. The old link block has eigenvalues
`(2,1/2,1/2)`, so its Schur complement proves inertia `(3 positive,1
negative)` without a symmetry assumption on the `q_i`.

At the hinge spanned by `(y,u_1)`, orthogonal projection of `u_2,u_3` gives

```text
P11=(4c1^2-4c1c2+4c2^2+3rho)/[4(c1^2+rho)],
P22=(4c1^2-4c1c3+4c3^2+3rho)/[4(c1^2+rho)],
P12=(2c1^2-2c1c2-2c1c3+4c2c3+rho)/[4(c1^2+rho)],

theta=acos[P12/sqrt(P11 P22)].
```

The timelike hinge area and its positive pole weight are

```text
A_u/a^2=(1/4)sqrt[4rho+(1-q_u-rho)^2],
w_u=d(A_u/a^2)/drho
   =(1+q_u+rho)/[4sqrt(4rho+(1-q_u-rho)^2)]>0.
```

All formulae were derived symbolically in the verifier. For all six ordered
shell-simplex types, an independent construction in `R^(1,3)` reconstructed
the edge Gram matrices and computed the facet normals directly. Its maximum
Gram and angle residuals were both `2.22e-16`.

## 3. The frozen witness

Select one vertex of the icosahedral link. Its graph-distance shells contain

```text
1,5,5,1
```

vertices. Freeze

```text
rho=1/4,
(q_0,q_1,q_2,q_3)=(x,3/2,4/5,3/2),
x in [11/25,9/20].
```

The 20 link triangles split into four types, five of each:

```text
(0,1,1), (1,1,2), (1,2,2), (2,2,3).
```

For every internal hinge `u`, define

```text
epsilon_u=2*pi-sum_(five incident simplices) theta,
E=sum_u epsilon_u w_u.
```

After removing the common nonzero conversion between `rho` and pole proper
time, `E=0` is the pole equation.

Arb interval arithmetic certifies

```text
E(11/25) = -0.0070255587018337244864... < 0,
E(9/20)  =  0.0140502821146620281889... > 0,
dE/dx > 0 throughout [11/25,9/20].
```

Therefore continuity proves a root and strict monotonicity proves it is the
only root in this bracket. High-precision bisection reports

```text
x*=0.443330898357481257452440941829425383965150474294...
```

The decimal is not promoted to an exact constant; the theorem is the
rational bracket plus interval signs and derivative.

On the entire closed bracket, not merely at the root, the verifier certifies:

- every final edge square is positive;
- every principal minor of all final tetrahedral Gram matrices is positive;
- every four-simplex is Lorentzian;
- every internal hinge is timelike;
- every projected angle stays on a real nonsingular branch;
- every area weight is strictly positive.

The smallest lower bound among the final-tetrahedron principal minors is
`0.3049`. Hence the witness is not sitting on a causal degeneracy.

## 4. How stationarity occurs

At the root, the four deficit orbits are

```text
( 2.071667979072602...,
 -0.569222065406002...,
  0.192909410405034...,
  0.423969940527738... ).
```

The corresponding positive weight orbits are approximately

```text
(0.404728756567238..., 0.55,
 0.511860573674895..., 0.55).
```

Thus the equation vanishes by weighted cancellation between positive and
negative curvature defects. None of the four orbits is locally flat. This
distinguishes the solution from the symmetric Euclidean golden construction,
where the common deficit itself vanished.

At the root, the largest angle is `84.7148718349...` degrees, while the
symmetric fivefold target is `72` degrees. Arb also certifies an angle above
`72` degrees uniformly on the complete root bracket. Therefore the exact
symmetric theorem `theta<72 degrees` cannot be generalized to arbitrary
asymmetric boundary data.

## 5. Canonicity attack

The combinatorial calculation finds the full icosahedral graph automorphism
group of order `120`, transitive on all twelve link vertices. Consequently a
scalar final-edge assignment invariant under the bare local geometry has
only one orbit. That is precisely the symmetric family already proved to
have no zero-volume Lorentzian stationary pole.

After selecting one link vertex, its order-ten stabilizer has four orbits of
sizes `1,5,5,1`; the witness is invariant only under this smaller group. The
data therefore contain:

```text
bare-symmetry length orbits:       1,
selected-vertex length orbits:     4,
pole equations in the frozen test: 1.
```

The certified nonzero derivative does prove that this particular root is
regular with respect to `x`. It does not select the chosen direction,
`rho=1/4`, or the three frozen ratios `3/2,4/5,3/2`. Variable-minus-equation
counting is not used as a proof of a solution-family dimension.

Accordingly:

- existence is **DERIVED** under the complete frozen hypotheses;
- the numerical boundary witness is **FITTED / STRUCTURAL**;
- a physical tick selected by the theory remains **OPEN**.

## 6. Physical meaning

This is the first real Lorentzian stationary local tent found on this
600-cell carrier without adding a volume term. It shows that the carrier is
not intrinsically incapable of local Lorentzian Regge dynamics.

It does not yet show that physical space is a 600-cell, that one root defines
a universal tick, or that repeated local tents compose into a causal global
evolution. Those claims require at least:

1. a target-independent rule selecting the asymmetric boundary state;
2. compatible pre/post canonical data for adjacent tent moves;
3. a constraint algebra or an equivalent gauge-redundancy construction;
4. a conflict-free global update schedule;
5. a refinement/continuum limit with effective Lorentz symmetry and GR.

The next falsifiable gate is selection, not another parameter search.

## 7. Status ledger

| Claim | Status |
|---|---|
| The general asymmetric simplex determinant is Lorentzian for `rho>0` | **DERIVED** |
| The displayed projected-angle and area formulas are exact | **DERIVED** |
| The frozen bracket stays strictly causal and nondegenerate | **DERIVED** |
| The frozen pole equation has a root in `[11/25,9/20]` | **DERIVED** |
| That root is unique within the frozen bracket | **DERIVED** |
| Stationarity is cancellation, not zero local curvature | **DERIVED** |
| The symmetric `theta<72 degrees` bound extends to all asymmetries | **REFUTED** |
| Bare icosahedral symmetry selects the witness | **REFUTED** |
| The witness was target-independent | **REFUTED; it was disclosed as target-found** |
| The theory currently selects a Lorentzian local tick | **OPEN** |
| Adjacent tents form a consistent global evolution | **OPEN** |
| This supplies `c`, `G`, the Planck time or Planck mass | **OPEN** |

## 8. Reproduction boundary

The witness, bracket, formulas, hypotheses and decision labels were committed
in `4a63f66` before the verifier existed. Only the targeted verifier was run;
it passed `21/21`. No full-suite or PDF build was performed.
