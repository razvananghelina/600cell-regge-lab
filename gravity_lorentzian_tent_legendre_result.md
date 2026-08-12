# The asymmetric Lorentzian tent generates a regular local Legendre map

Date: 2026-08-12

Parent protocol commit: `9b10ed4`

Mixed-causal angle correction commit: `bd4eaa5`

Registered verifier:
`reproducible/verify_gravity_lorentzian_tent_legendre.py`

Machine-readable result:
`reproducible/gravity_lorentzian_tent_legendre.json`

## Headline

> **DERIVED COMPUTATIONAL REGULAR LOCAL LEGENDRE MAP.** For the certified
> asymmetric Lorentzian zero-volume tent, using the fixed plus complex-angle
> branch, the complete local Regge action is real and its on-shell old/new
> mixed Hessian has stable numerical rank
>
> ```text
> rank W=12/12,
> s_min/s_max=0.0078960369... .
> ```
>
> After eliminating the one regular internal pole, the twelve old star
> lengths and their admissible pre-momenta locally determine the twelve new
> star lengths in this fixed-link sector.

> **OPEN GLOBAL CANONICAL DYNAMICS.** One regular local generating function
> is not yet a global evolution: adjacent-move constraint matching,
> overlapping schedules, variable link data, first-class closure and the
> continuum limit remain open.

The final targeted verifier passes `16/16`. No full suite was run.

The rank statement is a robust finite-difference computation under a
preregistered gap rule, not an exact symbolic determinant or interval proof.

## 1. Complete scope

The local four-ball is

```text
T_v=[v,v']*L_v,
```

with `L_v` the combinatorial icosahedron. In squared-length coordinates and
units `a=1`, vary

```text
p_u=s(v,u),       u=1,...,12,
q_u=s(v',u),      u=1,...,12,
s(v,v')=-rho.
```

All 30 link edges remain fixed at squared length one. The evaluation point is

```text
p_u=1,
rho=1/4,
q_u=(x*,3/2,4/5,3/2) on shells (1,5,5,1),
x*=0.443330898357481... .
```

Every four-simplex is Lorentzian, every final tetrahedron is spacelike and
every internal hinge is timelike. The action has no volume/cosmological,
matter or higher-curvature term.

The test covers the local 12-by-12 cone-edge star sector. It does not vary
the common link, the remaining 600-cell edges or a second tent move.

## 2. Complete complex action

Use the plus complex-angle branch of Borissova and Dittrich,
[arXiv:2303.07367](https://arxiv.org/abs/2303.07367). With signed squared
volumes,

```text
i S_Regge^+=sum_h sqrt(+V_h) epsilon_h^+,
S_L=-i(i S_Regge^+).
```

The verifier reconstructs the complete hinge census:

| hinge type | count | simplex incidence | boundary `k` |
|---|---:|---:|---:|
| internal `[v,v',u]` | 12 | 5 | 2 |
| old boundary `[v,u,w]` | 30 | 2 | 1 |
| new boundary `[v',u,w]` | 30 | 2 | 1 |
| common-link corner | 20 | 1 | 0 |

There are `92` hinges total. All 20 four-simplex Gram matrices have exactly
one negative eigenvalue.

At the witness:

```text
S_L=1.570646572792366...,
|Im S_L|=2.8e-16,
max |Im gradient S_L|=6.7e-16.
```

The twelve timelike-hinge curvatures are real. All 80 spacelike
boundary/corner curvatures are purely imaginary, so multiplication by
`-i` gives a real Lorentzian action. The closest logarithm argument stays
`0.271354...` away from zero.

Changing the common-link corner integer from `k=0` to `k=1` shifts the action
by the constant

```text
-27.2069904635... i
```

and changes every canonical derivative by exactly zero numerically. This is
because all link-triangle areas are frozen in the tested sector.

## 3. The sign error found by hostile verification

The first implementation used the tempting normalized-normal expression

```text
-N_ab/[sqrt(+N_aa)sqrt(+N_bb)]
```

as the complex cosine for every hinge. It passed the internal-angle and
reality checks but failed decisively:

```text
full-action/Schlaefli gradient error = 8.641,
relative Hessian antisymmetry        = 0.3811.
```

Therefore its preliminary rank was discarded.

The cause is exact and geometric. For adjacent facet normals of different
causal type, that expression has the wrong sign. The source formula A.38 is

```text
cos(theta_ab^+)
 =16 [partial V_sigma/partial s_ab]
  /[sqrt(+V_face_a)sqrt(+V_face_b)].
```

Across all 200 simplex-hinge incidences, the verifier finds:

- same causal type: A.38 equals the normalized-normal expression to
  `8.9e-16`;
- opposite causal type: A.38 equals its negative to `1.1e-15`.

With A.38 and the companion signed-volume sine formula A.39, 420 direct
directional checks at the witness and two nearby causal configurations obey
the complex Schlaefli identity with maximum residual `6.72e-10`.

This correction was source-forced, not chosen to obtain full rank. Because
the corrected full-rank result had already been seen, it was disclosed in
commit `bd4eaa5` before the independent action-Hessian confirmation was
added.

## 4. Pre/post momenta

The complex Schlaefli identity reduces first derivatives to

```text
partial_x S_L
 =-i sum_h partial_x[sqrt(+V_h)] epsilon_h^+.
```

Centered differences of the complete action reproduce all 25 components
with maximum relative error `2.81e-9`. The pole component is

```text
S_rho=-3.4e-15,
```

reproducing the independently certified stationary equation.

In squared-length coordinates,

```text
P_u^-=-partial S_L/partial p_u,
P_u^+=+partial S_L/partial q_u.
```

They are real and constant on the four stabilizer shells:

```text
P^- =(-1.09310442919, +0.08068234831,
      -0.48865353772, +0.31268558127),

P^+ =(-1.41763075108, +0.40589720792,
      -0.08299808958, +0.29964273756).
```

These are local canonical data of the frozen witness, not predicted particle
momenta or observable constants.

## 5. Pole elimination and mixed Hessian

Let

```text
b=(p_1,...,p_12,q_1,...,q_12).
```

The full `25x25` Hessian in `(b,rho)` was computed at three frozen steps:

```text
h=(2e-5,1e-5,5e-6).
```

Its maximum relative antisymmetry decreases to `6.27e-11`; the maximum
step discrepancy is `5.68e-9`. Its pole entry converges to

```text
S_rhorho=6.37870935921427...,
```

matching the independent Arb midpoint with maximum relative error `6.1e-9`.

Eliminate the pole using its nonzero Hessian:

```text
H_eff=S_bb-S_b,rho (S_rho,rho)^(-1) S_rho,b.
```

The old/new block

```text
W=partial^2 S_on-shell/(partial p partial q)
```

has singular values at the finest step

```text
1.24352168682, 1.24352168681,
0.72220977456, 0.72220977455,
0.66507655293,
0.63804809558, 0.63804809557,
0.63202943333,
0.06389733029,
0.03970923334, 0.03970923332,
0.00981889315.
```

All twelve exceed the relative rank threshold `1e-7`. The smallest/largest
ratio is

```text
0.007896036917...,
```

far from the preregistered inconclusive band. Singular values across the
three Hessian steps drift by at most `9.1e-8` relative.

The rank is unchanged by passing from positive squared lengths to ordinary
lengths, because this applies invertible diagonal Jacobians on the two sides
of `W`.

## 6. Independent full-action confirmation

After the sign correction was disclosed, the mixed block was rebuilt without
differentiating the Schlaefli gradient. Direct centered second differences of
the complete complex action supplied

```text
S_pq, S_p,rho, S_rho,q, S_rho,rho
```

and a second Schur complement. It gives

```text
rank=12,
s_min/s_max=0.00789712315...,
relative mixed-block error=7.56e-6,
maximum singular-value error=1.38e-4.
```

All are inside the correction protocol's frozen acceptance thresholds. The
rank conclusion therefore does not rest on the reduced-gradient code path.

## 7. Mathematical and physical consequence

The nonzero pole Hessian first gives `rho=rho(p,q)` locally. The full-rank
mixed block then makes the local Legendre relation invertible in the twelve
cone-edge directions. Under the stated finite-dimensional hypotheses, the
inverse-function theorem supplies a local canonical map between old and new
star data.

This is real progress toward dynamics:

```text
admissible Lorentzian solution
 -> regular internal pole
 -> real Hamilton principal function
 -> real pre/post momenta
 -> full-rank local star Legendre map.
```

It does **not** prove that 600-cell space is our universe, that the move has a
universal duration, or that many moves form GR. The most immediate missing
test is gluing: apply two adjacent tent moves and determine whether their
shared-edge post/pre momenta and constraints match.

## 8. Status ledger

| Claim | Status |
|---|---|
| Complete 92-hinge complex action on the frozen tent | **DERIVED COMPUTATIONAL** |
| The selected branch gives a real Lorentzian action and momenta | **DERIVED COMPUTATIONAL** |
| The naive normal cosine works for mixed-causal facets | **REFUTED** |
| A.38/A.39 satisfy the per-simplex complex Schlaefli identity | **DERIVED COMPUTATIONAL** |
| The pole entry matches the independent Arb derivative | **DERIVED COMPUTATIONAL CROSS-CHECK** |
| On-shell star mixed rank is `12/12` | **DERIVED COMPUTATIONAL, TWO ROUTES** |
| The local fixed-link Legendre relation is regular | **DERIVED COMPUTATIONAL** |
| Boundary state is uniquely selected by bare geometry | **REFUTED / NOT REQUIRED FOR A LAW** |
| Link-edge variations are included | **NOT CLAIMED** |
| Adjacent tent constraints and momenta match | **OPEN** |
| A global non-overlapping causal schedule exists | **OPEN** |
| The continuum theory has diffeomorphism symmetry and GR dynamics | **OPEN** |
| `c`, `G`, Planck time or Planck mass follows | **OPEN** |

## 9. Reproduction history

The original rank protocol was committed as `9b10ed4`. Its unqualified normal
formula failed `2/12` checks, so that output was rejected. The source-forced
mixed-causal correction and the already-seen preliminary rank were disclosed
in `bd4eaa5`. Only afterward was the independent direct-action Hessian route
added. The final targeted verifier passes `16/16`.

No full suite or PDF build was run.
