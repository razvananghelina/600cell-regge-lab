# The asymmetric Lorentzian tent has a regular local pole equation

Date: 2026-08-12

Preregistered protocol commit: `1035c54`

Witness/result commit inherited: `cc71574`

Registered verifier:
`reproducible/verify_gravity_lorentzian_tent_regular_evolution.py`

Machine-readable result:
`reproducible/gravity_lorentzian_tent_regular_evolution.json`

## Headline

> **DERIVED REGULAR LOCAL POLE.** At the certified asymmetric Lorentzian
> zero-volume tent solution, the internal Regge equation has
>
> ```text
> partial E/partial rho
> =6.3787093592142737476176222376... > 0.
> ```
>
> Arb separates this derivative from zero on the entire preregistered root
> bracket. Therefore the implicit-function theorem supplies a unique smooth
> local pole `rho=rho(q)` for every sufficiently nearby set of boundary
> lengths.

> **SUBSEQUENT ADVANCE.** The follow-up
> `gravity_lorentzian_tent_legendre_result.md` constructs the complete local
> boundary/corner action and finds a regular `12x12` pre/post star Legendre
> map. Constraint matching with adjacent moves remains unconstructed.

The targeted verifier passes `10/10`. No full suite was run.

## 1. Framing correction: laws do not select all initial data

The previous canonicity audit correctly found that bare icosahedral symmetry
does not choose the target-found four-shell boundary state. Treating that as
a necessary failure of dynamics was too strong.

In a canonical theory, a law relates boundary configurations and momenta; it
does not normally choose the unique state of the universe without initial
data. In canonical Regge calculus, new edge lengths may be free at one
Pachner/tent step and can be fixed later by pre-constraints. Hamilton's
principal function generates the canonical relation. This established
external formalism is summarized by Dittrich and Hoehn,
[arXiv:1108.1974](https://arxiv.org/abs/1108.1974).

Thus two statements must be kept separate:

- **DERIVED NON-PREDICTION:** the bare 600-cell does not predict the frozen
  `q_u` profile;
- **NOT REQUIRED FOR DYNAMICS:** a unique symmetry-selected boundary profile
  is not a prerequisite for an evolution law.

What the move must do first is solve its bulk equation for supplied boundary
data. That is the gate tested here.

## 2. Frozen equation and hypotheses

No geometry was changed from the preregistered witness. Use

```text
T_v=[v,v']*L_v,
rho=tau^2/a^2,
(q_0,q_1,q_2,q_3)=(x,3/2,4/5,3/2),
rho_0=1/4.
```

All old cone/link edges have squared length `a^2`; every final tetrahedron is
strictly spacelike; every four-simplex has signature `(-,+,+,+)`; internal
hinges are timelike; the volume coefficient, matter and higher-curvature
terms vanish. Only the pole is a bulk variable. The twelve `q_u` are boundary
data.

For each internal hinge,

```text
epsilon_u=2*pi-sum_(five simplices) theta,
w_u=d(A_u/a^2)/d rho>0,
E(rho,q)=sum_u epsilon_u w_u.
```

The root remains

```text
x*=0.443330898357481257452440941829425383965150474294...
```

defined rigorously by its rational bracket and interval signs, not by this
decimal.

## 3. Interval regularity theorem

The verifier rebuilds the combinatorial icosahedron and the complete angle,
area and deficit sum. It independently reproduces

```text
E(11/25)<0,
E(9/20)>0,
partial E/partial x>0.
```

Before this result, the protocol froze a one-box Arb test followed, only if
needed, by exactly 16 rational subintervals. The single complete bracket was
already sufficient:

```text
partial E/partial rho on x in [11/25,9/20]
    is enclosed by [6e+0 +/- 0.874],
```

whose lower endpoint is strictly positive. At the isolated root the tighter
ball is

```text
E_rho=6.3787093592142737476176222376236858093299...
```

The root is not on a causal boundary: the smallest certified lower bound
among all principal minors of the final tetrahedral Gram matrices is
`0.331695327`.

Since `E=0` and `E_rho!=0`, the implicit-function theorem proves a unique
local function

```text
rho=rho(q_1,...,q_12)
```

near this solution. This is a genuine local theorem, not dimension counting.

## 4. Boundary response

The twelve derivatives `E_(q_u)` are constant on the four stabilizer shells.
When one changes all lengths in a shell together, the implicit responses are

```text
d rho/d q_shell =
(-0.330540583197568...,
 +0.019559921904892...,
 -0.648166608783582...,
 +0.062697525643235...).
```

For an individual edge in shells of sizes `(1,5,5,1)`, they are

```text
(-0.330540583197568...,
 +0.003911984380978...,
 -0.129633321756716...,
 +0.062697525643235...).
```

The collective derivatives equal the sums of their individual edge
derivatives exactly within Arb enclosures. Independent 100-digit centered
finite differences reproduce every pole, individual-edge and shell
derivative with maximum declared relative error

```text
1.01e-16.
```

These numbers are response coefficients of the frozen local move, not
universal constants and not candidates for fitting to observations.

## 5. Pole Hessian

With `tau=sqrt(rho)`, the stationarity derivative is, up to the common action
normalization,

```text
dS/dtau=2 tau E(tau^2,q).
```

At `tau_0=1/2` and `E=0`,

```text
d^2S/dtau^2=E_rho=6.37870935921427... > 0.
```

Thus the frozen-sign action is nondegenerate and locally convex along the
pole direction. Reversing the overall action sign reverses “minimum” versus
“maximum” but cannot destroy regularity. No stability statement is made in
the full boundary/bulk configuration space.

## 6. What is now established

The mathematical chain is now

```text
canonical icosahedral tent carrier
 -> admissible asymmetric Lorentzian solution
 -> nonzero pole Hessian
 -> locally unique internal pole for nearby boundary data.
```

This is stronger than “one numerical zero exists”. It proves that the zero
is not a fragile tangent contact or a lapse-like flat direction at this
curved configuration.

It still does not establish:

1. the Lorentzian boundary/corner action on a fixed global branch;
2. pre- and post-momenta for the twelve old/new spatial edges;
3. the rank and null directions of the mixed Lagrangian Hessian;
4. constraint matching between adjacent tent moves;
5. a conflict-free global update schedule;
6. a continuum limit, universal light cone, `G` or Planck units.

The next correct gate is the first three items together: construct the full
Lorentzian Hamilton principal function and its Legendre map. Another search
for a prettier root would add no evidence.

## 7. Status ledger

| Claim | Status |
|---|---|
| The previous asymmetric root is independently reproduced | **DERIVED** |
| `E_rho` is nonzero at the root | **DERIVED WITH ARB** |
| Its sign stays positive on the complete root bracket | **DERIVED WITH ARB** |
| The pole is locally unique for nearby boundary data | **DERIVED BY IFT** |
| The pole-only Hessian is nondegenerate | **DERIVED** |
| The frozen-sign pole direction is locally convex | **DERIVED, POLE-ONLY** |
| Boundary response respects the four stabilizer shells | **DERIVED** |
| A unique boundary state must be selected for dynamics to exist | **FRAMING REFUTED** |
| Bare geometry predicts the target-found boundary state | **REFUTED** |
| The fixed-link 12-edge star pre/post Legendre map is regular | **DERIVED COMPUTATIONALLY SUBSEQUENTLY** |
| The global/adjacent-move canonical map is regular | **OPEN** |
| Adjacent moves satisfy a first-class constraint algebra | **OPEN** |
| A global physical clock has been derived | **OPEN** |

## 8. Reproduction boundary

The hypotheses, derivative test, 16-box fallback and decision labels were
committed in `1035c54` before the derivative sign was computed. The targeted
verifier passes `10/10`. No full suite or PDF build was run.
