# Result: generic duration freedom with two exceptional velocity obstructions

Date: 2026-08-21

## Provenance

```text
prior-art gate                                      62d92ab
frozen primary protocol                            c577419
registered primary implementation                  bd68c52
preserved direct-limit timeout                     85f6752
exact scaled-jet correction                        d8dde35
preserved radical-normalization timeout            cf909ce
frozen positive-radical inventory                  71b8312
first normalization implementation                 1046d01
preserved radical-recombination failure            10be749
recombination correction                           d6bc954
preserved first complete OPEN result               b181aa6
frozen exceptional/numeric adjudication            8acf139
corrected primary implementation                   98acd61
primary result and adversarial protocol             44a6ab2
frozen adversarial composition criterion           9f08aa0
registered adversarial implementation              bc7e1d0
preserved first adversarial disagreement           aae64d5
adversarial convergence adjudication               3eb07fb
accepted adversarial artifact                      58c6dff
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_generic_velocity_next_order.py
  13/13 PASS

reproducible/verify_gravity_600cell_generic_velocity_next_order_adversarial.py
  12/12 PASS
```

Accepted artifacts:

```text
reproducible/gravity_600cell_generic_velocity_next_order.json
SHA-256 4bc69490fc83a193b6ac2cbd8dbe291415a13b60e4dbcce4f499bf70152e5b18

reproducible/gravity_600cell_generic_velocity_next_order_adversarial.json
SHA-256 3ab16e6d19b527590b3dce6e8b3caa093efb6cc504a2a7824362ffc529a83a05
```

All four failed/timeout primary artifacts and the first adversarial
disagreement remain in the repository.  No full suite was run.

## Headline

```text
GENERIC_NEXT_ORDER_EXCEPTIONAL_BRANCHES_ADVERSARIALLY_CORROBORATED
```

> **DERIVED EXACT / STRUCTURAL, adversarially corroborated:** on the fixed
> homogeneous cellular 600-cell Regge-plus-conserved-dust action, the same
> generic incoming canonical state admits an arbitrary-small-duration
> one-slab jet through the first correction for every nonzero velocity except
> one time-reversal pair.  At that pair the equations have no endpoint jet at
> this order.  Neither case selects a nonzero duration.

The two exceptional points are obstructions, not ticks and not a derived
limiting speed.

## Exact theorem

For

```text
L_minus=1,
M=mu(v),
p0=p(v),
L_plus=exp(v h+a h^2),
rho=h^2,
v real and nonzero,
```

use the already-corroborated leading mass and momentum branch.  Set

```text
x=v^2,
r=sqrt(x+4),
q=sqrt(3x+8),
epsilon(x)=2*pi-5*acos((x+2)/(2(x+3))),

K(x)=10*r-(x+3)*q*epsilon(x),
B(x)=5*x*r+2*(x+3)*q*epsilon(x),

R(v)=1440*v/[r*q*(x+3)*(x+4)].
```

The first nonzero lapse residual factors exactly as

```text
C1(v,a)=R(v)[K(v^2)*a+B(v^2)].
```

The momentum residual satisfies the exact cross identity

```text
C1(v,0)*coefficient_a(P1)
 -P1(v,0)*coefficient_a(C1)=0.
```

Therefore, whenever `K!=0`, both equations have the same unique endpoint
coefficient

```text
a(v)=-B(v^2)/K(v^2).
```

This is a compatible formal local jet for arbitrary `h` through the
registered order; the equations determine the endpoint acceleration
coefficient but do not determine `h`.  Convergence to an exact analytic
finite-height family is not proved.

## Complete exceptional set

Define

```text
H(x)=(x+3)*sqrt((3x+8)/(x+4))*epsilon(x),
K(x)=sqrt(x+4)[10-H(x)].
```

The proof of completeness is analytic, not a sign scan:

- `epsilon(0)>0` follows from the exact value
  `cos(2*pi/5)=(sqrt(5)-1)/4<1/3`;
- `epsilon'(x)=5/[(x+3)sqrt(x+4)sqrt(3x+8)]>0`;
- `(3x+8)/(x+4)` has derivative `4/(x+4)^2>0`;
- hence `H` is strictly increasing on `x>=0`;
- `H(0)<pi*sqrt(2)<10`, while `H(x)->infinity`.

Consequently there is exactly one positive root `x_star` of `K`:

```text
x_star =5.6463444131338775223511342535853545173861284825292474951319...
v_star =2.3762037819037906586291844426115698345081682919472298042590...
```

The complete exceptional velocity set is

```text
v in {-v_star,+v_star}.
```

At either point `K=0`, but `B>0` and `R!=0`, so

```text
C1(v,a)=R(v)B(v^2) != 0
```

for every real `a`.  Thus the degree drop does not create additional
solutions; it removes the local branch.

## Mechanically independent replication

The primary route first rewrote the complete action in scaled variables and
then differentiated its exact jet.  The adversarial route reversed the
decisive order:

1. differentiate the unexpanded action in `rho`, `L_minus` and `L_plus`;
2. only then set `L_plus=1+tau*q`, `rho=tau^2`;
3. take partial derivatives at `(tau,q)=(0,v)`;
4. insert `q'=a+v^2/2`;
5. compare with the primary artifact only afterward.

It reproduced `C1`, `P1`, `K`, `B` and the cross identity exactly.  A new
bracket `5<x_star<6` agreed with the primary root beyond 80 decimals.  New
one-slab controls at

```text
v in {-6/5,2/5,11/5},
a in {-1/5,2/9},
h in {1/500,1/1000}
```

all converged at the preregistered first order.  Shifting the mass or momentum
by `1/10` produced exactly `-4*pi/5` or `-1/10`; shifting `K` moved the root.

## Conditional composition

On the generic domain `K!=0`, the primary exact route found that two
stationary half slabs have

```text
endpoint defect       =0,
final momentum defect =0,
action defect         =0
```

at the registered next order.  The adversarial route did not rederive these
three identities symbolically; it tested the unexpanded two-slab action at
new velocities.  Seventeen of eighteen initial convergence controls passed.
The action defect at `v=11/5`, near the exceptional point where
`a=-29.27279...`, was pre-asymptotic.  The preserved failure had order
`0.5678`; the registered smaller-height diagnostic gave

```text
0.92575, 0.96478, 0.98283,
```

converging to first order.  The exact composition statement is therefore
**DERIVED EXACT on the primary route, with independent numerical controls**,
not a second symbolic derivation.

## Physical status ledger

| Claim | Status |
|---|---|
| Generic duration remains free through the first correction | **DERIVED EXACT / ADVERSARIALLY CORROBORATED** |
| Exactly two exceptional velocities have no local endpoint jet | **DERIVED EXACT / ADVERSARIALLY CORROBORATED / STRUCTURAL** |
| Next-order stationary half-step composition on `K!=0` | **DERIVED EXACT primary; independently numerically controlled** |
| `v_star` is a limiting signal speed or `c` | **NOT DERIVED** |
| Exceptional pair is stable under carrier refinement | **OPEN** |
| Duration remains free at the following order | **OPEN** |
| Nonzero isolated relational tick | **NOT DERIVED** |
| Absolute classical tick | **DERIVED NEGATIVE under global scale covariance** |
| External novelty | **OPEN** |

## Physical interpretation

This result moves the tick boundary one order further without crossing it.
For almost every velocity the formal equations still allow any sufficiently
small proper height through this order and merely adjust the endpoint
coefficient `a(v)`.  At the two critical velocities even that jet is absent.
There is nowhere an isolated positive `h` selected by this calculation.

The number `v_star` is dimensionless and belongs to this homogeneous
finite-carrier action.  Calling it the speed of light would require a
propagation calculation for inhomogeneous modes, a conversion between spatial
and temporal units and refinement stability; none is present here.

## Post-result prior-art audit

The learned terms `singular discrete Legendre transform`, `degree drop` and
`variational-integrator resonance` recover known structural mechanisms:

- [Bahr--Dittrich](https://arxiv.org/abs/0905.1670) relate broken discrete
  gauge symmetry in curved Regge calculus to lapse-dependent
  pseudo-constraints.
- [Schmitt--Leok](https://arxiv.org/abs/1609.02309) relate variational-
  integrator resonance behavior to ill-posed boundary-value problems for
  exact discrete generating functions.
- [Marrero--Martin de Diego--Martinez](https://arxiv.org/abs/1608.01586)
  provide the exact-discrete-Lagrangian reference against which a finite
  variational integrator is compared.

These sources make a discretization singularity a plausible interpretation;
they do not identify the present `K(x)` or prove that `v_star` has physical
meaning.  No located primary source prints this 600-cell coefficient.
External novelty remains **OPEN**.

## Next load-bearing gate

Add the next endpoint coefficient without changing the state:

```text
L_plus=exp(v h+a(v)h^2+c h^3).
```

On the generic domain `K!=0`, extract the next lapse and momentum residuals
and ask whether one common `c(v)` exists.  This is target-free and has no
finite search bound.

- A common `c(v)` moves local duration freedom one more order and suggests
  looking for an all-orders identity.
- No common `c(v)` gives the first generic local pseudo-constraint
  obstruction.
- An isolated finite tick still requires a separate exact finite-root census
  and carrier/refinement selection; an asymptotic obstruction alone is not a
  tick.
