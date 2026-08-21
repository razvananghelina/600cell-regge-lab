# Result: generic-velocity leading reparametrization

Date: 2026-08-21

## Provenance

```text
prior-art gate                                  1fcab34
frozen primary protocol                        9472a15
frozen numerical convergence criterion         8af359d
registered primary implementation              dfa2688
preserved primary first failure                51bd4ce
corrected primary implementation               a5ebcea
primary result and adversarial protocol        7ed16ee
frozen adversarial radical normalizations      f424f31
registered adversarial implementation          74653d6
preserved adversarial first failure            4da1054
corrected adversarial implementation           0367d96
accepted adversarial artifact                  1fbf150
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_generic_velocity_composition.py
  11/11 PASS

reproducible/verify_gravity_600cell_generic_velocity_composition_adversarial.py
  9/9 PASS
```

Accepted artifacts:

```text
reproducible/gravity_600cell_generic_velocity_composition.json
SHA-256 8ded36f1fa00307fcb23369c25290c9f5bd701709762d6a865437c2507eabfc9

reproducible/gravity_600cell_generic_velocity_composition_adversarial.json
SHA-256 cd46c6c9d38e1b14fc32f09f8a2cf72039d28ed136ef4081fe0edba149a9b6b2
```

Both first-run failures and their correction protocols remain in the
repository.  No full suite was run.

## Headline

```text
GENERIC_VELOCITY_LEADING_REPARAMETRIZATION_ADVERSARIALLY_CORROBORATED
```

> **DERIVED EXACT / STRUCTURAL, adversarially corroborated:** on the exact
> homogeneous cellular 600-cell Regge-plus-conserved-dust action, the
> generic-velocity tangent has an interval-independent leading principal
> function, lapse constraint and incoming momentum.  One interval and two
> half intervals therefore agree at leading order for the same tangent state.

This does not select a tick.  It proves precisely that the interval factor is
still a reparametrization at leading order away from the turning point.

## Exact leading map

Use the dimensionless variables

```text
mu=M/L_minus,
tau=s e,
log(L_plus/L_minus)=s v e+O(e^2),
rho=s^2 e^2.
```

Define

```text
r(v)       =sqrt(v^2+4),
theta(v)   =acos((v^2+2)/(2(v^2+3))),
eta(v)     =asinh(v/sqrt(8(v^2+3))),
epsilon_v  =2*pi-5*theta(v).
```

Then

```text
S=s e L0(v,mu)+O(e^2),

L0(v,mu)
 =360 r(v) epsilon_v
  -1200 sqrt(3) v eta(v)
  -8 pi mu,

C(v,mu)
 =L0-v dL0/dv
 =1440 epsilon_v/r(v)-8 pi mu,

p_pre(v)
 =180 v epsilon_v/r(v)-600 sqrt(3) eta(v).
```

None of these expressions contains `s`.  The constraint is affine in `mu`
with nonzero coefficient and therefore has exactly one branch,

```text
mu(v)=180 epsilon_v/(pi r(v)).
```

The exact cosine bounds put `theta(v)` in `(pi/3,acos(1/3)]`, so this mass is
positive for every real `v`.  The mass is even, the momentum is odd and the
static control is recovered:

```text
mu(0)=90[2*pi-5*acos(1/3)]/pi,
p_pre(0)=0.
```

## Independent route and controls

The primary verifier first took the primitive action limit and derived its
Hamilton--Jacobi data.  The adversarial verifier instead:

1. differentiated the complete unexpanded action with respect to `rho` and
   `L_minus`;
2. substituted the linear tangent `L_plus=1+v*tau`, `rho=tau^2`;
3. took the exact direct limits;
4. normalized only seven preregistered positive radical factorizations;
5. read the primary formulas only after those expressions existed.

All action, constraint, momentum and mass formulas matched exactly.  New
100-decimal controls at

```text
v in {-7/10,3/10,13/10},
s in {3/4,1/3},
e in {1/300,1/600}
```

converged with first-order halving estimates inside the frozen interval
`[0.8,1.2]`.  Deleting the boundary contribution leaves the certified nonzero
defect

```text
-3600 v^2/[(v^2+3)sqrt(3v^2+8)],
```

and shifting the selected mass by `1/10` leaves exactly `-4*pi/5`.

The first adversarial run returned `7/9` because the code used structural
equality on `-8*pi*(v^2+4)/(v^2+4)`.  Exact simplification gave `-8*pi`; the
failure and correction were committed before rerunning.  No equation, point,
tolerance or outcome rule changed.

## Reconciliation with the turning-point no-go

There is no contradiction with the earlier absence of a same-state analytic
half-step at `v=0`.  The two calculations use different singular scalings:

```text
generic velocity:    log L=O(e),   p_pre=O(1),
turning point:        log L=O(e^2), p_pre=O(e).
```

The generic tangent restores the expected leading reparametrization.  The
turning-point expansion probes the next nonzero order immediately and finds
that the old scaled-lapse family changes the incoming momentum.  Whether the
generic family also develops a finite-step obstruction is therefore a
next-order question, not answered by the leading theorem.

## Physical status ledger

| Claim | Status |
|---|---|
| Generic-velocity leading map is independent of interval factor | **DERIVED EXACT / ADVERSARIALLY CORROBORATED** |
| Leading one-versus-two tangent composition | **DERIVED EXACT / STRUCTURAL** |
| Unique positive mass branch for each real velocity | **DERIVED EXACT** |
| Next-order same-state composition | **OPEN** |
| Nonzero isolated relational duration `tau/L` | **NOT DERIVED** |
| Absolute classical tick | **DERIVED NEGATIVE under global scale covariance** |
| Seconds or Planck time from dimensionless geometry alone | **IMPOSSIBLE under the present hypotheses** |
| External novelty of the explicit formulas | **OPEN** |

## Framing verdict

This is a consistency result for the discrete dynamics, not the derivation of
a physical clock.  At leading order the current action says that the same
tangent history may be parameterized by `s=1`, `s=1/2` or any other positive
factor.  Calling any one of those values the fundamental tick would add a
choice that the equations did not make.

The next load-bearing test must retain one more order in `e`, preserve the
same incoming canonical state and conserved mass, and compare one coarse slab
against two fine slabs including their intermediate momentum.  Only an
isolated nonzero solution for the dimensionless duration can begin to derive
a relational tick.  An absolute unit would still require a dimensionful
physical input or a separately derived dimensional-transmutation mechanism.

