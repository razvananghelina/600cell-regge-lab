# Scale-free asymptotic map and fifth-slab result

Date: 2026-08-22.

## Verdict

**DERIVED:** after the exact normalization

```text
x=m q,
y=(p_infinity-pi)/m^2,
p_infinity=60*pi-300*sqrt(3)*log(2),
```

the compactified `m -> 0+` boundary of the fixed homogeneous 600-cell dust
map is a continuous one-parameter family of fixed points.  It does not select
a universal fixed point.

**DERIVED COMPUTATIONAL, FIVE-SLAB SCOPED:** the previously accepted unique
four-slab history has three real algebraic fifth-slab roots and exactly one
physical continuation.  A primary equal-`p` census and an adversarial
equal-`mu` census agree beyond 55 digits; the latter nests at 110- and
180-digit target precision.

**HISTORICAL OPEN, NOW RESOLVED FOR THE REPRESENTATIVE BRANCH:** this note did
not prove an invariant region or infinite history.  The later rigorous
half-strip theorem now proves a unique successor at every finite step for the
accepted branch-B seed.  Convergence to a history-independent member of the
limiting family and external novelty remain **OPEN**.

**DERIVED NEGATIVE FOR THE STRONGER FRAMING:** the limit does not derive an
absolute tick, a universal scale ratio, `c`, `G`, Planck units or particle
masses.  The original hope that compactification might select one universal
asymptotic point is false for this map: every physical `x>60` is fixed at
leading order.

## Complete hypotheses

The result uses only:

- the fixed homogeneous tetrahedral-frustum 600-cell Regge action;
- zero cosmological constant;
- conserved global dust;
- the committed pre/post canonical-momentum convention;
- positive incoming mass parameter;
- positive proper slab height and positive endpoint scale;
- the unique accepted history through slab four for the finite forecast.

No measured constant, continuum target, fitted coefficient, adjustable branch
rule or finite plotting box enters the acceptance test.

## Provenance and look-ahead control

The normalized variables were motivated after observing the already accepted
`q2,q3,q4` sequence.  Therefore their discovery is not blind and is not
presented as such.  The first genuinely out-of-sample datum is the fifth slab.

The public commit order is:

```text
c39f9ca  preregister exact map, coefficients and fifth-slab forecast
251aaee  register the primary verifier before its first execution
46014f5  record the primary fifth-slab artifact
67d1f3f  register the adversarial verifier before its first execution
38c908c  record the adversarial artifact
```

The frozen fifth forecast preceded both fifth-root censuses:

```text
x5-x4 = 0.0111863718238540698850889120790053991...,
x5    = 125.342979632764278093271017319287842440...,
q5    = 1006.534925564110861553259465876776082....
```

Commit ordering guards against changing that forecast after seeing the root.
It does not make the post-hoc choice of `x,y` blind and it does not replace an
independent human replication.

## Exact map

Define

```text
epsilon(q)=2*pi-5*acos[(q^2+2)/(2(q^2+3))],

mu(q)=180*epsilon(q)/(pi*sqrt(q^2+4)),

p(q)=180*q*epsilon(q)/sqrt(q^2+4)
     -600*sqrt(3)*asinh[q/sqrt(8(q^2+3))].
```

For

```text
U(m,x)=mu(x/m)/m,
V(m,x)=[p(x/m)-p_infinity]/m^2,
```

the exact root equation and update are

```text
4*pi*(U-1)+x*(V+y)=0,
r=2/U-1,
y_plus=-r*((r+1)*V+y).
```

Direct algebra gives the exact drift identity

```text
y_plus-y = 4*(U-1)/U^2 * (V+2*pi*U/x).
```

Two mechanically different calculations give

```text
U = 60/x
    +(-120-300*sqrt(3)/pi)*m^2/x^3
    +O(m^4),

V = -120*pi/x^2
    +(360*pi+900*sqrt(3))*m^2/x^4
    +O(m^4).
```

The primary used a symbolic `t=1/q` series.  The adversarial verifier instead
differentiated expressions regular at `t=0` and divided by the relevant
factorials; it did not call the primary series construction.

With

```text
A=120*pi+300*sqrt(3),
f(x)=4*pi*(x-30)/x^2,
```

the root curve and drift are

```text
y=f(x)+A*m^2/x^4+O(m^4),

y_plus-y=A*(60-x)*m^2/(900*x^3)+O(m^4).
```

At `m=0`, the physical expanding branch is

```text
x>60,
0<y<pi/30,
h/m=(x-60)/(30*x)>0,
r=(x-30)/30>1.
```

For each `y` in that interval, the limiting quadratic has an inner root in
`(30,60)` with negative height and one physical outer root in `(60,infinity)`.
On every outer root, `y_plus=y`.  This proves a continuous fixed family and
refutes selection of a single limiting point.

Implicit coefficient matching gives

```text
x_plus-x =
[(2*pi+5*sqrt(3))/(60*pi)]*(1-r0^(-2))*m^2+O(m^4),

r0=(x-30)/30.
```

Thus at every fixed `x>60`, sufficiently small positive `m` gives `x_plus>x`
and `y_plus<y`.  This is a local asymptotic result, not an induction theorem.

## Out-of-sample fifth slab

The primary result is

```text
all real roots     3
physical roots     1

q5                 1006.53493784425818414891995223361831716...
h5                 0.002163977529932147003985768834462395102...
L5/L4              3.17811898858662493889113408589035492510...
x5                 125.342981162001085132821642505768069828...
```

The other two roots are rejected by the preregistered inequalities:

```text
q5=-0.2295717047360...   L5/L4=-0.9390846881931...  negative endpoint
q5=316.6988622583965...  h5=-0.0021639775284...     negative height
```

The forecast error is

```text
x5-x5_forecast = 0.000001529236807039550625186480227387879...,
error/m3^4      = 0.000062347025773433066192126058734500857....
```

The preregistered signs also hold:

```text
x5-x4   = +0.0111879010606611094357140985592327870...,
y4-y3   = -0.00000517803715601510895242253121276261...,
y5-y4   = -0.000000512630025996080418460201870036236....
```

The numerical agreement is supporting evidence for the already derived
coefficient, not a fitted theorem.

## Independent root census

The primary verifier partitions the complete real line at equal-`p`
stationary points of the un-divided root equation.  The adversarial verifier
uses instead

```text
R(q)=p(q)-pi+4*pi*(mu(q)-m)/q,
R'(q)=4*pi*(m-mu(q))/q^2,
```

so its stationary points are the independently enumerated solutions of
`mu(q)=m`.  It separately excludes `q=0`, reconstructs four slabs by
redifferentiating the complete action, and validates the physical fifth root
against that action.  Wrong momentum rescaling, reversed momentum sign, dust
mass reset and boost omission all change the state or coefficients and fail.

The primary returned `11/11`; the adversarial verifier returned `12/12`.
Only these targeted verifiers were run.  No full-suite result is claimed.

Accepted artifacts:

```text
gravity_600cell_finite_height_asymptotic_map.json
a93837d2bbec340ddbac528c0be4da52aefe45c8f0d4310496eb1aef6a7b19b6

gravity_600cell_finite_height_asymptotic_map_adversarial.json
5215b2f07140be44f9e864b2688afa5e8e522b310a33ee5f7efa6cfccebc7405
```

## Status ledger

| Claim | Status |
|---|---|
| Exact normalized root equation and drift identity | **DERIVED** |
| Asymptotic coefficients through `O(m^2)` | **DERIVED** by two routes |
| Continuous physical fixed family at `m=0` | **DERIVED** |
| One universal asymptotic fixed point | **DERIVED NEGATIVE** |
| Signs of the first finite-`m` drift at fixed `x>60` | **DERIVED local asymptotic** |
| Unique physical continuation through slab five | **DERIVED COMPUTATIONAL, FIVE-SLAB SCOPED** |
| Forecast agreement | **PATTERN / out-of-sample numerical control** |
| Invariant physical region for every later slab | **DERIVED in the later invariant-region theorem** |
| Infinite deterministic history | **DERIVED for the representative homogeneous branch; generic `v` OPEN** |
| History-independent limiting member of the family | **OPEN** |
| Fundamental or absolute tick | **DERIVED NEGATIVE from scale covariance / not selected here** |
| External novelty of the exact coefficient theorem | **OPEN** |

## Next falsifiable gate

The invariant-region gate proposed here was later passed by two mechanically
distinct rigorous certificates.  The next gate is the global basin over the
original incoming state `v`, followed by nonhomogeneous perturbations.  The
induction establishes a dimensionless relational history, not an absolute
time unit.
