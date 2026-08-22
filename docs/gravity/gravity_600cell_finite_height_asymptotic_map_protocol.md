# Protocol: scale-free asymptotic map and fifth-slab forecast

Date: 2026-08-22.

Status: frozen before any fifth-slab solve or artifact.

Accepted input:

```text
fourth-slab consolidation commit  2dca534
public-tree cleanup commit        6333518
```

## 1. Complete hypotheses

Use only:

- the fixed homogeneous tetrahedral-frustum 600-cell Regge action;
- zero cosmological constant;
- conserved global dust;
- the committed pre/post momentum convention;
- positive incoming mass parameter `m`;
- positive proper slab height;
- positive endpoint scale;
- the unique accepted history through slab four when making the
  out-of-sample forecast.

No measured constant, continuum target, fitted coefficient, adjustable
branch rule or finite plotting box is admissible.

## 2. Exact canonical relation

Define

```text
epsilon(q)=2*pi-5*acos[(q^2+2)/(2(q^2+3))],

mu(q)=180*epsilon(q)/(pi*sqrt(q^2+4)),

p(q)=180*q*epsilon(q)/sqrt(q^2+4)
     -600*sqrt(3)*asinh[q/sqrt(8(q^2+3))],

p_infinity=60*pi-300*sqrt(3)*log(2).
```

For an incoming state `(m,pi)`, every slab root obeys

```text
4*pi[mu(q)-m]+q[p(q)-pi]=0,

h=[p(q)-pi]/[2*pi*mu(q)],
r=1+h*q=2*m/mu(q)-1,

m_plus=m/r,
pi_plus=p(q)+[p(q)-pi]/r.
```

## 3. Frozen scale-free variables

For `m>0`, set

```text
x=m*q,
y=(p_infinity-pi)/m^2,

U(m,x)=mu(x/m)/m,
V(m,x)=[p(x/m)-p_infinity]/m^2.
```

The exact root equation and map must be derived as

```text
4*pi(U-1)+x(V+y)=0,

r=2/U-1,

y_plus=-r[(r+1)V+y].
```

The decisive exact drift identity to prove is

```text
y_plus-y =
4*(U-1)/U^2 * [V+2*pi*U/x].
```

A verifier that merely substitutes stored decimal roots does not pass this
gate.

## 4. Frozen asymptotic claims

Derive the `m->0+` series at fixed `x>0` symbolically:

```text
U = 60/x
    +(-120-300*sqrt(3)/pi)*m^2/x^3
    +O(m^4),

V = -120*pi/x^2
    +(360*pi+900*sqrt(3))*m^2/x^4
    +O(m^4).
```

Define

```text
A=120*pi+300*sqrt(3).
```

The preregistered consequences are:

```text
y = f(x)+A*m^2/x^4+O(m^4),
f(x)=4*pi*(x-30)/x^2,

y_plus-y =
A*(60-x)*m^2/(900*x^3)+O(m^4).
```

At the compactified boundary `m=0`, the physical expanding branch is

```text
x>60,
0<y<pi/30,

h/m=(x-60)/(30*x),
r=(x-30)/30>1.
```

For each `y in (0,pi/30)`, the limiting quadratic has two positive roots.
The root in `(30,60)` has negative height and the unique limiting physical
root lies in `(60,infinity)`.

The limiting map must satisfy `y_plus=y` for every such physical root. If
so, the asymptotic boundary contains a **continuous fixed family**, not a
selected universal fixed point.

## 5. Next-order physical-root drift

Let `x_plus` denote the outer physical root for the next incoming state.
Writing

```text
r0=(x-30)/30,
```

derive without fitting

```text
x_plus-x =
[(2*pi+5*sqrt(3))/(60*pi)]*(1-r0^(-2))*m^2
+O(m^4).
```

Thus for every fixed `x>60`, sufficiently small positive `m` predicts

```text
y_plus<y,
x_plus>x.
```

This is a local asymptotic statement. It does not by itself prove an
invariant region, an infinite history or a history-independent limit.

## 6. Out-of-sample fifth-slab forecast

Use the accepted fourth incoming state and physical fourth root only after
the symbolic formulas above are implemented. The frozen input is

```text
m3 =0.395744374852478801317435698526508699892262815172057749108055761...
q4 =316.698862258396521252257613203538772441458276208574268200419575...
r4 =3.177924571855581867540887867308593202015971928861904193867865...
```

Hence

```text
x4=m3*q4
  =125.331793260940424023385928407208837040957261007383495059962...,

r0,4=x4/30-1
    =3.1777264420313474674461976135736279013652420335794....
```

Before solving the fifth equations, freeze the next-order forecast

```text
x5-x4 approximately 0.0111863718238540698850889120790053991...,
x5     approximately 125.342979632764278093271017319287842440...,
q5     approximately 1006.534925564110861553259465876776082....
```

The approximation has no fitted tolerance. Report its raw error and the
error divided by `m3^4`; do not convert agreement into a theorem.

Independently enumerate every real fifth-slab root on the complete real line.
Classify physical roots only by `h5>0` and `r5>0`.

## 7. Required controls and replication

Primary route:

- obtain the series by a symbolic `t=1/q` expansion;
- prove the exact drift identity algebraically;
- enumerate fifth roots through the equal-`p` stationary partition.

Adversarial route:

- redifferentiate the complete action independently;
- obtain asymptotic coefficients from derivatives at `t=0`, not by calling
  the primary series;
- enumerate fifth roots through equal-`mu` points and
  `R'(q)=4*pi[m-mu(q)]/q^2`;
- read the primary artifact only after its own root census is complete.

Hostile controls must change at least one derived coefficient or fifth
incoming state: wrong momentum rescaling, reset dust mass, reversed momentum
sign and omission of the boost term.

## 8. Frozen outcomes

### CONTINUOUS_ASYMPTOTIC_FIXED_FAMILY_AND_UNIQUE_FIFTH_SLAB

All exact identities and coefficients pass; the fifth slab has exactly one
physical root; its drift has the preregistered signs. Label the symbolic
family **DERIVED**, the fifth continuation **DERIVED COMPUTATIONAL,
five-slab scoped**, and any infinite-history reading **OPEN**.

### ASYMPTOTIC_FAMILY_DERIVED_BUT_HISTORY_STOPS

The symbolic family passes but the accepted history has no physical fifth
root. The finite history stops before slab five; asymptotic local existence
does not govern this state.

### FIFTH_SLAB_BRANCHES

The symbolic family passes but more than one physical fifth root exists.
Finite-horizon uniqueness is lost.

### ASYMPTOTIC_MAP_REFUTED

Any exact identity or frozen coefficient fails under either derivation.

### ASYMPTOTIC_MAP_OPEN

Incomplete global root census, unresolved branch convention, inadequate
precision nesting or disagreement between the two mechanical routes.

No outcome derives an absolute tick, `c`, `G`, Planck units or particle
masses.

