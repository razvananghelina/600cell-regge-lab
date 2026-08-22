# Protocol: invariant half-strip for the exact normalized canonical map

Date: 2026-08-22.

Prior-art gate commit: `f1ae36a`.

Status: frozen before evaluating any new invariant-region inequality.

## 1. Scope and complete hypotheses

Use only:

- the fixed homogeneous tetrahedral-frustum 600-cell Regge action;
- zero cosmological constant;
- conserved global dust;
- the already committed canonical pre/post momentum convention;
- positive proper slab height and positive endpoint scale;
- the exact finite-height root classification and the accepted history through
  slab five.

No measured constant, continuum target, fitted coefficient, alternate branch
rule or finite numerical trajectory may enter the proof.

The result, if positive, is a theorem for a homogeneous discrete cosmology.
It is not local gravity, continuum general relativity or a derivation of an
absolute time unit.

## 2. Provenance disclosure and frozen region

The thresholds below were chosen after the accepted history was known:

```text
D = {(m,x): 0<m<=2/5, x>=125}.
```

They are rational padding around the accepted fourth incoming pair

```text
m3=0.3957443748524788013...,
x4=125.3317932609404240....
```

Therefore membership of that history is not out-of-sample evidence and the
region is not claimed to be dynamically selected.

An earlier conversational draft proposed the bounded wedge

```text
125<=x<=126-m^2.
```

Before any invariant-region evaluation, exact inspection showed that an upper
boundary is artificial: proving `x_plus>x` makes the half-strip invariant and
is stronger.  This protocol replaces the unexecuted wedge by `D`.  If `D`
fails, the wedge or another smaller region will not be substituted in this
mission.  The failure is the result.

## 3. Exact regularized map

Set

```text
z=1/x,
t=m*z,
```

so the unbounded half-strip compactifies to

```text
0<m<=2/5,
0<z<=1/125.
```

Define expressions analytic at `t=0`:

```text
epsilon(t)=2*pi-5*acos[(1+2*t^2)/(2*(1+3*t^2))],

M(t)=180*epsilon(t)/(pi*sqrt(1+4*t^2)),

P(t)=180*epsilon(t)/sqrt(1+4*t^2)
     -600*sqrt(3)*asinh[1/sqrt(8*(1+3*t^2))],

p_infinity=P(0)=60*pi-300*sqrt(3)*log(2),

W(t)=[P(t)-P(0)]/t^2,
W(0)=-120*pi.
```

Then

```text
U(m,z)=z*M(m*z),
V(m,z)=z^2*W(m*z),

Y(m,z)=-V-4*pi*z*(U-1).
```

`Y` is the exact root curve: an incoming state with root `x=1/z`
has `y=Y(m,z)`.  Its exact outgoing state is

```text
r=2/U-1,
m_plus=m/r,
y_plus=-r*((r+1)*V+Y).
```

No decimal root may be substituted into these formulas in the domain proof.

## 4. Frozen invariant-region gates

For every `(m,z)` in the compactified domain, prove all of the following with
strict outward-rounded inequalities for `m,z>0` and continuous boundary
limits at `m=0` or `z=0`:

1. `0<U<1`, hence the represented current root has positive height,
   `r>1`, and `0<m_plus<m`.
2. `y_plus>0`.
3. `partial Y(m,z)/partial z>0` on the larger closed rectangle
   `0<=m<=2/5`, `0<=z<=1/125`.
4. The same-`x` bracket is strict:

   ```text
   Y(m_plus,z)-y_plus>0.
   ```

5. `Y(m_plus,0)=0<y_plus`.

Items 2--5 imply that the next root equation

```text
Y(m_plus,z_plus)=y_plus
```

has exactly one solution in `0<z_plus<z`, hence

```text
x_plus=1/z_plus>x>=125.
```

Together with `m_plus<m<=2/5`, this maps `D` strictly into itself.

The decisive same-`x` gap vanishes on the compactification boundary.  The
certificate must prove the regular quotient

```text
[Y(m_plus,z)-y_plus]/(m^2*z^2)>0
```

through its continuous extension; simply dropping the axes from a finite
grid is forbidden.

## 5. Complete physical-root uniqueness

An invariant root inside `D` is not enough if another physical root exists.
The verifier must also prove, without a finite root box:

```text
R(q)=p(q)-pi+4*pi*(mu(q)-m)/q,
R'(q)=4*pi*(m-mu(q))/q^2.
```

Use the already certified unimodality of `mu` and independently check the
sign consequences needed here:

- on `q>0`, `R` has at most two roots; at most its outer root can have
  positive height, and the invariant-region root is that outer root;
- `q=0` is excluded because `m<=2/5<mu(0)=30`;
- on `q<0`, any physical root would require `m<mu(q)<2m`;
- prove `p(q)>0` for `q<0`, `pi_plus<p_infinity<0`, and a strict bound showing
  that the negative correction `4*pi*(mu-m)/q` cannot cancel
  `p(q)-pi_plus` on that putative physical interval.

Only after these global facts pass may the successor be called the unique
physical root.

## 6. Certified computation

Primary route:

1. derive every algebraic identity in Sections 3--5 with SymPy;
2. represent `M`, `W`, `Y`, `y_plus` and the normalized gap by Taylor models
   about `t=0` with outward-rounded Arb remainder bounds;
3. use `|t|<=1/250`, which strictly contains the full domain
   `|m*z|<=1/312.5`;
4. use even Taylor terms through degree 12 and bound the remainder from the
   corresponding next derivative on the full rational interval;
5. certify the compact `(m,z)` rectangle by deterministic bisection at 192
   decimal digits, splitting the normalized widest coordinate first and the
   left coordinate on ties;
6. freeze maximum depth 28 and preserve the complete unresolved leaf list if
   that depth is reached;
7. certify the accepted `(m3,x4)` membership using its frozen artifact hash,
   but do not count membership as evidence for invariance.

Every interval endpoint must be rational and every transcendental enclosure
must come from Arb outward rounding.  Floating-point sampling, fitted
tolerances and agreement at finitely many trajectory points do not pass.

## 7. Controls

Known-pass controls:

- at `m=0`, recover

  ```text
  Y(0,z)=4*pi*z-120*pi*z^2,
  partial_z Y(0,z)=4*pi*(1-60*z)>0
  ```

  on `0<=z<=1/125`;
- reproduce the accepted fifth root from the current pair only after the
  domain theorem has been evaluated.

Known-fail controls:

- `x=60` must fail strict positive-height invariance at the compactified
  boundary;
- omission of the boost term must change at least one Taylor coefficient and
  the normalized same-`x` gap;
- wrong outgoing momentum rescaling and dust-mass reset must change the next
  incoming state.

## 8. Outcome hierarchy

### `INVARIANT_HALF_STRIP_CERTIFIED`

Every exact, interval and global-root gate passes.  Label the half-strip
**DERIVED COMPUTATIONAL WITH RIGOROUS INTERVAL CERTIFICATE** and infer by
induction that the accepted branch-B history has a unique physical successor
for every later finite step.  Label this an **infinite relational history in
the frozen homogeneous model**, not a fundamental time flow.

The use of infinite extendibility to prefer branch B remains **STRUCTURAL**:
it is a global-in-time selector, was motivated after the finite bifurcation,
and is not promoted to a local physical law.

### `CANDIDATE_HALF_STRIP_REFUTED`

A rigorously certified point or box violates any required inequality.  Record
the first lexicographic witness and stop.  Do not replace `D` by the earlier
wedge or tune either threshold.

### `INVARIANT_HALF_STRIP_OPEN`

Use this for an unresolved interval leaf, failed remainder enclosure,
incomplete global-root proof or resource limit.  Numerical scans cannot
upgrade `OPEN`.

No outcome establishes genericity over the original incoming parameter `v`,
nonhomogeneous stability, local gravitational degrees of freedom, continuum
GR, an absolute tick, `c`, `G` or Planck units.

## 9. Required adversarial gate

A positive primary result is not accepted until a separately registered
verifier:

- avoids the primary Taylor-model construction for the decisive gap;
- uses a different compactification or derivative-integral enclosure;
- rederives the complete action and outgoing state;
- attacks both coordinate axes, threshold perturbations and legitimate sign
  conventions;
- reads the primary artifact only after its own domain certificate is built.

Disagreement leaves the theorem **OPEN**.
