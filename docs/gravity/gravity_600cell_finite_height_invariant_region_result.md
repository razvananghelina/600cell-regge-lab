# Invariant half-strip and infinite homogeneous history

Date: 2026-08-22.

## Verdict

**DERIVED COMPUTATIONAL WITH TWO MECHANICALLY DISTINCT RIGOROUS INTERVAL
CERTIFICATES:** for the fixed homogeneous 600-cell dust canonical relation,
every state in

```text
D={(m,x): 0<m<=2/5, x>=125}
```

has exactly one physical successor, and that successor satisfies

```text
0<m_plus<m,
x_plus>x>=125.
```

Thus `D` maps strictly into itself.

**DERIVED BY INDUCTION, REPRESENTATIVE-SEED SCOPED:** the accepted branch-B
state

```text
m3=0.3957443748524788013...,
x4=125.3317932609404240...
```

lies inside `D`.  It therefore has one physical successor at every later
finite step.  Together with the earlier result that branch A has no third
slab, exactly one of the two second-slab branches at the frozen incoming state
`v=3/2` admits a complete forward history in this homogeneous model.

**STRUCTURAL:** using complete forward extendibility to prefer branch B is a
global admissibility condition.  It is not a local equation of motion, was
motivated after the bifurcation was known, and does not establish a physical
arrow of time.

**OPEN:** the basin in the original incoming state `v`, nonhomogeneous
stability, convergence to a member of the limiting fixed family, total proper
duration, continuum general relativity and external novelty.

**NOT DERIVED:** an absolute tick, `c`, `G`, Planck units, particles or local
gravitational degrees of freedom.

## Complete hypotheses

The theorem uses only:

- the fixed homogeneous tetrahedral-frustum 600-cell action;
- zero cosmological constant;
- conserved global dust;
- the committed pre/post canonical-momentum convention;
- positive proper slab height and positive endpoint scale;
- the exact all-real finite-height root classification;
- the accepted branch-B history through slab five.

No measured constant, continuum target, fitted coefficient, alternate branch
rule or finite root box enters the domain proof.

The rational thresholds were chosen after the accepted fourth state was
known.  The region is therefore not blind and is not dynamically selected.
The nonacceptance diagnostic `x>=124` also passes, so `125` is demonstrably
not a sharp physical boundary.  This does not invalidate invariance, but it
prevents interpreting the particular half-strip as a selected object.

## Exact normalized map

Set

```text
z=1/x,
u=m^2*z^2,

epsilon(u)=2*pi-5*acos[(1+2*u)/(2*(1+3*u))],

M(u)=180*epsilon(u)/(pi*sqrt(1+4*u)),

P(u)=180*epsilon(u)/sqrt(1+4*u)
     -600*sqrt(3)*asinh[1/sqrt(8*(1+3*u))],

P0=60*pi-300*sqrt(3)*log(2),
W(u)=[P(u)-P0]/u,
U=z*M,
V=z^2*W,
Y=-V-4*pi*z*(U-1).
```

The current state is on the exact root curve `y=Y(m,z)`.  Its outgoing data
are

```text
r=2/U-1,
m_plus=m/r,
y_plus=-r*((r+1)*V+Y).
```

With

```text
N=P-P0+2*pi*u*M,
Bbar=N/u^2,
C=W+4*pi*M,
```

the same-`x` gap has the exact continuous factorization

```text
[Y(m_plus,z)-y_plus]/(m^2*z^2)
 =4*(1-U)*Bbar/M^2
  -z^2*(1-r^(-2))*Cbar,
```

where `Cbar` is a mean value of `-C'` between the two mass arguments.  This
includes both compactification axes and is not inferred from a finite grid.

## Rigorous signs

The independent global Lagrange certificate encloses the complete interval

```text
0<=u<=4/390625
```

at both 160 and 256 decimal digits.  Representative strict lower endpoints
from the 256-digit record are

```text
M                         > 59.9936880469
Bbar                      > 896.554263117
C                         > 376.911798755
-C'                       > 896.035723639
1-U                       > 0.519972881258
y_plus/z                  > 9.54860504162
partial_z Y               > 6.53339044810
normalized same-x gap     > 0.460493701291
negative-root margin      > 166.647732695.
```

The bounds imply

```text
r>1,
0<m_plus<m,
0<y_plus<Y(m_plus,z),
Y(m_plus,0)=0,
partial_z Y(m_plus,z)>0.
```

Hence the next-root equation has one solution with `0<z_plus<z`, so
`x_plus>x`.

For complete real-root uniqueness use

```text
R(q)=p(q)-pi+4*pi*(mu(q)-m)/q,
R'(q)=4*pi*(m-mu(q))/q^2.
```

The already certified unimodality of `mu` leaves only the outer positive root
as a possible positive-height solution.  `q=0` is excluded by
`m<=2/5<mu(0)=30`.  On the negative axis a physical root would require
`m<mu<2m`; the rigorous tail margin above prevents its momentum correction
from cancelling `p(q)-pi_plus`.  Therefore the successor found inside `D` is
the only physical real root.

## Primary certificate and preserved failures

The primary proof used exact integral Taylor identities for `W`, `Bbar` and
`-C'`.  Its first run returned `14/14`, but Arb's compact pretty printer hid
some strictly positive lower endpoints.  That first artifact was preserved.
After a separately frozen reporting-only correction, the primary again
returned `14/14` with explicit endpoints.

The first adversarial method partitioned `u` into 64 rational leaves and used
direct quotient interval evaluation away from zero.  It returned `5/12` and

```text
INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN.
```

All 63 direct-quotient leaves lost the cancelling first-order dependency in
`Bbar` and `-C'` at both 160 and 256 digits.  Their bounds were essentially
unchanged with precision, so this was not roundoff and not a wrong-sign
witness.  Its 968 KiB artifact remains committed with SHA-256

```text
f7d1f36e5ed679c39d1c38dbc21509ae52211f6735b38a2da46046fb798f54d5.
```

The resolution route was frozen only after preserving that failure.  It used
one degree-six Maclaurin polynomial with global Lagrange derivative
remainders, no derivative integrals, no direct quotient interval evaluation,
no subdivision and no grid.  It returned `12/12` and

```text
INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED.
```

Its artifact SHA-256 is

```text
813e05bd66b47cc3ae1cd35d0a2eddb9c645a850d84abeaad37d15b14a6a380f.
```

The primary and resolution bounds overlap, their signs agree, and both
reproduce the fifth successor from a complete redifferentiation of the action.
Reversed momentum sign, boost omission, the naively divided zero quotient and
the `x=60` boundary all fail as required.

## Public provenance

```text
f1ae36a  complete the pre-result prior-art gate
559ca7a  preregister the invariant half-strip
8e30b1c  freeze the primary integral-remainder method
8dd1d60  register the primary verifier before its first run
7a40539  preserve the first primary artifact with lossy display
8f770cd  freeze the reporting-only correction
9d5b479  implement explicit endpoint reporting
f5d8d0f  record the corrected 14/14 primary certificate
26ef9c3  preregister the direct-quotient adversarial route
87d9aac  register it before its first run
bfa4db4  preserve its 5/12 OPEN artifact
a851018  freeze the global-Lagrange resolution protocol
aa1da9b  register the resolution before its first run
ec83f69  record the 12/12 adversarial corroboration.
```

Commit ordering proves that the domains, degrees, precisions, signs and
controls preceded their corresponding executions.  It does not make the
post-hoc half-strip blind and does not replace outside replication.

## Post-result literature audit

The search was repeated using the resulting terms `normalized canonical
map`, `invariant region`, `complete forward history`, `600-cell dust` and
`causal endpoint`.

Dittrich and Höhn, *Canonical simplicial gravity*,
[arXiv:1108.1974](https://arxiv.org/abs/1108.1974), explicitly explain that
data free at one discrete move can become fixed a posteriori by pre-constraints
from subsequent moves.  Therefore dependence of admissibility on later steps
is a **KNOWN structural feature** of canonical simplicial dynamics, not by
itself proof of a malformed boundary problem.  It still does not turn the
present global selector into a local physical law.

De Felice and Fabri,
[arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093), report a causal
stopping point in a Sorkin evolution of a dust 600-cell.  Their carrier
evolution and equations differ from the frozen frustum relation here, so the
results are a warning/control rather than a contradiction.

Jercher and Steinhaus,
[arXiv:2312.11639](https://arxiv.org/abs/2312.11639), derive matter-dependent
existence inequalities and monotone branches for a spatially flat cuboidal
frustum model with a massless scalar.  This is structurally close but not the
same cells, matter or normalized map.

No checked primary source states the present exact half-strip or its
coefficient bounds.  Search failure is not proof of novelty; external novelty
remains **OPEN** pending expert review.

## What the theorem does and does not repair

It repairs the finite-horizon weakness: the statement is no longer that one
branch survived five slabs.  The accepted branch has a unique continuation at
every later finite step.

It does not repair the representative-state weakness.  `v=3/2` was chosen as
a convenient preregistered incoming state, not derived.  Nor does it introduce
local degrees of freedom: the model remains homogeneous and contains only a
global scale mode.

Accordingly the result is a real theorem about this reduced canonical
relation, not yet local gravitational physics and not a theory of nature.

## Status ledger

| Claim | Status |
|---|---|
| Exact invariant half-strip `D` | **DERIVED COMPUTATIONAL / RIGOROUS** |
| Unique physical successor for every point of `D` | **DERIVED** |
| Accepted branch-B history continues at every finite step | **DERIVED BY INDUCTION, REPRESENTATIVE-SCOPED** |
| Branch B is the only complete branch among the two at `v=3/2` | **DERIVED / STRUCTURAL** |
| The threshold `x=125` is sharp or selected | **DERIVED NEGATIVE / post-hoc** |
| Infinite extendibility is a local evolution law | **OPEN / not established** |
| Basin over the original incoming state `v` | **OPEN** |
| Limit point and total proper duration | **OPEN** |
| Nonhomogeneous stability and local gravity | **OPEN** |
| Absolute tick or limiting speed | **NOT DERIVED** |
| External novelty | **OPEN** |

## Next falsifiable gates

1. Classify the number of physical branches and entry into the invariant
   region as a function of the original incoming state `v`.  This attacks the
   arbitrariness of `v=3/2`; another isolated state does not suffice.
2. Then break homogeneity and evolve preregistered perturbation modes.  Only
   their stable linear propagation can test local gravitational physics and
   an effective space-time wave operator.

Only the targeted invariant-region verifiers and documentation/static guards
are claimed at this checkpoint.  No full-suite result is claimed.
