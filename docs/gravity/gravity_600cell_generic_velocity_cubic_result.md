# Result: cubic fixed-state obstruction, not a derived tick

Date: 2026-08-21

## Provenance

```text
prior-art gate                                      70e7ca2
frozen primary protocol                            d2efdf4
registered primary before first execution          628b769
preserved first exact OPEN census                   fcd5af5
frozen global-classification proof                 3de74e7
first classification implementation                8c29f61
preserved structural-equality failure              ea1023d
algebraic equality correction                       50bbc56
accepted primary artifact                          08bdde5
frozen adversarial protocol                        4523385
registered adversarial before first execution      ab63090
accepted adversarial artifact                      150a406
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_generic_velocity_cubic.py
  9/9 PASS

reproducible/verify_gravity_600cell_generic_velocity_cubic_adversarial.py
  10/10 PASS
```

Accepted artifacts:

```text
reproducible/gravity_600cell_generic_velocity_cubic.json
SHA-256 1d35b46cd4db20df0af3ed3e6b5de676d69753cf5059e0eb607d1eec949b9103

reproducible/gravity_600cell_generic_velocity_cubic_adversarial.json
SHA-256 b5167d597a927f8b441a096c31034aa04efa435883284dc2d9bfbd3b9cb3ff0d
```

The first complete `OPEN` artifact and the subsequent `8/9` structural-
equality failure are preserved in git history.  No full suite was run.

## Headline

```text
GENERIC_CUBIC_FIXED_STATE_OBSTRUCTION_ADVERSARIALLY_CORROBORATED
```

> **DERIVED NEGATIVE, scoped and adversarially corroborated:** for the fixed
> homogeneous cellular 600-cell Regge-plus-conserved-dust action, no real
> nonzero generic incoming velocity admits a same-state `C^3` endpoint jet.
> The lapse and incoming-momentum equations require incompatible cubic
> coefficients.  Their cross-resultant is strictly positive on the complete
> registered domain.

This kills the arbitrary-small-duration branch for this finite action.  It
does not select an isolated positive duration and therefore does not derive a
tick.

## Complete hypotheses

Use the certified positive-Lorentzian homogeneous cellular slab at zero
cosmological constant, with conserved dust and

```text
L_minus=1,
M=mu(v),
p0=p(v),
rho=h^2,
h>0,
```

where `mu(v)` and `p(v)` are the accepted leading state.  On

```text
v real,
v!=0,
K(v^2)!=0,
```

freeze the unique accepted lower-order coefficient

```text
a(v)=-B(v^2)/K(v^2)
```

and ask for a three-times differentiable endpoint history with jet

```text
L_plus=exp(vh+a(v)h^2+c h^3+o(h^3)).
```

The theorem does not cover a changed mass, changed incoming momentum,
nonhomogeneous modes, a refined/perfect action, the separately classified
turning point `v=0`, or the already-obstructed pair `K=0`.

## Exact obstruction

Define

```text
x=v^2,
r=sqrt(x+4),
q=sqrt(3x+8),
theta(x)=acos((x+2)/(2(x+3))),
epsilon(x)=2*pi-5*theta(x),

K(x)=10*r-(x+3)q epsilon(x),
B(x)=5*x*r+2(x+3)q epsilon(x).
```

After the lower orders vanish, write the next lapse and momentum
coefficients as affine functions of `c`:

```text
C2(v,c)=C_c(v)c+C_0(v),
P2(v,c)=P_c(v)c+P_0(v).
```

Both implementations independently prove the recursive slope identities

```text
C_c(v)=coefficient_a(C1),
P_c(v)=coefficient_a(P1).
```

In particular,

```text
C_c(v)=1440 v K(v^2)/[
  sqrt(v^2+4)sqrt(3v^2+8)(v^2+3)(v^2+4)
].
```

It is nonzero on the complete registered domain.  Therefore the lapse
equation has exactly one candidate

```text
c_lapse=-C_0/C_c.
```

The exact cross-resultant is

```text
Delta(v)=C_0(v)P_c(v)-P_0(v)C_c(v)
        =129600 epsilon(v^2)^2/(v^2+4).
```

At the lapse root the remaining momentum residual is

```text
P2(v,c_lapse)=-Delta(v)/C_c(v),
```

Thus a common `c` exists if and only if `Delta(v)=0`.

## Exact global positivity

For `x>=0`, set

```text
z(x)=(x+2)/(2(x+3)).
```

The adversarial proof uses the direct exact inequalities

```text
z(x)-1/3=x/[6(x+3)]>=0,
1/3-cos(2*pi/5)=(7-3sqrt(5))/12>0.
```

The final positivity follows from `7^2>9*5`.  Since the real principal
`acos` is strictly decreasing,

```text
theta(x)<2*pi/5,
epsilon(x)>0.
```

Therefore

```text
Delta(v)>0
```

for every real `v`, hence on every registered component separated by `v=0`
and the already-excluded roots of `K`.  There is no exceptional cubic velocity and no
degree-drop rescue.  At `K=0` the quadratic jet already fails; at `v=0` the
generic tangent expansion is not the applicable stratum.

## Mechanically independent replication

The primary route used the exact scaled action in five independent variables
`(lm,lp,q,w,tau)` and its path Hessian.  The adversarial route instead:

1. differentiated the complete unscaled action in `rho` and `L_minus` first;
2. only then substituted `L_plus=1+tau*q`, `rho=tau^2`;
3. reduced to the two-variable functions `(tau,q)`;
4. obtained the second path coefficient directly from

   ```text
   (1/2)f_tt+q1 f_tq+(1/2)q1^2 f_qq+q2 f_q;
   ```

5. formed `Delta` before reading any primary coefficient;
6. compared the two complete expressions only afterward.

It reproduced both full coefficients and `Delta` exactly.  Its disjoint
110-decimal controls used

```text
v in {-7/5,2/3,5/2},
c in {-2/7,1/5},
h in {1/1200,1/2400,1/4800,1/9600},
```

and every direct full-action quotient converged at the preregistered first
order.  Shifting the mass by `1/13` gave exactly `-8*pi/13`; deleting the
exponential contribution to the cubic quotient changed it by
`va+v^3/6`; the lapse root left an exactly nonzero momentum residual.

## What the result does and does not mean

| Claim | Status |
|---|---|
| Same-state generic endpoint jet exists through quadratic order | **DERIVED EXACT / ADVERSARIALLY CORROBORATED** |
| Same-state generic `C^3` endpoint jet exists for the fixed action | **DERIVED NEGATIVE / ADVERSARIALLY CORROBORATED** |
| Any exceptional generic velocity restores the cubic jet | **DERIVED NEGATIVE** |
| The cubic mismatch is a local pseudo-constraint obstruction | **STRUCTURAL, strongly supported by exact form and prior art** |
| The mismatch selects a nonzero interval | **DERIVED NEGATIVE for this asymptotic gate** |
| An isolated finite-height solution exists away from `h=0` | **OPEN** |
| The obstruction survives action improvement or carrier refinement | **OPEN** |
| Absolute classical tick | **DERIVED NEGATIVE under global scale covariance** |
| Fundamental relational tick | **NOT DERIVED** |
| External novelty of the explicit 600-cell coefficient | **OPEN** |

The result is stronger than merely failing to find a branch numerically: it
proves that a smooth same-state branch cannot pass through `h=0` within the
registered homogeneous finite action.  It is narrower than a no-dynamics
theorem.  Isolated finite roots, nonanalytic histories, extra degrees of
freedom and improved actions remain outside its hypotheses.

## Post-result prior-art audit

The learned term is precisely a higher-order consistency condition or
pseudo-constraint:

- [Dittrich--Hoehn, arXiv:0912.1817](https://arxiv.org/abs/0912.1817)
  state that higher-order Regge dynamics break the linearized symmetry and
  generate consistency conditions on lower-order background gauge
  parameters; their quadratic constraints depend on those parameters.
- [Bahr--Dittrich, arXiv:0905.1670](https://arxiv.org/abs/0905.1670)
  identify broken diffeomorphism symmetry in curved Regge calculus with
  pseudo-constraints rather than exact gauge constraints.
- [Bahr--Dittrich, arXiv:0907.4323](https://arxiv.org/abs/0907.4323)
  develop improved and perfect actions as the route for recovering continuum
  dynamics and gauge symmetry on a discrete carrier.
- [Vermeeren, arXiv:1505.05411](https://arxiv.org/abs/1505.05411)
  constructs modified Lagrangians for variational integrators, providing the
  relevant comparison framework for formal interpolating dynamics.

The first two sources predict the qualitative mechanism.  They do not print
the present positive `Delta(v)`, use this conserved-dust 600-cell action or
prove the present global no-jet theorem.  The search therefore prevents a
claim that the mechanism is new; novelty of the explicit coefficient remains
**OPEN**.

Static registry audit after adding both targeted verifiers: `402/402`
distinct registrations, two deliberate exclusions, zero duplicates, zero
unregistered files and zero missing files.  This was a registry audit, not a
full-suite execution.

## Next load-bearing gate

The direct tick question is now an exact finite-height problem.  For each
fixed incoming state `(M,p0)=(mu(v),p(v))`, solve the two unexpanded equations

```text
F(1,L_plus,h^2;mu(v))=0,
p_pre(1,L_plus,h^2;mu(v))=p(v)
```

for the two unknowns `(L_plus,h)` with `h>0`, without imposing the failed
Taylor endpoint.  A complete bounded root census must be preregistered before
evaluation.

- No positive roots would close the finite homogeneous tick route for this
  action.
- Isolated positive roots would be candidates only; composition, stability,
  refinement and scale selection would still be required before calling one
  a physical tick.
- A family of roots would retain reparametrization rather than select time.

Because global scale covariance remains exact, even a dimensionless
`h/L_minus` root cannot by itself produce seconds or Planck time.
