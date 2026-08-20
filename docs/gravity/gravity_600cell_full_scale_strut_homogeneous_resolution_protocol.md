# Protocol: exact homogeneous weak-pole canonical line

Date: 2026-08-20  
Status: **preregistered, target-disclosed symbolic/computational test**

## Frozen provenance and scope

This protocol follows the prior-art gate and exploratory disclosure committed in
`16e28ec`. It tests the mathematical operator underlying the unique homogeneous
sector of the already frozen carrier/action intersection. It does not reopen any
nonhomogeneous sector and does not test the omitted pole equation or physical
interpretation.

Required inputs are pinned by SHA-256 in the verifier:

- prior-art gate and disclosure;
- closed homogeneous frustum action source/artifact (`16/16`,
  `HOMOTHETIC_FRUSTUM_ACTION_INVARIANT`);
- canonical-lapse source/artifact (`7/7`,
  `HOMOTHETIC_CANONICAL_LAPSE_SELECTED`);
- complete carrier source/artifact and exact symbolic carrier result;
- binary intersection artifact;
- primary multiprecision artifact (`17/17`, homogeneous `OPEN`);
- adversarial nonhomogeneous artifact (`7/7`).

No coefficient, candidate ratio, physical target or tolerance may be fitted.

## Independent closed-action construction

Using new SymPy symbols `L_minus,L_plus,rho`, reconstruct rather than import the
frozen real homogeneous gravitational action

```text
S = 360 (L_minus+L_plus) h [2 pi-5 acos(cosine)]
    +600 sqrt(3) (L_minus^2-L_plus^2) asinh(boost),

h       = sqrt(rho+(L_plus-L_minus)^2/4),
cosine  = ((L_plus-L_minus)^2+2 rho)
          /(2*((L_plus-L_minus)^2+3 rho)),
boost   = (L_plus-L_minus)
          /sqrt(8*((L_plus-L_minus)^2+3 rho)).
```

The dust term is independent of the lower spatial length and therefore does not
enter the old spatial momentum used below. The verifier must check the formula
text against the frozen action artifact and retain its action-invariance controls.

Define the old collective momentum, up to an irrelevant common sign and
normalization, by

```text
p_minus = (L_minus/2) * partial S/partial L_minus.
```

At fixed `L_minus`, define

```text
p_s = L_plus * partial p_minus/partial L_plus,
p_z = rho    * partial p_minus/partial rho,
lambda = L_plus/L_minus.
```

These are derivatives with respect to `s=log lambda` and `z=log rho`.

## Exact generator and convention control

Construct, without numerical null-vector input,

```text
sigma = -lambda*p_z,
c     =  p_s.
```

Here the carrier coordinate is `sigma=delta lambda`, while the homogeneous
action coordinate is `delta s=sigma/lambda`. SymPy must prove identically

```text
p_s*(sigma/lambda) + p_z*c = 0.
```

It must also prove the carrier length responses are precisely the differentials
of the homothetic geometry:

```text
pole:       c,
diagonal:   (L_minus^2*sigma-rho*c)
            /(lambda*L_minus^2-rho),
upper edge: 2*sigma/lambda.
```

Corruption control: replace the correct conversion by `delta s=sigma`. The
resulting old-momentum differential must simplify to

```text
p_s*p_z*(1-lambda),
```

up to the disclosed sign convention, and must be nonzero at the frozen
nonstatic background. A sign-flipped generator must also fail.

## Background nondegeneracy and independent numerical bridge

Evaluate `p_s,p_z` and `-lambda*p_z/p_s` independently from the closed action at
160 and 220 decimal digits using automatic differentiation, never a finite
difference. Require:

- `abs(p_s)>1e5` and `abs(p_z)>1e-1` at both levels;
- relative P160/P220 changes below `1e-140`;
- even/odd background states agree at their already frozen parity bound;
- the analytic ratio agrees with the independently stored canonical-lapse
  endpoint-Jacobian ratio inside its frozen derivative error budget;
- after undoing the disclosed column scalings, the mean P100 candidate agrees at
  absolute error below `1e-30`; both five-component within-group spreads are
  nonzero and below `1e-30`, and remain diagnostics rather than inputs;
- the even/odd normalized analytic-generator projector distance is below `1e-70`.

The numerical bridge identifies the previously observed vector; it is not the
proof of the exact line.

## Nullity closure

The exact generator gives the D vector

```text
(sigma/scale_column_norm) repeated 5 times,
(c/strut_column_norm)     repeated 5 times,
```

and the K common vector with the strut component repeated in both the carrier and
canonical-input blocks. The exact generating-function and stationary-diagonal
identities imply these vectors lie in the mathematical D and K kernels.

Require the frozen P160 exhaustive drop-one-minor records to certify rank at
least 9 for D and 14 for K in each parity. Therefore:

```text
rank(D)=9,  nullity(D)=1,
rank(K)=14, nullity(K)=1.
```

The two parity generator projectors must agree within the frozen parity error.

## Outcome hierarchy

1. `HOMOGENEOUS_LINE_CONTROL_FAILED`: provenance, action identity, branch,
   parity, corruption or background nondegeneracy fails.
2. `HOMOGENEOUS_LINE_SYMBOLIC_DISAGREEMENT`: an exact differential identity
   fails.
3. `HOMOGENEOUS_LINE_NUMERICAL_BRIDGE_DISAGREEMENT`: exact identities pass but
   the independent stored constructions do not identify the same line.
4. `HOMOGENEOUS_LINE_NULLITY_OPEN`: the line is derived but the frozen lower-rank
   certificates do not close uniqueness.
5. `HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE`: all checks pass.

Outcome 5 is **DERIVED STRUCTURAL/COMPUTATIONAL** and still requires a mechanically
different adversarial verifier before consolidation. It means exactly one weak-pole
canonical response, not gauge, propagation, a solution tangent, a clock or a tick.

Only the new targeted verifier and static registry checks may run. No full suite.
