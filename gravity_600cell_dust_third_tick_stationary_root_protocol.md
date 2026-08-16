# Preregistration: target-independent stationary roots for tick three

Date: 2026-08-16

Prior-art gate: `7b9a676`.

Status: frozen before evaluating the third-slab stationary equation and
before parsing or comparing the second tick's post-momentum.

## 1. Fixed carrier and target firewall

Use both derived order-24 staircase schedules and the complete 100-decimal
Lorentzian Regge plus fixed-mass dust action.  Hold fixed

```text
A1 = -3.11605957669450169173470644419863944122165192557277135128791e-6,
B2 = -9.34818705890582713633822299265753373027428194008991504419612e-6,
R2 = -1.42370275520098029961300545242474815338378370661665379256974e-5.
```

Hash the accepted second-tick artifact and require SHA-256

```text
936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70.
```

Do **not** JSON-parse that artifact, do not load any `post_momentum`, target,
junction residual or previous failed solver endpoint, and serialize
`target_parsed=false`.

For absolute upper log `C=log(L3/L0)` at the inherited lapse `R2`, use

```text
q_old    = exp(2*B2)*L0^2,
q_new    = exp(2*C)*L0^2,
rho      = exp(R2)*rho0,
diagonal = exp(B2+C)*L0^2-rho.
```

Define `G(C,R2)` as the mean of the five pole equations and `P(C,R2)` as
the mean of the thirty pre-momenta.  No momentum value affects root discovery.

## 2. Frozen main grid and diagnostic sentinels

Use exactly 257 main-grid nodes

```text
x=C/abs(A1)=-8+k/16,  k=0,...,256.
```

Also evaluate the 18 diagnostic points

```text
x=+-2^n,  n=4,...,12.
```

Sentinels may diagnose outer signs and branch survival but do not extend the
finite root-exhaustion claim.

## 3. Deterministic candidate extraction

Mark every main-grid node satisfying

```text
abs(G)<1e-25
```

and merge adjacent marked nodes into one node cluster.  Independently record
every adjacent pair with strict sign change.  If a sign bracket touches a node
cluster, merge them as one candidate; otherwise retain both.

For every strict sign bracket perform exactly 80 bisections.  Retain the
opposite-sign half and use the final midpoint.  A node cluster representative
is the marked grid node having minimal `abs(G)`, with lower grid index as the
tie-breaker.  Sort final roots by increasing `C`.

Record before any later target comparison:

- exact candidate count `N`;
- source node/bracket indices;
- `(C,x,G,P)` and complete 30-component pre/post momenta;
- the complete 35 residuals and Lorentzian branch diagnostics.

## 4. Full and derivative gates

At each reported root require

```text
max abs(30 diagonal residuals) < 1e-60,
max abs(5 pole residuals)      < 1e-25,
max abs(all 35 residuals)      < 1e-25,
diagonal spread                < 1e-60,
pole spread                    < 1e-60.
```

Independently recompute the `(G,P)` Jacobian with respect to absolute `(C,R)`
using

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Use the already frozen entrywise error test.  Require `G_C` and

```text
D=G_C*P_R-G_R*P_C
```

to exceed `100` times their calibrated error envelopes.  Repeat the complete
root and derivative evaluation for the odd schedule and require

```text
|C_even-C_odd| < 1e-27,
max abs(local_even-local_odd) < 1e-24,
max abs(pre_even-pre_odd)     < 1e-22,
all determinant signs equal.
```

## 5. Predictions recorded before evaluation

- **STRUCTURAL:** exact time reversal of tick two predicts one stationary root
  at `C=A1`, i.e. `x=-1`.
- **PATTERN:** continuation of the two-tick sequence predicts a contracting
  stationary root near `C=6*A1`, i.e. `x=-6`.
- **PATTERN:** the predicted candidate count on the main grid is `N=2`.

These are diagnostics, not discovery or acceptance rules.  A different count
or position must be reported without changing the grid.

## 6. Mechanical outcomes and scope

Assign exactly one:

1. `THIRD_TICK_ROOT_CONTROL_FAILED`;
2. `THIRD_TICK_ROOT_BRANCH_FAILED`;
3. `THIRD_TICK_ROOT_ENUMERATION_OPEN` if any candidate fails root/full/
   derivative/parity gates;
4. `THIRD_TICK_STATIONARY_ROOTS_ENUMERATED` otherwise.

Outcome 4 exhausts sign-changing and node roots on the frozen grid only.
Tangential even-multiplicity roots between nodes and roots outside
`x in [-8,8]` remain **OPEN**.

No target hit or physical branch is selected here.  The artifact and result
note must be committed before a later verifier may parse `p_post,2`.

Only the new targeted verifier will be run.  The full suite will not be run.
