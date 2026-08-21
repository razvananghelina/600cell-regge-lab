# Preregistration: same-state one-versus-two weak-lapse composition

Date: 2026-08-21

Prior-art gate commit: `4f4b5c2`.

Status: frozen before constructing the coefficient equations, before solving
for any half-step branch, and before evaluating any coarse/fine discrepancy.

## 1. Frozen action and state

Use only the exact cellular action and conventions in
`gravity_600cell_relational_step_composition_prior_art.md`.  Put

```text
epsilon = 2*pi-5*acos(1/3),
M       = (90/pi)*epsilon,
L0      = 1,
x       = e^2,
p0(e)   = 180*epsilon*e.
```

Here `p0` is the canonical post-momentum of the exact static coarse slab.
For each fixed `e`, it is held literally identical in the coarse and fine
histories.  Across the asymptotic sequence `e -> 0`, this state approaches the
time-symmetric turning point.  Consequently this protocol tests local
composition on that one state family only; it is not a generic nonzero-
velocity semigroup theorem.

Define exact logarithmic derivatives

```text
F(Lm,Lp,rho) = rho * partial S/partial rho,
Pminus       = (Lm/2) * partial S/partial Lm,
Pplus        = (Lp/2) * partial S/partial Lp.
```

The canonical equations are `F=0`, `Pminus+p_in=0`, and at a shared
boundary `Pplus(previous)+Pminus(next)=0`.  The final physical canonical
momentum is `Pplus`.

## 2. Coarse branch

The one-slab coarse history uses

```text
Lc    = exp(Ac*x+Bc*x^2+O(x^3)),
rhoc  = e^2 exp(Rc*x+O(x^2)),
```

and solves

```text
F(1,Lc,rhoc)=0,
Pminus(1,Lc,rhoc)+p0(e)=0.
```

The exact coefficient branch must reproduce the already blind-derived
contracting `n=1` cellular coefficients.  This is a provenance control, not
a fitted target.

## 3. Fine branch enumeration

The two-slab fine history uses the same `(L0,p0,M)` and the ansatz

```text
L1   = exp(A1*x+B1*x^2+O(x^3)),
L2   = exp(A2*x+B2*x^2+O(x^3)),
rho1 = (e^2/4) exp(R1*x+O(x^2)),
rho2 = (e^2/4) exp(R2*x+O(x^2)).
```

It must solve all four equations

```text
F(1,L1,rho1)=0,
Pminus(1,L1,rho1)+p0(e)=0,
F(L1,L2,rho2)=0,
Pplus(1,L1,rho1)+Pminus(L1,L2,rho2)=0.
```

Construct the exact weak-lapse coefficient equations through the first order
that determines all six coefficients.  Enumerate every real solution in
`Q(sqrt(2),epsilon)` satisfying the stated half-lapse ansatz.  Record

```text
N_fine = exact number of admissible coefficient branches.
```

No numerical seed, sign or closeness to the coarse endpoint may discard a
branch.  If symbolic elimination cannot establish completeness, the result
is `OPEN`, not a uniqueness claim.

## 4. Exact leading composition gate

For every admissible fine branch derive, before nonlinear evaluation,

```text
scale coefficient:
  Cq = A2-Ac,

momentum coefficient:
  Cp = coefficient of e in [Pplus_fine-Pplus_coarse],

proper-time leading ratios:
  sqrt(rho1)/e -> 1/2,
  sqrt(rho2)/e -> 1/2.
```

A branch is `LEADING_COMPOSITION_COMPATIBLE` only if `Cq=Cp=0` exactly and
both proper-time limits are exactly `1/2`.  Report the hit fraction over all
`N_fine` branches.  One compatible branch among many is a look-elsewhere
pattern and cannot establish a selected relational map.

## 5. Full nonlinear controls

Only exact-leading-compatible branches proceed.  At

```text
e in {1/100,1/200,1/400},
```

solve the unexpanded action equations at 100 decimal digits, initialized by
the branch's disclosed series.  No restart or alternate root search is
allowed.  Require:

```text
max equation residual       < 1e-60,
positive L and rho,
finite real action/derivatives,
actual tau_i/e -> 1/2 on the registered branch.
```

Define normalized discrepancies

```text
Eq(e) = abs(log(L2/Lc))/e^2,
Ep(e) = abs(Pplus_fine-Pplus_coarse)/abs(p0(e)),
Et(e) = abs[(sqrt(rho1)+sqrt(rho2))-sqrt(rhoc)]/e.
```

If a discrepancy is resolved above `1e-70`, require strict decrease under
both halvings and observed halving orders at least `1.8`.  There is no upper
bound: exact composition or superconvergence is permitted.  If it remains
below `1e-70` at all three points, label it `NUMERICALLY_EXACT_AT_100_DPS`.

This is a second-order-or-better local consistency gate because the leading
scale, momentum and total-time sizes are respectively `O(e^2)`, `O(e)` and
`O(e)`.

## 6. Calibration and hostile control

Before the Regge comparison:

1. use the exact constant-force flow

   ```text
   q(h)=q0+h*p0+f*h^2/2,
   p(h)=p0+f*h
   ```

   and prove algebraically that one step `h` equals two steps `h/2` from the
   same `(q0,p0)`;
2. replace the first fine incoming momentum by its own static value
   `p0/2` and require the same-state gate to reject it exactly.

The second control is the old lambda-family mistake made mechanically
falsifiable.

## 7. Registered outcomes

Assign exactly one headline outcome.

### `SAME_STATE_HALF_STEP_BRANCH_ABSENT`

**DERIVED NEGATIVE, scoped** if complete exact elimination gives
`N_fine=0`.  The current homogeneous action then has no half-lapse branch of
the registered turning-point form from the same canonical state.

### `TEMPORAL_REFINEMENT_LEADING_MISMATCH`

**DERIVED NEGATIVE, scoped** if `N_fine>0` but no branch passes both exact
leading composition identities.

### `TEMPORAL_REFINEMENT_NONUNIQUE`

**STRUCTURAL / OPEN selection** if more than one exact-leading-compatible
branch survives.  Print all branches and nonlinear outcomes; do not select
one post hoc.

### `TEMPORAL_REFINEMENT_NUMERICALLY_UNRESOLVED`

**OPEN** if exact branch enumeration succeeds but any required nonlinear
solve, residual, reality or order gate is unresolved.

### `RELATIONAL_STEP_COMPOSITION_QUADRATIC`

**DERIVED COMPUTATIONAL LOCAL** only if there is exactly one compatible fine
branch and all three normalized discrepancies converge at second order or
better.

### `RELATIONAL_STEP_COMPOSITION_EXACT`

**DERIVED COMPUTATIONAL LOCAL** only if the unique branch is exact at the
registered symbolic orders and every nonlinear discrepancy remains below
the numerical floor.

## 8. Interpretation firewall

No outcome derives an absolute tick.  A positive composition result supports
a relational refinement limit and weakens, rather than strengthens, the
claim that one lattice slab is fundamental.  A negative result establishes a
finite-carrier pseudo-constraint/refinement failure under the complete
hypotheses; it does not turn the surviving coarse lapse into a fundamental
time quantum.

This protocol does not test generic nonzero initial velocity, spatial
refinement, anisotropic modes, a continuum Einstein limit, `c`, `G`,
`hbar`, Planck units, inflation or particle masses.

Only the new targeted verifier will be run.  The full suite will not be run.

