# Protocol: target-blind analytic weak-lapse jet of the cellular 600-cell map

Date: 2026-08-17

Prior-art commit: `b77856a`

Status: frozen before differentiating or expanding the closed cellular action
and before parsing any committed tick artifact.

## 1. Provenance firewall

This mission has two separately committed stages.

### Stage A: blind derivation

The blind verifier may read only:

```text
reproducible/gravity_600cell_homothetic_frustum_action.json
SHA-256 c0226a47607113930a31259d0cbee8ea33df2f7b0ba9416f9dbe5d647cede52d
```

and the formula frozen below.  It must not open, import or parse any file
whose name contains

```text
weak_lapse, canonical_lapse, second_tick, third_tick, fourth_tick.
```

Its artifact must contain

```text
tick_artifacts_parsed = false
```

and print the complete coefficient sequence before any numerical target
comparison.  Commit the Stage-A artifact and result note first.

### Stage B: disclosed comparison

Only after the Stage-A commit may a second registered verifier compare its
frozen coefficients with the existing `n<=4` artifacts.  Its input must be
the exact Stage-A artifact hash.  This ordering makes the integer sequence a
prediction of the action calculation rather than a re-fit to the ticks.

## 2. Frozen action and variables

Set `L0=1` during symbolic work.  Define

```text
epsilon3 = 2*pi-5*acos(1/3),
mu       = M/L0 = (90/pi)*epsilon3,
e        = tau_static/L0 > 0.
```

For one slab use exactly

```text
Delta = L_plus-L_minus,
h = sqrt(rho+Delta^2/4),
c = (Delta^2+2*rho)/(2*(Delta^2+3*rho)),
b = Delta/sqrt(8*(Delta^2+3*rho)),

S(L_minus,L_plus,rho)
 = 360(L_minus+L_plus)h[2*pi-5*acos(c)]
 + 600*sqrt(3)(L_minus^2-L_plus^2)asinh(b)
 - 8*pi*mu*sqrt(rho).
```

The static slab is indexed by zero:

```text
L_-1=L_0=1,
rho_0=e^2.
```

For `n=1,2,3,4`, put

```text
L_n   = exp(A_n e^2+B_n e^4+O(e^6)),
rho_n = e^2 exp(R_n e^2+O(e^4)).
```

No integer form for `A_n`, `B_n` or `R_n` is assumed.  `B_n` is a nuisance
coefficient required to decide the leading lapse correction honestly; it is
not a comparison target.

## 3. Equations to expand

Define the exact logarithmic lapse equation

```text
F_n = partial S(L_(n-1),L_n,rho_n)/partial log(rho_n).
```

Define the exact seam equation at `L_(n-1)` by

```text
G_n = partial/∂log(L_(n-1)^2)
      [S(L_(n-2),L_(n-1),rho_(n-1))
       +S(L_(n-1),L_n,rho_n)].
```

For `n=1`, the first term is the frozen static slab
`S(1,1,e^2)`.  This is equivalent to matching its exact post-momentum to the
pre-momentum of slab one; no stored momentum target is read.

Symbolically substitute the ansatz and extract the first nonzero
coefficients.  The expected static cancellations may remove the `e` term;
the verifier must discover, not assume, the first surviving order.  It must
expand through at least `e^5` so that a claimed `e^3` leading system has an
explicit next-order remainder.

Solve the coefficient equations recursively over the exact field generated
by rational numbers, `sqrt(3)`, `pi` and `acos(1/3)`.  At every step record:

- the two coefficient equations before substitution of the solution;
- the rank of the leading system for `A_n`;
- the next `2x2` Jacobian with respect to `(B_n,R_n)`;
- its exact determinant and rank;
- every algebraic solution;
- the unique real solution continuous with the contracting branch.

If more than one admissible real solution survives, do not choose by
proximity to the stored ticks: report nonuniqueness.

### 3.1 Preserved first blind-run correction

The registered implementation `80f3164` was run once without reading any
tick artifact.  It established exactly

```text
F_n/e starts at x=e^2,
G_n/e starts at x^0,
```

but then stopped before producing a coefficient artifact because the leading
`F` and `G` equations both fixed only `A_n`; `R_n` was absent.  At `n=1`, the
two equations share the same nonzero factor.  Thus the protocol's proposed
`2x2` leading solve was false: the leading system has rank one.

This is a scientific finding, not merely a software exception.  It is the
symbolic counterpart of the very soft lapse direction in the numerical
Jacobian.  The corrected frozen procedure is:

1. solve and cross-check the rank-one leading equations for `A_n`;
2. retain `B_n e^4` in `log L_n`;
3. solve the next lapse and seam coefficients jointly for `(B_n,R_n)`;
4. require this next `2x2` system to have exact rank two.

No target coefficient, stored tick state or integer ratio was read before
this correction.  The two-stage provenance firewall and all comparison
criteria are unchanged.

## 4. Blind checks

Stage A must pass all of the following.

1. Reconstruct the closed action from the frozen artifact and verify the
   static identity and static collective momentum independently.
2. Symbolic differentiation and independent arbitrary-precision centered
   differences agree at three generic nonstatic rational points below
   `1e-60`.
3. The first nonzero weak-lapse coefficient order is mechanically identified.
4. The recursive rank-one `A_n` systems and next rank-two `(B_n,R_n)` systems
   for `n=1..4` are exact and have their ranks stated.
5. Substitution of the selected coefficients makes every registered leading
   equation exactly zero.
6. At `e in {1/100,1/200,1/400}`, substitute the truncated `n<=4` series into
   the exact equations.  Residuals must decrease with the order predicted by
   the first omitted term; record observed halving orders but do not fit a
   coefficient.
7. Derive the sequences of increments, cumulative positions, lapse
   corrections and post-momenta from the coefficient solution.  Print all
   ratios even when they are nonintegers.

Numerical controls cannot rescue a failed exact coefficient equation.

## 5. Continuum control fixed in advance

The volume-matching map already used by the repository is

```text
L = zeta R,
zeta^3 = pi^2 sqrt(2)/50.
```

For closed continuum Friedmann dust at its maximum radius,

```text
ddot R/R = -1/(2R^2).
```

The static Regge slab has equal endpoints with opposite canonical momenta;
the conservative continuum interpretation is that its boundaries are the
two half-step points around a time-symmetric maximum.  Under that explicitly
stated half-step convention, the first post-static log-scale increment is

```text
log(L_1/L_0) = -(zeta^2/2)e^2+O(e^4).
```

Stage A must report the exact ratio of its derived discrete coefficient
`A_1` to `-zeta^2/2`.  Agreement is a control, not an acceptance condition.
No alternate time offset may be chosen after seeing the ratio.

## 6. Stage-A outcomes

### `CELLULAR_WEAK_LAPSE_JET_DERIVED`

Report **DERIVED** only if the exact series, rank, substitution and
independent residual-order controls all pass and there is one admissible
contracting coefficient branch through `n=4`.

### `CELLULAR_WEAK_LAPSE_JET_NONUNIQUE`

Report **DERIVED STRUCTURAL** if more than one admissible real coefficient
branch survives the frozen equations.  Print all branches; do not compare
them with ticks in Stage A.

### `CELLULAR_WEAK_LAPSE_JET_REFUTED`

Report **DERIVED NEGATIVE** if the action has no coefficient branch satisfying
the registered canonical equations through `n=4`.

### `CELLULAR_WEAK_LAPSE_JET_OPEN`

Use **OPEN** for unresolved symbolic limits, branches or numerical controls.

## 7. Stage-B comparisons and outcomes

Freeze these disclosed inputs only after the Stage-A artifact commit:

```text
weak-lapse artifact SHA-256
  500be1c4e2d7ec4104b9773bc1cfc71065c9d930607119eb616367d18fa5d8f9,

fourth-tick artifact SHA-256
  4d8d03957675a6f454c1ad05102ffd1711f48c2e5a19f09b2898a60c9f07020d.
```

The comparison verifier must:

1. compare every blind normalized coefficient through `n=4` with the
   committed Richardson intercept and its already frozen uncertainty band;
2. convert the blind absolute coefficient using the committed `tau0/L0` and
   compare with the stored `u1/lambda^2` and `v1/lambda^2` limits;
3. feed every committed finite-lapse `n<=4` state into the cellular lapse and
   seam equations and report residuals without re-solving;
4. confirm that no staircase carrier is needed for the comparison.

Report `CELLULAR_JET_EXPLAINS_FOUR_TICKS` only if all disclosed comparisons
pass.  Otherwise report `CELLULAR_JET_DOES_NOT_EXPLAIN_FOUR_TICKS` and list
the first failed observable.  A match proves an internal analytic
explanation of the committed weak-lapse law, not continuum convergence.

## 8. Scope fixed before evaluation

- Exact leading homogeneous canonical jet: tested.
- Existing four-tick integer law: compared only in Stage B.
- Continuum Friedmann coefficient: one fixed-resolution control.
- Full finite-lapse closed-form solution: not required.
- Spatial refinement: **OPEN**.
- Anisotropic modes, gauge separation and gravitational waves: **OPEN**.
- Fundamental lapse, `c`, Planck units and particle masses: **OPEN**.
- External novelty: **OPEN**.
