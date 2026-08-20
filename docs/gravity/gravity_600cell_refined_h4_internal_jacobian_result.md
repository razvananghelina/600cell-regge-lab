# Result: the refined H4 internal Jacobian is full rank but schedule dependent

Date: 2026-08-20

Prior-art commit: `4ea4430`  
Protocol commit: `9f1721c`  
Registered implementation commit: `84aa1dc`

## Complete hypotheses

Use `K0=P(sd K_600)`, every one of its 24 rank-colour staircase schedules,
the exact rank-derived geometry, the corrected complex Lorentzian Regge
action with boundary terms, equal fixed spatial boundaries, supplied
`tau0=0.0102`, the previously selected total dust mass and the conditional
local `P1` weights.

At the already certified off-shell induced flat fill, keep the twelve
boundary edge-square types fixed. Differentiate the ten total internal orbit
equations with respect to the logarithms of all six cross-diagonal and four
positive lapse-square coordinates. No Newton update, root, gauge condition,
physical target or effective boundary operator is used here.

## Reproducibility and controls

The targeted verifier passed `9/9` twice. Both executions wrote the same
artifact with SHA-256

```text
b900021c21df67c1de1ae18929be302b0d47d2f267c4a919388711a0a0bf5eaa
```

The maximum step-refinement difference was `2.01e-31` and the maximum
100-versus-140-decimal difference was `4.40e-86`. Raw imaginary parts and raw
Hessian antisymmetry remained inside the frozen entrywise envelopes. Six
independent centred second differences of the action reproduced their matrix
quadratic forms with maximum relative discrepancy `3.03e-37`. Every schedule
and its reverse agreed to at most `1.31e-125` entrywise.

These are reproducibility and internal-consistency controls, not the
mechanically different adversarial replication required for a later physical
claim.

## DERIVED COMPUTATIONAL result

Every one of the 24 internal Jacobians has

```text
rank = 10
inertia (positive, zero-compatible, negative) = (8, 0, 2).
```

For the lexicographically first schedule the eigenvalues are approximately

```text
-0.0183689761, -0.0182118419,
+0.0176033899, +0.0192099588,
+5.93917719, +15.3290252, +21.7883165,
+669.158924, +9232.50105, +12167.1137.
```

Thus no eigenvalue lies remotely inside the numerical uncertainty envelope,
and there is no local lapse/gauge null direction at this off-shell point. The
exact induced common-lapse tangent has image norm
`0.0183305189938364...`, not zero.

The 24 matrices form exactly 12 classes. Each class contains one schedule and
its time reverse; no two other schedules share a matrix. This is the first
schedule dependence in the internal quadratic response. Since it is measured
off shell, it is **STRUCTURAL**, not yet physical schedule dependence of an
on-shell canonical map.

The frozen outcome is

```text
REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_MULTIPLE_CLASSES
```

## Unapplied Newton diagnostic

For the first schedule, the frozen linear Newton proposal in coordinate order

```text
(cross_01,cross_02,cross_03,cross_12,cross_13,cross_23,
 rho_0,rho_1,rho_2,rho_3)
```

is approximately

```text
(+0.00208522,+0.00154997,+0.00137153,
 +0.00575932,+0.00381842,+0.01119168,
 -2,-2,-2,-2).
```

Its norm is `4.000022704...`; it was not applied. Because the coordinates
are logarithmic, the linear proposal would shrink all four positive lapse
squares by `exp(-2)`. This is a **PATTERN** suggesting attraction toward the
degenerate zero-lapse boundary. It is not proof that every positive
Lorentzian root is absent, and a failed local Newton run could not supply such
a proof.

## Framing attack and consequence

Full rank licenses a square, ungauged ten-equation stationary solve. It does
not establish existence, uniqueness or physicality of a root. The two
negative Hessian directions are also not instabilities of a solution because
the expansion point is off shell.

The 12 classes rule out solving one arbitrary colour order and calling its
answer canonical. At minimum one representative of every class must be
solved, with the time-reversed partner used as a control. Equality of the
initial residuals does not imply equality of the nonlinear roots.

## Next kill test

Preregister a bounded, high-precision nonlinear solve of all ten internal
equations for one representative of each of the 12 classes. Preserve positive
lapse squares, positive cross diagonals and the common Lorentzian branch, and
monitor approach to degeneracy explicitly.

- A certified finite positive root in every class advances the route.
- Inequivalent on-shell roots leave temporal canonicity unresolved unless an
  independent rule selects or sums schedules.
- Evidence that local solvers approach only lapse zero remains **OPEN** unless
  supplemented by a mathematical exclusion certificate.
- A proof that no common-branch positive root exists closes this refined
  static-background route.

Nothing here derives a tick, evolution law, `c`, `G`, Planck scale, graviton
or particle mass.
