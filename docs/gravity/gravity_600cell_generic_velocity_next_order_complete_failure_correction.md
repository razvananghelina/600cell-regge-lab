# Correction protocol: first complete next-order execution

Date: 2026-08-21

Registered implementation: `bd68c52`.  Scaled-jet correction: `d8dde35`.
Frozen radical inventory: `71b8312`.  Recombination correction: `d6bc954`.

Preserved first complete artifact:

```text
reproducible/gravity_600cell_generic_velocity_next_order_first_complete_failure.json
SHA-256 eac5353f2ac93090a4c3a94f64fbc6b711519eed473ae14a4937670fae0e9c30
```

The execution returned `8/10` and the honest outcome
`GENERIC_NEXT_ORDER_OPEN`.  It derived equal generic one-slab roots and zero
conditional composition defects, but this is not yet accepted.

## 1. Hostile structural reductions

The mass-shift defect was left as

```text
-4*pi*(v^2+4)/[5(v^2+4)],
```

and the momentum-shift defect retained equivalent inverse-function and
composite-radical spellings.  Apply the already-frozen `stable_normalize`
function before testing both defects.  Require exact residuals `-4*pi/5` and
`-1/10`; no hostile target changes.

## 2. Frozen numerical failure

Eleven of the twelve `(v,a,observable)` controls met the preregistered
first-order interval.  The sole failure was

```text
v=3/2, a=0, momentum:
errors=(2.256955964e-6,1.459798247e-6),
observed order=0.6286092769.
```

Do not change or delete this failed gate.  Before any adjudication, derive the
next error coefficient exactly or evaluate the same unexpanded residual at
the separately registered diagnostic heights

```text
h in {1/1600,1/3200,1/6400,1/12800}
```

at 100 decimals.  These points diagnose whether the frozen pair was
pre-asymptotic because of coefficient cancellation; they do not replace it.
If consecutive orders do not approach one, the scaled-jet result is
disagreed with and remains `OPEN`.

## 3. Missing exceptional-velocity proof

The verifier incorrectly set `classification_complete=True` merely because
both expressions have generic polynomial degree one in `a` and their generic
roots simplify to the same expression.  A symbolic degree does not prove that
the two leading coefficients are nonzero at every real `v!=0`.

Before a global verdict:

1. record both coefficients of `a` after frozen normalization;
2. prove their real nonzero domains, or enumerate every real nonzero zero;
3. at every degree-drop velocity, substitute directly into the unsolved
   `C1` and `P1` equations and classify all roots;
4. only then set `classification_complete`.

Numerical sign scans are controls, not a proof of the exceptional set.

## 4. Composition scope

The recorded zero endpoint, momentum and action defects remain conditional on
the global one-slab classification.  They may not be promoted while the
exceptional-velocity proof or numerical diagnostic is open.

No equation, state, sample already used, branch, physical label or outcome
hierarchy changes.

