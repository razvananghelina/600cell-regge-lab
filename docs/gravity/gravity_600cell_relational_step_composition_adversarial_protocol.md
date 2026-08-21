# Adversarial protocol: same-state half-step obstruction

Date: 2026-08-21

Primary prior-art commit: `4f4b5c2`  
Primary protocol commit: `a025951`  
Primary registered implementation: `55cb3a2`

Status: frozen after the primary verifier returned
`SAME_STATE_HALF_STEP_BRANCH_ABSENT`, before constructing or running the
adversarial verifier.

## 1. Claim under attack

On the exact homogeneous cellular Regge-plus-conserved-dust action, use

```text
L_minus=1,
p_initial=180*epsilon*e,
rho_fine=(e^2/4) exp(O(e^2)),
log L_plus=O(e^2).
```

The primary result claims that the first nominal half-step cannot satisfy
both its complete lapse equation and the same incoming canonical momentum at
leading weak-lapse order.  Therefore a two-half-step history never reaches
its second slab within the registered branch class.

This is not a claim about generic nonzero initial velocity, nonanalytic
branches, or other temporal carriers.

## 2. Independent calculation order

Do not reuse the primary `S/e` series or either stored leading polynomial.

1. Reconstruct the unexpanded cellular action.
2. Differentiate that exact action first to obtain `F=rho*dS/drho` and
   `Pminus=(Lminus/2)*dS/dLminus`.
3. Only then substitute

   ```text
   Lminus=1,
   Lplus=exp(A*e^2),
   rho=e^2/4.
   ```

4. Independently evaluate

   ```text
   lim_(e->0+) F/e^3,
   lim_(e->0+) [Pminus+180*epsilon*e]/e.
   ```

5. Eliminate `A` exactly after substituting the physical
   `epsilon=2*pi-5*acos(1/3)`.  The primary verdict is corroborated only if
   the two limit equations have no common real root and their resultant is
   rigorously nonzero.

The use of exact-action-first differentiation is the load-bearing mechanical
difference from the primary action-series-first derivation.

## 3. Arbitrary-precision falsification

At the two separate candidate roots selected by the lapse and momentum limit
equations, evaluate the full unexpanded derivatives at 100 decimals for

```text
e in {1/100,1/200,1/400}.
```

Require:

- the equation used to define each candidate has the predicted vanishing
  leading limit;
- the other equation converges to the predicted nonzero scaled obstruction;
- the nonzero obstruction exceeds the Richardson drift by a factor of at
  least `100` at the finest pair.

No binary64 root decision is admissible.

## 4. Controls

- Reproduce the exact static action and momentum from the unexpanded
  derivatives.
- With the fine incoming momentum changed to `90*epsilon*e`, require the
  nonzero lapse root to satisfy both leading equations.  This proves the
  verifier can recover the old changed-state branch.
- Perturb the same-state coefficient from `180` to `90` only in the hostile
  control; no other coefficient may be adjusted.

## 5. Outcomes

- `SAME_STATE_HALF_STEP_ABSENCE_ADVERSARIALLY_CORROBORATED` only if every
  exact limit, elimination, high-precision and control gate passes.
- `PRIMARY_HALF_STEP_NO_GO_REFUTED` if the exact-action-first route finds a
  common admissible root.
- `HALF_STEP_ADVERSARIAL_DISAGREEMENT` for any unresolved symbolic or
  precision disagreement.  The consolidated claim then remains **OPEN**.

Only this targeted adversarial verifier will be run.  The full suite will
not be run.

