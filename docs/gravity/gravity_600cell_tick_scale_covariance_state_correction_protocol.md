# State-correction protocol for the adversarial tick covariance audit

Date: 2026-08-21

Original scale protocol: `fd1f8a8`.

Exact flag-incidence/state diagnosis: commit `0fa8947`.

Status: frozen before evaluating any direct action on the corrected state.

All failed artifacts remain historical evidence.  Only the new targeted
verifier and static registry audit may be run.

## 1. Exact implementation error to correct

The frozen state has three separately indexed perturbations:

```text
old[0:30]      *= exp(1e-6*((i mod 7)-3)),
internal[0:35] *= exp(1e-6*((i mod 5)-2)),
new[0:30]      *= exp(1e-6*((i mod 11)-5)).
```

The first adversarial implementation concatenated `internal+new` and applied
the internal modulus-five rule to all 65 entries.  This violates the original
protocol.  No action, physical coefficient, scale factor or tolerance is being
changed here.

## 2. Mandatory state and carrier gates

Load definitions only from the primary and direct carrier sources.  Require
exact equality of all slab, stabilizer, edge, triangle, coordinate and sign
maps as in the accepted portions of the incidence audit.

Construct `old`, `internal`, and `new` independently and concatenate only after
their separate perturbations.  Expand both carrier states to every labelled
edge.  Require exact equality between carriers and exact invariance under all
24 stabilizer elements.

Parse the frozen source files structurally and require the primary source to
contain the separate `(7,3)`, `(5,2)`, `(11,5)` construction while the failed
adversarial source contains the recorded combined modulus-five construction.
This makes the repair cause mechanically checkable rather than narrative.

Load the flag-incidence artifact only as a combinatorial control.  Require 260
triangle, 100 simplex and 1000 flag orbits per parity, 24000 flags, and zero
shortcut/exact coefficient mismatches.

## 3. Corrected binary64 direct audit

Using the direct 2400-simplex evaluator and independently added dust term,
evaluate the corrected base state for both parities.  For
`alpha in {3/5,7/4}` evaluate both

```text
(alpha^2 q, alpha M),
(alpha^2 q, M)          [hostile fixed-mass control].
```

For the simultaneous scaling require

```text
S_scaled = alpha^2 S_base,
dS_scaled/dq = dS_base/dq
```

for all 95 raw squared-length derivatives within

```text
128*eps_binary64*2400*max(1,1/minimum_Gram_modulus).
```

Require that envelope below `1e-3`.  Every fixed-mass action and pole-gradient
defect must exceed `max(1e-8,100*envelope)`.  Retain all original branch and
nonzero-support gates.

## 4. Corrected arbitrary-precision direct audit

At 80 decimals, using the already disclosed literal 2400-simplex action
implementation (not the primary orbit action), evaluate for each parity

```text
alpha in {1,3/5,7/4}.
```

Require all six states on the Lorentzian branch and

```text
relative_error(S_alpha,alpha^2 S_1) < 1e-55.
```

Only after all direct values exist, load the primary artifact.  Require each
high-precision direct base action to agree with its own stored primary action
within `1e-45`.  Also evaluate the primary action freshly at the corrected
base and require it to reproduce the stored primary value within `1e-45`.

The corrected binary64 base action must agree with the high-precision direct
value inside its propagated envelope.  No cross-parity equality is required.

## 5. Outcome hierarchy

Assign exactly one:

1. `TICK_SCALE_STATE_CORRECTION_CONTROL_FAILED` for provenance, state,
   incidence, branch, support, hostile-control or envelope failure;
2. `TICK_SCALE_CORRECTED_IMPLEMENTATIONS_DISAGREE` if either direct scaling
   path fails or corrected direct and primary actions disagree;
3. `ABSOLUTE_CLASSICAL_TICK_NO_GO_ADVERSARIALLY_CORROBORATED` only if every
   gate passes.

Outcome 3 establishes, under the complete stated hypotheses:

> **DERIVED EXACT / ADVERSARIALLY CORROBORATED:** the zero-cosmological-
> constant classical Regge-plus-dust equations are globally scale covariant
> when all geometrized masses scale with the geometry, so they cannot select
> an absolute nonzero tick.

It leaves dimensionless `tau/L`, relative lapse, relational dust time and any
independently justified scale-breaking physics open.

