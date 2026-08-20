# Result: the inherited refined H4 fill is off shell

Date: 2026-08-20

Prior-art commit: `5518fa7`  
Protocol commit: `3d36c54`  
First failure record: `f64809b`  
Control-correction protocol: `f463784`  
Corrected implementation: `3ddd66b`

## Complete hypotheses

Use `K0=P(sd K_600)`, its rank colouring, every one of the 24 standard
colour-ordered staircase triangulations of `K0 x I`, the corrected complex
Lorentzian Regge action including boundary terms, the exact rank-derived
projected chamber geometry, supplied `tau0=0.0102`, the previously selected
total dust mass and the conditional local `P1` mass weights.

Restrict the edge metric only by the full spatial `H4` symmetry shared by all
schedules.  The invariant slab has 22 positive coordinates: six old and six
new boundary edge-square types, six cross-diagonal types and four vertical
positive lapse squares whose Lorentzian edge squares are negative.

The tested point has equal old/new boundaries, one common lapse square and
cross diagonals induced by the flat static product.  No coordinate is solved
or fitted in this census.

## Reproducibility

The corrected targeted verifier passed `12/12` twice.  Both executions wrote
the identical 518,506-byte artifact with SHA-256

```text
283be37bc7530a3cc4fce9e279272359f107f09fb7b1b0eaff141059bfb4e018
```

The first execution stopped before its first action evaluation and remains
documented separately.  Its two bookkeeping corrections were preregistered
before the successful run.

Independent controls include 100- versus 140-decimal arithmetic and centred
finite differences at `1e-15` and `5e-16`.  The maximum analytic/Richardson
relative disagreement over 44 checked derivatives is
`3.37e-54`.  The certified per-edge zero envelope is `1.00e-60`.

## DERIVED exact combinatorial facts

Every schedule has

```text
57,600 pentachora,
149,280 distinct triangles,
28 H4 triangle types,
0 mixed incidence signatures.
```

Every triangle of one state type has the same complete abstract
simplex/hinge incidence signature.  This directly validates the 22-variable
reduction rather than inferring it from a dimension count.

The exact rank formula reconstructs

```text
V_chord = 19.1479329183128449930044943534...
s0      = 1.01018895214456343488076005472...
M       = 2.36580263014636331515849057251...
```

and the `P1` weights place exactly `M/4` in each rank class.  The differences
from the older rounded artifacts are `2.01e-15` in volume and `1.27e-13` in
mass.

## DERIVED variational result

All 24 schedules give one and the same internal residual vector.  The action
spread is `1.84e-136`, and the maximum time-reversal-pair difference is
`3.36e-140`.

All `24*6=144` cross-diagonal residuals are zero-compatible.  Their displayed
magnitudes lie between approximately `5e-142` and `2e-140`, far inside the
`1e-60` envelope.

All `24*4=96` vertical residuals are certified nonzero.  For any schedule the
per-edge log residuals are

| rank | vertices | per-edge log residual |
|---:|---:|---:|
| 0 | 120 | `-3.06316391703898833459649578502e-4` |
| 1 | 720 | `+4.87292463215819812173890650252e-5` |
| 2 | 1,200 | `+3.18367351901373692665778942454e-5` |
| 3 | 600 | `-6.08852876253933493020927508206e-5` |

The corresponding total rank-orbit log residuals are

```text
(-0.0367579670044678600151579494203,
 +0.0350850573515390264765201268182,
 +0.0382040822281648431198934730945,
 -0.0365311725752360095812556504924).
```

They cancel in the common induced-lapse direction: the maximum residual over
all schedules is `2.55e-138`.  Thus the old collective lapse equation is
reproduced to high precision while each of the four equations hidden inside
it fails decisively.

The frozen outcome is

```text
REFINED_H4_INDUCED_FILL_OFF_SHELL
```

## Meaning and framing attack

**DERIVED NEGATIVE:** the inherited static product with conditional `P1` dust
is not a stationary point of the refined invariant Regge system.  Therefore
its action Hessian cannot be used as the action-generated canonical boundary
map.  The earlier acceleration result remains a valid homothetic
minisuperspace projection, but it is not a solution of the enlarged equations.

**DERIVED POSITIVE CONTROL:** the diagonal equations vanish and the complete
residual vector is schedule independent at this point.  This is stronger than
a numerical accident but is not schedule independence of the dynamics: it is
only equality of actions and first derivatives at one off-shell flat fill.

Changing the four rank masses after seeing these residuals could cancel the
four equations directly.  That would be fitted matter and is forbidden.  The
`P1` distribution may be replaced only by an independently derived matter
principle preregistered before these values are compared.

## Next falsification boundary

Keep both equal spatial boundaries, the selected mass and the conditional
`P1` weights fixed.  For every schedule solve the full ten internal equations
for six cross diagonals and four vertical lapse squares, starting from the
induced fill with a frozen continuation/trust rule.  The full ten variables
must remain free even though the initial cross residuals vanish; Hessian
coupling can move them.

- If no common-branch Lorentzian stationary fill exists, the refined static
  background route closes.
- If stationary fills exist but differ physically between schedules, temporal
  canonicity fails unless an independent rule selects or sums schedules.
- If all fills are related and their effective boundary Hessians agree, the
  route advances to non-invariant sectors and only later to refined mode
  transport.

This result derives no evolution, graviton, tick, `c`, `G` or Planck scale.

