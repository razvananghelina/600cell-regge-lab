# Protocol: full pole equation on the unique homogeneous weak line

Date: 2026-08-20  
Status: **preregistered deterministic classification**

## Frozen inputs

Pin by SHA-256:

- this protocol and its prior-art/framing gate;
- the canonical-lapse source and accepted `7/7` artifact;
- the primary exact homogeneous-line `10/10` artifact;
- the repaired adversarial homogeneous-line `7/7` artifact and consolidation
  note.

Require both staircase parities to retain the same accepted endpoint and
calibrated rank-two Jacobian.  No new root, finite difference, fitted direction,
precision level or tolerance is allowed.

## Deterministic calculation

For each parity load the operational endpoint matrix

```text
J=[[a,b],[c,d]]
```

and its already calibrated entrywise error budget `epsilon`.  Set

```text
v=(-d,c).
```

Check exactly in the loaded decimal arithmetic

```text
[c,d] dot v = 0,
[a,b] dot v = -det(J).
```

For entrywise perturbations bounded by `epsilon`, use the determinant error
bound

```text
B = epsilon*(abs(a)+abs(b)+abs(c)+abs(d)) + 2*epsilon^2.
```

Require

```text
abs(det(J)) > 100*B,
abs(det(J)-stored_determinant) < 1e-45.
```

Construct the independently stored primary line in `(delta s,delta z)` as

```text
v_primary=(-p_z,p_s)
```

from the exact-action bridge.  Require its normalized projector to agree with
that of `v` below `1e-30`.  Require the two parity projectors below `1e-70` and
their certified determinant signs to agree.

Controls:

- a planted rank-one `2 x 2` matrix must give zero pole derivative on its
  momentum-row kernel;
- the accepted matrix must give a pole derivative separated from zero by the
  determinant interval;
- all upstream outcome and parity gates must remain true.

## Mechanical outcomes

1. `HOMOGENEOUS_POLE_TRANSVERSALITY_CONTROL_FAILED`;
2. `HOMOGENEOUS_POLE_TRANSVERSALITY_CONVENTION_DISAGREEMENT`;
3. `HOMOGENEOUS_WEAK_LINE_FULLY_NULL_OPEN` if the determinant interval contains
   zero;
4. `HOMOGENEOUS_WEAK_LINE_TRANSVERSE_TO_POLE_EQUATION` if every check passes.

Outcome 4 is **DERIVED COMPUTATIONAL from already certified endpoint data**.  It
classifies the weak line as off-shell for the full fixed-input equations and the
nonstatic endpoint as locally isolated.  It does not refute the endpoint or the
existence of a discrete evolution map, and it does not derive an absolute tick,
`c`, `G` or Planck units.

Only the targeted verifier and static registry check may run.  No full suite.

