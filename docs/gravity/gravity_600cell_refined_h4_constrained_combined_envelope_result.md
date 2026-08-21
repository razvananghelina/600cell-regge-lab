# Combined constrained action/Hessian envelope: corroborated

Date: 2026-08-21

## Provenance

| stage | commit |
|---|---|
| frozen algebraic protocol | `5cefff6` |
| verifier registered before execution | `711105b` |

The verifier passed `10/10` twice with byte-identical artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_combined_envelope.json
SHA-256 34e2d598a6f608c9436217024138b32f0095e5df8e32d3ff91df2b182843aa0d.
```

It evaluated no geometry, action, Hessian, solve or class census.  No full
suite or deferred nonlinear census was run.

## Result

For every one of the twelve frozen direct-action directions, the stored
action ladder has target-independent convergence ratios

```text
R differences: 16,16,
X differences: 64.
```

The exact entrywise propagation

```text
e_total=e_action+||y||_1^2 e_K
```

accepts `12/12` action/Hessian comparisons.  The maximum
error-to-combined-envelope ratio is

```text
0.00066391091.
```

Removing the Hessian uncertainty rejects all `12/12`, while every corrupted
quadratic lies at least `2.85e23` combined envelopes away.

The frozen outcome is

```text
REFINED_H4_COMBINED_ENVELOPE_CORROBORATED.
```

## Meaning

- **DERIVED EXACT / COMPUTATIONAL:** the combined bound follows from
  `|y^T Delta K y|<=||y||_1^2 e_K` and accepts every independently extrapolated
  action value.
- **DERIVED NEGATIVE:** the former fixed relative `1e-28` control treated the
  finite-difference Hessian as much more accurate than its own frozen
  envelope allowed.
- **PROCEDURAL:** this licenses a preregistered repair and complete rerun of
  the primary verifier.  It does not retroactively accept its tentative
  single schedule class.
- **OPEN:** the primary rerun and mechanically independent replication.

