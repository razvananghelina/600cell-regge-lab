# Result: the shifted sign gate is blocked by binary serialization

Date: 2026-08-18

## Headline

The enlarged uncertainty at the shifted tick is not presently a structural
failure of the negative-shape candidate.  In all `16/16` disclosed cells the
certified error is dominated by re-enclosing a binary64 tangent midpoint
before reconstructing the principal-function block.

The preregistered outcome is

```text
SHIFTED_PRECISION_BINARY_SERIALIZATION_DOMINANT.
```

This result authorizes a separate direct high-precision reconstruction from
the action Hessian.  It does **not** certify the `30` negative modes, and it
does not permit deleting the binary half-ULP enclosure from an archived
midpoint.

## Provenance ledger

| stage | commit |
|---|---|
| post-result literature and framing gate | `52f337b` |
| frozen precision protocol | `86a9129` |
| registered verifier before execution | `c028aa0` |
| deterministic artifact | `67e99fb` |

The verifier was executed twice.  Its JSON artifact was byte-identical, with
SHA-256

```text
409e428ca05b4b6c6e380d7af6d84fc3834afc1fdb95065d6f0a05618e1d2cee.
```

Only the targeted verifier was run; the full suite was not run.

## What was measured

For sectors `4,5`, both parities and all four frozen derivative variants, the
committed `T_2` tangent was converted to 80-decimal Flint balls in two ways:

1. a rigorous binary mode that includes the mandatory half-ULP radius of each
   stored binary64 midpoint;
2. a diagnostic stored-ball-only mode that omits that half-ULP.

The second object is explicitly non-rigorous.  It is used only to attribute
the error, never to infer a sign.

All `32` boundary-twist determinant balls exclude zero.  All `16` rigorous
reconstructions reproduce the committed `Kminus` balls, and all `16`
downstream restricted-error sums reproduce exactly within the frozen binary
roundoff envelope.

## Quantitative result

Across the complete `16`-cell census,

```text
rigorous/counterfactual S_2,10 radius ratio
    1.4893576448e50 ... 1.6667652926e50,

Kminus radius / (Kzero radius + Kplus radius)
    45.0510610413 ... 45.0510610877.
```

Thus every cell crosses both preregistered dominance thresholds:

```text
16/16  SERIALIZATION_DOMINANT_RESOLVED.
```

For a representative cell, the reconstructed restricted error is about
`6.5669e-6`.  Roughly `58.38%` is already present in the restricted source
ball and `41.62%` is the certified carrier-lift contribution; the final
binary arithmetic floor is only about `1.08e-5` of the total.  The upstream
`Kminus` block is therefore the load-bearing source.

## Interpretation

- **DERIVED COMPUTATIONAL:** the shifted OPEN result is caused by a specific
  information-losing numerical interface: high-precision tangent balls are
  serialized through binary64 midpoints and later inverted again.
- **DERIVED COMPUTATIONAL:** this interface dominates the other shifted
  principal-function coefficient radii by a factor of about `45`.
- **NOT DERIVED:** removing that interface will resolve the signs.  A direct
  high-precision calculation must still be run with rigorous enclosures.
- **OPEN:** persistence, rotation or loss of the old rank-`15` fibers;
  their projectors; any reduced propagator; physical instability; inertia of
  matter; graviton interpretation; continuum limit; and `c`.

Here “inertia” in the certified calculations means the signature of a
Hermitian form—the counts of positive, negative and zero directions—not the
observed inertial mass of a particle.

## Next load-bearing calculation

Reconstruct the `T_2` and `T_3` Hamilton principal-function blocks directly
from the high-precision action Hessians, without an intermediate binary64
tangent archive.  Freeze the derivative schedules, inputs, error propagation
and sign thresholds before inspecting the resulting stiffness signs.  Only
that calculation can decide whether the `30` negative-shape modes persist.
