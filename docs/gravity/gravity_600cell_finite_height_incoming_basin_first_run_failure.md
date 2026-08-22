# First incoming-basin discovery run: implementation failure

Date: 2026-08-22.

Verifier commit: `7be5a04`.

Status: **PRESERVED IMPLEMENTATION FAILURE; NO SCIENTIFIC VERDICT**.

The first execution of
`reproducible/verify_gravity_600cell_finite_height_incoming_basin_discovery.py`
was made only after the discovery protocol and verifier had been committed and
pushed.  It evaluated all 3072 frozen incoming nodes, then reported

```text
[FAIL] all frozen Chebyshev nodes have complete all-real branch trees
       nodes=3072; depth=4
[FAIL] known representative signatures are reproduced only after discovery
       v=3/2 labels=[]; second counts=[None, None, None]
RuntimeError: control state lacks its unique first slab
```

No JSON artifact was written because the post-discovery control raised before
the final serialization block.

## Cause

At an initial state

```text
m=mu(v), pi=p(v),
```

the diagonal zero-height solution `q=v` obeys both

```text
E(mu(v),p(v),v)=0,
E_q(mu(v),p(v),v)=p(v)-p(v)=0.
```

The frozen protocol explicitly requires this real root to be recorded and
then excluded as the known diagonal zero-height slab.  The first
implementation instead sent every stationary zero of `E` to the generic
multiple-root ambiguity gate before applying the initial diagonal exclusion.
Consequently every initial node was labelled unresolved.  This is a mismatch
between the implementation and the preregistered protocol, not evidence about
the physical branch diagram.

## Permitted correction

The correction must do only the following:

1. pass the frozen initial `v` into the complete root census;
2. identify exactly one stationary root matching `q=v`;
3. verify the two displayed identities at the accepted residual scale;
4. serialize that root in the all-real census;
5. exclude it only when applying the positive-height physical gate;
6. retain the generic `UNRESOLVED` verdict for every other stationary zero.

The protocol commit, input hashes, 3072 nodes, precision, depth, tree budget,
physical gates, candidate rules and outcome criteria remain unchanged.  The
corrected implementation must be committed before its first execution.

## Second execution: tail-anchor implementation failure

Corrected-verifier commit: `d2c791d`.

The second execution passed provenance, thresholds and the first 128 complete
incoming trees. It then raised

```text
TypeError: bad operand type for abs(): 'NoneType'
```

inside the analytic-tail bracketing helper. A later canonical state had no
finite stationary point. In that case `E` is monotone on the whole real axis;
the two analytic tail signs still provide a complete root bracket, but the
implementation tried to size the tail probe from a nonexistent finite
boundary. The protocol neither requires nor permits this state to be
discarded.

The conforming correction is to use `q=0` solely as the finite magnitude
anchor when an interval has no finite endpoint. The root census already
checks that `E(0)` is noncritical before tail bracketing. No numerical box,
grid, branch rule, precision or outcome criterion changes. No JSON artifact
was produced by this execution.
