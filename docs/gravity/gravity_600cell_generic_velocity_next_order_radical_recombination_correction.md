# Correction protocol: prevent radical recombination

Date: 2026-08-21

Frozen radical inventory commit: `71b8312`.  First normalization
implementation commit: `1046d01`.

Preserved third-execution record:

```text
reproducible/gravity_600cell_generic_velocity_next_order_third_timeout.json
SHA-256 f6ceb8496ef87d58a04d13b678224854d51f19348a64f2becb25edc3caa92147
```

The eight exact factorizations passed internally, but the implementation then
called

```text
powsimp(...,force=True)
```

after replacing composite positive square roots by the two primitive
radicals.  That call recombined the primitives into new expanded polynomial
square roots.  The frozen post-normalization inventory gate rejected them,
and root simplification was interrupted.  No root or physical outcome was
computed.

Remove `powsimp(force=True)` after radical normalization.  Use only
`together`, `cancel` or collection in the endpoint coefficient, none of which
may join positive square-root factors.  Apply the frozen normalizer once more
after any generic simplification and require the inventory gate to pass
before constructing the two polynomials in `a`.

Do not add the recombined products to the frozen inventory: they were created
by the simplifier rather than by the derivative.  No equation, branch,
sample, threshold or outcome rule changes.
