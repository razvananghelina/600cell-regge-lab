# Local branch stability: primary certificate

Date: 2026-08-22.

Status: **PRIMARY CERTIFIED / PENDING ADVERSARIAL REPLICATION**.

## Provenance

- prior-art gate: `5bace00`;
- frozen protocol: `6512791`;
- verifier registered before first run: `241a7d4`;
- first symbolic/control failure preserved: `0a276b4`;
- explicit radical-positivity correction: `a9c68dd`;
- second Arb-export failure and final narrow correction: `676eb5d`.

The failed runs are documented in
`gravity_600cell_finite_height_local_signature_first_run_failure.md`. Neither
correction changed a root bracket, precision, physical gate, terminal label or
acceptance criterion.

## Targeted result

The corrected verifier returned

```text
RESULT: 10/10 checks passed
OUTCOME: LOCAL_SIGNATURE_PRIMARY_CERTIFIED
```

Its artifact SHA-256 is

```text
9f524cc22df8cfb5083f372481b3efd19868252b85551d56378327eea7a6d613.
```

Only this targeted verifier was run. No full-suite claim is made.

## Primary theorem

Under the fixed homogeneous tetrahedral-frustum 600-cell action, zero
cosmological constant, conserved global dust, committed canonical-momentum
convention, positive height and positive endpoint scale, there exists an
unspecified `epsilon>0` such that every incoming state on

```text
(m,pi)=(mu(v),p(v)), |v-3/2|<epsilon,
```

has the same ordered physical tree through slab four:

```text
one physical first child;
two physical second children;
lower-q branch: DEAD at the next update;
higher-q branch: one third and one fourth child, then ENTERED_D.
```

The certificate proves this locally because:

- the diagonal `q=v` zero-height solution is an exact persistent,
  nondegenerate tangency with `p'(3/2)<0`;
- every other real root in all five visited states is simple;
- every stationary value, origin value and analytic-tail coefficient relevant
  to the complete real-root count is strictly separated from zero;
- every non-diagonal `h` and `r` gate is strict;
- every physical branch is disjointly ordered and every outgoing denominator
  excludes zero;
- the branch-B entry satisfies strict inequalities

  ```text
  2/5-m > 0.0042556251475,
  m*q-125 > 0.3317932609404;
  ```

- all five physical edges satisfy the independently redifferentiated complete
  action, with the largest printed residual below `6e-116`.

The exact radical identity needed for `partial_q E=p-pi` was proved by
factorization and positivity, not by numerical comparison.

## Interpretation

- **DERIVED PRIMARY:** the complete strict hypotheses for local constancy at
  `v=3/2`.
- **STRUCTURAL:** this removes the literal isolated-point objection. The
  accepted history is stable on some nonzero open neighbourhood of the
  one-parameter incoming curve.
- **OPEN:** adversarial corroboration, any explicit or maximal radius, and the
  global incoming basin.
- **NOT DERIVED:** that the neighbourhood is large or physically selected,
  that complete extendibility is a local law, or that the result covers the
  full `(m,pi)` plane or nonhomogeneous physics.

The margin `2/5-m` is small, and the invariant half-strip itself was selected
post-hoc. Therefore `generic` and `large basin` remain inadmissible even if the
primary certificate is corroborated.

## Required next gate

Before consolidation, preregister a mechanically different replication that
does not use interval Newton or the discovery roots as its decisive step. A
direct outward-rounded sign-bisection proof on fixed rational brackets, with
monotonicity and analytic tails, is the frozen candidate method.

