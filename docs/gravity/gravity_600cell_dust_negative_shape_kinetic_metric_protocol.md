# Protocol: canonical kinetic-metric Rouche certificate

Date: 2026-08-18

Prior-art gate commit: `cdcaf8d`

This protocol is **TARGET-DISCLOSED**.  The Euclidean midpoint and literal
coefficient-ball count `15/0/15`, the Euclidean minimum sampled ratio `8.17`,
and the failure of all Euclidean `100x` covers are already known.  No output of
the kinetic-metric calculation may change the metric, carrier, error transport,
contour or outcome ordering below.

## 1. Frozen cells and carrier

Reconstruct the same `16` cells used by the committed root-count verifier:

```text
parity  = even, odd
sector  = 4, 5
variant = operational_primary, operational_shadow,
          validation_primary, validation_shadow.
```

Use the exact frozen hashes, conformal carrier, shape carrier, negative
stiffness eigenspace, Gamma/Omega restrictions, invariance residuals and error
formulas of that verifier.  The committed root-count JSON must have outcome
`NEGATIVE_SHAPE_ROOT_COUNT_SAFETY_OPEN`, `13/13` passing checks, all `16`
literal covers passing, all `16` Euclidean `100x` covers open, and midpoint
counts `15/0/15` everywhere.

## 2. Kinetic metric fixed before evaluation

On the shape carrier form

```text
B_S = - W^* M W,
B_- = E_-^* B_S E_-.
```

Here `M` is the symmetrized centered kinetic midpoint, `W` is the inherited
shape basis and `E_-` is the inherited `15`-column negative-stiffness basis.

Propagate the kinetic error first to `B_S` with the frozen shape-subspace
formula and then to `B_-` with the frozen negative-subspace formula.  Classify
`B_-` as positive-resolved only if

```text
lambda_min(B_-) > 100 epsilon_B.
```

No regularization or eigenvalue clipping is permitted.  If positivity is not
resolved, the mission stops `KINETIC_OPEN`.

For a positive-resolved cell compute by Hermitian eigendecomposition the unique
positive midpoint square root and inverse

```text
S     = B_-^(1/2),
S_inv = B_-^(-1/2).
```

Record `||S^2-B_-||`, `||S S_inv-I||`, `kappa_2(B_-)`, and
`kappa_2(S)=sqrt(kappa_2(B_-))`.  Both reconstruction residuals must be below

```text
1000 eps_machine m max(1, ||B_-||, ||S||^2)
```

with `m=15`.

## 3. Mandatory error transport

Use

```text
Q_B(z) = S Q(z) S_inv,
A2_B = S (I+Gamma) S_inv,
A1_B = S (-2I+Omega) S_inv,
A0_B = S (I-Gamma) S_inv.
```

The only admitted transported coefficient errors are

```text
epsilon_Gamma_B = kappa_2(S) epsilon_Gamma,
epsilon_Omega_B = kappa_2(S) epsilon_Omega.
```

The contour evaluation floor remains the preregistered

```text
1000 eps_machine m max(1, 2||A2_B||+||A1_B||+||A0_B||).
```

It covers floating construction/evaluation of the fixed midpoint transform.
Omitting `kappa_2(S)`, transporting only favorable perturbation directions, or
optimizing another similarity is a control failure.

## 4. Continuous contour calculation

Run exactly the previous deterministic unit-circle cover on `Q_B`:

```text
initial intervals = 256
maximum depth     = 32
maximum evaluated = 2,000,000
safeties           = 1 and 100
```

The pointwise perturbation is

```text
delta_B(theta)
 = epsilon_Gamma_B |exp(2 i theta)-1| + epsilon_Omega_B.
```

No interval, parity, sector or derivative variant may be removed.  Record the
same complete cover ledger as before, including the minimum sampled ratio and
the first failed point.

The exact roots and root counts are similarity-invariant.  Do not diagonalize
the companion again as a new search.  Transfer the already certified midpoint
count only if the transformed literal cover passes and the square-root controls
pass.

## 5. Controls

For every cell:

1. multiply `B_-` by the fixed scalar `7`; the transformed coefficient
   matrices must agree with those from `B_-` to within `100` times the
   transformation floor, the transported condition factor must agree to the
   same relative floor, and both cover verdicts must be identical;
2. conjugate `Q` by the fixed reversal permutation (unitary); its singular
   values at the `256` preregistered interval centres must agree with the
   unweighted values within `100` evaluation floors, confirming invariance
   under a unitary coordinate relabelling;
3. verify explicitly that the transported errors equal
   `kappa_2(S)` times the Euclidean errors stored in the ledger;
4. verify that no matrix other than the positive square root of `B_-` enters
   the transform.

## 6. Outcomes

Outcome order:

1. `NEGATIVE_SHAPE_KINETIC_METRIC_CONTROL_FAILED` if provenance, carrier,
   square-root, scalar-rescaling, unitary or error-transport controls fail;
2. `NEGATIVE_SHAPE_KINETIC_METRIC_OPEN` if any `B_-` is not
   positive-resolved;
3. `NEGATIVE_SHAPE_KINETIC_LITERAL_OPEN` if any transformed literal cover
   fails;
4. `NEGATIVE_SHAPE_KINETIC_SCHEDULE_DEPENDENT` if transformed `100x` verdicts
   differ among parity/sector/variant cells;
5. `NEGATIVE_SHAPE_KINETIC_SAFETY_OPEN` if all transformed literal covers pass
   and all transformed `100x` covers remain open;
6. `NEGATIVE_SHAPE_LOCAL_HYPERBOLIC_KINETIC_RESOLVED` only if all transformed
   `100x` covers pass and all inherited counts are `15/0/15`;
7. `NEGATIVE_SHAPE_KINETIC_MIXED_OPEN` for any complete but uniformly
   nonbinary pattern not covered above.

If outcome 5 holds, the numerical-margin branch is closed: no further
coordinate norm may be tried.  The next physical mission must be evolution
across independently solved, nonidentical slabs.  Outcome 6 would still be a
local frozen statement and would not alter that physical next step.

