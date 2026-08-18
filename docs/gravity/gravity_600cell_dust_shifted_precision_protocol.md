# Preregistered protocol: shifted principal-function error budget

Date: 2026-08-18

Prior-art gate commit: `52f337b`.

Status: **TARGET-DISCLOSED, PREREGISTERED BEFORE ANY PRINCIPAL-BLOCK
SERIALIZATION RATIO OR DOWNSTREAM ERROR FRACTION IS COMPUTED.**

The precision-open shifted result and its two sectors are disclosed.  This is
not a new sign test.  The verifier is forbidden to load an eigenvalue or sign
label.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_two_step_full_tangent.json` | `f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc` |
| `gravity_600cell_dust_two_step_full_tangent.npz` | `ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d` |
| `verify_gravity_600cell_dust_two_step_full_tangent.py` | `c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717` |
| `gravity_600cell_dust_shifted_jacobi.json` | `63b37b6000146d5d53dbbc01da5c9aba9a5e3373b8bc3830a404ef0f681ecf31` |
| `gravity_600cell_dust_shifted_jacobi.npz` | `d2f507c4a2fa11c5d7a808849c199a986278516f422cf43654f6de153ab170d0` |
| `verify_gravity_600cell_dust_shifted_jacobi.py` | `fc070de1a5f89524f119fd09ae25611559cb0bb64defbd34486b7980f058470d` |
| `gravity_600cell_dust_shifted_centered.json` | `265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47` |
| `gravity_600cell_dust_shifted_centered.npz` | `c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8` |
| `verify_gravity_600cell_dust_shifted_centered.py` | `a3c45e3e636057d83a663d3248dd023f7d04ec6e544c698f9116307822be337a` |
| `gravity_600cell_dust_shifted_shape_stiffness.json` | `14fe5bc91e3ae4712c6ea19b8120785e2facd364e1ceb194009123fa353a4315` |
| `verify_gravity_600cell_dust_shifted_shape_stiffness.py` | `031d0dd1cab45d0093015fcab7ce7b56e098a5742895eed71f5a531aee31c2a6` |

Require the inherited outcomes and exact array counts.  Only the `T_2`
tangent, shifted coefficient and shifted centered archives are numerical
inputs.  The stiffness JSON may supply `eta_S` and the final restricted error,
but the verifier must reject source code containing access to
`A_eigenvalues`, `A_sign_labels`, `sign_counts` or a target sign.

## Complete 16-cell census

Use both parities, disclosed sectors `4,5` and all four derivative variants:

```text
2 * 2 * 4 = 16 cells, each with position dimension 30.
```

For the committed `T_2` block, independently reconstruct

```text
T_2=[A B; C D],
S_2,10=C-D B^-1 A.
```

Use 80-decimal Flint balls in two declared modes.

### Certified binary mode

For each real and imaginary binary64 component use

```text
stored radius + half ULP.
```

This must reproduce the committed shifted `Kminus` midpoint and enclose its
stored radius componentwise.  Failure is a control contradiction.

### Stored-ball-only counterfactual

Use the same binary midpoint but omit the half-ULP, retaining only the stored
Flint radius.  This object is explicitly **NOT a valid enclosure of the
pre-serialization high-precision value**.  It is used only to measure how much
of the propagated radius is caused by the mandatory binary midpoint
re-enclosure.  It cannot be used for a sign, projector or physics claim.

For each mode record input component-radius Frobenius norms and the output
`S_2,10` radius Frobenius norm.  Define

```text
R_serial = radius_F(S10 certified binary)
           / radius_F(S10 stored-ball-only).
```

## Coefficient and downstream budget

From the shifted Jacobi archive record the radii of `Kminus,Kzero,Kplus` and

```text
R_K = radius_F(Kminus)/(radius_F(Kzero)+radius_F(Kplus)).
```

From the shifted centered `V` archive and committed shape-carrier `eta_S`,
reconstruct exactly the existing restricted-error sum

```text
epsilon_HV = ||radius_HV||_F
             + 1000 eps n max(1,||H_V||_2),

carrier_lift = 2 eta_S (||H_V||_2+epsilon_HV),

arithmetic = 1000 eps n max(1,||H_V||_2),

epsilon_VS = epsilon_HV + carrier_lift + arithmetic.
```

Require agreement with the stored `restricted_A_error` up to a binary
roundoff envelope.  Record all three absolute contributions and fractions.

## Frozen classifications

Per cell:

```text
SERIALIZATION_DOMINANT_RESOLVED  R_serial > 100 and R_K > 10,
SERIALIZATION_MIXED              neither ratio contradicts, but at least one
                                 does not cross its resolved threshold,
SERIALIZATION_NOT_DOMINANT       R_serial < 10 or R_K < 1.
```

Values between the `10/100` or `1/10` bands are mixed.  The complete outcome
uses the first applicable branch:

1. `SHIFTED_PRECISION_AUDIT_CONTROL_FAILED`;
2. `SHIFTED_PRECISION_SERIALIZATION_NOT_DOMINANT` if any cell is resolved not
   dominant;
3. `SHIFTED_PRECISION_ATTRIBUTION_MIXED` if none is not-dominant but any is
   mixed;
4. `SHIFTED_PRECISION_BINARY_SERIALIZATION_DOMINANT` only if all `16` cells
   cross both dominance thresholds and all downstream sums reproduce.

Outcome 4 authorizes a separate preregistered direct high-precision
principal-function reconstruction from the action Hessian.  It does not
authorize deleting half-ULPs from the current archive or choosing a smaller
error by hand.

## Exclusions

- no stiffness eigenvalue, sign, target rank or projector;
- no comparison against a desired negative result;
- no change to derivative steps or action;
- no physical instability, wave, speed or continuum claim;
- no full-suite run.
