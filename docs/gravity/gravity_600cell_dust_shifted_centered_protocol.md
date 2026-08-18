# Preregistered blind protocol: shifted centered Jacobi coefficients

Date: 2026-08-18

The blind shifted Jacobi operator was committed in `ee57bcc`.  No shifted
centered coefficient, inertia or spectrum has been evaluated.  This is the
unique algebraic centered decomposition of that committed stencil; no target
or fitted normalization is admitted.

## Frozen input

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_shifted_jacobi.json` | `63b37b6000146d5d53dbbc01da5c9aba9a5e3373b8bc3830a404ef0f681ecf31` |
| `gravity_600cell_dust_shifted_jacobi.npz` | `d2f507c4a2fa11c5d7a808849c199a986278516f422cf43654f6de153ab170d0` |
| `verify_gravity_600cell_dust_shifted_jacobi.py` | `fc070de1a5f89524f119fd09ae25611559cb0bb64defbd34486b7980f058470d` |

Require outcome `SHIFTED_JACOBI_CERTIFIED`, pass count `8/8`, exactly `560`
arrays and the recorded NPZ hash.  No earlier centered, conformal, shape or
negative-mode artifact may be loaded.

## Unique centered decomposition

For each of the two schedules, seven sectors and four derivative variants use

```text
M_2 = (K_-^(2)+K_+^(2))/2,
N_2 = (K_+^(2)-K_-^(2))/2,
V_2 =  K_-^(2)+K_0^(2)+K_+^(2).
```

All `56` Flint determinant balls of `M_2` must exclude zero before forming

```text
Gamma_2 = M_2^-1 N_2,
Omega_2 = M_2^-1 V_2.
```

Require entrywise ball identities reconstructing all three original
coefficients and their normalized forms.  Store midpoint/radius arrays for
exactly `M_2,N_2,V_2,Gamma_2,Omega_2`, hence `560` arrays.

Record, without a target:

- Hermitian inertia of `(M_2+M_2^*)/2`;
- adjoint-defect ratios of `M_2,N_2,V_2`;
- singular and eigenvalue diagnostics for `Gamma_2,Omega_2`;
- reality classification of the `Omega_2` spectrum;
- fourteen even/odd schedule comparisons using the inherited `10/100` bands.

## Outcome hierarchy

1. `SHIFTED_CENTERED_CONTROL_FAILED` for provenance, carrier or archive
   failure.
2. `SHIFTED_CENTERED_MASS_SINGULAR` if any determinant ball contains zero.
3. `SHIFTED_CENTERED_IDENTITY_FAILED` if a centered identity fails.
4. `SHIFTED_CENTERED_SCHEDULE_DEPENDENT` if a comparison is resolved
   dependent.
5. `SHIFTED_CENTERED_SCHEDULE_OPEN` if none is dependent but one is open.
6. `SHIFTED_CENTERED_CERTIFIED` otherwise.

The inertia and spectrum do not affect the outcome.  No physical wave,
polarization, speed or growth label is allowed at this stage.  Only the
targeted verifier is run; the full suite is excluded.
