# Weak-scale collective path: frozen result

Date: 2026-08-13

Prior-art and coordinate correction: `ad5f0ad`

Frozen protocol: `8380f0d`

Implementation commit: `caf5caa`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_weak_scale_path.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_weak_scale_path.json`

Targeted run: **11/11 implementation checks passed**.  The full suite was not
run.

## 1. Result

Both even and odd schedules receive

```text
ALL_13_PATH_POINTS_STATIONARY_WITHIN_ERROR.
```

All thirteen parity comparisons receive

```text
PARITY_AGREES_WITHIN_ACTION_ERROR.
```

The frozen verdict is therefore:

> **DERIVED COMPUTATIONAL ON THE FROZEN GRID:** all thirteen collective
> points in both schedules are stationary within the smaller-step action
> error and the four-soft-mode scale.  A continuous stationary family remains
> **PATTERN**, not an analytic interval theorem.

## 2. Numerical scale

Across each schedule and all thirteen values `t in [-0.1,0.1]`:

```text
sixth-order log-gradient norm        4.766e-19 ... 7.859e-19
empirical gradient error             2.272e-14 ... 3.067e-14
quotient correction proxy            1.750e-17 ... 2.886e-17
correction-proxy error                1.469e-12 ... 1.982e-12
collective scalar                     1.731e-23 ... 2.337e-23
maximum imaginary contamination      below 1.14e-92
```

The reported gradients are only `2.56e-5` of their conservative empirical
error at worst.  Thus the correct reading is “consistent with zero extremely
deep inside the error envelope,” not a measurement of a nonzero
`1e-19` gradient.

All 5,460 displaced action geometries remained Lorentzian and away from
branch boundaries.

## 3. The precision correction behaved as predicted

The old derivative steps were `2e-4,1e-4,5e-5`; the new steps were ten times
smaller.  The median old/new empirical-error ratio is

```text
9.970054e3,
```

with range `9.9669e3 ... 9.9729e3`.  This is the expected fourth-order
Richardson improvement of approximately `1e4`, obtained without changing a
scientific threshold.

The even/odd gradient-row differences are only

```text
1.80e-20 ... 2.97e-20,
```

or about `4e-7` of their combined empirical envelopes.  The apparent odd
schedule failure in the preceding binary64 transverse solve was therefore a
numerical artifact, not schedule-dependent physics.

## 4. Physical interpretation

- **DERIVED COMPUTATIONAL ON GRID:** the regular-boundary dust sandwich has
  no selected collective duration at any of the thirteen tested lapse values.
- **PATTERN:** the whole displayed interval is a continuous stationary lapse
  family.
- **STRUCTURAL / KNOWN:** lapse freedom at a symmetric Regge slice is not a
  new gravitational mechanism.
- **REFUTED:** the previous claimed odd correction of order `1e-2`; it came
  from projecting raw residuals onto a logarithmic Hessian.
- **NOT DERIVED:** a physical tick, speed limit, Planck time or Planck mass.

This strengthens the negative statement about time selection: `tau=0.0102`
is input coordinate/lapse data on the regular sandwich, not a duration chosen
by the geometry.

## 5. What remains open

Thirteen samples, however precise, do not prove an identity for every real
`t`.  The next exact question is whether the local equations can be reduced
symbolically on

```text
q+rho=l0^2
```

so that their `rho` dependence cancels after the dust balance and 600-cell
incidence counts are inserted.  A positive symbolic identity would upgrade
the continuous family from **PATTERN** to **DERIVED**.

Even that would cover only the regular boundary.  The physically decisive
downstream test remains nonlinear continuation under a preregistered
zero-sum boundary deformation, where the lapse may cease to be gauge and be
fixed by a pseudo-constraint.

## 6. Status ledger

| Claim | Status |
|---|---|
| All 13 points stationary within weak-scale action error | **DERIVED COMPUTATIONAL, both parities** |
| Even and odd schedules agree on those points | **DERIVED COMPUTATIONAL within error** |
| Every point of the continuous interval is stationary | **PATTERN / analytic proof open** |
| `tau` is selected on the regular sandwich | **REFUTED** |
| Precision quotient and 29 zero-sum tangent responses remain valid | **DERIVED COMPUTATIONAL LINEAR** |
| Those tangent responses integrate to finite solutions | **OPEN** |
| Collective lapse remains gauge after shape deformation | **OPEN** |
| Full 840-edge carrier or multiple slabs | **NOT TESTED** |
