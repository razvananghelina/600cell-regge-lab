# Preserved second run failure: audited converter namespace

Date: 2026-08-22

History-precision correction commit: `00236cb`.

Status: **IMPLEMENTATION FAILURE PRESERVED; NO TANGENT RESULT ACCEPTED.**

## Controls reached before failure

The second targeted execution passed:

- every frozen hash and the unique registry entry;
- the accepted input outcomes;
- the corrected branch-B history control, with maximum high-precision
  committed error `1.77807e-68`;
- the exact conformal-symplectic scale identities;
- all 43 imported geometry certificates;
- literal equality of the `r1` and `r2` representation carriers in both
  parities, with zero basis distance at 180 digits;
- the even normalized and direct-physical Lorentzian derivative controls;
- the preliminary even kernel scale classification
  `SCALE_LIFT_CONFIRMED`, with maximum imaginary remainder about
  `4.06445e-153`.

These values are incomplete intermediate controls.  They are not an accepted
second-tangent result.

## Failure

At the first `d=3` sector, after the pre-Legendre solve, the wrapper
`ball_records_from_built` attempted to call `acb_midpoint_and_radii` as a
module global.  The audited function had deliberately been loaded into the
isolated `core` namespace, so Python raised

```text
NameError: name 'acb_midpoint_and_radii' is not defined
```

No JSON or NPZ artifact was produced.

## Frozen mechanical correction

Change the wrapper signature to accept the converter explicitly and call it
as

```text
ball_records_from_built(built_by_level, core["acb_midpoint_and_radii"])
```

No equation, matrix, basis, derivative, precision, tolerance, classifier or
outcome condition changes.
