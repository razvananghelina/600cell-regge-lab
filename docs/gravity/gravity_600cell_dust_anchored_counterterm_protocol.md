# Preregistration: anchored counterterm corollary on pure-momentum rays

Date: 2026-08-17

Prior-art commit: `20e9a26`.

Status: frozen before parsing configuration-only defects of the 16 committed
pure-momentum cases.

## 1. Frozen inputs

Require

```text
nonlinear case seed
2104c69ba6b21d3a3d92c7071d7f2702cb7d33f7f0e3ff17954f64c469f0c01d

committed nonlinear result
a1e00071fa41f986dfaee84ea6e7689a14c50823f6c87d76889e6cb9346a7e3f.
```

Require the upstream result to report `8/8`, 32 cases and
`BROKEN: 32`.  Read only its already committed outputs, solver corrections,
case metadata and the seed's unique physical-edge permutation.  Do not
re-evaluate the action or alter a nonlinear solve.

## 2. Frozen census

Select mechanically every case with

```text
sector = MOMENTUM.
```

Require exactly

```text
4 directions x 2 signs x 2 levels = 16 cases,
levels = {0.5,1.0}.
```

Require each selected input ray to have zero configuration component by its
frozen metadata.  No case may be selected or removed by output size.

## 3. Configuration-only defect

For each parity and calibration variant read the first 30 components of

```text
output=(log q_new[30],p_post[30]/p_star).
```

Map the even configuration to odd coordinates with the already selected
physical-edge permutation.  Define

```text
d_q = ||q_new,odd,operational - P q_new,even,operational||_2.
```

Define a conservative empirical uncertainty

```text
u_q = ||q_even,operational-q_even,validation||_2
    + ||q_odd,operational-q_odd,validation||_2
    + sum of all four stored full-output correction norms
    + 1e-70.
```

Using a full-output correction norm in a configuration-only comparison can
only enlarge this proxy.

Classify each case:

- `ANCHORED_CONSISTENT` if `d_q <= 10*u_q`;
- `ANCHORED_REFUTED` if `d_q > 100*u_q`;
- `OPEN` otherwise;
- `OPEN_SOLVE` if any of its four upstream solves was not successful.

## 4. Scaling diagnostic

For each direction and sign, if both levels are `ANCHORED_REFUTED`, report

```text
order = log2(d_q,full/d_q,half).
```

Label `[1.5,2.5]` as `QUADRATIC_COMPATIBLE`; otherwise label the resolved
order without interpreting it as a force law.  Scaling does not enter the
main outcome.

## 5. Mechanical outcome

- any `ANCHORED_REFUTED`:
  `ANCHORED_ENDPOINT_COUNTERTERM_REFUTED_ON_FROZEN_RAYS`;
- no refutation but at least one `OPEN` or `OPEN_SOLVE`:
  `ANCHORED_ENDPOINT_COUNTERTERM_OPEN`;
- all 16 `ANCHORED_CONSISTENT`:
  `ANCHORED_ENDPOINT_COUNTERTERM_CONSISTENT_ON_FROZEN_RAYS`;
- provenance/census/implementation failure:
  `ANCHORED_ENDPOINT_COUNTERTERM_CONTROL_FAILED`.

A refutation is a **DERIVED COMPUTATIONAL COROLLARY** under the complete
anchoring hypotheses.  It is not a refutation of unanchored endpoint shifts,
general canonical transformations, improved/perfect actions or the full
carrier.  Finite consistency is only a **PATTERN**.

Only the new lightweight verifier will run.  The full suite and the expensive
action verifier will not run.
