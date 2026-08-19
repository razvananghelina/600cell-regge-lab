# Universal local canonical-data lift

Date: 2026-08-19

## Frozen provenance

- prior-art gate: commit `fae9940`;
- target-disclosed protocol: commit `dd302d8`;
- first registered implementation: commit `ac5f80f`;
- first reconciliation-control failure: commit `0ed3880`;
- preregistered coordinate correction: commit `592266a`;
- first complete artifact: commit `6967ed7`;
- artifact SHA-256:
  `0a569e48189c56bc081efcee33f7826fedd52afb93b6135ddb2fec385b56fbdf`.

The first execution wrote no artifact and printed no block coefficient.  It
stopped because the control confused normal-displacement amplitude `nu` with
raw strut squared-length variation `s`.  The exact basis change

```text
s_v  = 6 (lambda - 1) sigma_v - 2 tau nu_v,
nu_v = 3 (lambda - 1) sigma_v / tau - s_v / (2 tau)
```

was recorded before the corrected rerun.  The 48-unknown face constraints,
support target, corruption, and outcome hierarchy were unchanged.

The corrected targeted verifier passed 13/13.  No full-suite run was
performed.

## Independent construction

For each legitimate exact construction, one unknown 6-by-8 block `X` was
used identically on all 600 tetrahedral cells:

```text
cell flex = X (sigma_v0,...,sigma_v3, s_v0,...,s_v3)^T.
```

Substitution into every global residual coefficient produced 51,320 exact
affine constraints on only 48 unknowns.  This did not call or copy the earlier
3600-variable global elimination.

Every construction gave:

```text
affine rank                    48 / 48
candidate obstructions             0
nonzero entries of X             48 / 48
nonzero global residual rows        0 / 6000
corrupted-image obstructions        > 0
```

The complete rational blocks and their canonical hashes are stored in the
JSON artifact.  Reversing all face orientations leaves the baseline block
exactly unchanged.  Odd local relabelling and overall metric-sign reversal
produce the corresponding exact sign/permutation changes.

## Exact locality

Because all 48 local coefficients are nonzero, direct membership checks give

```text
support(each cell-flex row)
    = the four sigma and four strut coordinates at that cell's vertices;

support(each sigma_v or s_v column)
    = exactly the 20 tetrahedral cells containing v.
```

Each data column has exactly `20 * 6 = 120` nonzero flex coefficients.  These
are exact set equalities, not support-count coincidences.

**DERIVED.** The unique rational global lift is the repetition, in local
vertex order, of a universal exact 6-by-8 cell block.  It is strictly local
on 600-cell vertex stars.

## Reconciliation with the failed old formula

After the exact `(sigma,nu)` to `(sigma,s)` basis change, both the old and new
physical displacement blocks induce the same ten local squared-length data.
For every construction:

```text
new flex block != old flex block;
rank(new flex block - old flex block) = 3;
all 48 entries of that difference are nonzero;
Jacobian * (new physical response - old physical response) = 0.
```

Thus their difference is an exact three-dimensional family inside the local
six-dimensional Poincare kernel.  The old hand-chosen representative fails
global face gluing; the complete geometry selects a different local Poincare
representative.  This explains the old negative without retracting it.

For both rational representatives, baseline and alternate-right-inverse
coordinate blocks differ, while the reconstructed physical displacement
responses agree exactly.  Hence the result is not an artifact of the chosen
right-inverse graph.

## Scientific status

**DERIVED.** The first-order kinematic carrier and its local response are now
exact over `Q`, independently reproduced, convention-stable, and equipped
with a discriminating corrupted-image control.

**STRUCTURAL.** Local vertex-star response is compatible with established
discrete-conformal and flat-background vertex-displacement mechanisms.  The
external novelty of this explicit 600-cell block is **OPEN**.

**OPEN.** The result supplies no action, canonical momenta, symplectic form,
constraint/gauge split, propagating tensor mode, clock, tick, `c`, `G`, or
Planck scale.  In particular, arbitrary kinematically admissible strut data
are not yet physical time.

## Next falsifiable calculation

Use the local block to pull the exact Regge boundary Hessian back to the 240
canonical data coordinates.  Before computing its spectrum, preregister:

1. the complete action and boundary-term convention;
2. which variables are varied and which background equations are required;
3. the expected gauge/null directions without using observed eigenvalues;
4. rank and small-eigenvalue decisions in exact or certified arithmetic;
5. the separation between the 120 conformal scale coordinates and the 120
   strut/lapse coordinates.

That is the first calculation capable of deciding whether the arbitrary
struts remain gauge, become constraints, or acquire dynamical coupling.

