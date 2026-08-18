# Protocol: can the licensed spectral data select the gravitational metric?

Date: 2026-08-12

## Provenance

This is a post-recognition hostile audit.  Before this protocol was written it
was already known that symmetry and exact Whitney refinement leave the
continuum

```text
g_u=(1-u)g_R+u g_0,   0<=u<=1,
```

and that the repository contains a finite heat trace, exact finite moments and
the continuum round de Rham curvature coefficient.  It was also recognized
that a heat trace without its time and a spectral action without its cutoff
function do not by themselves define a variational metric selector.  The
tests below freeze the precise obstruction and actively test the possible
counterclaim that the exact ratio `c1/(2 c0)=31/11` supplies the missing scale.
This is not a blind discovery protocol.

No measured gravitational target, preferred endpoint, fitted cutoff function
or fitted heat time may enter the audit.

## Complete scope and hypotheses

1. The metric configuration space under consideration is exactly the already
   certified round--Regge family above.  Here `g_R` is the pushed-forward flat
   facet metric and `g_0` is the unit round metric on `S3`.

2. A **licensed spectral functional** is one whose operator, trace, function,
   scale and metric dependence are already defined in an active authoritative
   source.  Historical `exp*.py` scans and formulas explicitly labelled open,
   pattern or continuum template do not license new dynamics.

3. The source census is restricted to the active construction:

   - the 2640-state incidence Kähler--Dirac matrix with Euclidean cochain
     adjoint and its moments `Tr(I)`, `Tr(D^2)`, `Tr(D^4)/2`;
   - its heat family `Tr exp(-t D^2)`;
   - the exact Whitney operators with mass matrices on the fixed Regge metric;
   - the smooth continuum de Rham heat coefficient at the round metric.

4. A metric is **dynamically selected** only if one already licensed scalar
   functional `S_n(u)` is defined on the complete exterior-cochain carrier at
   every registered refinement level, including both endpoints, and has the
   same unique minimizer over its entire already allowed cutoff/heat-scale
   range.  Choosing a time, sign, moment combination or endpoint after seeing
   its value is forbidden.

5. Scalar `P1` spectra on the projected barycentric carrier are an admissible
   diagnostic but cannot certify the complete Kähler--Dirac spectral action:
   degrees 1, 2 and 3, the exact rank-edgewise tower and exact curved mass
   matrices must be present for a positive selection result.

6. Smooth and Regge/conical heat coefficients may not be identified by
   analogy.  A common singular heat theorem for the full exterior operator is
   required before their coefficients can be compared as one action.

## Frozen exact tests

1. Audit the active sources and record, without inference, which of the
   operator, function, scale and variable-metric domain each one supplies.

2. Prove the metric-scale/heat-time covariance.  For

   ```text
   g' = c^2 g,
   ```

   the Hodge Laplacian obeys `Delta_g'=c^-2 Delta_g`, and hence

   ```text
   K_g'(t)=K_g(t/c^2).
   ```

   Therefore an unnormalized heat time cannot independently fix an overall
   metric scale.

3. Test the proposed finite ratio

   ```text
   r(D)=Tr(D^2)/(2 Tr(I))=31/11.
   ```

   Under `D^2 -> c^-2 D^2`, prove `r -> c^-2 r`.  Determine whether either
   `r` or its reciprocal defines a scale independently of the same operator
   normalization it is meant to vary.

4. Give an exact finite spectral control showing that heat-trace ordering can
   reverse with heat time.  Use the preregistered spectra

   ```text
   A={0,1,10},   B={0,2,3}.
   ```

   The comparison must use exact rational values `x=exp(-t)` at two fixed
   points, not an optimized numerical scan.

5. Verify that positive even moments and a positive heat cutoff do not supply
   an interior overall scale without an additional normalization.  Distinguish
   this scaling statement from shape selection at fixed volume.

6. Audit whether the existing scalar projected-refinement code supplies a
   complete exterior operator and whether the existing exact Whitney tower
   supplies the round endpoint.  A search space unable to represent both
   claims may not be cited as evidence for either.

7. Audit the mathematical category at the endpoints: `g_0` is smooth with
   distributed scalar curvature, while `g_R` is piecewise flat with curvature
   on a skeleton.  Record whether the repository contains the singular
   all-form heat coefficients needed for a common comparison.

## Decision boundary

- **SELECTED METRIC:** an already licensed complete functional has a unique,
  cutoff-independent and refinement-stable minimizer in `u`.
- **DERIVED CURRENT-CONSTRUCTION NO-GO:** every licensed complete functional
  is either fixed-metric, lacks a metric-dependent operator, or retains a
  free function/scale whose allowed variation can alter the preference.
- **ILL-POSED NUMERICAL GATE:** the two endpoints do not yet share one complete
  metric-dependent Kähler--Dirac family or one applicable heat expansion.  In
  that case no numerical `u` will be reported, because it would compare
  different constructions rather than evaluate one action.

The latter two outcomes do not prove that emergent gravity is impossible.
They prove that the present repository has not yet supplied the dynamical
metric selector.  A future full variable-metric operator and independently
fixed action remain admissible new constructions.

