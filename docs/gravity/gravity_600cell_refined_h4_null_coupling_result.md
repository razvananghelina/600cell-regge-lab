# Refined H4 product-lapse null coupling: consolidated result

Date: 2026-08-21

## Verdict

The one-dimensional kernel of the refined internal Hessian is the tangent to
the exact product-duration family, but it is **not** a null direction of the
full boundary-plus-internal Hessian.  It couples to the boundary through the
nonzero covector

```text
c = H_bi n = (1/2) g_boundary.
```

Consequently the linear internal equation

```text
H_ib delta_b + H_ii delta_i = 0
```

is solvable only on the eleven-dimensional boundary hyperplane

```text
c^T delta_b = 0.                                      (1)
```

This closes the interpretation of the singular direction as an unconstrained
gauge freedom.  It does not yet establish a continuum Hamiltonian constraint,
a propagator, a physical tick, `c`, `G`, Planck units or particle physics.

## Frozen provenance

| stage | commit |
|---|---|
| prior-art gate | `32854ba` |
| primary protocol | `752891b` |
| primary verifier registered | `fab5f4d` |
| primary result recorded | `fd0e5f3` |
| adversarial protocol | `9df1247` |
| adversarial verifier registered | `96c7119` |
| first formal disagreement preserved | `f89d2d0` |
| threshold correction preregistered | `a92ab1f` |
| first incomplete implementation | `24560c8` |
| incomplete correction preserved | `8c79422` |
| final one-line correction | `df85637` |

The primary verifier passed `16/16` twice with byte-identical artifact

```text
reproducible/gravity_600cell_refined_h4_null_coupling.json
SHA-256 6b6fbd95b07f365b3fcac332fa3546021e8d756a510af0184bc974e52d5efa79.
```

The corrected adversarial verifier passed `11/11` twice with byte-identical
artifact

```text
reproducible/gravity_600cell_refined_h4_null_coupling_adversarial.json
SHA-256 5c1f596958f9d878c8d9d3ccb6ecc8359f72164e8f36dd9930fb71ddc1351ce9.
```

Only the two targeted verifiers and a static registry audit were used.  The
full suite and the deferred nested root census were not run.

## Primary calculation

Use `u=log(tau^2)`.  The geometry-derived, unnormalised internal tangent is

```text
n_cross,rs = -tau0^2/q_cross,rs,
n_rho,r    = 1.
```

At all `24` staircase schedules and all `72` tested product-family points
`tau/tau0 in {1/2,1,2}`, the maximum internal residual is `1.182e-76`.
Directional differentiation gives `H_ii n=0`; the maximum computed internal
image is `4.108e-44`, inside the maximum propagated numerical envelope
`9.740e-41`.  Together with the independently frozen internal inertia
`(9,1,0)`, this identifies the complete one-dimensional internal null line.

The boundary coupling is nonzero.  Its common twelve-component row is

```text
(0.0182753971831513703513...,
 0.00121804711711642586920...,
 0.0000327518273255248437897...,
 0.0366773216171774596419965...,
 0.000494989505267934316830...,
 0.0191118520096157871206302...,
 repeated on the new boundary layer).
```

All `24` rows agree within their frozen envelopes and have rank one.  The
maximum schedule difference is `9.173e-44`; the maximum old/new reversal
difference is `2.184e-125`.  Componentwise the result agrees with
`(1/2)g_boundary` to `6.493e-44`.  The smallest nonzero-coupling resolution
ratio is `3.766e38`, so the decision is not threshold-marginal.

## Mechanically independent reconstruction

The adversarial verifier does not import or execute the primary verifier, the
full-Hessian verifier or the Lorentzian action evaluator.  Starting from the
actual spatial incidences, it rebuilds each rank-pair curvature `C_rs` and
derives locally

```text
A^2                           = -l^2 tau^2/4,
dA/d log(l^2)                 = i l tau/4,
g_boundary,rs                 = tau C_rs/4,
c_rs=d g/d log(tau^2)         = tau C_rs/8.
```

It also verifies the vertical product-null coefficient exactly from
`m_r=K_r/(8*pi)`.  The independently reconstructed coupling differs from the
primary finite-difference rows by at most `6.493e-44`, only
`0.000666667` of the corresponding frozen primary envelope.  Boundary
reconstruction error is `3.5e-77`; the vertical coefficient and layer-swap
errors are exactly zero.  Wrong factor, sign, dropped-component and swapped
pair controls all fail by resolved amounts.

**Scope limitation:** the adversarial route independently reconstructs the
coupling covector and the vertical curvature-balance identity.  It does not
independently rebuild the six cross-coordinate components of `H_ii n=0`;
those remain supported by the primary finite-family and directional tests.

## Threshold-correction ledger

The first adversarial protocol demanded component agreement below `1e-68`
although the already frozen primary finite-difference rows were known only to
roughly `6e-41`--`10e-41`.  Its formal outcome was therefore disagreement even
though the observed error, `6.493e-44`, lay about three orders of magnitude
inside every relevant envelope.  That result is preserved in `f89d2d0`.

The correction was preregistered before rerunning: compare each row with its
stored primary envelope and require error/envelope `<1`, while retaining the
`1e-68` gate for the independently high-precision boundary reconstruction.
The first code correction missed the old absolute threshold in the auxiliary
rank-one predicate; that failed run is preserved in `8c79422`.  Commit
`df85637` changes only that predicate.  The final double run is the accepted
adversarial result above.

## Status ledger

- **DERIVED COMPUTATIONAL / STRUCTURAL, adversarially corroborated:** the
  product-lapse tangent is the unique internal null line, and its boundary
  coupling is the schedule-independent rank-one row
  `c=tau*C/8=(1/2)g_boundary`.
- **DERIVED STRUCTURAL:** equation (1) is a necessary boundary compatibility
  condition in the frozen invariant log-squared-edge coordinates.
- **DERIVED NEGATIVE:** the singular direction is not an unconstrained
  full-Hessian gauge direction.  An ordinary Schur complement or an arbitrary
  Moore--Penrose replacement is not justified.
- **OPEN:** the constrained effective quadratic form on `ker(c^T)`, after
  quotienting the internal null line, and whether that form is identical for
  all 24 schedules.
- **OPEN:** whether the compatibility row has the continuum interpretation of
  a Hamiltonian constraint.
- **OPEN:** external novelty; no post-result literature claim is made here.
- **NOT COMPUTED:** nonhomogeneous propagation, a dispersion relation, a tick,
  `c`, `G`, Planck units and particle masses.

## Next falsifiable gate

Construct the effective boundary bilinear form only on `ker(c^T)`, eliminating
the nine nonsingular internal directions and quotienting the one-dimensional
internal null line.  The comparison across the 24 schedules must be
basis-independent on that hyperplane.  A positive result would establish only
a schedule-independent constrained `H4` quadratic form; because this is still
the invariant sector, it cannot by itself determine a propagation speed.

