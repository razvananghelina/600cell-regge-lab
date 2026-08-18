# Preregistration: gauge-quotient boundary response

Date: 2026-08-13

Prior-art gate commit: `ff8f404`

Exact collective lapse result: `790fc7f`

Five-direction Schur result: `dc927a5`

Original mixed-Jacobian record: `7d5e9fc`

Status: **frozen before contracting the exact lapse tangent with the
internal--final-boundary mixed block**.

No mixed compatibility row, constrained boundary direction, quotient response
or response spectrum has been inspected before this protocol.

## 1. Frozen linear problem

At the published stationary dust sandwich, use logarithmic coordinates

```text
u in R^35  : 30 staircase diagonals + 5 positive pole magnitudes,
v in R^30  : final-boundary squared lengths.
```

The internal logarithmic equations are

```text
E = (x/24) partial S_total/partial x.
```

Their linearization is

```text
H delta_u + B delta_v = 0,

H = partial E/partial u,   shape 35 x 35,
B = partial E/partial v,   shape 35 x 30.
```

The exact collective path gives the unfitted raw tangent

```text
w_raw = (-tau^2/(l0^2-tau^2) repeated 30 times,
          1                         repeated  5 times),
w = w_raw / norm(w_raw).
```

Use this analytic vector, not the weakest numerical eigenvector.

## 2. Frozen quotient basis and regularity audit

Construct a deterministic Householder matrix that maps `w` to the first
coordinate vector.  Fix its sign by the standard vector `w+sign(w_0)e_0` and
take columns 1 through 34 as the orthonormal quotient basis `Q`.  Require

```text
Q^T Q = I_34,   Q^T w = 0.
```

Load the recorded Richardson `H`, symmetrize it, and remove measured gauge
leakage by the already justified projection

```text
P = I-w w^T,
H_perp = P H P,
H_Q = Q^T H_perp Q.
```

This projection is permitted only because commit `790fc7f` independently
established the exact null path.  Report rather than hide `norm(H w)` before
projection.

The quotient is classified `QUOTIENT_REGULAR` only if all of the following
upstream and current gates agree:

- the recorded `30 x 30` staircase block has rank 30 at relative thresholds
  `1e-7,1e-9,1e-11`;
- all four already frozen relative Schur eigenvalues exceed
  `100*epsilon_5`;
- `H_Q` has 34 singular values above the absolute threshold `1e-9`;
- its four smallest singular values agree with the frozen four relative
  eigenvalues to normalized error below `3e-3`.

Otherwise label it `QUOTIENT_UNRESOLVED`.  Full quotient rank is a scientific
outcome, not a verifier PASS target.

## 3. Frozen high-precision boundary compatibility row

The solvability condition is

```text
c delta_v = 0,    c = w^T B.
```

Compute `c` independently from the complete action at 90 decimals.  For each
of the thirty printed final-boundary orbit coordinates `j`, evaluate the
mixed rectangle

```text
M_h[j] = [S(t,+h e_j)-S(t,-h e_j)
          -S(-t,+h e_j)+S(-t,-h e_j)]/(4 h^2),
```

where

```text
rho(t)=tau^2 exp(t),   q(t)=l0^2-rho(t),
f_j(v)=l0^2 exp(v),
```

at exactly `h=5e-4` and `h=2.5e-4`.  Then

```text
c_fine = M_fine/(24 norm(w_raw)),
c_R    = (4 M_fine-M_coarse)/(3*24 norm(w_raw)),
epsilon_c = norm_2(c_R-c_fine),
epsilon_floor = max(epsilon_c,1e-35).
```

All 240 displaced action geometries per parity must remain Lorentzian and off
branch boundaries under the certified binary64 evaluator.  Report all action
imaginary parts.

Also contract the previously recorded analytic `B` with `w` and compare it
with `c_R`, normalized by `max(1,norm(B))`; frozen agreement tolerance is
`3e-6`.  The high-precision action row controls the scientific label.

Assign exactly one compatibility label:

- `ALL_30_BOUNDARY_DIRECTIONS_COMPATIBLE` if
  `norm(c_R) <= 10*epsilon_floor`;
- `ONE_LINEAR_BOUNDARY_CONSTRAINT` if
  `norm(c_R) > 100*epsilon_floor`;
- `BOUNDARY_COMPATIBILITY_UNRESOLVED` otherwise.

No component or boundary combination is selected after seeing the row.

## 4. Frozen response construction

If the quotient is regular:

- for all-thirty compatibility, use `Y=I_30`;
- for one constraint, construct a deterministic Householder basis `Y` for
  `ker(c_R)`, with the same sign convention, giving 29 columns;
- for unresolved compatibility, do not construct a response.

Solve

```text
H_Q A = -Q^T P B Y,
delta_u = Q A delta_z.
```

Report:

- the quotient solve residual;
- singular spectrum and relative ranks of the `35 x dim(Y)` response;
- response norm and condition;
- unprojected linear residual and its component along `w`;
- both phase-parity differences.

The solve is an implementation PASS only when its quotient residual is below
`1e-5`; rank or amplification of the response is a scientific result.

## 5. Frozen outcome labels

Combine the two independent classifications:

- `REGULAR_QUOTIENT_ALL_30_LINEAR_RESPONSES`;
- `REGULAR_QUOTIENT_29_RESPONSES_PLUS_ONE_BOUNDARY_CONSTRAINT`;
- `QUOTIENT_OR_COMPATIBILITY_UNRESOLVED`.

The first two advance the route only at linear order.

## 6. Claim boundary

Even the strongest outcome does not yet invoke the nonlinear
implicit-function theorem for all thirty boundary directions.  We know an
exact gauge family only on the regular boundary.  Away from it, the remaining
equation may become a higher-order pseudo-constraint.  This audit establishes
only the tangent response modulo collective lapse.

It does not identify gauge-invariant gravitons, select a physical boundary
perturbation, derive a clock, cover the unreduced 840-edge space or propagate
multiple slabs.  A nonlinear displaced-root test must be separately
preregistered after this result, using deterministic boundary directions and
checking all 35 equations rather than only 34 selected ones.
