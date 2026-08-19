# Preregistration: corrected 120-column strut carrier, before target comparison

Date: 2026-08-19

Status: freeze this protocol before constructing the first accepted-slab
carrier matrix.  This phase is forbidden from loading any Hessian, Schur,
tangent, singular-vector or eigenvector artifact.

## Frozen provenance

| input | SHA-256 |
|---|---|
| prior-art gate | `e0064e73d161f7ba64b5a5c0c14ace4276cd0ce21d22fcb04b329506132064ca` |
| exact local-lift source | `4389861a4b64d043325e0661ae9c2340f61e5c8eb50399c9fd2083a360dadbc1` |
| exact local-lift artifact | `0a569e48189c56bc081efcee33f7826fedd52afb93b6135ddb2fec385b56fbdf` |
| exact local-lift result | `646972a19450f1734ef522cb0b9693cc809b19d7895eb21823b20332a958d56d` |
| frustum-equivalence source | `99f47f0cfc70d2c0784d002cc08898e29f28a53e51930e6683c95629af128587` |
| frustum-equivalence artifact | `7e7c23efaf24a2c99a68f3b302b9ef575e0f777ef46f73ccaea9f99e1ecd58dc` |
| frustum-equivalence result | `b63808c260f12711ab25bdc72414f36c6c0f89f9420619c448c875d2dac7b093` |
| complete one-slab geometry source | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| accepted homothetic tick | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |

The verifier may read exactly these inputs.  In particular, it may not read
the full-lapse Schur or full-boundary tangent artifacts.

## 1. Independent symbolic theorem

Use independent symbols `lambda`, `a_parallel`, `a_square` and an arbitrary
positive bottom-edge scale.  Construct directly

```text
s_u = a^2,
s_v = (a+(lambda-1)e)^2,
d_uv = (a+lambda e)^2.
```

Differentiate at fixed boundary scales and eliminate `delta(a.e)`.  Require
SymPy to reduce identically to

```text
(lambda-1) delta d_uv + delta s_u-lambda delta s_v = 0.
```

The calculation must reject `lambda=1`; no pseudoinverse or limiting value is
allowed.

## 2. Mechanically different rational controls

Read the two committed baseline local blocks at

```text
(lambda,tau)=(2,5),(3,11).
```

Independently rebuild the `10 x 16` upper-edge-plus-strut Jacobian, its
deterministic pivot right inverse and six-dimensional exact kernel.  Do not
call the local-lift verifier.  Reconstruct the physical `16 x 8` response
from each serialized `6 x 8` block and require:

1. all ten declared natural-length derivatives are exact;
2. for all twelve ordered cross diagonals and all four pure-strut columns,
   the response equals

   ```text
   (-delta s_u+lambda delta s_v)/(lambda-1);
   ```

3. all four pole responses reproduce the identity.

This is an adversarial control: the symbolic trapezoid proof does not use the
face-gluing block, while the block reconstruction does not assume the
trapezoid result.

## 3. Accepted-slab carrier

Reconstruct the fixed 600-cell slab directly and require, for both schedule
parities,

```text
720 old boundary edges,
720 new boundary edges,
720 oriented cross diagonals,
120 poles,
24 stabilizer elements.
```

Use the accepted state `(s,r)` only through

```text
lambda = exp(s),
rho = rho0 exp(r),
q_diag = lambda L0^2-rho.
```

Require `lambda>0`, `lambda!=1`, `rho>0`, `q_diag>0`.  In lexicographically
sorted internal-edge order followed by sorted final-boundary-edge order,
construct `G` with columns indexed by the 120 logical vertices:

```text
pole (v,v+120):                 G[row,v] = 1,
diagonal (u,v+120):             G[row,u] =  kappa,
                                     G[row,v] = -lambda kappa,
new-boundary edge:              zero,

kappa = rho/((lambda-1) q_diag).
```

No coefficient may be chosen from a dynamic operator.

## 4. Exact incidence and equivariance gates

Before any floating spectrum, require by edge labels and coefficient roles:

1. every pole row has support one, every diagonal row support two, and every
   final-boundary row support zero;
2. every column has one pole plus the twelve oriented diagonals incident on
   that vertex;
3. `rank(G)=120`, certified by the literal pole identity block rather than a
   tolerance;
4. every one of the 24 schedule-stabilizer permutations intertwines rows and
   vertex columns, including the source/target coefficient roles;
5. the sum of all 120 columns is exactly the analytic collective-lapse
   column: one on poles, `-rho/q_diag` on every diagonal, zero on the final
   boundary;
6. the even and odd carriers use their own oriented diagonal sets; no
   post-hoc orbit matching is allowed.

As a negative control, delete the source-endpoint coefficient from the
lexicographically first diagonal.  The corrupted matrix must fail both the
collective-column identity and at least one stabilizer intertwining check.

## 5. Target-blind numerical census

At 100 decimal digits, and independently after conversion to binary64,
record for each parity:

- all 120 singular values of `G` in descending order;
- numerical rank under explicitly printed precision-scaled thresholds;
- condition number;
- the gain `||G u||` on the normalized uniform vertex vector `u`;
- the 119 eigenvalues/square-root gains of the quadratic form
  `Q^T G^T G Q` on its Euclidean orthogonal complement, where `Q` is a
  deterministic Householder basis fixed before evaluation;
- the coupling norm `||Q^T G^T G u||`, so the uniform line and its complement
  are not silently treated as invariant when they are not;
- even/odd distances between the two ordered singular multisets.

These are intrinsic carrier diagnostics, not a comparison with any dynamic
target.  No desired count or target value is an acceptance condition.

## 6. Look-elsewhere and artifact discipline

The artifact must serialize the complete sparse row supports and coefficients
of both carriers, their full singular multisets and every control above.  Its
classification is committed before a second protocol may load a dynamic
target.

The candidate count is exactly one carrier per schedule parity: equation (2)
of the prior-art gate fixes it.  Rescaling individual columns, mixing them,
selecting a graph orientation, or adding scale columns after seeing a target
is forbidden.

## 7. Mechanical outcomes

Assign in this order:

1. `CORRECTED_STRUT_CARRIER_CONTROL_FAILED` if provenance, geometry,
   symbolic identity, rational-block reconstruction, incidence,
   equivariance, collective or corruption controls fail;
2. `CORRECTED_STRUT_CARRIER_NUMERICALLY_OPEN` if the exact controls pass but
   the 100-digit and binary64 intrinsic singular censuses disagree beyond
   their declared roundoff gate;
3. `CORRECTED_STRUT_CARRIER_FROZEN` if all controls pass and the complete
   target-blind artifact is written.

No outcome in this phase identifies a Schur sector, a hyperbolic mode, gauge,
a pseudo-constraint, a graviton, a physical instability, a tick, `c`, `G` or
Planck units.

## 8. Execution discipline

Register and commit the verifier before its first execution.  Run only that
targeted verifier and the static registry guard.  Freeze its first artifact
in a separate commit.  Only then write and commit the target-comparison
protocol.
