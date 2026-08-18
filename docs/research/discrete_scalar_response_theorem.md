# Discrete Scalar-Response Theorem

Working note for the `gamma_PPN` repair before editing the paper text.

## Exact theorem

On the 600-cell simplicial complex, let

- `d_0 : R^V -> R^E`
- `d_1 : R^E -> R^F`
- `Delta_0 = d_0^T d_0`
- `B = d_0 d_0^T`.

Then the Moore-Penrose pseudoinverses satisfy

- `B^+ d_0 = d_0 Delta_0^+`.

For a point source at vertex `v`, define

- `delta_v in R^V`
- `Phi_v = Delta_0^+ delta_v`
- `h_v = B^+ d_0 delta_v`.

Then exactly:

- `h_v = d_0 Phi_v`
- `d_1 h_v = 0`
- `d_0^+ h_v = Phi_v`
- `Delta_0 Phi_v = P_0 delta_v`

where `P_0 = I - (1/N) J` is the projection to the zero-mean subspace.

## Meaning

This proves that the static edge response produced by a vertex source is fully
encoded by a single scalar potential. The response lives in the exact sector
and has no independent coexact contribution.

So the safe statement is:

- the tested discrete model has no independent static coexact scalar response;
- the potential recovered from the edge field is exactly the same scalar
  potential that generated it.

## What this does not prove

This is not yet a derivation of the continuum PPN parameter `gamma` in the
standard 4D post-Newtonian sense. In particular, it does not by itself provide:

- a metric ansatz with continuum `Psi` and `Phi`
- the standard PPN field equations
- a derivation of `gamma = Psi/Phi` for a 4D relativistic metric theory

So the paper should talk about:

- `gamma_disc = 1` as a discrete scalar-response analogue,

not about a proved continuum `gamma_PPN = 1`.
