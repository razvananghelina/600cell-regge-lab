# Minimal Edge-Lift Theorem

## Main improvement

The previous bounded-search route for the `Z[phi]` quantum numbers was:

`C1 + C2 => unique assignment`  (paper version)

then:

`C1 + Q + C2 => unique assignment`

then:

`C1 + S => unique assignment`

where `S` is the simple-edge condition

`dz in {0} union {+/-phi^r} union {+/-2phi^r}`.

The new point is that, for the neutral branch ordering used in the
bounded-search verifiers, `S` is not an independent extra axiom.

It is exactly the canonical **minimal edge lift** of the edge exponent jumps.


## Setup

Use the neutral node order

`[e, u, d, s, mu, c, tau, brA, brB]`

with bare exponents

`n = [0, 3, 5, 11, 11, 16, 17, 19, 26]`

on the McKay tree edges

`(e,u), (u,d), (d,s), (s,mu), (mu,c), (c,tau), (tau,brA), (c,brB)`.

For each edge `(i,j)`, define the bare exponent jump

`Delta n_ij = n_j - n_i`.

Write the edge lift in `Z[phi]` as

`dz_ij = da_ij + db_ij phi`

with the consistency condition

`5 da_ij + 6 db_ij = Delta n_ij`.


## Minimal-L1 principle

For each edge, choose the unique integer solution `(da,db)` minimizing

`|da| + |db|`.

This gives a canonical section of the quotient map

`pi : Z^2 -> Z,   pi(a,b) = 5a + 6b`

restricted to the actual edge-jump set

`{0, 1, 2, 3, 5, 6, 10}`.


## Explicit edge lifts

The unique minimal-L1 lifts are:

- `0  -> (0,0)   = 0`
- `1  -> (-1,1)  = phi^-1`
- `2  -> (-2,2)  = 2phi^-1`
- `3  -> (3,-2)  = -phi^-3`
- `5  -> (1,0)   = 1`
- `6  -> (0,1)   = phi`
- `10 -> (2,0)   = 2`

So the simple-edge set is recovered automatically:

`dz in {0, 1, 2, phi, phi^-1, 2phi^-1, -phi^-3}`

which is a subset of

`{0} union {+/-phi^r} union {+/-2phi^r}`.


## Reconstruction of node values

Because the McKay graph is a tree, once the root is fixed to

`z_e = 0`

the edge lifts integrate uniquely to node values:

- `z_e   = (0,0)`
- `z_u   = (3,-2)`
- `z_d   = (1,0)`
- `z_s   = (1,1)`
- `z_mu  = (1,1)`
- `z_c   = (2,1)`
- `z_tau = (1,2)`
- `z_brA = (-1,4)`
- `z_brB = (4,1)`

This is exactly the bounded-search uniqueness solution.


## Emergent consequences

Once the node values are reconstructed in this way, the earlier constraints
become consequences:

- `C1`: each node norm is `0`, `+/-1`, `+/-5`, or `+/-19`
- `Q`: every edge satisfies `dz / sigma(dz) = +/- phi^k`
- `C2`: `sum_edges N(dz) = -6`
- node sum: `sum_nodes N(z) = +6`
- flatness: `sum_nodes N(z) + sum_edges N(dz) = 0`


## Status

This is a genuine theoretical improvement over the original search-based
statement:

- it is constructive, not just filtering;
- it explains why the simple-edge condition appears;
- it reduces the uniqueness problem to a local canonical lift rule on edges.

What it does **not** yet solve:

- it still uses the neutral branch ordering `brA/brB`;
- it does not by itself identify which neutral endpoint should be interpreted
  physically as `t` or `b`;
- it does not yet prove that the minimal-L1 rule is forced by deeper dynamics.


## Verifier

Implemented in:

`reproducible/verify_minimal_edge_lifts.py`

Run:

```powershell
python reproducible\verify_minimal_edge_lifts.py
```
