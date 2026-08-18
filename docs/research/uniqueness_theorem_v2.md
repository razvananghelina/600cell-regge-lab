# Uniqueness Theorem v2

## Goal

Replace the current uniqueness statement

`C1 (irreducibility) + C2 (edge flatness) => unique assignment`

with a stronger and more structural formulation:

`C1 + Q + C2 => unique assignment`

where `Q` is an edgewise Galois-quantization condition.


## Definitions

Let the nine McKay-tree nodes carry elements

`z_f = a_f + b_f phi in Z[phi]`

with fixed bare exponents

`n_f = 5 a_f + 6 b_f in {0,3,5,11,11,16,17,19,26}`.

Write the McKay tree as the 8 edges

`(e,u), (u,d), (d,s), (s,mu), (mu,c), (c,tau), (tau,brA), (c,brB)`.

For an edge `e = (i,j)`, define

`dz_ij = z_j - z_i`.


## Constraints

### C1. Irreducibility

Each `z_f` is zero, a unit, or an irreducible element of `Z[phi]`.

Equivalently:

`|N(z_f)| in {0,1} union {p prime : p mod 5 in {0,1,4}}`.


### Q. Edge Galois Quantization

For every McKay edge `(i,j)`, the Galois ratio is quantized:

`dz_ij / sigma(dz_ij) = +/- phi^{k_ij}`

for some integer `k_ij`.

Equivalently, the logarithmic Galois phase

`log_phi |dz_ij / sigma(dz_ij)|`

is an integer on every edge.


### C2. Flatness

The total signed edge norm is fixed:

`sum_(i,j in E) N(dz_ij) = -b1 = -6`.


## Computational Theorem (proved by verifier)

Within the bounded search window `|t_f| <= 15` for the parametrization

`(a,b) = (-n - 6t, n + 5t)`,

the triple of constraints `C1 + Q + C2` selects a unique assignment:

`[(0,0), (3,-2), (1,0), (1,1), (1,1), (2,1), (1,2), (-1,4), (4,1)]`.

This is verified by

`reproducible/verify_uniqueness_quantized.py`

with counts:

- `543,129,600` assignments satisfy `C1` in `|t|<=15`
- `20` assignments satisfy `C1 + Q`
- `1` assignment satisfies `C1 + Q + C2`


## Physical Edge Data of the Selected Solution

For the unique selected assignment, the edge differences are:

- `z_u - z_e = -phi^{-3}`
- `z_d - z_u = 2 phi^{-1}`
- `z_s - z_d = phi`
- `z_mu - z_s = 0`
- `z_c - z_mu = 1`
- `z_tau - z_c = phi^{-1}`
- `z_brA - z_tau = 2 phi^{-1}`
- `z_brB - z_c = 2`

Hence the edge ratios are:

- `-phi^{-6}`
- `-phi^{-2}`
- `-phi^{2}`
- `+phi^{0}`
- `+phi^{0}`
- `-phi^{-2}`
- `-phi^{-2}`
- `+phi^{0}`

So the selected solution is not only flat; its edge transport is quantized in
integral powers of the fundamental unit `phi`.


## Why This Version Is Better

The original theorem uses `C1 + C2` only. That leaves a very large search
space, and the present proof is only a bounded search.

The strengthened route is better because:

1. `Q` is local on edges, not global on the whole assignment.
2. `Q` is Galois-theoretic, so it has a direct arithmetic interpretation.
3. `Q` shrinks the candidate set drastically before `C2` is applied.
4. `C2` becomes the final selector of a small quantized set, instead of doing
   all the work by itself.


## Suggested New Theorem Statement

> **Theorem (Quantized Norm Selection, bounded-search form).**
> Let `z_f = a_f + b_f phi in Z[phi]` be assigned to the nine McKay-tree nodes,
> with fixed bare exponents `n_f = 5 a_f + 6 b_f` equal to
> `{0,3,5,11,11,16,17,19,26}`.
> Assume:
> 1. each `z_f` is zero, a unit, or an irreducible element of `Z[phi]`;
> 2. for every McKay edge `(i,j)`, the ratio `dz_ij / sigma(dz_ij)` is of the
>    form `+/- phi^k` with `k in Z`;
> 3. `sum_(i,j in E) N(dz_ij) = -6`.
>
> Then, within the search domain `|t_f| <= 15`, there is a unique assignment,
> namely
> `[(0,0), (3,-2), (1,0), (1,1), (1,1), (2,1), (1,2), (-1,4), (4,1)]`.


## What Remains Open

This is not yet a global analytic theorem. Two steps remain:

1. Prove an a priori bound on the allowed `t_f`, or replace bounded search by
   an infinite-family exclusion argument.
2. Derive `Q` from a more primitive framework principle, rather than taking it
   as an axiom.


## Best Next Move

The next mathematical target should be:

**derive `Q` from the arithmetic of edge differences**

for example by showing that allowed edge Wilson lines must lie in the subgroup

`{0} union { +/- m phi^r : m in {1,2}, r in Z }`

or an equivalent quantized transport class.
