# Adversarial protocol: direct sign-bisection local branch stability

Date: 2026-08-22.

Primary result commit: `c5a202b`.

Status: frozen before the adversarial verifier exists or any listed rational
bracket is evaluated.

## 1. Purpose and independence boundary

Attempt to falsify the primary local-signature theorem at `v=3/2` with a
mechanically different real-root certificate.

This route must not use:

- interval Newton or a Krawczyk operator;
- decimal root or stationary-point seeds from the discovery artifact;
- the primary root balls to construct or refine a bracket;
- an adaptive scan in `v` or an explicit neighbourhood radius;
- the primary full-action residual as its decisive check.

Use only exact rational brackets, outward-rounded endpoint signs, derivative
signs on complete brackets, analytic tails and recursive state balls. The
primary artifact may be compared only after the adversarial tree has been
assembled in memory.

## 2. Fixed equations and hypotheses

Use exactly the hypotheses and functions in
`gravity_600cell_finite_height_local_signature_protocol.md`, including

```text
E=4*pi*(mu(q)-m)+q*(p(q)-pi),
E_q=p(q)-pi,
h=(p(q)-pi)/(2*pi*mu(q)),
r=2*m/mu(q)-1,
m_plus=m/r,
pi_plus=p(q)+(p(q)-pi)/r.
```

The symbolic radical factorization from the primary failure resolution is a
shared exact identity, not the decisive numerical method. Recheck it directly.

## 3. Frozen rational brackets

The complete ordered brackets are frozen as follows. `diag` denotes the exact
persistent root and stationary point `q=3/2`.

| State path | Stationary brackets for `p(q)=pi` | Root brackets for `E=0` |
|---|---|---|
| `root` | `diag`, `[5,6]` | `diag`, `[9,10]` |
| `root/c0` | `[1,2]`, `[17,18]` | `[0,1]`, `[9,10]`, `[31,32]` |
| `root/c0/c0` | `[-1,0]` | `[-3,-2]`, `[0,1]` |
| `root/c0/c1` | `[1,2]`, `[55,56]` | `[-1,0]`, `[31,32]`, `[99,100]` |
| `root/c0/c1/c0` | `[1,2]`, `[177,178]` | `[-1,0]`, `[99,100]`, `[316,317]` |

These integer brackets are disclosed after the primary result. Their widths or
locations carry no physical meaning. Every bracket must pass the internal sign
and monotonicity proof; a stored midpoint match cannot rescue a failed bracket.

## 4. Direct certified bisection

Represent every endpoint as an exact `Fraction` and every point evaluation as
an Arb ball at 220 decimal digits or more.

For a bracket `[a,b]` of a scalar function `f`:

1. require strict opposite signs at `a` and `b`;
2. require `f'` to have one strict sign on the complete Arb interval `[a,b]`;
3. bisect at exact dyadic midpoints for at most 420 steps;
4. retain the sign-changing half whenever the midpoint sign is strict;
5. if the midpoint ball first contains zero after at least 240 steps, stop and
   retain the current sign-changing bracket;
6. require the final endpoints still to have strict opposite signs and the
   derivative to exclude zero on the final interval.

This proves existence by the intermediate-value theorem and uniqueness by
strict monotonicity. A midpoint ambiguity is not treated as a root value.

For every stationary bracket use `f=p(q)-pi`, `f'=p'(q)`. For every non-diagonal
root bracket use `f=E`, `f'=p(q)-pi`.

## 5. Completeness and recursion

At every visited state independently certify:

- the position of `pi` relative to `p_star`, `p_infinity`, and zero;
- the resulting complete stationary count on both axes;
- strict `E` signs at all stationary enclosures other than the exact diagonal;
- nonzero `E(0)` and both analytic-tail coefficients;
- that the number and ordering of sign-reversal intervals equals the frozen
  root-bracket list;
- pairwise disjoint final root intervals;
- strict `h`, `r`, `mu(q)-2*m`, and `E_q` signs for every non-diagonal root;
- the exact recursive counts `(2,1)`, `(3,2)`, `(2,0)`, `(3,1)`, `(3,1)`;
- terminals `DEAD` and `ENTERED_D` on the same ordered paths;
- strict entry inequalities `0<m<2/5`, `q>0`, `m*q>125`.

Propagate the state using the entire bisection root interval. Do not replace it
with a decimal midpoint.

## 6. Frozen positive and negative controls

Before the gravity tree, run the same direct bisection engine on

```text
f(q)=q^2-2.
```

- positive control: `[1,2]` must certify exactly one root;
- negative control: `[2,3]` must be rejected because its endpoint signs agree.

After the adversarial tree is assembled in memory:

- require `r-(1+h*q)` to contain zero with enclosure width below `1e-150` on
  every physical edge;
- compare the ordered root intervals and strict-sign pattern with the primary
  artifact, without using the comparison to repair any result;
- require the hostile entry inequality `m*q>126` to fail.

## 7. Outcomes

### `LOCAL_SIGNATURE_ADVERSARIALLY_CORROBORATED`

Both polynomial controls behave as frozen; every rational stationary and root
bracket is certified directly; the complete recursive tree and strict terminal
gates agree with the primary result; all post-construction controls pass.

This permits consolidation of the existence of an unspecified open
neighbourhood. It does not certify an explicit radius or global basin.

### `LOCAL_SIGNATURE_ADVERSARIAL_OPEN`

Use for any failed bracket sign, derivative dependency, early midpoint
ambiguity, root-count disagreement, recursive-state disagreement or control
failure. Preserve the artifact and keep the local theorem **OPEN** until the
two methods are reconciled.

## 8. Interpretation boundary

Even corroboration proves only local structural stability on the special
one-parameter incoming curve. It neither derives `v=3/2` nor makes the
extendibility criterion local, generic, physically selected or
nonhomogeneous.

