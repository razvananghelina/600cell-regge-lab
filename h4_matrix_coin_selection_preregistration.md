# Preregistration: does H4 chamber geometry select a matrix coin?

Date: 2026-08-11

## Motivation and provenance boundary

Commit `a837d06` closed positive scalar reweighting of the four fixed local
directions.  A matrix coin can evade that scalar obstruction, but it also
introduces enough coefficients to fit a desired dispersion unless the
geometry sharply restricts it.

This protocol therefore counts the allowed matrix family **without testing
isotropy or a continuum target**.  The local Gram matrices and their scalar
failure are already known; this is target-blind only with respect to the coin
parameter census, not a claim of historically blind geometry.

## Frozen carrier and symmetry class

Use the active chamber carrier

\[
\mathcal H=\mathbb C^{14400}\otimes\mathbb C^4.
\]

The H4 action is transitive on the complete flag chambers and preserves the
four right-Coxeter colour labels.  A pointwise local coin is

\[
C_{\rm loc}=\bigoplus_k C_k.
\]

First prove the standard transitivity statement: equivariance forces all
`C_k` to be the same matrix.  Then record the full unitary freedom of this
constant matrix.  Do not call constancy uniqueness: `U(4)` has 16 real
parameters.

## Two preregistered metric tensors

At chamber 0 construct exactly the two already frozen local matrices

\[
G_d=(d_i\mathbin\cdot d_j)_{ij},
\qquad
G_u=(u_i\mathbin\cdot u_j)_{ij},
\]

for literal geodesic steps and unit directions.  Verify their labelled Gram
matrices are chamber-independent before using the representative.

For each `G`, compute:

1. all eigenvalues and their multiplicities;
2. the complex dimension of the commutant
   `{X in M_4(C): [X,G]=0}`, both from multiplicities
   `sum m_lambda^2` and from an explicit 16-column commutator matrix;
3. the real dimension and factorization of its unitary subgroup
   `product_lambda U(m_lambda)`;
4. the dimension of the intersection with the published robust block class
   `X=I_2 tensor Y`, `Y in M_2(C)`.

No free phase or block coefficient may be optimized.

## Complete coefficient-free spectral census

For a Gram matrix with `r` distinct eigenspaces, enumerate every Hermitian
unitary spectral function

\[
C_\epsilon=\sum_{j=1}^r\epsilon_jP_j,
\qquad \epsilon_j\in\{-1,+1\}.
\]

Identify `C` and `-C`, which differ only by a global phase for a single coin
application.  The complete count must therefore be

\[
N_{\rm spectral}=2^{r-1}.
\]

Record the full multiset of ranks of the `+1` eigenspaces.  Do not compare
these coins with isotropy in this step.

One member is designated before the calculation:

\[
C_0=2P_{\ker G}-I_4.
\]

The zero eigenspace is distinguished algebraically rather than by a desired
dispersion.  Verify whether it is one-dimensional, whether `C_0` is real,
symmetric and unitary, and whether it belongs to the paper's restricted
`I_2 tensor U(2)` class.

The reflection about `sqrt(p)`, where `p` is the unique positive zero-drift
probability vector, is **not** part of this spectral census.  It is a second,
structural Szegedy-style construction and must be preregistered separately
before any dynamical comparison.

## Decision boundaries

- **DERIVED NONSELECTION:** H4 equivariance leaves a positive-dimensional
  non-global unitary family.
- **DERIVED METRIC NONSELECTION:** adding `[C,G]=0` still leaves more than a
  global `U(1)`.
- **DERIVED BLOCK COLLAPSE:** imposing both Gram commutation and the paper's
  block form leaves only scalar matrices; this would select no directional
  mixing.
- **DERIVED CANONICAL CANDIDATE:** a one-dimensional kernel makes `C_0` a
  coefficient-free nontrivial reflection, even though the broader axioms do
  not force it.

These outcomes are compatible: the geometry can provide one canonical
candidate without uniquely selecting it from all equivariant coins.

## Scope

This protocol counts local matrix freedom only.  It does not test the
one-period distribution, Floquet spectrum, Dirac cone, refinement, mass or
light speed.  Any later dynamical test must commit the candidate list first
and report a full hit fraction.
