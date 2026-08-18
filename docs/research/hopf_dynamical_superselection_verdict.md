# Existing symmetry and Hessian dynamics do not induce Hopf superselection

Date: 2026-08-10

Protocol commit: `0807c44` (mechanisms and decision boundary frozen before
their fixed algebras were computed).

Verifier: `reproducible/verify_hopf_dynamical_superselection.py`.  Targeted
result: `14/14`.

## Complete scope

On the already derived six-fibration carrier `H_F=C^6`, this audit uses only:

1. the exact effective `A5` permutation action;
2. the complete Hessian

   ```text
   H_X(i,j)=3 Tr(X(Box_i Box_j+Box_j Box_i))
   ```

   of the previously audited cubic on the normalized five-dimensional
   Hopf--Box field space.

No bath, noise model, state, new Hamiltonian, regulator, label weights or
diagonal projection is added.  Exact superselection means that all six
minimal projections `E_ii` are fixed observables.

## Symmetry averaging selects the wrong algebra

The exact Reynolds map

```text
E_A5(M)=(1/60) sum_g P_g M P_g^-1
```

has superoperator rank two.  Its image is

```text
span{I,J} ~= C+C,
```

as expected from the multiplicity-free permutation module `1+5`.  In
particular,

```text
E_A5(E_00)=I/6,
E_A5(E_01)=(J-I)/30.
```

Thus no individual fibration projection is fixed, and a collective
off-diagonal transition survives.  On the order-parameter module `W`, the
same group average is zero: averaging the orbit erases every nontrivial `X`
rather than selecting one of the six directions.

**DERIVED NEGATIVE:** `A5` covariance cannot be reinterpreted as six-sector
superselection.  Its canonical expectation forgets the individual labels.

## Hessian-generated dephasing also selects the wrong algebra

For self-adjoint `H_X`, infinite-time Heisenberg averaging projects onto the
commutant `{H_X}'`.  All six label projections can be fixed only if `H_X` is
diagonal in the label basis.

The exact linear map

```text
W -> offdiag(H_X)
```

has rank five, equal to `dim W`.  Its kernel is zero.  Consequently

```text
H_X label-diagonal iff X=0.
```

The normalized field condition `Tr(X^2)=7200` excludes this sole case.

At every desired point `X=Box_i`:

- all 36 commutators `[H_Box_i,E_jj]` are nonzero;
- the graph of nonzero off-diagonal Hessian entries is connected;
- `C(F) intersect {H_Box_i}' = C I`;
- the exact Hessian eigenvalue multiplicities are `1,1,2,2`, so the full
  commutant has dimension `1^2+1^2+2^2+2^2=10`.

The degeneracies therefore leave a ten-dimensional noncommutative fixed
algebra, while only one diagonal dimension survives.  Neither number is the
six-dimensional label algebra.

**DERIVED NEGATIVE:** the full cubic Hessian dynamically mixes labels at
every nonzero normalized field, including all six hoped-for vacua.

## Framing attack

Calling either construction a “canonical conditional expectation” is not
enough.  Canonicity identifies the averaging operation; superselection is a
claim about its image.  Here both images are computed exactly and neither is
`C(F)`.

The diagonal `D_aux` does preserve the label projections, but it was built
from the diagonal response `Phi`.  Using its evolution to justify that same
diagonal response would assume the conclusion.

## Status ledger

- **DERIVED:** the `A5` twirl fixes `span{I,J} ~= C+C`.
- **DERIVED NEGATIVE:** it fixes no individual fibration projection and
  retains the collective transition `(J-I)/30`.
- **DERIVED:** orbit averaging kills the whole nontrivial order-parameter
  module `W`.
- **DERIVED:** `X -> offdiag(H_X)` is injective on `W`.
- **DERIVED NEGATIVE:** no normalized `X` gives label-diagonal Hessian
  dynamics.
- **DERIVED NEGATIVE:** at every `Box_i`, Hessian dephasing preserves only
  scalar diagonal observables.
- **DERIVED NO-GO IN THE STATED ARENA:** neither exact, parameter-free
  averaging mechanism already present in the repository derives six-label
  superselection.
- **OPEN:** a separately derived physical environment or locality law.  It
  would be new input and must not be presented as implicit in the current
  construction.

## Next admissible bypass

The failure of superselection does not yet exclude dynamics on the complete
label Hessian.  A narrower question remains: can a basis-independent spectral
functional of the full `H_X`, with every off-diagonal channel retained,
select the Hopf directions directly?

The non-fitted audit must enumerate the complete characteristic data of
`H_X` on its five-dimensional physical subspace before comparing extrema
with `Box_i`.  A fitted combination of moments, or a heat-kernel regulator
chosen after viewing the extrema, would have no evidential value.
