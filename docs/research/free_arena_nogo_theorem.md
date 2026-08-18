# Free-arena no-go claim: refutation and corrected theorem

Date: 2026-07-27 (third session)

## Decision

The proposed free-arena no-go is **REFUTED**.

The failure occurs in its first step.  On

`H=C[G] tensor C^m`

with `A=R(C[G]) tensor I_m`,

`A' = L(C[G]) tensor M_m(C)`,

not merely `L(C[G])`.  Order zero only puts `JAJ^-1` inside this larger
commutant.  It does not force a factor-preserving right-to-left exchange.

An exact counterexample is constructed for `G=Q8`, `m=16`.  It has:

- the free left regular `Q8` action;
- a faithful right group-algebra action;
- nonzero self-adjoint equivariant `D`;
- form-like even doubling with signs `(J^2,JD,Jgamma)=(+,+,-)`;
- exact order zero and first order on the full group-algebra spanning set;
- nonzero inner one-forms.

This directly contradicts the claimed mutual incompatibility of order zero,
KO reality, and nonzero fluctuations on every free arena.

A corrected **factor-preserving no-go theorem** does hold.  It explains the
primal and primal--dual failures, but its extra hypothesis excludes the Q8
counterexample.  All finite counterexample checks are in
`reproducible/verify_free_orbifold_arenas.py`; the earlier 2I instances remain
certified by `verify_inner_fluctuations.py` and
`verify_primal_dual_triple.py`.

## 1. What order zero actually forces

Let

`A=R(C[G]) tensor I_m subset End(C[G] tensor C^m)`.

The double-centralizer theorem for the regular representation gives

`A'=L(C[G]) tensor M_m(C)`.

If `J` is antiunitary and order zero holds, then

`JAJ^-1 subset A'`.

The image has the same complex dimension `|G|` as `A`, but `A'` has dimension
`|G|m^2`.  Equality with `L(C[G]) tensor I_m` follows only for `m=1`, or
after an additional factor-preservation hypothesis.  For `m>1`, the image
may use the multiplicity matrices and may be entangled with the left regular
factor.

For a proper star-subalgebra `B subset R(C[G])`, order zero gives only

`JBJ^-1 subset B'`.

Here `B'` is larger still.  Dimension counting cannot promote this inclusion
to a side exchange.  In particular a three-block
`C+H+M3(C)` corner has no theorem forcing its opposite action to occupy the
corresponding left corner.  That requires a specified bimodule
representation.  **DERIVED correction.**

## 2. The corrected factor-preserving theorem

**Theorem (DERIVED).**  Let

`H=C[G] tensor C^m`,

`A=R(C[G]) tensor I_m`,

and let `D in A_left' = R(C[G]) tensor M_m(C)` be `G`-equivariant.  Suppose
an antiunitary `J` satisfies:

1. `JAJ^-1=L(theta(C[G])) tensor I_m` for a star anti-automorphism `theta`;
2. conjugation by `J` sends
   `R(C[G]) tensor M_m(C)` into
   `L(C[G]) tensor M_m(C)` without exchanging the regular factor with a
   regular representation hidden in `C^m`;
3. `JD=epsilon' DJ`, `epsilon' in {+1,-1}`.

Then

`D in Z(C[G]) tensor M_m(C)`.

Consequently `[D,A]=0` and `Omega_D^1(A)=0`.

### Proof

Equivariance gives

`D in R(C[G]) tensor M_m(C)`.

By hypotheses 2 and 3,

`D=epsilon' J D J^-1
   in L(C[G]) tensor M_m(C)`.

Therefore

`D in (R(C[G]) intersection L(C[G])) tensor M_m(C)`.

For the regular representation,

`R(C[G]) intersection L(C[G])=Z(C[G])`.

Hence `D` has central group-algebra coefficients.  Central left and right
multiplications agree, and every such coefficient commutes with the right
algebra.  Thus `[D,a]=0` for all `a in A`, and every represented inner
one-form vanishes.  QED.

This proof permits any fixed permutation or antiunitary transformation of
the multiplicity indices: it uses the total operator-space intersection, not
individual coefficient labels.  What it forbids is exchanging the regular
factor with a second regular module inside the multiplicity space.

### Proper subalgebras

If `B subset A` is proper, factor-preserving reality still implies the same
central conclusion only when hypotheses 1--2 are imposed for the full
`C[G]` action.  Order zero for `B` alone does not.  The correct general
statement is merely

`JBJ^-1 subset B'`;

neither side exchange nor centrality follows without more bimodule data.
This is the exact scope for an SM-type corner.

## 3. Exact free Q8 counterexample

Write

`H0=C[Q8]_x tensor C[Q8]_y`.

The physical symmetry acts by `L_g` on the `x` factor, and

`A=R_x(C[Q8])`.

Define

`J0(delta_x tensor delta_y)
 = delta_(y^-1) tensor delta_(x^-1)`

anti-linearly.  Then

`J0^2=1`,

`J0 R_g^x J0^-1=L_(g^-1)^y`.

The opposite action lives in the multiplicity factor and commutes with
`A`, so order zero holds.

Let `i,j` be quaternion generators and set

`X=i_complex (R_i^x-R_(-i)^x)`.

`X` is self-adjoint and noncentral.  Put

`D0=X+J0 X J0^-1`.

Then:

- `D0=D0*`;
- `D0` commutes with every free left `Q8` action;
- `J0D0=D0J0`;
- `[[D0,R_a^x],J0R_b^xJ0^-1]=0`, because the first commutator acts on `x`
  and the opposite algebra acts on `y`;
- `[D0,R_j^x]!=0`.

Thus order zero, first order, reality, equivariance, and nonzero fluctuations
already coexist.

For an even KO6 version, take `H=H0 direct-sum H0` and

`Gamma=diag(+1,-1)`,

`D=[[0,D0],[D0,0]]`,

`J=[[0,J0],[J0,0]]`.

Exact `128 by 128` matrices give

`J^2=+1`, `JD=DJ`, `JGamma=-Gamma J`, `{D,Gamma}=0`,

while order zero, first order, and `[D,A]!=0` remain true.  Since

`H = C[Q8] tensor C^16`

as a free left module, this is a direct counterexample to the proposed
theorem on its stated class of arenas.

**Post-audit correction:** these are finite real-even data satisfying the
listed axioms, but the manufactured even double demonstrably fails
orientability and Poincare duality.  For every group, a trivial double with
identical sheet representations has

`pi(a) J pi(b) J^-1 = diag(S,S)`,

so its intersection pairing is

`Tr(diag(1,-1) diag(S,S))=0`.

The Q8 control therefore proves compatibility only of the displayed
bilinear/KO axioms and nonzero fluctuations.  It is not a Poincare-dual,
manifold-like finite spectral triple.  The actual `128 x 128` matrices and
all `64` group-element spanning pairs are checked by the verifier.
**DERIVED negative.**

## 4. Finite 2I instances

The earlier 2I results are consistent with the corrected theorem:

1. **Primal `22 Reg`.**  Orbitwise inversion is factor-preserving and makes
   order zero/first order hold for the right algebra, but the noncentral
   incidence coefficients give neither `JD=DJ` nor `JD=-DJ`.  The canonical
   left placement is real-compatible but has zero fluctuations.

2. **Primal plus dual `44 Reg`.**  Pure cellular star has KO6 signs but
   preserves right multiplication, so order zero fails.  Adding
   factor-preserving inversion repairs the opposite side but loses the
   `JD` sign.

3. **Synthetic Q8.**  The multiplicity factor is itself a regular module and
   `J` exchanges the two regular coordinates.  This violates the
   factor-preserving hypothesis and evades the center intersection.

The finite results therefore identify two different escape mechanisms:

- non-free/stabilizer modules;
- free modules whose multiplicity space carries additional regular
  bimodule structure.

## Status ledger

### Strengthened

- **DERIVED:** the exact commutant is
  `L(C[G]) tensor M_m(C)`.
- **DERIVED:** corrected factor-preserving no-go theorem.
- **DERIVED scoped counterexample:** free Q8 KO6 data with order zero, first
  order, and nonzero fluctuations, but with failed orientability and
  identically zero intersection form.
- **DERIVED:** the previous primal/primal--dual split obstruction is a
  factor-preserving result, not a universal free-arena theorem.

### Refuted / downgraded

- **REFUTED:** freeness alone makes order zero, KO reality, and nonzero
  fluctuations incompatible.
- **REFUTED:** order zero forces a right-to-left exchange when `m>1`.
- **DOWNGRADED:** non-free arenas are not the only unbroken-symmetry escape.

### Open

- classification of entangled opposite representations inside
  `L(C[G]) tensor M_m`;
- whether the 2I multiplicity `m=22` or `44` has a geometrically selected
  second regular/bimodule organization;
- the remaining finite-triple axioms and physical content of such an
  entangled construction.
