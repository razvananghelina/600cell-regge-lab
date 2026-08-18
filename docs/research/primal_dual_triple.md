# Primal--dual 600/120-cell arena and the real-structure gate

Date: 2026-07-27 (second session)

## Decision

The doubled cellular arena and its Hodge-star KO signs are **DERIVED**.  The
real spectral-triple gate is **DERIVED negative for every enumerated
star/sign variant**.

The oriented cellular dual of the boundary of the 600-cell is the boundary of
the 120-cell.  Its cells are indexed by complementary primal cells, so its
f-vector is

`(600,1200,720,120)`.

The left `2I` action is free on every dual layer.  Thus

`C_primal = C[2I] tensor C^22`,

`C_dual   = C[2I] tensor C^22`,

`H = C_primal direct-sum C_dual = C[2I] tensor C^44`.

There is an exact pure-Hodge variant with signs

`(J^2, JD/DJ, J gamma/gamma J)=(+,+,-)`,

the KO6 table.  This verifies the predicted grading and reality signs.
However, pure Hodge star sends the canonical right algebra to another right
action.  Since `C[2I]` is noncommutative, order zero fails, and the displayed
first-order witness fails as well.

Composing Hodge star with the previously derived orbitwise inversion sends
right multiplication to left multiplication.  It repairs order zero and
first order exactly, but then `JD` has neither sign:

- for the degree-even star convention,
  `nnz(JD-DJ)=22000`, `nnz(JD+DJ)=40760`;
- for the degree-alternating convention these two numbers exchange.

Changing the star orientation, sheet-return sign, or antiunitary phase does
not evade this split obstruction.  Therefore no enumerated candidate
simultaneously has a valid KO reality table, order zero, first order, and
nonzero fluctuations.  A real spectral triple has **not** been constructed.

All finite certificates are in
`reproducible/verify_primal_dual_triple.py`.

## 1. Explicit oriented cellular dual

Let `P^k` be primal `k`-cochains in the increasing-simplex orientation used
by the existing Kähler--Dirac verifier.  Its coboundaries are

`d_k:P^k -> P^(k+1)`.

Define the dual cellular layer

`Q^j = { dual(sigma) : sigma a primal (3-j)-cell }`.

This is an explicit basis bijection, not a claim that the nonsimplicial dual
cells are primal simplices.  Orient the dual cells by the convention

`q_j = d_(2-j)^T`, for `j=0,1,2`.

Thus:

- dual vertices are primal tetrahedra: 600;
- dual edges are primal triangles: 1200;
- dual faces are primal edges: 720;
- dual solids are primal vertices: 120.

The verifier constructs all six integer incidence matrices and checks

`d_1 d_0=d_2 d_1=0`, `q_1 q_0=q_2 q_1=0`,

and the three transpose relations exactly over `Z`.

### Exact homology certificate

The old primal ranks were verified numerically.  This verifier strengthens
them to exact ranks.  Gaussian elimination over `F2` gives the lower
certificates

`rank(d_0,d_1,d_2)=(119,601,599)`.

Connectedness gives `rank(d_0)<=119`, while `d^2=0` successively gives

`rank(d_1)<=720-119=601`,

`rank(d_2)<=1200-601=599`.

The lower and upper bounds coincide, so the ranks over `Q` are exact.  The
dual ranks are exactly `(599,601,119)`.  Both complexes have

`b=(1,0,0,1)`.

**DERIVED exact**, with no diagonalization.

## 2. Free action and module chart

The primal verifier again enumerates all 120 signed cell actions.  Every
orbit representative has 120 distinct targets.  The orbit counts are

`P: (1,6,10,5)`.

Cellular duality transports the signed action from a primal `k`-cell to its
dual `(3-k)`-cell, giving

`Q: (5,10,6,1)`.

Freeness is preserved because a group element stabilizes a dual cell exactly
when it stabilizes its indexed primal cell.  The explicit primal
representatives are those printed by
`verify_inner_fluctuations.py`; the explicit dual representatives are their
formal duals in reverse degree.  This proves

`P=22 Reg`, `Q=22 Reg`, `H=44 Reg`.

**DERIVED.**  The construction is canonical as an oriented cellular dual
after the displayed dual-orientation convention.  A Euclidean metric volume
normalization is not used: this is the unit-weight combinatorial star.

## 3. Star, Dirac, and the finite sign family

Set

`D_P=d+d*`, `D_Q=q+q*`, `D_tot=diag(D_P,D_Q)`.

Adjoints here use the explicitly specified orthonormal cellular metric:
every oriented primal and dual cell basis vector has norm one.  This is the
same Hilbert convention used for the established primal `D`; it is a
**STRUCTURAL metric convention**, not a computed circumcentric volume
Hodge star.  With it, `D_tot` is exactly self-adjoint, form-odd, and has
four harmonic modes.

For `lambda in {+1,-1}`, define

`S_lambda:P^k -> Q^(3-k)`,

`S_lambda(sigma)=lambda^k dual(sigma)`.

The incidence convention gives the exact identity

`S_lambda D_P = lambda D_Q S_lambda`.

Let `sigma in {+1,-1}` be the return-sheet sign and define the antiunitary
linear part

`U_(lambda,sigma)
 = [[0, sigma S_lambda^T],[S_lambda,0]]`,

with coefficient conjugation understood.  Exact sparse identities give:

| `lambda` | `sigma` | `J^2` | `JD` | `J gamma_form` | KO reading |
|---:|---:|---:|---:|---:|---|
| `+1` | `+1` | `+1` | `+DJ` | `-gamma J` | KO6 |
| `+1` | `-1` | `-1` | `+DJ` | `-gamma J` | KO2 |
| `-1` | `+1` | `+1` | `-DJ` | `-gamma J` | no standard even KO table |
| `-1` | `-1` | `-1` | `-DJ` | `-gamma J` | no standard even KO table |

The grading sign is always negative because `k` and `3-k` have opposite
parity.  Multiplication of an antiunitary by `e^(i theta)` changes none of
these signs:

`(e^(i theta)J)^2=e^(i theta)e^(-i theta)J^2=J^2`.

Thus factors `+1,-1,+i,-i`, or any other unit phase do not add sign variants.
**DERIVED finite enumeration.**

## 4. Algebra, order conditions, and the obstruction

Let

`A_R=R(C[2I])`

act by the same right convolution on all 44 orbit coordinates.  The
Galois-flipped realization is obtained by the simultaneous block exchange

`rho1<->rho7`, `rho2<->rho8`.

For the full group algebra it is the same algebra with relabelled Wedderburn
factors, so the following order-condition verdict is unchanged.

### Pure Hodge star

Because `S_lambda` is equivariant and does not invert the group coordinate,

`J R_b J^-1 = R_(theta(b))`,

where `theta` is identity for the displayed action and may be the Galois
automorphism for the relabelled realization.  Either way the opposite action
is still right multiplication.  Exact noncommuting group elements `s,t`
give

`[R_s,J R_t J^-1] != 0`.

Hence order zero fails for all four pure-star variants.  The verifier also
computes a nonzero

`[[D,R_s],J R_t J^-1]`,

so the same witness fails first order.  This includes the KO6 variant.
**DERIVED negative.**

The obstruction is intrinsic to every equivariant pure cellular star and
every degree/sheet/phase sign: none changes right into left multiplication.
It also applies to either Galois block choice.  The noncommutative `H` and
`M3(C)` factors of either SM-type corner therefore fail order zero as well.

### Hodge star composed with orbitwise inversion

Let `I(delta_g tensor e_alpha)=delta_(g^-1) tensor e_alpha` on both sheets
and use `J=S_lambda I K` with either sheet-return sign.  Now

`J R_b J^-1=L_(theta(b)^*)`,

so left/right commutation proves order zero.  Since `D` and every `R_a` are
left-equivariant,

`[[D,R_a],L_c]=0`,

and first order also holds exactly.

But inversion is not compatible with the derived incidence coefficients.
All four composed variants have neither `JD=DJ` nor `JD=-DJ`, with exact
residuals:

| `lambda` | either `sigma` | `nnz(JD-DJ)` | `nnz(JD+DJ)` |
|---:|---|---:|---:|
| `+1` | `+/-1` | 22000 | 40760 |
| `-1` | `+/-1` | 40760 | 22000 |

The return-sheet sign and antiunitary phase cannot turn a nonzero residual
into zero.  **DERIVED negative for this entire composed family.**

## 5. Exact nonzero one-form census

Although no candidate passes the full gate, the nonzero-fluctuation statement
and its dimension can be computed exactly.

The extracted `D` contains 124 group coefficients in 112 multiplicity
blocks, supported on 13 group elements.  In a Wedderburn block
`A_i=M_(n_i)(C)`, let `U_i` be the span of the identity and all projected
coefficient blocks of `D`.  Exact matrices over `Q(sqrt(5),i)` give

`n_i = (1,2,3,4,5,6,4,2,3)`,

`r_i=dim U_i=(1,4,9,12,12,12,12,4,9)`.

The dimension formula is exact.  Since

`A_i tensor A_i^op = End(A_i)`,

a universal one-form acts as a linear map on `A_i` that vanishes on the
identity.  Restriction to `U_i` therefore has dimension

`n_i^2(r_i-1)`.

Summing gives

`dim_C Omega_D^1(A_R)=sum_i n_i^2(r_i-1)=1191`.

The doubled diagonal representation does not double this number: its primal
component determines the dual component through the same algebra element and
the dual incidence relation.

The represented calculus is star-stable, so its self-adjoint part has real
dimension 1191.  Every one-form is off-diagonal in form degree and hence has
zero ordinary Hilbert-space trace.  The standard trace-unimodularity
condition removes no further one-form direction:

`candidate self-adjoint unimodular field dimension = 1191`.

Separately, the unimodular gauge Lie algebra of the full Wedderburn algebra
has real dimension

`sum_i n_i^2-1=119`.

These are **DERIVED candidate-space dimensions**, not dimensions of physical
gauge fields in a real spectral triple, because the gate fails.

## 6. Axiom-by-axiom verdict

| property | pure KO6 star | star plus inversion |
|---|---:|---:|
| finite Hilbert space `H=44 Reg` | holds | holds |
| `D=D*`, nonzero | holds | holds |
| `{D,gamma}=0` | holds | holds |
| `J^2=+1` choice | holds | holds |
| `J gamma=-gamma J` | holds | holds |
| `JD=DJ` | holds | **fails: neither sign** |
| order zero for `A_R` | **fails** | holds |
| first order for `A_R` | **fails** | holds |
| nonzero inner one-forms | holds, dimension 1191 | holds, dimension 1191 |
| real spectral triple | **not constructed** | **not constructed** |

The failure is intrinsic within the enumerated derived family:

- every pure equivariant star preserves the right/right placement and fails
  order zero for a noncommutative algebra;
- every star--inversion variant inherits the exact nonzero inversion/Dirac
  residual, merely exchanging its plus/minus counts when `lambda` changes.

This is not a theorem against every imaginable primal--dual antiunitary.
An additional geometrically derived unitary mixing the 44 orbit labels could
in principle alter both facts.  No such unitary is currently supplied by the
cellular duality, so that larger classification remains **OPEN**.

## 7. Physics gate

Task 4 is not activated.  In particular:

- neither SM-type block choice defines a verified real spectral triple;
- no physical commutant generator `Y` is licensed;
- Route C anomaly forcing has no valid module/generator input;
- the `M15/M16` census and generation count do not start;
- the 1191 candidate one-form directions are not Yukawa parameters;
- no fluctuated-`D` mass spectrum is physically licensed.

Therefore the registered `Z[phi]` mass-exponent comparison is not run.  This
is a **protocol skip**, not a null match, and creates zero new look-elsewhere
trials.

## Status ledger

### Strengthened

- **DERIVED:** explicit oriented 120-cell dual incidence complex.
- **DERIVED exact:** primal and dual ranks and Betti numbers.
- **DERIVED:** free dual action and `H=44 Reg`.
- **DERIVED:** the complete four-row cellular-star sign table, including a
  genuine KO6 sign solution before order conditions.
- **DERIVED:** exact order-zero/first-order failures for pure star.
- **DERIVED:** exact order-zero/first-order success but `JD` failure for all
  star--inversion variants.
- **DERIVED:** nonzero one-form dimension 1191 and candidate unimodular
  self-adjoint dimension 1191.

### Downgraded / closed

- **DERIVED negative:** the primal--dual star does not by itself complete the
  real spectral triple.
- **DERIVED negative:** the KO6 star and right-convolution algebra occupy the
  wrong mutual slot: `JAJ^-1` remains right multiplication.
- **DERIVED negative:** adding inversion fixes the algebra slot but destroys
  the `JD` sign.
- The earlier phrase “the only missing axiom is J-compatibility” is too
  narrow: pure Hodge `J` has the correct signs but exposes an independent
  order-zero obstruction.

### Open

- a geometrically derived orbit-mixing unitary that could simultaneously
  implement the opposite algebra and intertwine `D`;
- classification beyond the finite star/orientation/inversion/Galois family;
- a verified real spectral triple with nonzero fluctuations;
- `Y`, anomaly forcing, multiplets, generations, Yukawa blocks, and the
  frozen mass-spectrum comparison.
