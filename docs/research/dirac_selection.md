# C3 Dirac-selection audit

Date: 2026-07-23

## Decision

The derived cuts do **not** select a unique Dirac.

The exact cut sequence is

| stage | independent real dimension | status |
|---|---:|---|
| all self-adjoint, parity-odd, `C3`-equivariant operators on `W` | 148 | DERIVED foundation |
| Krajewski legality and first order on the explicit `16+14` witness | 132 | DERIVED (`d0`) |
| KO6 doubling and `JD=DJ` | 132 | DERIVED conditional on the established doubling |
| arithmetic Galois/sigma condition | 132 | no further derived cut |
| generic gauge quotient and overall scale | 122 | DERIVED generic stratum |
| polynomial spectral-action criticality | not finite | DERIVED negative for the natural quadratic/quartic tests |

Thus the selection principle remains **OPEN**, with a generic
122-real-dimensional projective gauge moduli space before a spectral
functional is chosen.  Moreover, the first nonconvex polynomial spectral
action already has a circle of gauge-inequivalent critical points, so
criticality alone is not a finite selector.

All finite claims are checked by
`reproducible/verify_dirac_selection.py`.

## 0. Precise arena

Use `2I = SL(2,5)` and choose any subgroup `H` of order three.  Exact
enumeration finds ten such subgroups.  For one `H`, its normalizer has order
12, hence its conjugacy orbit has size `120/12=10`; direct conjugation gives
all ten subgroups.  Therefore all order-three subgroups are conjugate.

**DERIVED:** the concrete `C3` choice is canonical up to `2I` conjugacy.

Keep yesterday's Krajewski witness, with labels
`(1,1bar,2,3,3bar)`:

`H+ = A + 2B`, with `A=(2,2)` and `B=(3bar,2)`,

`H- = C + 2B`, with `C=(2,1bar)`.

Its `C3` weight vectors are

`A=(2,1,1)`, `B=(2,2,2)`, `C=(0,1,1)`.

They sum to `(6,5,5)` and `(4,5,5)`, of dimensions 16 and 14.
For the upper odd block `T:H+ -> H-`, equivariant Hom dimensions are dot
products of weight vectors.  First order permits two Krajewski cells to
couple exactly when they share a left or right label.  The complete census is

| source to target | copies | legal? | complex dimension |
|---|---:|---|---:|
| `A -> C` | 1 | shared left `2` | 2 |
| `A -> B` | 2 | shared right `2` | 16 |
| `B -> C` | 2 | no shared label | 0 (8 excluded) |
| `B -> B` | `2 x 2` | identical cell | 48 |

Hence `dim_C T_legal=66`.  Self-adjointness fixes the reverse block to
`T*`, so

`d0 = dim_R D_legal,self-adjoint = 2*66 = 132`.

Without first order, `dim_C T=74`, or 148 real self-adjoint dimensions.
The legal space is therefore a real codimension-16 linear subspace of the
148-dimensional arena.  **DERIVED.**

## 1. Ordered cuts

### (a) Self-adjointness, legality, and first order

The result is `d0=132`, as above.  No exhibited block was silently promoted
to the full answer: all four source/target cell types and both `B`
multiplicities were counted.

### (b) KO6 reality

A single `16+14` sheet cannot carry the required parity-reversing real
structure.  Use the established doubled space

`Htilde = H + sigma(H)`

with opposite sheet chirality and the sheet-swap antiunitary `J`.  Before
reality the two legal sheets give 264 real parameters.  The equation
`JD=DJ` says that the second sheet is the anti-linear outer-twisted conjugate
of the first.  It is determined, not independently chosen.  Therefore the
fixed space has 132 real parameters.

**DERIVED conditional:** given the previously established doubled/Galois
KO6 construction, reality neither kills nor further selects the single-sheet
legal Dirac data.

### (c) Galois/sigma covariance

There are two distinct notions which must not be conflated.

1. The **DERIVED** outer/Galois sheet covariance permutes `2 <-> 2'` and
   `3 <-> 3'`.  Together with anti-linear sheet swap, this is already the
   relation imposed by `JD=DJ`; it leaves 132 real parameters.
2. An arithmetic equation on the entries of `D` is **not derivable**.
   General complex Dirac coefficients do not carry a specified
   `Q(sqrt(5))` lattice on which `sigma` acts.  Declaring all entries to lie
   in that field, or declaring a preferred basis, would be an additional
   axiom rather than a real-linear dimension cut.

Optionally imposing `sigma(D)=D` after choosing a
`Q(sqrt(5))` coefficient form restricts coefficients to `Q`, but this is a
countable arithmetic subset, not a canonically defined real submanifold, so
no honest real dimension is assigned.  A bare equation
`sigma(D)=phi D` is actually fatal: applying sigma twice gives
`D=sigma(phi)phi D=-D`, hence `D=0`.  The successful scalar relation
`z_b=phi sigma(z_t)` is not an involutive fixed-point law on one variable and
cannot be transferred to `D` without additional paired data.

**OPEN:** derive an integral coefficient lattice and a unique sigma law.

### (d) Gauge quotient and scale

The grading-preserving bimodule commutant is

`C_A + M2(C)_B+ + C_C + M2(C)_B-`,

whose unitary group has real dimension `1+4+1+4=10`.  On a generic legal
Dirac, its stabilizer is only the common scalar `U(1)`.  This is exact: the
twelve internal `B -> B` coefficient matrices can generically span `M2`,
forcing both `U(2)` factors to the same scalar; nonzero `A -> B` and
`A -> C` blocks then tie both remaining phases to it.  Thus the generic
conjugation orbit has dimension nine.

After quotienting by that orbit and by positive overall scale, the generic
moduli dimension is

`132 - 9 - 1 = 122`.

**DERIVED:** 122 is the dimension of the generic (principal) stratum.
Special Diracs have larger stabilizers and form lower-dimensional singular
strata; the quotient is not globally a manifold.

### (e) Spectral-action criticality

For `S_f(D)=Tr f(D^2)`, unconstrained variation gives

`D f'(D^2)=0`;

inside the legal linear space its orthogonal projection must vanish.  No
polynomial `f` is selected by the finite theory.

- If `f(x)=a x+b x^2` with `a>=0`, `b>0`, pairing the critical equation
  with `D` gives `a Tr(D^2)+2b Tr(D^4)=0`, so `D=0`.  This is the unique
  critical point and is removed by the nonzero scale quotient.
- For the first symmetry-breaking choice
  `f(x)=b x^2-a x`, `a,b>0`, the two legal `A <-> C` character channels
  may both have modulus `sqrt(a/(2b))`.  They solve the full critical
  equation.  Bimodule gauge removes their common phase but not their
  relative phase, leaving an exact `S1` of critical gauge classes.

**DERIVED negative:** even the natural quartic polynomial does not pin `D`
to finitely many points.  **OPEN:** a derived `f`, normalization constraint,
or additional geometric term capable of reducing the 122-dimensional
moduli.

## 2. Prize checks

### Hypercharge

There is no uniquely surviving `D`.  On the generic stratum, the stabilizer
of `D` in the bimodule commutant is only the common scalar `u(1)`.  Its
generator acts with one constant charge on the whole witness.  The mixed
and gravitational anomaly equations force that common charge to zero; it
does not yield the nontrivial Standard Model tuple.  Special nongeneric
Diracs have enhanced, choice-dependent commutants, so choosing one of their
abelian generators would insert the answer.

**DERIVED negative on the generic stratum:** Route C does not produce SM
hypercharge.  **OPEN:** a selected nongeneric `D` and a derived
generation-blind commutant generator.

### Generations

On one sheet,

`W|C3 = (10,10,10) = 10 Reg(C3)`.

After KO6 doubling it is `(20,20,20)=20 Reg(C3)`, and each doubled chirality
separately has `(10,10,10)=10 Reg(C3)`.

This equality permits a basis in which `C3` cyclically permutes three
10-dimensional sectors.  It does not canonically provide such sectors:
the character decomposition instead gives three 10-dimensional isotypic
spaces on which `C3` acts by three different scalars, and identifying their
multiplicity spaces requires an arbitrary `U(10)` matching.  Nor is a
10-dimensional sector the derived `M15`/`M16` Standard Model generation
module.

**DERIVED:** the exact isotypic census above.  **PATTERN:** interpreting the
three regular positions as generations.  It is not yet SM-shaped.

### Theorem-chain status

The strengthened chain is

`600-cell -> 2I -> central-parity segregation -> maximal C3/C5 escapes`

`-> C5 seed obstruction and C3 witness -> d0=132`

`-> KO6 J-fixed dimension 132 -> generic gauge/scale moduli dimension 122`.

The chain stops there.  No unique `D`, mass matrix, Yukawa matrix, or
hypercharge tuple is predicted by these cuts.

## Status ledger

### Strengthened

- **DERIVED:** all ten `C3` subgroups are conjugate.
- **DERIVED:** the full legal witness space has `d0=132`, codimension 16 in
  the 148-dimensional self-adjoint odd equivariant space.
- **DERIVED conditional:** KO6 reality on the established double leaves 132
  independent real parameters.
- **DERIVED:** the generic gauge-and-scale quotient has dimension 122.
- **DERIVED negative:** the quartic spectral action has a critical `S1`.
- **DERIVED:** `W|C3=10 Reg(C3)`; the double is `20 Reg(C3)`.

### Downgraded or delimited

- Sheet covariance is derived; arithmetic sigma covariance of arbitrary
  Dirac coefficients is not.
- The moduli dimension 122 describes the generic stratum, not singular
  enhanced-symmetry strata.
- `C3` regularity is a generation **PATTERN**, not a generation derivation.
- A generic surviving `u(1)` is only the common scalar and does not give SM
  hypercharge.

### Open

- a derived arithmetic coefficient lattice and sigma law;
- a derived finite spectral polynomial or other variational functional;
- selection of a nongeneric gauge orbit of `D`;
- hypercharge, a canonical generation splitting, and all physical
  mass/Yukawa data.
