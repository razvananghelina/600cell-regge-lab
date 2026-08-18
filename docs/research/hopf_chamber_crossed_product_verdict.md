# Canonical chamber lift of the six-axis crossed product

Date: 2026-08-11

Protocol commit: `a5c485d`.

Registered verifier:
`reproducible/verify_hopf_chamber_crossed_product.py`.
Targeted result: `25/25` in 2.8 seconds.

No matter character, mass, coupling or Standard-Model target was used.  The
full verifier suite was not run, following the user's instruction.

## Headline

There is a canonical faithful representation of

```text
B_R = R(A5/D5) crossed_product A5
    = M6(R)+M6(R)+M12(R)+M12(R)
```

on each free 60-chamber orientation sheet.  This is a genuine positive
construction missing from the earlier repository-state audit.  It is not,
however, a real spectral triple: order zero is impossible for **every** real
structure on the fixed doubled representation, and the derived chamber
adjacency also fails connectedness.

The result is therefore a **DERIVED CANONICAL-LIFT NO-GO**, not a physical
discovery or a Standard-Model construction.

## Canonical representation

The exact icosahedron has 12 vertices, 30 edges, 20 faces and 120 complete
flags.  Its 60 rotations split the flags into two free orbits of size 60.
Sending a chamber to the antipodal axis containing its vertex gives the
equivariant projection

```text
p : A5/1 -> A5/D5.
```

On either sheet the representation is fixed without a basepoint:

```text
pi(delta_x) e_c = 1_(p(c)=x) e_c,
pi(u_g) e_c     = e_(g c).
```

The 360 matrices `pi(delta_x u_g)` have pairwise disjoint ten-entry supports
and are exactly linearly independent.  The image therefore has the full
dimension 360 and the representation is faithful.  Its real Wedderburn
module multiplicities are `(1,1,2,2)` for simple-module dimensions
`(6,6,12,12)`.

This corrects, without contradicting, the scoped statement in
`hopf_six_existing_operator_lift_verdict.md`: that audit proved that none of
the seven operators already committed at that time supplied a lift.  The
chamber-to-axis covariant action tested here is a new, explicitly
preregistered construction.

## Universal order-zero obstruction on this representation

The exact commutant of the represented algebra has real dimension 10 on one
sheet and 40 on the doubled 120-state carrier.  If an invertible antilinear
real structure `J` obeyed order zero, then

```text
J pi(B_R) J^-1  subset  pi(B_R)'.
```

Conjugation by `J` preserves faithfulness and real vector-space dimension.
The left-hand side would therefore be a 360-dimensional algebra inside a
40-dimensional vector space.  This is impossible.

Thus:

- **DERIVED:** no real structure whatsoever can make this fixed faithful
  doubled chamber representation satisfy order zero;
- the conclusion does not depend on restricting `J` to a geometric
  permutation;
- the exhaustive geometric census is an independent consistency check, not
  the logical basis of the universal conclusion.

For that census, all 60 improper icosahedral symmetries reverse `gamma`,
preserve the chamber adjacency `D`, and normalize the represented algebra.
Exactly 16 square to `+1` (identity times inversion and 15 order-two
rotations times inversion).  Order zero passes for `0/60`, hence for `0/16`
of the `J^2=+1` candidates.  The 360-by-360 basis census contains 39,360
ordered noncommuting pairs.

First order was deliberately not evaluated: the preregistered protocol only
reaches it after order zero survives, and none does.

## Independent connectedness obstruction

The derived 3-regular chamber adjacency has nonzero inner one-forms, but it
is `A5`-invariant.  Consequently the entire represented group algebra
`R[A5]` commutes with it.

The exact commutator map on the 360-dimensional crossed-product image has
rank 300.  A modular rank lower bound meets the explicit 60-dimensional
group-algebra kernel, certifying over the rationals that

```text
dim ker(a -> [D,pi(a)]) = 60.
```

Connectedness requires dimension one, so it fails independently of order
zero.

## Hostile framing audit

1. **Scope.**  This is not a no-go for every representation of `B_R`, every
   enlarged carrier, or every alternative algebra.  It closes the natural
   faithful representation selected by the chamber-to-axis projection on
   the fixed 120-state chamber carrier.
2. **No fitted basis.**  No base chamber, right-regular identification or
   unitary conjugation was selected.  Such choices would introduce new
   symmetry-breaking data.
3. **Why freeness was insufficient.**  The free chamber orbit removes the
   stabilizer obstruction and successfully gives faithfulness, but
   faithfulness makes the represented algebra much larger than its
   commutant.  Freeness solves the lift and simultaneously exposes the
   order-zero obstruction.
4. **Not novelty in the literature.**  The dimension obstruction is an
   elementary finite-dimensional representation-theoretic consequence.  I
   do not know that the particular chamber construction has prior literature
   novelty, and no novelty search was used as evidence.
5. **No physics gate opened.**  The positive result is an exact algebraic
   bridge.  Since it cannot satisfy order zero or connectedness, it does not
   yet define the noncommutative geometry needed for matter physics.

## Status ledger

- **DERIVED:** the two chamber sheets are free `A5` orbits of size 60.
- **DERIVED:** chamber-to-antipodal-axis incidence canonically defines a
  faithful 360-dimensional crossed-product image on each sheet.
- **DERIVED:** its commutant dimensions are 10 per sheet and 40 doubled.
- **DERIVED UNIVERSAL WITHIN THE FIXED REPRESENTATION:** every possible `J`
  fails order zero by the dimension inequality `360 > 40`.
- **DERIVED:** all 60 geometric improper symmetries normalize the same image;
  all fail order zero.
- **DERIVED:** nonzero one-forms exist, but the `D`-commutant has dimension
  60, so connectedness fails.
- **OPEN:** whether a different, independently selected carrier or a smaller
  algebra has a faithful bimodule with a sufficiently large opposite
  commutant and a geometrically selected connected Dirac operator.
- **NO TARGET COMPARISON:** the result has no fitted physical comparison.

## Programme consequence

The missing link is now sharply located.  The six-axis crossed product can be
put on the chamber geometry canonically and faithfully; the failure is not a
missing lift.  The failure is bimodule balance: the same 120 states do not
have enough commutant room for a faithful opposite action, and the invariant
adjacency retains the full group-algebra symmetry.

Any continuation must change at least one of the carrier, the algebra, or the
representation, and must derive that change independently.  Merely choosing
a different `J` or changing coefficients in the existing chamber adjacency
cannot evade either exact obstruction.
