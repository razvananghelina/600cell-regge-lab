# Adversarial protocol: generic-velocity leading map

Date: 2026-08-21

Primary registered implementation: `dfa2688`  
Primary correction protocol: `51bd4ce`  
Corrected primary implementation: `a5ebcea`

Status: frozen after the corrected primary verifier returned
`GENERIC_VELOCITY_LEADING_REPARAMETRIZATION` (`11/11`), before constructing
or running an adversarial verifier.

## 1. Claim under attack

For symbolic real velocity `v`, positive interval factor `s`, positive mass
ratio `mu` and `e->0+`, the primary result claims:

```text
S/(s e) -> L0(v,mu),
2F/(s e) -> C(v,mu),
p_pre    -> p(v),
```

with no `s` dependence, one positive constraint branch `mu(v)`, even mass,
odd momentum and exact leading one-versus-two composition.

## 2. Mechanically independent derivation

Do not reuse the primary primitive-limit or Hamilton--Jacobi derivation.

1. Reconstruct the complete unexpanded cellular action.
2. Differentiate it first with respect to `rho` and `L_minus`.
3. Substitute the exact linear tangent path

   ```text
   L_minus=1,
   L_plus=1+s*v*e,
   rho=s^2*e^2.
   ```

   Prove separately that it has the same leading tangent as the primary
   exponential path.
4. Ask SymPy for the direct limits of `2F/(s e)` and `p_pre` before reading
   the primary formulas.
5. Normalize only the following exact positive-radical factorizations, with
   `u=v^2>=0`:

   ```text
   9u^4+120u^3+592u^2+1280u+1024
     =(u+4)^2(3u+8)^2,

   9u^3+84u^2+256u+256
     =(u+4)(3u+8)^2,

   3u^3+32u^2+112u+128
     =(u+4)^2(3u+8),

   27u^4+360u^3+1776u^2+3840u+3072
     =3(u+4)^2(3u+8)^2,

   3u^5+50u^4+331u^3+1088u^2+1776u+1152
     =(u+3)^2(u+4)^2(3u+8).
   ```

6. Only after the normalized direct expressions are frozen, compare them
   with the primary artifact hash
   `8ded36f1fa00307fcb23369c25290c9f5bd701709762d6a865437c2507eabfc9`.

The direct derivative/limit plus radical normalization is the decisive
mechanical difference from the primary action-limit/HJ route.

## 3. Branch and symmetry checks

Independently solve the direct constraint for `mu`.  Require:

- exactly one branch because the coefficient of `mu` is nonzero;
- positivity for every real `v` from exact cosine bounds;
- mass even and momentum odd;
- the `v->0` static mass and zero momentum controls;
- exact `s` independence before substituting a numerical interval factor.

## 4. New arbitrary-precision controls

Use 100 decimals and points not used by the primary verifier:

```text
v in {-7/10,3/10,13/10},
s in {3/4,1/3},
e in {1/300,1/600}.
```

Require the same frozen first-order error test as the primary protocol:
resolved errors decrease with halving order in `[0.8,1.2]`; errors below
`1e-70` count as exact at precision; mixed resolution is `OPEN`.

## 5. Hostile controls

1. Delete the lateral-boundary `rho` derivative contribution while retaining
   the lateral one.  The leftover velocity-dependent term must be nonzero for
   every registered nonzero control velocity.
2. Replace the constraint-selected mass by `mu(v)+1/10`.  The direct leading
   constraint must shift by exactly `-4*pi/5`.

These controls show that the cancellation and mass branch are active rather
than tautological.

## 6. Outcomes

- `GENERIC_VELOCITY_LEADING_REPARAMETRIZATION_ADVERSARIALLY_CORROBORATED`
  only if every exact, branch, symmetry, new precision and hostile gate
  passes.
- `PRIMARY_GENERIC_VELOCITY_RESULT_REFUTED` if the normalized exact direct
  limits have a certified nonzero difference from the primary formulas.
- `GENERIC_VELOCITY_ADVERSARIAL_DISAGREEMENT` for any unresolved limit,
  radical branch or numerical gate.  The primary result then remains
  unaccepted.

No outcome establishes next-order composition, a fundamental tick or an
absolute time unit.  Only the targeted adversarial verifier will be run.

