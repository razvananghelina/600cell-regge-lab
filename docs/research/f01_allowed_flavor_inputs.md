# F01 Closure: Allowed Flavor Inputs

## Exact claim

Determine the maximal flavor-relevant input set that remains admissible after
the collapse of the gauge route.

## Formal setting

Audit of the exact core and of the already-accepted conditional arithmetic
results, under the rule that no gauge-derived input may be used.

## Inputs used

- [gauge_route_damage_inventory.md](D:\infinity\ToE\science\gauge_route_damage_inventory.md)
- [one_integer_paper_exact_core.tex](D:\infinity\ToE\science\one_integer_paper_exact_core.tex)
- [flavor_first_program.md](D:\infinity\ToE\science\flavor_first_program.md)

## Output status

- `Derived lemma`

## Partition of Inputs

### A. Exact admissible flavor inputs

These can be used directly in the flavor-first program:

1. the Fibonacci seed and bootstrap `a_1 = 5`;
2. the arithmetic field `Z[phi]`;
3. the generation-count theorem:
   exactly three stable unit sectors on the chiral line;
4. the McKay `\widetilde{E}_8` shadow as exact representation-theoretic
   structure;
5. the exact exponent lattice form
   \[
   n = 5a + 6b;
   \]
6. the exact arithmetic norm
   \[
   N(a+b\phi)=a^2+ab-b^2;
   \]
7. the exact constructive `(a,b)` theorem, but only in its declared
   conditional form on the exponent set.

### B. Conditional flavor inputs

These may be used only with explicit labels:

1. the exponent set
   \[
   \{0,3,5,11,11,16,17,19,26\}
   \]
   as an accepted conditional input;
2. the reading of the three stable unit sectors as the three physical families;
3. the reading of lattice exponents as effective suppression charges;
4. the use of `A_5`-golden-ratio flavor literature as analogy, not as a
   derived symmetry of the exact chain.

### C. Forbidden inputs

These may not be used in the flavor-first program:

1. any derived gauge group claim;
2. any use of `alpha`, `alpha_s`, `sin^2 theta_W` as foundational flavor data;
3. any electroweak selector or Higgs/W/Z mass relation as input;
4. any argument that presupposes the failed `fiber -> A_5 -> gauge` route.

## Proof / derivation

The generation theorem, lattice arithmetic, and `(a,b)` construction remain
valid independently of the gauge no-go, because they depend on:

- `a_1 = 5`;
- `Z[phi]`;
- McKay/arithmetic structure;
- explicit integer-lattice constructions.

By contrast, any quantity whose physical reading required a derived gauge
sector is removed from the admissible input list.

Therefore the surviving flavor program is not empty: it has a genuine exact
arithmetic core, but it must be formulated as a flavor-structure program rather
than as a full SM derivation.

## What is still not proved

This step does not yet prove that the admissible data form a compelling flavor
theory. It proves only that the flavor-first program has a nonempty exact input
set that is independent of the failed gauge route.

## Decision

`F01` is closed as `Derived lemma`.
