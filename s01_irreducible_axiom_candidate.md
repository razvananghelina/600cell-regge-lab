# S01 Irreducible Axiom Candidate

Scop: sa formulam onest candidatul de axiomă ireductibilă pentru `S01`,
în cazul în care nu găsim un principiu structural mai slab care să derive
`no-branching`.

## Context

Până acum au căzut următoarele direcții:

1. `minimality` naivă;
2. `productive = non-invertible`;
3. `one generator + non-invertible + self-return`;
4. `one generator + self-dual + trivial pointed subcategory`.

În schimb, a rămas în picioare un theorem condițional foarte curat:

- dacă `X \otimes X` nu produce niciun simplu netrivial nou și `1` apare în
  `X \otimes X`, atunci Fibonacci urmează riguros.

Deci problema `S01` s-a redus efectiv la întrebarea:

- este `no-branching` derivabil,
  sau este chiar axioma fondațională minimă?

## Candidate Irreducible Axiom

> **Axiom S01*.**
> There exists a nontrivial simple self-dual object `X` in a rigid semisimple
> monoidal category such that:
> 1. `1` occurs in `X \otimes X`;
> 2. every simple summand of `X \otimes X` is isomorphic either to `1` or to
>    `X`;
> 3. `X` is not invertible.

În formă scurtă:

\[
X \otimes X \in \mathrm{span}_{\mathbb{Z}_{\ge 0}}\{1,X\},
\qquad
1 \subset X \otimes X,
\qquad
X \not\text{ invertible}.
\]

Din asta urmează:

\[
X \otimes X \cong 1 \oplus X.
\]

## Why This Is a Legitimate Minimal Seed

Această axiomă nu spune direct:

- `X \otimes X = 1 \oplus X`

ci spune ceva structural mai puțin specific:

- auto-interacțiunea lui `X` nu ramifică în noi tipuri simple;
- unitatea reapare;
- `X` nu este invertibil.

Totuși, împreună cu rigiditatea, acestea sunt suficiente pentru Fibonacci.

Deci, dacă este acceptată ca axiomă, seed-ul devine extrem de mic.

## Physical Reading

Lectura fizico-structurală naturală este:

- auto-referința pură nu produce complexitate nouă de tip;
- ea se reproduce pe sine și își reintroduce originea;
- nu generează „ramuri” ontologice suplimentare la prima auto-compoziție.

În limbaj scurt:

- `pure self-reference is reproductive, not branching`.

## Why This Might Be Irreducible

Există un motiv serios să credem că `no-branching` ar putea fi ireductibil:

- orice încercare de a-l înlocui cu un criteriu de „minimalitate” riscă să
  reformuleze exact același conținut;
- contraexemplele standard arată că ipoteze mai slabe lasă loc ramificării;
- până acum, tot ce este strict mai slab a eșuat.

Asta nu dovedește ireductibilitatea absolută.

Dar justifică poziția:

- `no-branching` este candidatul principal la axioma primitivă minimă.

## Decision Criterion

Vom accepta `Axiom S01*` ca axiomă fondațională doar dacă una dintre următoarele
devine clară:

1. găsim un `no-go theorem` convingător că niciun principiu substanțial mai
   slab nu poate exclude ramificarea;
2. după o căutare disciplinată, toate candidatele naturale mai slabe cad;
3. orice reformulare mai slabă se dovedește echivalentă, în esență, cu
   `no-branching`.

## Strategic Consequence

Dacă `Axiom S01*` este acceptată, atunci seed-ul nou al întregului framework
devine:

- existența unui obiect simplu auto-dual a cărui auto-interacțiune este
  non-branching și non-invertibilă.

Din acel punct:

- Fibonacci
- `\phi`
- `a_1 = 5`
- restul lanțului

pot curge exact ca înainte.

## Honest Status

Status curent:

- `candidate irreducible axiom`, nu încă decizie finală.

Este însă, în acest moment, cel mai bun candidat la seed fondațional minim.
