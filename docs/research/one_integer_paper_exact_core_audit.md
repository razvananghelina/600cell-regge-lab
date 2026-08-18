# Exact Core Audit

Scop: audit pas-cu-pas al versiunii `one_integer_paper_exact_core.tex`.

Regula:

- Nu pastram formulare mai tari decat justifica demonstratia.
- Pentru fiecare pas trebuie sa existe un status unic:
  - `Theorem`
  - `Conditional theorem`
  - `Computational fact`
  - `Axiom`
  - `Open`
  - `Too strong / remove`

Checklist standard pentru fiecare pas:

1. Care este afirmatia exacta?
2. Ce intrari foloseste?
3. Unde este partea fortata?
4. O putem transforma in teorema / lema conditionala / no-go?
5. Daca nu, o retrogradam la axiom sau o scoatem.

## Step 01

- `ID`: S01
- `Location`: Seed assumption before bootstrap
- `Claim`: Fibonacci fusion rule `\tau \otimes \tau = 1 \oplus \tau` is the starting seed.
- `Current status`: `Axiom`
- `Current status`: `Irreducible axiom`
- `Why forced`:
  - nu este derivat dintr-un principiu fizic anterior in paper;
  - este alegerea de pornire a intregului cadru.
- `Target form`:
  - sa ramana axiom explicit, nu teorema falsa;
  - dar sa fie constrans cat mai tare de teoreme de clasificare externe.
- `Needed work`:
  - folosim clasificarea lui Ostrik pentru rank 2: singurele cazuri sunt pointed si Fibonacci;
  - asta da o minimalitate restrictionata reala, chiar daca nu elimina axioma.
- `Decision`: keep as explicit axiom, but now supported by a rank-2 uniqueness theorem in the restricted class.
  Final closure:
  - after systematic elimination of weaker candidates, the seed is accepted as
    an irreducible axiom in the stronger no-branching form recorded in
    [s01_closure_decision.md](D:\infinity\ToE\science\s01_closure_decision.md)
  - this is stronger and cleaner than merely postulating Fibonacci by name
  - the theorem layer then recovers Fibonacci from rigidity plus classification
  Active foundational subprogram:
  - [zero_postulate_program.md](D:\infinity\ToE\science\zero_postulate_program.md)
  - [rank2_self_reference_theorem.md](D:\infinity\ToE\science\rank2_self_reference_theorem.md)
  - [productive_self_reference_axioms.md](D:\infinity\ToE\science\productive_self_reference_axioms.md)
  Current sharp target:
  - do not try to re-derive Fibonacci coefficients directly;
  - instead try to derive `rank 2 + non-pointed` from a minimal productive self-reference principle, then use classification.

## Step 02

- `ID`: S02
- `Location`: Bootstrap selection of `a_1 = 5`
- `Claim`: The equation `d_1(a_1) = \phi(a_1)` has unique positive-integer solution `a_1 = 5`.
- `Current status`: `Theorem`
- `Why forced`:
  - partea matematica nu pare fortata;
  - singura vulnerabilitate este interpretarea ei ca vacuum selection in sens fizic.
- `Target form`:
  - theorem matematic exact;
  - physical reading optional si separata.
- `Needed work`:
  - resolved by reformulating the theorem without hidden self-reference:
    `phi` now comes from `S01`, and the theorem asks for the unique integer
    `n` with `d_1(n)=phi`.
- `Decision`: closed as theorem; see
  [s02_bootstrap_closure.md](D:\infinity\ToE\science\s02_bootstrap_closure.md)

## Step 03

- `ID`: S03
- `Location`: Passage from `a_1 = 5` to 600-cell / `2I`
- `Claim`: `a_1 = 5` fixes the binary icosahedral / 600-cell realization.
- `Current status`: `Derived lemma` / `Scoped no-go theorem`
- `Why forced`:
  - exista un pas de selectie a realizarii geometrice;
  - este plauzibil si natural, dar nu inca teoremizat in core.
- `Target form`:
  - ori theorem de selectie intr-o clasa bine definita;
  - ori conditional statement: "if one chooses the `2I` / 600-cell realization".
- `Needed work`:
  - formulare precisa a clasei de alternative;
  - separare intre criterii intrinseci si criterii importate;
  - posibila tinta rezonabila:
    - `within the class of regular 4D polytopes compatible with Z[phi], 2I / McKay-E8, and a 12-dimensional local candidate sector, the 600-cell is unique`.
- `Decision`: closed in the strongest honest form; see
  [s03_realization_closure.md](D:\infinity\ToE\science\s03_realization_closure.md)
  - `Derived lemma`:
    if one seeks a regular convex 4-polytope realization carrying the
    bootstrap field and the intrinsic `H4/E8` package, then the 600-cell is
    uniquely selected by `|V|=120=|2I|` and local degree `12`.
  - `Scoped no-go theorem`:
    full absolute selection from `a_1 = 5` alone is not proved and, inside the
    regular-convex class, is blocked by the `H4` dual ambiguity until the
    vertex-count and local-degree selectors are imposed.

## Step 04

- `ID`: S04
- `Location`: Spectral data of the 600-cell
- `Claim`: exact adjacency/Laplacian spectrum with multiplicities in `Z[\phi]`.
- `Current status`: `Computational fact`
- `Why forced`:
  - nu pare fortat, dar trebuie sa fie clar ca este calcul finit exact.
- `Target form`:
  - proposition with reproducible exact computation.
- `Needed work`:
  - optional: enrich with a concise proof strategy or exact characteristic polynomial factorization.
- `Decision`: closed as exact computational fact; see
  [s04_spectrum_closure.md](D:\infinity\ToE\science\s04_spectrum_closure.md)

## Step 05

- `ID`: S05
- `Location`: McKay shadow
- `Claim`: McKay graph of `2I` is affine `E_8`.
- `Current status`: `Theorem`
- `Why forced`:
  - nu este fortat, este standard.
- `Target form`:
  - exact standard theorem.
- `Needed work`:
  - resolved: treat as standard theorem plus local computational confirmation.
- `Decision`: closed; see
  [s05_mckay_closure.md](D:\infinity\ToE\science\s05_mckay_closure.md)

## Step 06

- `ID`: S06
- `Location`: Hopf-fiber operator input
- `Claim`: chosen Hopf fibration and fiber adjacency define `A_fiber`.
- `Current status`: `Computational fact` / `Derived uniform statement`
- `Why forced`:
  - exista o alegere de fibratie;
  - trebuie clarificat daca rezultatele sunt independente de alegere.
- `Target form`:
  - "for any Hopf fibration ..." if true;
  - otherwise conditional on chosen fibration.
- `Needed work`:
  - test / prove fibration-independence for kernel results used in core.
- `Decision`: closed; see
  [s06_hopf_closure.md](D:\infinity\ToE\science\s06_hopf_closure.md)
  - exhaustive finite enumeration over the order-10-subgroup Hopf fibrations yields exactly 6 distinct fibrations;
  - across all 6, the verified stable invariants are:
    - unique nontrivial coefficient `b_1 = 6`,
    - `dim ker(Box) = 9`,
    - `ker(Box)` lies in the same adjacency-spectral sector `E_A(12) ⊕ E_A(6phi) ⊕ E_A(6phi')`, equivalently `rho_0 ⊕ rho_1 ⊕ rho_8`,
    - spectral-gap ratio `lambda_1(L_cross)/lambda_1(L_fiber) = a_1 = 5`;
  - no wider claim of full fibration-invariance is made.

## Step 07

- `ID`: S07
- `Location`: exact wave operator coefficient
- `Claim`: unique coefficient `b_1 = 6`.
- `Current status`: `Computational fact`
- `Why forced`:
  - the selection mechanism is strongly constrained by the commutative
    decomposition, but the full proof of uniqueness is not yet recorded
    cleanly in the exact-core manuscript;
  - at present, uniqueness is verified exactly on the six discrete fibrations,
    but not yet promoted to a general theorem in the paper.
- `Target form`:
  - exact theorem from joint spectral analysis of `A` and `A_fiber`.
- `Needed work`:
  - optional future upgrade: prove cleanly why `c=6` is the only nontrivial
    integer compatible with the joint spectrum, not just the only one found in
    verification.
- `Decision`: closed as computational fact on the verified six-fibration class;
  see
  [s07_wave_coefficient_closure.md](D:\infinity\ToE\science\s07_wave_coefficient_closure.md)

## Step 08

- `ID`: S08
- `Location`: edge kernel and `A_5` decomposition
- `Claim`: exact 13-d kernel, 12-d nontrivial sector, factorization through `A_5`.
- `Current status`: `Computational fact` / `Theorem`
- `Why forced`:
  - the old strong reading risked conflating two different 12-dimensional
    spaces:
    - the nontrivial edge-kernel sector inside `R^{720}`;
    - the fiber permutation module inside `R^{12}`;
  - this identification is not automatic and needed to be kept open.
- `Target form`:
  - weak uniform closure across all 6 fibrations;
  - no identification of the two 12-dimensional spaces unless an explicit
    equivariant map is constructed.
- `Needed work`:
  - completed by exact verification across all six fibrations with full
    character decomposition on the edge kernel and on the fiber permutation
    module.
- `Decision`: closed in weak uniform form; see
  [s08_edge_closure.md](D:\infinity\ToE\science\s08_edge_closure.md)
  - `Computational fact`:
    for every one of the six fibrations,
    `ker(Box_1)=rho_0 \oplus 2 rho_5` and has dimension 13;
  - `Theorem`:
    the action on fiber labels factors through `A_5`;
  - `Computational fact`:
    for every one of the six fibrations, the fiber permutation module is
    `1 \oplus 3 \oplus 3' \oplus 5`;
  - `Open`:
    no canonical identification is yet proved between the 12-dimensional
    edge-kernel sector and the 12-dimensional fiber permutation module.

## Step 09

- `ID`: S09
- `Location`: unique compact Lie-algebra candidate
- `Claim`: from the 12-d adjoint decomposition one gets `u(1) + su(2) + su(3)`.
- `Current status`: `Conditional theorem` / `Blocked by S08b no-go`
- `Why forced`:
  - requires the existence of a compact connected Lie group whose adjoint restricts that way;
  - theorem is if-then, not emergent gauge theory;
  - moreover, the step-by-step chain now has a negative result:
    the quotient-compatible `A_5` bridge from the fiber permutation module to
    the nontrivial edge-space sector is obstructed.
- `Target form`:
  - keep as conditional theorem in the manuscript;
  - do not use in the active derivation chain before the bridge step is closed.
- `Needed work`:
  - if `S09` is to survive in the main chain, it now needs a different input
    than the fiber permutation module.
- `Decision`: keep as conditional theorem in the manuscript, but block it in
  the exact step-by-step chain unless a different source of the
  12-dimensional local sector is found.

## Step 09b

- `ID`: S09b
- `Location`: direct edge-endomorphism gauge route
- `Claim`: the exact nontrivial edge sector `2 rho_5` might itself generate the
  target gauge algebra through `End_{2I}(2 rho_5)`.
- `Current status`: `No-go theorem`
- `Why forced`:
  - this was the only clean gauge-recovery route independent of the failed
    fiber-to-`A_5` bridge.
- `Target form`:
  - compute the Schur type of `rho_5` and read off the canonical compact Lie
    algebra on the multiplicity space.
- `Needed work`:
  - completed; see
    [s09_edge_endomorphism_nogo.md](D:\infinity\ToE\science\s09_edge_endomorphism_nogo.md)
- `Decision`:
  - the Frobenius-Schur indicator of `rho_5` is `-1`, so
    `End_{2I}(rho_5)=H` and `End_{2I}(2 rho_5)=M_2(H)`;
  - the canonical compact Lie algebra is `sp(2)`, not the SM gauge algebra;
  - therefore this direct edge-sector route fails as well.

## Step 10

- `ID`: S10
- `Location`: generation-count theorem
- `Claim`: three unit sectors on `a=1`.
- `Current status`: `Theorem`
- `Why forced`:
  - theorem itself not forced;
  - physical reading as fermion generations would be forced.
- `Target form`:
  - exact arithmetic theorem only.
- `Needed work`:
  - clarify chirality input if retained.
- `Decision`: keep theorem, avoid stronger reading.

## Step 11

- `ID`: S11
- `Location`: constructive `(a,b)` uniqueness
- `Claim`: unique assignment for a given exponent set.
- `Current status`: `Conditional theorem`
- `Why forced`:
  - depends on importing a chosen exponent set.
- `Target form`:
  - explicit if-then theorem.
- `Needed work`:
  - either derive exponent set internally in core or keep clearly conditional.
- `Decision`: keep as conditional theorem for now.

## Step 12

- `ID`: S12
- `Location`: exact scalar response
- `Claim`: `B^+ d_0 = d_0 \Delta_0^+` and consequences.
- `Current status`: `Theorem`
- `Why forced`:
  - none if kept discrete.
- `Target form`:
  - exact theorem.
- `Needed work`:
  - none urgent.
- `Decision`: keep.

## Step 13

- `ID`: S13
- `Location`: simplicial spectral-action coefficients
- `Claim`: exact triple `(2640,14880,55920)` and reduced `(11,62,233)`.
- `Current status`: `Computational fact`
- `Why forced`:
  - no issue if treated as exact finite data.
- `Target form`:
  - exact proposition.
- `Needed work`:
  - optional derivation details.
- `Decision`: keep.

## Step 14

- `ID`: S14
- `Location`: continuum interpretations excluded from core
- `Claim`: couplings, physical masses, gravity, cosmology are not in exact core.
- `Current status`: `Decision`
- `Why forced`:
  - these would be forced if promoted prematurely.
- `Target form`:
  - explicit exclusion list.
- `Needed work`:
  - maintain discipline.
- `Decision`: keep excluded until upgraded honestly.
