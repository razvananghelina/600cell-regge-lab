# Codex mission: a non-central real structure from the context geometry

Date prepared: 2026-07-28
Paste the whole of section "PROMPT" into a fresh Codex session with
`cwd = /home/razvan/science`.

---

## PROMPT

HUNT MISSION. Hostile standard: attack the framing, do not just execute it.
You have refuted my proposals before (the free-arena dichotomy, the
`A (x) A^op = End(A)` identification, the Hecke/groupoid route, and two of my
over-broad scope claims). Do it again if this one is wrong. A clean refutation
is a better deliverable than a forced positive.

### Rules (binding, from CLAUDE.md)

- RULE ZERO: invent nothing. If you do not know, say "I don't know" and check.
- Label every statement DERIVED / STRUCTURAL / PATTERN / OPEN.
- Honest negatives are valid deliverables.
- No PDF builds.
- Mandatory verifier in `reproducible/`, registered in `run_all.py`. That file
  now has a coverage guard that FAILS with exit code 2 if any `verify_*.py` on
  disk is unregistered. Do not bypass it; register your verifier.
- Interpreter: `/home/razvan/science/.venv/bin/python`
  (numpy 2.5.1, scipy 1.18.0, sympy 1.14.0, z3, ortools 9.15, networkx).
  The old `/tmp/science-python-deps` path is dead — it vanished on reboot and
  silently made the whole suite unrunnable.

### What is already settled — do NOT redo

- **Trivial doubling lemma.** Identical sheets (`pi(a)=diag(pi0,pi0)`,
  `gamma=diag(+1,-1)`, `J` off-diagonal) force the intersection form to vanish,
  for any group and any projections. Confirmed on the real 128x128 Q8 control,
  which therefore fails Poincare duality demonstrably, not merely "untested".
- **The design filter, with three confirmed instances** (trivial double, Galois
  `C^9` double, `M1 (+) M3`): a doubling can rescue Poincare duality only if the
  swap does NOT preserve what the graded trace sees — node dimensions and McKay
  parity. **Apply this filter FIRST to every new candidate and discard early if
  it fails.** It is cheap and it has been right three times.
- **Chamber carrier.** A complete finite real-even triple with unimodular PD via
  `INTEGRAL_SURVIVOR_BITS` (`det = 1`; the `det = 4` PALPABLE certificate would
  have FAILED strict PD, since `det = Pf^2`). But its algebra is forced
  commutative.
- **The forcing mechanism** (independently reverified): the geometric `J` there
  is multiplication by the CENTRAL element, so `J`-conjugation is trivial — all
  60 rotations are fixed by it. Hence `J A J^-1 = A` for any one-sided
  convolution algebra, order zero reduces to `[A,A] = 0`, and `A` must be
  commutative. Therefore `U(A)` is a torus and no non-abelian gauge sector can
  exist on that algebra.

### The one door left open — this mission is about it

The forcing argument requires `J A J^-1 = A`. It says **nothing** about algebras
that are not invariant under `J`-conjugation, and **nothing** about carriers
whose `J` is not central. Today's structural results hand you exactly such
candidates.

Established today, and independently reverified:

- `C10` has exactly **six** conjugates in `2I`. Any two distinct ones intersect
  in exactly the centre `{+-1}`, and any two distinct ones **generate all of
  `2I`** (verified for all 15 pairs).
- `N_2I(C10) = Dic5`, order 20, non-abelian, `C10` of index 2 with all ten
  outside elements of order 4. So the six contexts are the transitive `2I`-set
  `2I/Dic5`, and there are 12 oriented cosets.
- `Ind_{C10}^{2I} chi_1 = 2 + 4 + 6` and `Ind chi_3 = 2' + 4 + 6`, both
  multiplicity free. This is forced: `C10` contains the centre (`c^5 = z`), so
  `chi` is odd there, so only odd irreps can occur; the odd irreps have
  dimensions 2, 2, 4, 6, and the only multiplicity-free decompositions of
  dimension 12 are those two. The amalgamated union is the full odd McKay
  sector, dimension 14, and it does not depend on which of the six conjugate
  `C10` is chosen.
- **Uniqueness worth recording if it is not yet in a note:** among ALL finite
  subgroups of `SU(2)`, `2I` is the only one in which `C10` is non-normal. In
  cyclic groups every subgroup is normal; in binary dihedral groups elements of
  order 10 lie in the cyclic part `<a>` and every subgroup of `<a>` is closed
  under inversion hence normal; `2T` and `2O` have no elements of order 10.
  Verify this and record it — it is stronger than the "exceptional ladder"
  restriction and it survives the `Dic5` counterexample.

### TASK A — main hunt: a NON-CENTRAL real structure

Build carriers from the context geometry above and give them a `J` that is
**not** central, so the commutativity forcing cannot apply.

- **A1.** `C[2I/Dic5] = C^6`, and the oriented `C^12`, with `J` induced by a
  permutation of the six contexts — in particular an ODD permutation — rather
  than by the central element.
- **A2.** `M1 (+) M3` with `J` the Galois swap `2 <-> 2'`.
  **I predict this fails PD by the design filter**, because Galois preserves
  both dimension (`dim 2 = dim 2'`) and parity (both odd), so the graded traces
  cancel in pairs — exactly as they did for the `C^9` Galois double. Confirm or
  refute that prediction explicitly *before* investing further in A2.
  Also check separately: does any one-form distinguish the two sheets? If `D`
  only connects the shared `4 (+) 6`, the sole non-trivial feature of the
  doubling is invisible to the dynamics, and the carrier is non-trivial only on
  paper.
- **A3.** Any carrier where `J` exchanges two distinct contexts. Because two
  distinct contexts meet only in the centre and generate `2I`, such a `J` is
  manifestly non-central, so the forcing theorem does not cover it.

For each candidate the decisive question is single and blunt:

> **Is there a NONCOMMUTATIVE algebra satisfying order zero?**

If no, say so and give the structural reason. If yes, run the full gate list:
first order; KO6 signs against a derived `gamma` and `J`; orientability in
metric dimension zero; connectedness; nonzero inner one-forms; and the
intersection form's rank AND determinant over `Z`, reporting the Pfaffian when
antisymmetric (`det = Pf^2`, so `det = 4` means `Pf = 2`, non-degenerate over
`Q` but NOT unimodular over `Z`, hence strict PD fails).

Do not fit `D` or `gamma`. Do not enlarge a carrier merely to create room. Do
not claim a Standard Model.

### TASK B — the non-`J`-invariant subalgebra door

Independently of A, attack the forcing theorem's hypothesis directly. Order zero
requires `[A, J A^* J^-1] = 0`; the forcing needs `J A J^-1 = A`. Search for
star-subalgebras of `M_120(C)` (or of the relevant carrier's endomorphisms) that
satisfy order zero **without** being `J`-invariant, and ask whether any is
non-commutative. A full classification is likely out of reach — scope it
honestly, for example to algebras generated by small sets of derived operators,
and state the scope. If you can prove that order zero plus non-commutativity
forces `J`-invariance after all, that is a stronger theorem than anything in
Task A and it would close the door permanently. Either outcome is publishable.

### TASK C — cheap SAT query, solvers are now installed

Decision problems are far easier for SAT than the global `dim Omega`
minimisation. The `C^36` chamber partition has trivial `A5` stabiliser (orbit
60), i.e. the founding symmetry is completely broken. Ask z3 or CP-SAT:

> Does ANY legal partition have a non-trivial `A5` stabiliser (`C2`, `C3`, `C5`,
> or larger)?

If the answer is no for every non-trivial subgroup, that is a clean DERIVED
no-go — it would show that breaking `A5` is not a choice but a necessity on that
carrier, and it would correctly scope a claim I earlier made too broadly.

### Deliverables

- One new markdown note per task actually attempted, each with a correct status
  ledger and explicit evasion boundaries.
- Verifiers in `reproducible/`, registered (the coverage guard will enforce it).
- Re-run the FULL suite with the venv interpreter and report the real pass
  count. It was **67/67 in 319.8s** after the ten previously unregistered
  verifiers were added. Do not trust any count quoted to you, including that
  one — re-run it.

### Two warnings about presentation

- Numbers that look like bridges and are not. "14" now means two unrelated
  things: the sum of dimensions `2+2'+4+6`, and the regular rank
  `1^2+2^2+3^2 = 14` of the `C+H+M3` corner. Never place them side by side as
  if related.
- `K_6` is not a derivation. Six mutually equivalent objects give `K_6`
  automatically; degree `= 6-1 = 5` and Laplacian gap `= 6` are restatements of
  the cardinality, which is itself `|2I|/|Dic5| = 120/20`. Say "bootstrap
  closure, not an independent derivation of `a_1 = 5` or `b_1 = 6`", and keep
  that caveat in the final wording — it was written once already and then
  dropped from a polished headline.
