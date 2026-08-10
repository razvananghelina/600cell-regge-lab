# Codex mission: from the local rigidity no-go to a real theorem

> **Resolved 2026-08-10:** B1 was refuted by an exact noncommutative all-gate
> witness on the fixed carrier. See `chamber_b1_refutation.md`. The prompt
> below is retained as the historical mission specification.

Prepared after `chamber_noncomm_no_go_refutation.md` (18/18) refuted the
proposed commutativity theorem.  Paste section "PROMPT" into a fresh Codex
session with `cwd = /home/razvan/science`.

---

## PROMPT

PROOF MISSION. Hostile standard: attack the framing, do not just execute it.
You have refuted my last two proposals and you were right both times. Do it
again if this one is wrong; a clean refutation beats a forced positive.

### Rules (binding, from CLAUDE.md)

- RULE ZERO: invent nothing; say "I don't know" and check.
- Label everything DERIVED / STRUCTURAL / PATTERN / OPEN.
- Honest negatives are valid deliverables.
- No PDF builds.
- Verifier in `reproducible/`, registered in `run_all.py` (the coverage guard
  exits 2 on any unregistered `verify_*.py`; deliberate exclusions go in
  `DELIBERATELY_SKIPPED` with a reason).
- Interpreter `/home/razvan/science/.venv/bin/python`
  (numpy, scipy, sympy, z3, ortools, networkx).

### Two errors of mine you already corrected -- do not let me repeat them

1. My commutativity theorem omitted **orientability** and **connectedness**.
   `A = M2(C) + C^3` lives exactly in that gap. Any statement I propose from
   here on must list its hypotheses completely before you spend time on it.
2. I claimed the 84 unimodular partition witnesses were strong evidence for a
   general commutativity theorem. They are not: the search parametrised only
   partition algebras, which are commutative by construction. That is circular
   and carries zero evidential weight. **Never cite a search as evidence for a
   claim that the search space could not have falsified.**

### What is established and must not be redone

- `chamber_noncomm_no_go_refutation.md`, 18/18: on the fixed chamber carrier
  with derived `D`, `gamma`, `J`, the algebra `A = M2(C) + C^3` is faithful,
  unital, noncommutative, and satisfies order zero, first order,
  `[gamma, A] = 0`, nonzero one-forms and a unimodular intersection form
  (`Pf = det = 1`). It FAILS metric-dimension-zero orientability and FAILS
  connectedness.
- The valid local no-go, for the exact 30-block `C5` witness: with cell
  dimensions `c_ij = dim(P_i J(P_j) H_+)`, the census over the 90 nonzero
  cells is `70 x 1, 10 x 2, 10 x 3`; the graph joining `i ~ j` whenever
  `c_ij = 1` is connected on all 30 nodes; order zero forces every nonzero
  cell dimension to be a multiple of `n_i n_j`; a unit cell therefore forces
  `n_i = n_j = 1`; connectivity propagates that to all 30 blocks. Hence that
  witness admits no noncommutative amplification retaining the same central
  supports. Scope: says nothing about nonlocal central projectors or different
  Krajewski multiplicities -- which is exactly where the counterexample lives.

### TARGET A -- does the rigidity generalise? (tractable, do this first)

The local no-go rests on one contingent fact: for that one witness, the
**unit-cell graph** (`i ~ j` iff `c_ij = 1`) happens to be connected. Ask
whether that is a coincidence or a consequence.

> **A1 (empirical, cheap, do it before any proof attempt).** For every
> unimodular witness already in hand -- the 84 from the `C5`-invariant
> contraction enumeration, plus the `C^36` PALPABLE and INTEGRAL certificates
> -- compute the cell census and test whether the unit-cell graph is connected.
> Report the distribution. If even one legal witness has a disconnected
> unit-cell graph, A2 is false as stated and you should say so immediately
> rather than trying to prove it.

> **A2 (the theorem, if A1 survives).** For a partition algebra on the chamber
> carrier satisfying first order, orientability and connectedness, is the
> unit-cell graph necessarily connected? If yes, then no such witness admits a
> noncommutative amplification with the same central supports -- a genuine
> generalisation of the local no-go from one witness to the whole class.

Note honestly what A2 would and would not give: it covers **amplifications of
partition algebras only**. It leaves the counterexample's regime untouched.
Do not let the abstract get ahead of that.

### TARGET B -- the regime where the counterexample lives

The counterexample uses nonlocal central projectors and different Krajewski
multiplicities. It also fails orientability and connectedness. So those two
axioms are the untested levers, and the honest question is:

> **B1.** On the fixed chamber carrier with derived `D`, `gamma`, `J`: does a
> unital faithful `*`-algebra satisfying order zero, first order,
> `[gamma, A] = 0`, nonzero one-forms, a nondegenerate intersection form,
> **metric-dimension-zero orientability and connectedness** have to be
> commutative?

That is my previous statement with the two missing hypotheses restored. I do
NOT claim it is true. Treat it as a conjecture with an explicit hypothesis
list, and attack it from both sides:

- try to repair the counterexample -- can `M2(C) + C^3`, or another nonlocal
  algebra, be modified to satisfy orientability and connectedness while
  staying noncommutative? If yes, B1 is refuted and that is the deliverable.
- try to prove it -- the natural route is to show that orientability plus
  connectedness force the central projectors to be local (indicator functions
  of chamber blocks), which would reduce B1 to A2. Whether that reduction is
  even true is the crux; test it before assuming it.

Useful structural facts, verified: `gamma` has eigenvalues `+-1` with sixty
each, so `A` lies in the commutant `M60 + M60`; writing `A = A_+ + A_-` on the
two orientation sheets, `J` exchanges them, so order zero reads
`[A_+, theta(A_-)] = 0`, a **mutual**-commutant condition, not
self-commutation. For `A = M_k` that only gives `k^2 <= 60`, i.e. `k <= 7`.
So order zero alone cannot force commutativity -- consistent with the
counterexample -- and any proof must use first order, orientability or
connectedness. `k <= 7` also bounds the finite case analysis if a direct
proof stalls.

### Constraints

Do not fit `D` or `gamma`; do not enlarge the carrier; do not claim a Standard
Model. Every witness that has ever passed all gates is a partition algebra,
hence commutative, hence `U(1)^k` -- so even a positive answer to B1 closes a
question rather than opening the physics gate. Say that plainly in the note.

### Deliverables

- One markdown note with a complete hypothesis list for every statement and a
  correct status ledger.
- A verifier in `reproducible/`, registered.
- Re-run the FULL suite and report the real pass count. Last independently
  confirmed: **72/72 in 566.9s** before your refutation verifier was added; the
  registry is now 73. Do not trust any count quoted to you, including that one.
