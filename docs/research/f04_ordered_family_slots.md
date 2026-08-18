# F04 Closure: Ordered Family Slots

## Exact claim

Determine the strongest honest reading of the three exact unit sectors as
candidate generations.

## Formal setting

Three-slot theorem on the chiral line
\[
a=1,\qquad b\in\{0,1,2\},
\]
together with the admissible grading
\[
w(a,b)=5a+6b.
\]

## Inputs used

- [f03_suppression_charge_scaffold.md](D:\infinity\ToE\science\f03_suppression_charge_scaffold.md)
- [f02_minimal_flavor_object.md](D:\infinity\ToE\science\f02_minimal_flavor_object.md)
- [one_integer_paper_exact_core.tex](D:\infinity\ToE\science\one_integer_paper_exact_core.tex)

## Output status

- `Derived lemma`
- `No-go warning`

## Exact structural result

On the exact chiral line `a=1`, the three unit sectors are
\[
(1,0),\qquad (1,1),\qquad (1,2).
\]
Applying the admissible grading gives
\[
w(1,0)=5,\qquad w(1,1)=11,\qquad w(1,2)=17.
\]

Therefore the three exact slots are not just counted; they are canonically
ordered by the exact grading:
\[
5 < 11 < 17,
\]
with uniform spacing
\[
11-5=17-11=6.
\]

## Consequence

The exact theorem of three slots plus the admissible grading produce:

- three family slots;
- a canonical internal ordering;
- an exact hierarchical separation between them.

So the slots are structurally differentiated, not merely counted.

## Strongest honest reading

The strongest honest reading at this stage is only:

- `ordered family slots`

and not yet:

- `ordered candidate generations`;
- `physical generations`.

Reason: the ordering
\[
5<11<17
\]
is an arithmetic consequence of the linear grading on the line `a=1`. By
itself, it does **not** supply a physical hierarchy of masses.

To get a mass hierarchy one still needs an additional monotone map such as
\[
w \mapsto m(w),
\]
for example of suppression type `m \sim \epsilon^{w}` or `m \sim \phi^{-w}`.
No such map is derived at `F04`.

## Failure criterion check

The failure criterion for `F04` was:

- if the three slots have no internal structural differentiation, they remain
  only family slots.

That failure does **not** occur:

- the slots are exactly ordered by `w(1,b)=5+6b`.

However, a stronger physical reading fails at this step:

- ordering in `w` is not yet ordering in mass;
- therefore the step cannot promote the slots beyond family-slot language.

## Decision

`F04` is closed in split form:

- `Derived lemma`:
  the exact three slots are canonically ordered family slots with hierarchy
  `5 < 11 < 17`;
- `No-go warning`:
  this ordering does not by itself imply a physical mass hierarchy, because no
  monotone mass map `w \mapsto m(w)` has yet been derived.
