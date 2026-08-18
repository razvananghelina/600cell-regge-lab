# Prior-art gate: target-independent stationary-root enumeration

Date: 2026-08-16

Nonmonotone base diagnosis: `f9d7ada`.

Status: **completed before evaluating a new root grid**.

## Exact object

At the fixed second lower geometry `a1` and fixed lapse log `r1`, enumerate
all resolved zeros on a preregistered finite interval of

```text
G(b,r1) = mean(five complete pole equations).
```

No canonical momentum target enters the enumeration.  Every root's canonical
pre-momentum and derivative data are recorded only as target-independent
observables.

## KNOWN

Symmetry-reduced Lorentzian Regge cosmologies can have distinct contracting,
expanding and static branches.  A recent frustum model explicitly finds
monotone contracting/expanding branches and matter-dependent existence
conditions:

- Jercher and Steinhaus, *Cosmology in Lorentzian Regge calculus: causality
  violations, massless scalar field and discrete dynamics*,
  <https://arxiv.org/abs/2312.11639>.

Canonical pre/post evolution and time reversal are standard properties of a
discrete action:

- Dittrich and Hoehn, <https://arxiv.org/abs/1108.1974>.

Closed Regge FLRW models have a distinguished time-symmetric slice and may
terminate at null struts:

- Liu and Williams, <https://arxiv.org/abs/1501.07614>.

Multiple numerical roots or time-reversed branches are therefore not by
themselves a new physical discovery.

## CONTROL / OPEN

The committed bracket attempts already show one sign change on the negative
`b` side and a residual below `1e-25` at `b=0`.  They do not enumerate all
roots, because nested symmetric endpoints are not a root-isolation grid.

The new calculation must record the exact grid size, every sampled value,
every isolated sign bracket, every near-zero grid node and the distinct root
count before any branch is selected.

A finite interval/grid cannot prove absence of tangential even-multiplicity
roots or roots outside its domain.  Its result is **DERIVED COMPUTATIONAL ON
THE FROZEN DOMAIN**; analytic global completeness remains **OPEN**.

No located source gives this exact root multiset for the present 600-cell dust
action.  External novelty remains **OPEN**.
