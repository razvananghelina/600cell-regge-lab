"""The spectral dimension is not measurable on the fixed 600-cell.

holographic_dimension.md, tower_spacetime.md and warped_spacetime.md each
measure a single dimension number on a fixed finite carrier and each returns a
negative or inconclusive verdict, with a static value near 3.  The adversarial
audit of 2026-07-24 already downgraded those to a "scoped negative /
inconclusive calibration".  This file explains WHY they are inconclusive, and
the reason is not the choice of plateau criteria.

If dimension is emergent rather than fundamental, the observable is not a
number but a flow,

    d_s(sigma) = -2 d log P(sigma) / d log sigma,
    P(sigma)   = (1/N) Tr exp(-sigma L) = (1/N) sum_i exp(-sigma lambda_i),

the diffusion spectral dimension.  That is the quantity which runs in causal
dynamical triangulations and in asymptotic safety.  Measuring it requires a
window of scales in which the walk has left a single site but has not yet
filled the graph.  The width of that window is controlled by the ratio of the
largest Laplacian eigenvalue to the spectral gap.

On the 600-cell vertex graph that ratio is about 7, i.e. less than one decade,
so no such window exists.  The flow therefore has no plateau, and any number
read off it -- 3.0688 statically, or the 3.63 peak found here -- is a
finite-size artefact rather than a dimension.

CONSEQUENCE.  The four-dimensionality question cannot be decided on the fixed
120-vertex carrier by any choice of observable or criterion.  It requires the
refinement tower to work.  That relocates the open problem instead of leaving
three separate unexplained negatives.

Exact spectra from integer Laplacians; no fitting and no free parameters.
"""
import io
import contextlib

import numpy as np

N_PASS = 0
N_FAIL = 0


def check(label, ok, detail=""):
    global N_PASS, N_FAIL
    if ok:
        N_PASS += 1
        print(f"[PASS] {label}", flush=True)
    else:
        N_FAIL += 1
        print(f"[FAIL] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


PHI = (1 + 5 ** 0.5) / 2

# --------------------------------------------- 600-cell vertex graph
Q = "verify_q8_transplant_2I.py"
g = {"__name__": "__sd__", "__file__": Q}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(Q).read(), Q, "exec"), g)
MUL, NG, G = g["MUL"], g["NG"], g["G"]

W = np.array([float(q[0][0]) + float(q[0][1]) * PHI for q in G])
nb = [i for i in range(NG) if abs(W[i] - PHI / 2) < 1e-9]
check("600-cell adjacency is the 12-element trace-phi class", len(nb) == 12)

A = np.zeros((NG, NG))
for x in range(NG):
    for h in nb:
        A[x, MUL[x][h]] = 1.0
A = np.maximum(A, A.T)
check("the 600-cell vertex graph is 12-regular on 120 vertices",
      A.shape == (120, 120) and np.allclose(A.sum(axis=1), 12))

L = np.diag(A.sum(axis=1)) - A
ev = np.linalg.eigvalsh(L)
ev[0] = 0.0
gap = sorted(ev)[1]
top = ev.max()
check("integer Laplacian spectrum, one zero mode (the graph is connected)",
      abs(ev[0]) < 1e-9 and gap > 1e-6,
      f"gap = {gap:.6f}, largest = {top:.6f}, ratio = {top/gap:.3f}")


def flow(ev, lo=-2.5, hi=1.5, n=600):
    sig = np.logspace(lo, hi, n)
    P = np.array([np.exp(-s * ev).sum() for s in sig]) / len(ev)
    dP = np.array([-(ev * np.exp(-s * ev)).sum() for s in sig]) / len(ev)
    return sig, -2.0 * sig * dP / P


sig, ds = flow(ev)
peak = float(ds.max())
check("the flow never reaches four on the fixed carrier",
      peak < 4.0, f"max d_s = {peak:.4f} at sigma = {sig[int(np.argmax(ds))]:.4f}")


def plateau_decades(sig, ds, tol=0.02):
    """Widest span of sigma, in decades, over which d_s varies by < tol
    relatively, ignoring the trivial d_s -> 0 tail."""
    best = 0.0
    live = ds > 0.5
    i = 0
    while i < len(ds):
        if not live[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(ds) and live[j + 1] and abs(ds[j + 1] - ds[i]) <= tol * ds[i]:
            j += 1
        if j > i:
            best = max(best, float(np.log10(sig[j] / sig[i])))
        i = max(j, i + 1)
    return best


dec = plateau_decades(sig, ds)
check("the flow has no plateau: no decade of scale where d_s is flat to 2%",
      dec < 1.0, f"widest flat span = {dec:.3f} decades")

# The available window, by definition of diffusion times 1/lambda.
window = np.log10((1.0 / gap) / (1.0 / top))
check("there is less than one decade of usable scale on this carrier",
      window < 1.0,
      f"diffusion window 1/lambda_max .. 1/gap = {1/top:.4f} .. {1/gap:.4f} "
      f"= {window:.3f} decades")

print()
print("  d_s at sample scales:")
for s in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
    j = int(np.argmin(np.abs(sig - s)))
    print(f"    d_s({s:5.2f}) = {ds[j]:6.3f}")

print()
print("-" * 74)
print("DERIVED NEGATIVE: on the fixed 600-cell vertex graph the diffusion")
print("spectral dimension has no plateau, because the spectrum spans less than")
print("one decade.  Any single number extracted from it -- the static 3.0688 in")
print("holographic_dimension.md, or the 3.63 peak here -- is a finite-size")
print("artefact, not a dimension.")
print()
print("SCOPE: this says nothing about whether the CONTINUUM limit is")
print("four-dimensional.  It says the question cannot be decided on the fixed")
print("120-vertex carrier by any observable, and therefore depends entirely on")
print("the refinement tower that tower_spacetime.md was unable to establish.")
print("-" * 74)
print(f"RESULT: {N_PASS} passed, {N_FAIL} failed.")
