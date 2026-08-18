#!/usr/bin/env python3
"""Guard the GitHub-facing Markdown layout and its curated entry links."""

from pathlib import Path
import re
import sys
from urllib.parse import unquote


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ALLOWED_ROOT_MARKDOWN = {"README.md", "CLAUDE.md"}
CURATED_INDEXES = (ROOT / "README.md", ROOT / "docs" / "README.md")
LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


root_markdown = {path.name for path in ROOT.glob("*.md")}
check(
    "only README.md and the binding CLAUDE.md remain at repository root",
    root_markdown == ALLOWED_ROOT_MARKDOWN,
    str(sorted(root_markdown)),
)

gravity = sorted((ROOT / "docs" / "gravity").glob("*.md"))
research = sorted((ROOT / "docs" / "research").glob("*.md"))
check(
    "both documentation collections are populated",
    len(gravity) >= 1 and len(research) >= 1,
    f"gravity={len(gravity)}, research={len(research)}",
)

broken = []
for source in CURATED_INDEXES:
    for match in LINK.finditer(source.read_text()):
        raw = match.group(1).strip()
        target = raw[1:-1] if raw.startswith("<") and raw.endswith(">") else raw
        target = target.split("#", 1)[0]
        if not target or target.startswith(("http://", "https://", "mailto:", "#")):
            continue
        resolved = (source.parent / unquote(target)).resolve()
        if not resolved.exists():
            broken.append(f"{source.relative_to(ROOT)} -> {raw}")
check(
    "all curated local Markdown links resolve",
    not broken,
    "; ".join(broken),
)

required_runtime_notes = (
    ROOT / "docs" / "gravity" / "gravity_hamiltonian_constraint_gap_protocol.md",
    ROOT / "docs" / "gravity" / "gravity_metric_phase_space_canonicity_protocol.md",
    ROOT / "docs" / "gravity" / "gravity_tent_move_regge_protocol.md",
    ROOT / "docs" / "gravity" / "gravity_time_slab_canonicity_protocol.md",
    ROOT / "docs" / "research" / "dimension_reconciliation.md",
    ROOT / "docs" / "research" / "spectral_action_discrete_theorem.md",
    ROOT / "docs" / "research" / "hopf_whitney_metric_selection_result.md",
    ROOT / "docs" / "research" / "prediction_provenance_ledger.md",
    ROOT / "docs" / "research" / "round_a2_transverse_hessian_protocol.md",
    ROOT / "docs" / "research" / "round_regge_spectral_action_sign_protocol.md",
    ROOT / "docs" / "research" / "smooth_derham_a4_stabilization_protocol.md",
)
missing_runtime = [str(path.relative_to(ROOT)) for path in required_runtime_notes
                   if not path.is_file()]
check(
    "all verifier-consumed Markdown sources remain present",
    not missing_runtime,
    str(missing_runtime),
)

print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)
