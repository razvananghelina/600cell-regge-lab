#!/usr/bin/env python3
"""Guard the focused GitHub-facing layout of the 600-cell Regge lab."""

from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ALLOWED_ROOT_MARKDOWN = {"README.md", "CLAUDE.md"}
ALLOWED_ROOT_ENTRIES = {
    ".gitignore",
    "CLAUDE.md",
    "README.md",
    "commons",
    "docs",
    "reproducible",
}
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

tracked = subprocess.check_output(
    ["git", "ls-files"], cwd=ROOT, text=True
).splitlines()
public_entries = {path.split("/", 1)[0] for path in tracked}
check(
    "the tracked-facing root contains only the focused Regge layout",
    public_entries <= ALLOWED_ROOT_ENTRIES,
    str(sorted(public_entries - ALLOWED_ROOT_ENTRIES)),
)

gravity = sorted((ROOT / "docs" / "gravity").glob("*.md"))
check(
    "the gravity documentation collection is populated",
    len(gravity) >= 1,
    f"gravity={len(gravity)}",
)

legacy_paths = (
    ROOT / "legacy",
    ROOT / "docs" / "research",
    ROOT / "commons" / "constants.py",
)
check(
    "legacy fitted-theory trees are absent from the current layout",
    not any(path.exists() for path in legacy_paths),
    str([str(path.relative_to(ROOT)) for path in legacy_paths if path.exists()]),
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
    ROOT / "docs" / "gravity" / "CURRENT_STATUS.md",
    ROOT / "docs" / "gravity" / "gravity_600cell_finite_height_classification_result.md",
    ROOT / "docs" / "gravity" / "gravity_600cell_finite_height_composition_result.md",
    ROOT / "docs" / "gravity" / "gravity_600cell_finite_height_selector_result.md",
    ROOT / "docs" / "gravity" / "gravity_600cell_finite_height_third_slab_result.md",
    ROOT / "docs" / "gravity" / "gravity_600cell_finite_height_fourth_slab_result.md",
    ROOT / "docs" / "gravity" / "gravity_600cell_tick_scale_covariance_result.md",
    ROOT / "docs" / "gravity" / "gravity_600cell_projected_refinement_acceleration_comparison_result.md",
)
missing_runtime = [
    str(path.relative_to(ROOT)) for path in required_runtime_notes
    if not path.is_file()
]
check(
    "all current result-chain notes remain present",
    not missing_runtime,
    str(missing_runtime),
)

retained_verifiers = sorted(HERE.glob("verify_*.py"))
non_gravity = [
    path.name for path in retained_verifiers
    if path.name != "verify_documentation_layout.py"
    and not path.name.startswith("verify_gravity_")
]
check(
    "the public verifier inventory is gravity-only",
    not non_gravity,
    str(non_gravity),
)

print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)
