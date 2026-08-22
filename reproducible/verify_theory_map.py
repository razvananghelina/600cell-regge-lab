#!/usr/bin/env python3
"""Verify the authoritative theory map and machine-readable route ledger."""

from __future__ import annotations

import ast
import json
from pathlib import Path
import re
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
JSON_PATH = ROOT / "docs" / "theory_map.json"
MARKDOWN_PATH = ROOT / "docs" / "THEORY_MAP.md"
RUN_ALL_PATH = HERE / "run_all.py"
EXPECTED_LABELS = {"DERIVED", "STRUCTURAL", "PATTERN", "OPEN"}
EXPECTED_STATES = {
    "FOUNDATION",
    "ACCEPTED",
    "BOUNDED_NO_GO",
    "OPEN_GATE",
    "ACTIVE_GATE",
    "METHOD_CONTROL",
    "PATTERN_CONTROL",
}
ID_PATTERN = re.compile(r"^[A-Z][A-Z0-9]*(?:-[A-Z0-9]+)+$")
tests = 0
passed = 0


def check(label: str, condition: object, detail: str = "") -> None:
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def registered_scripts(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "scripts" for target in node.targets):
            continue
        value = ast.literal_eval(node.value)
        if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
            raise TypeError("run_all.py scripts is not a literal string list")
        return value
    raise ValueError("run_all.py has no literal scripts registry")


def has_dependency_cycle(routes_by_id: dict[str, dict]) -> bool:
    temporary: set[str] = set()
    permanent: set[str] = set()

    def visit(route_id: str) -> bool:
        if route_id in permanent:
            return False
        if route_id in temporary:
            return True
        temporary.add(route_id)
        for dependency in routes_by_id[route_id]["depends_on"]:
            if visit(dependency):
                return True
        temporary.remove(route_id)
        permanent.add(route_id)
        return False

    return any(visit(route_id) for route_id in routes_by_id)


data = json.loads(JSON_PATH.read_text())
routes = data.get("routes", [])
route_ids = [route.get("id") for route in routes]
routes_by_id = {route["id"]: route for route in routes if isinstance(route.get("id"), str)}

check(
    "schema and authority are frozen",
    data.get("schema_version") == 1
    and data.get("authority") == "docs/THEORY_MAP.md"
    and data.get("updated") == "2026-08-22",
)

check(
    "route identifiers are nonempty, valid and unique",
    bool(routes)
    and len(route_ids) == len(set(route_ids))
    and all(isinstance(route_id, str) and ID_PATTERN.fullmatch(route_id) for route_id in route_ids),
    f"routes={len(routes)}, distinct={len(set(route_ids))}",
)

bad_labels = [route["id"] for route in routes if route.get("status") not in EXPECTED_LABELS]
bad_states = [route["id"] for route in routes if route.get("route_state") not in EXPECTED_STATES]
check(
    "only binding evidence labels and declared route states occur",
    set(data.get("status_labels", [])) == EXPECTED_LABELS
    and set(data.get("route_states", [])) == EXPECTED_STATES
    and not bad_labels
    and not bad_states,
    f"bad_labels={bad_labels}, bad_states={bad_states}",
)

active = [route["id"] for route in routes if route.get("route_state") == "ACTIVE_GATE"]
check(
    "exactly one active gate agrees with the top-level pointer",
    len(active) == 1
    and data.get("active_gate") == active[0]
    and routes_by_id[active[0]].get("status") == "OPEN",
    f"active={active}, pointer={data.get('active_gate')}",
)

required_list_fields = ("hypotheses", "evidence", "depends_on", "search_terms")
missing_fields = []
for route in routes:
    scalar_ok = all(
        isinstance(route.get(field), str) and route[field].strip()
        for field in ("id", "title", "branch", "status", "route_state", "claim", "scope")
    )
    lists_ok = all(
        isinstance(route.get(field), list)
        and (field == "depends_on" or bool(route[field]))
        and all(isinstance(item, str) and item.strip() for item in route[field])
        for field in required_list_fields
    )
    if not scalar_ok or not lists_ok:
        missing_fields.append(route.get("id", "<missing-id>"))
check(
    "every route states its claim, scope, hypotheses, evidence and aliases",
    not missing_fields,
    str(missing_fields),
)

missing_evidence = []
for route in routes:
    for relative in route["evidence"]:
        if not (ROOT / relative).is_file():
            missing_evidence.append(f"{route['id']}:{relative}")
check(
    "every evidence path exists",
    not missing_evidence,
    str(missing_evidence),
)

bad_dependencies = [
    f"{route['id']}:{dependency}"
    for route in routes
    for dependency in route["depends_on"]
    if dependency not in routes_by_id or dependency == route["id"]
]
cycle = False if bad_dependencies else has_dependency_cycle(routes_by_id)
check(
    "all dependency IDs resolve and the route graph is acyclic",
    not bad_dependencies and not cycle,
    f"bad={bad_dependencies}, cycle={cycle}",
)

bad_aliases = []
for route in routes:
    aliases = route["search_terms"]
    normalized = [alias.casefold().strip() for alias in aliases]
    if len(aliases) < 3 or len(normalized) != len(set(normalized)):
        bad_aliases.append(route["id"])
check(
    "every route has at least three distinct duplicate-search aliases",
    not bad_aliases,
    str(bad_aliases),
)

bad_no_go = [
    route["id"]
    for route in routes
    if route["route_state"] == "BOUNDED_NO_GO"
    and not (
        route["status"] == "DERIVED"
        and isinstance(route.get("kill_scope"), str)
        and route["kill_scope"].strip()
        and isinstance(route.get("reopen_condition"), str)
        and route["reopen_condition"].strip()
    )
]
check(
    "every bounded no-go has a derived scope and reopening condition",
    not bad_no_go,
    str(bad_no_go),
)

bad_open = [
    route["id"]
    for route in routes
    if route["route_state"] in {"OPEN_GATE", "ACTIVE_GATE"}
    and not (
        route["status"] == "OPEN"
        and isinstance(route.get("next_test"), str)
        and route["next_test"].strip()
    )
]
check(
    "every open or active gate names its next falsifiable test",
    not bad_open,
    str(bad_open),
)

markdown = MARKDOWN_PATH.read_text()
check(
    "the human map contains a Mermaid graph and the bounded no-go ledger",
    "```mermaid" in markdown
    and "## No-go ledger and reopening rules" in markdown
    and "## Current decision path" in markdown,
)

missing_markdown_ids = [route_id for route_id in route_ids if f"`{route_id}`" not in markdown]
check(
    "every machine route ID occurs in the human map",
    not missing_markdown_ids,
    str(missing_markdown_ids),
)

required_references = {
    ROOT / "README.md": "docs/THEORY_MAP.md",
    ROOT / "docs" / "README.md": "THEORY_MAP.md",
    ROOT / "docs" / "gravity" / "CURRENT_STATUS.md": "docs/THEORY_MAP.md",
    ROOT / "reproducible" / "README.md": "verify_theory_map.py",
    ROOT / "CLAUDE.md": "docs/THEORY_MAP.md",
}
missing_references = [
    str(path.relative_to(ROOT))
    for path, needle in required_references.items()
    if needle not in path.read_text()
]
check(
    "curated indexes, current status and binding rules reference the map",
    not missing_references,
    str(missing_references),
)

scripts = registered_scripts(RUN_ALL_PATH)
duplicates = sorted({name for name in scripts if scripts.count(name) > 1})
check(
    "the map verifier is registered exactly once and the registry has no duplicates",
    scripts.count("verify_theory_map.py") == 1 and not duplicates,
    f"map_count={scripts.count('verify_theory_map.py')}, duplicates={duplicates}",
)

print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)
