#!/usr/bin/env python3
"""Seal every verification artifact into a deterministic SHA-256 manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUTPUT = HERE / "audit_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def included(path: Path) -> bool:
    return (
        path.is_file()
        and path != OUTPUT
        and "__pycache__" not in path.parts
        and "build" not in path.parts
        and path.suffix not in {".o", ".pyc"}
    )


def main() -> None:
    paths = sorted((path for path in HERE.rglob("*") if included(path)), key=lambda path: path.as_posix())
    files = {
        path.relative_to(REPO).as_posix(): {"bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in paths
    }
    manifest = {
        "schema": "etzion-silberstein-e6-f2-release-v2",
        "claim": (
            "The maximum dimension of a binary rank-distance-three linear code "
            "on the Ferrers diagram with column heights (5,5,5,5,1,1) is 11; "
            "row-cone propagation gives bound 12 and exact optimum 11 at every "
            "minimum rank distance at least three."
        ),
        "field": "F2",
        "support_column_heights": [5, 5, 5, 5, 1, 1],
        "minimum_rank_distance": 3,
        "etzion_silberstein_bound": 12,
        "exact_maximum_dimension": 11,
        "upper_bound_audits": {
            "primary_orbit_representatives": 186_562,
            "independent_orbit_representatives": 186_562,
            "raw_kernels": 8_382_465,
            "raw_unsat": 8_382_465,
            "raw_sat": 0,
            "raw_errors": 0,
            "structural_surviving_orbits": 4,
            "structural_surviving_mass": 70,
        },
        "lower_bound_certificate": {
            "dimension": 11,
            "words": 2_048,
            "minimum_nonzero_rank": 3,
            "rank_distribution": {"0": 1, "3": 605, "4": 1098, "5": 344},
        },
        "row_cone_family": {
            "field": "F2",
            "minimum_rank_distances": "all integers d >= 3",
            "singleton_bound": 12,
            "exact_maximum_dimension": 11,
            "finite_certificate_stages_checked": 11,
            "finite_check_distance_range": [3, 13],
        },
        "files": files,
    }
    OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": str(OUTPUT), "files": len(files), "sha256": sha256(OUTPUT)}, sort_keys=True))


if __name__ == "__main__":
    main()
