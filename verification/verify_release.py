#!/usr/bin/env python3
"""Verify the exact data shipped with the E6/F2 counterexample.

The checks use only Python's standard library.  They validate file integrity,
the explicit 11-dimensional code, all primary orbit records, the independent
classification and solver summaries, the raw all-kernel summaries, and the
four residual CNF certificates.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MANIFEST = HERE / "audit_manifest.json"
KERNELS_PER_CLASS = 2_794_155


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def binary_rank(rows: list[int], width: int) -> int:
    rows = list(rows)
    rank = 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(rows)) if (rows[index] >> column) & 1),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and ((rows[index] >> column) & 1):
                rows[index] ^= rows[rank]
        rank += 1
    return rank


def matrix_rank_from_columns(columns: list[int]) -> int:
    rows = [
        sum(((column >> row) & 1) << index for index, column in enumerate(columns))
        for row in range(5)
    ]
    return binary_rank(rows, 6)


def verify_dimension_11() -> dict:
    certificate = load_json(HERE / "dim11" / "e6_dim11_certificate.json")
    basis = certificate["basis_column_masks"]
    require(len(basis) == 11, "dimension-11 certificate does not have 11 generators")
    require(all(len(generator) == 6 for generator in basis), "generator width is not six")
    for generator in basis:
        require(all(0 <= value < 32 for value in generator[:4]), "full-column mask out of range")
        require(all(value in (0, 1) for value in generator[4:]), "short-column support violation")

    packed = [sum(column << (5 * index) for index, column in enumerate(row)) for row in basis]
    require(binary_rank(packed, 30) == 11, "certificate generators are dependent")
    distribution: Counter[int] = Counter()
    words: set[int] = set()
    for coefficients in range(1 << 11):
        columns = [0] * 6
        for index, generator in enumerate(basis):
            if (coefficients >> index) & 1:
                columns = [left ^ right for left, right in zip(columns, generator)]
        packed_word = sum(column << (5 * index) for index, column in enumerate(columns))
        words.add(packed_word)
        distribution[matrix_rank_from_columns(columns)] += 1

    expected = {0: 1, 3: 605, 4: 1098, 5: 344}
    require(len(words) == 2048, "certificate does not contain 2,048 distinct words")
    require(dict(sorted(distribution.items())) == expected, "unexpected rank distribution")
    return {
        "dimension": 11,
        "distinct_words": len(words),
        "minimum_nonzero_rank": min(rank for rank in distribution if rank),
        "rank_distribution": {str(key): value for key, value in expected.items()},
    }


def rank_from_column_masks(columns: list[int]) -> int:
    pivots: dict[int, int] = {}
    for value in columns:
        x = value
        while x:
            pivot = x.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = x
                break
            x ^= pivots[pivot]
    return len(pivots)


def singleton_cuts(heights: tuple[int, ...], distance: int) -> tuple[int, ...]:
    return tuple(
        sum(max(0, height - (distance - 1 - j)) for height in heights[j:])
        for j in range(distance)
    )


def verify_small_examples() -> dict:
    specifications = {
        "E_4": ("e4_dim2_certificate.json", (3, 3, 1, 1), 2, {0: 1, 3: 3}),
        "E_5": ("e5_dim6_certificate.json", (4, 4, 4, 1, 1), 6,
                {0: 1, 3: 57, 4: 6}),
    }
    summaries: dict[str, dict] = {}
    for label, (filename, expected_heights, expected_dimension,
                expected_distribution) in specifications.items():
        certificate = load_json(HERE / "small_examples" / filename)
        basis = certificate["basis_column_masks"]
        heights = tuple(certificate["support_column_heights"])
        distance = certificate["minimum_nonzero_rank"]
        require(certificate["field"] == "F2", f"{label}: wrong field")
        require(heights == expected_heights, f"{label}: wrong support")
        require(len(basis) == certificate["dimension"] == expected_dimension,
                f"{label}: wrong certificate dimension")
        require(min(singleton_cuts(heights, distance)) == expected_dimension,
                f"{label}: Singleton bound mismatch")
        require(all(len(generator) == len(heights) for generator in basis),
                f"{label}: generator width mismatch")
        require(all(0 <= value < (1 << height)
                    for generator in basis
                    for value, height in zip(generator, heights)),
                f"{label}: support violation")

        words: set[tuple[int, ...]] = set()
        distribution: Counter[int] = Counter()
        for coefficients in range(1 << expected_dimension):
            columns = [0] * len(heights)
            for index, generator in enumerate(basis):
                if (coefficients >> index) & 1:
                    columns = [left ^ right for left, right in zip(columns, generator)]
            words.add(tuple(columns))
            distribution[rank_from_column_masks(columns)] += 1
        require(len(words) == 1 << expected_dimension,
                f"{label}: dependent generators")
        require(dict(sorted(distribution.items())) == expected_distribution,
                f"{label}: rank distribution mismatch")
        summaries[label] = {
            "dimension": expected_dimension,
            "minimum_nonzero_rank": min(rank for rank in distribution if rank),
            "rank_distribution": {str(rank): count
                                  for rank, count in expected_distribution.items()},
            "singleton_cuts": list(singleton_cuts(heights, distance)),
        }
    return summaries


def verify_row_cone_family() -> dict:
    certificate = load_json(HERE / "dim11" / "e6_dim11_certificate.json")
    basis = certificate["basis_column_masks"]
    heights = tuple(certificate["support_column_heights"])

    for t in range(11):
        distance = 3 + t
        require(min(singleton_cuts(heights, distance)) == 12,
                f"row-cone stage {t}: Singleton bound mismatch")
        require(all(mask < (1 << height) for row in basis for mask, height in zip(row, heights)),
                f"row-cone stage {t}: support violation")

        distribution: Counter[int] = Counter()
        for coefficients in range(1 << len(basis)):
            columns = [0] * len(heights)
            for index, generator in enumerate(basis):
                if (coefficients >> index) & 1:
                    columns = [left ^ right for left, right in zip(columns, generator)]
            distribution[rank_from_column_masks(columns)] += 1
        expected = {0: 1, 3 + t: 605, 4 + t: 1098, 5 + t: 344}
        require(dict(sorted(distribution.items())) == expected,
                f"row-cone stage {t}: rank distribution mismatch")

        coned_basis: list[list[int]] = []
        for index, generator in enumerate(basis):
            tag = [0] * 12
            tag[index] = 1
            coned_basis.append([value << 1 for value in generator] + tag)
        basis = coned_basis
        heights = tuple(height + 1 for height in heights) + (1,) * 12

    return {
        "stages": 11,
        "distance_range": [3, 13],
        "singleton_bound": 12,
        "exact_dimension_lower_certificates": 11,
        "words_per_stage": 2048,
    }


def verify_primary() -> dict:
    specifications = [
        ("field", "standard_d_kernel_orbits.jsonl", "standard_d_all_orbits_cpp.jsonl", 3_240, 105, 189),
        ("II", "class2_kernel_orbits.jsonl", "class2_all_orbits_cpp.jsonl", 26_643, 117, 189),
        ("III", "class3_kernel_orbits.jsonl", "class3_all_orbits_cpp.jsonl", 156_679, 113, 189),
    ]
    summaries: dict[str, dict] = {}
    for label, orbit_name, result_name, expected_count, expected_min, expected_max in specifications:
        orbit_path = HERE / "primary" / orbit_name
        result_path = HERE / "primary" / result_name
        count = mass = nodes = conflicts = 0
        rank_min, rank_max = 10**9, -1
        with orbit_path.open(encoding="utf-8") as orbit_stream, result_path.open(encoding="utf-8") as result_stream:
            for count, pair in enumerate(zip(orbit_stream, result_stream, strict=True), start=1):
                orbit = json.loads(pair[0])
                result = json.loads(pair[1])
                require(result["rep_index"] == count - 1, f"{label}: nonconsecutive result index")
                for key in ("rank2_count", "orbit_size", "span"):
                    require(result[key] == orbit[key], f"{label}: orbit/result mismatch in {key}")
                require(result["status"] == "UNSAT", f"{label}: non-UNSAT primary record")
                mass += orbit["orbit_size"]
                rank_min = min(rank_min, orbit["rank2_count"])
                rank_max = max(rank_max, orbit["rank2_count"])
                nodes += result["nodes"]
                conflicts += result["conflicts"]
        require(count == expected_count, f"{label}: wrong primary orbit count")
        require(mass == KERNELS_PER_CLASS, f"{label}: orbit masses do not cover every kernel")
        require((rank_min, rank_max) == (expected_min, expected_max), f"{label}: rank-two range mismatch")
        summaries[label] = {
            "orbit_representatives": count,
            "coverage_mass": mass,
            "rank2_count_range": [rank_min, rank_max],
            "UNSAT": count,
            "SAT": 0,
            "nodes": nodes,
            "conflicts": conflicts,
        }
    require(sum(item["orbit_representatives"] for item in summaries.values()) == 186_562,
            "wrong aggregate primary orbit count")
    return summaries


def verify_classification() -> dict:
    data = load_json(HERE / "independent" / "spread_classification.json")
    require(data["good_quotient_vectors"] == 2_912, "classification quotient-vector count mismatch")
    require(data["normalized_spreads"] == 19_936, "classification spread count mismatch")
    require(data["lr_orbits"] == 3, "classification did not produce three LR orbits")
    masses = sorted(item["mass"] for item in data["orbits"])
    require(masses == [336, 2_800, 16_800], "classification orbit masses mismatch")
    require(sum(masses) == 19_936, "classification masses do not cover all normalized spreads")
    require(all(item["transpose_orbit"] == index for index, item in enumerate(data["orbits"])),
            "transpose does not fix every classification orbit")
    require(len({item["orbit"] for item in data["named_classes"]}) == 3,
            "named MRD classes do not occupy distinct LR orbits")
    return {"normalized_spreads": 19_936, "lr_orbit_masses": masses}


def verify_independent() -> dict:
    files = sorted(
        path
        for path in (HERE / "independent").glob("orbit_v3_*class*.json")
        if ".aggregate." not in path.name
    )
    require(files, "independent v3 terminal result files are missing")
    by_class: dict[int, list[dict]] = {}
    for path in files:
        data = load_json(path)
        require(data["checker"] == "affine-equation-dpll-48var-v3", f"wrong checker in {path.name}")
        require(data["mode"] == "orbits", f"wrong mode in {path.name}")
        require(data["sat"] == 0 and data["unsat"] == data["processed"], f"non-UNSAT result in {path.name}")
        by_class.setdefault(data["class"], []).append(data)
    require(set(by_class) == {1, 2, 3}, "independent results do not cover all three classes")

    expected = {1: (3_240, 900), 2: (26_643, 108), 3: (156_679, 18)}
    summaries: dict[str, dict] = {}
    for class_number, (orbit_count, automorphisms) in expected.items():
        parts = sorted(by_class[class_number], key=lambda item: item["range_start"])
        cursor = 0
        processed = coverage = unsat = nodes = conflicts = forced = 0
        for part in parts:
            require(part["total_kernels"] == KERNELS_PER_CLASS, "independent kernel total mismatch")
            require(part["automorphism_count"] == automorphisms, "independent automorphism count mismatch")
            require(part["total_orbits"] == orbit_count, "independent orbit count mismatch")
            require(part["range_start"] == cursor, "independent shard ranges have a gap or overlap")
            cursor = part["range_end"]
            processed += part["processed"]
            coverage += part["coverage_mass"]
            unsat += part["unsat"]
            nodes += part["search_nodes"]
            conflicts += part["conflicts"]
            forced += part["forced_equations"]
        require(cursor == orbit_count and processed == orbit_count and unsat == orbit_count,
                "independent shards do not cover every orbit")
        require(coverage == KERNELS_PER_CLASS, "independent shard masses do not cover every kernel")
        summaries[str(class_number)] = {
            "orbit_representatives": processed,
            "coverage_mass": coverage,
            "UNSAT": unsat,
            "SAT": 0,
            "nodes": nodes,
            "conflicts": conflicts,
            "forced_equations": forced,
            "shards": len(parts),
        }
    return summaries


def verify_raw_cuda() -> dict:
    specifications = [
        ("field", "class_field.json", 0, "229d914e2388c642", 15_361_512_744, [105, 189]),
        ("II", "class_II.json", 2, "f53eb945006316ff", 8_639_635_778, [117, 189]),
        ("III", "class_III.json", 3, "f03d7599f5c810b2", 8_182_534_029, [113, 189]),
    ]
    summaries: dict[str, dict] = {}
    for label, name, class_number, digest, nodes, rank_range in specifications:
        data = load_json(HERE / "raw_cuda" / name)
        require(data["class"] == class_number and data["start"] == 0, f"raw {label}: range/class mismatch")
        require(data["count"] == KERNELS_PER_CLASS, f"raw {label}: wrong kernel count")
        status = data["status_counts"]
        require(status == {"UNSAT": KERNELS_PER_CLASS, "SAT": 0, "constraint_overflow": 0, "solver_overflow": 0},
                f"raw {label}: bad terminal status")
        require(data["result_digest_fnv1a64"] == digest, f"raw {label}: digest mismatch")
        require(data["sum_nodes"] == nodes, f"raw {label}: node sum mismatch")
        require(data["constraint_count_range"] == rank_range, f"raw {label}: constraint range mismatch")
        summaries[label] = {
            "raw_kernels": data["count"],
            "UNSAT": status["UNSAT"],
            "SAT": 0,
            "errors": 0,
            "nodes": nodes,
            "digest": digest,
        }

    repeats = [
        load_json(HERE / "raw_cuda" / name)
        for name in ("determinism_repeat_a.json", "determinism_repeat_b.json", "cross_build_sm80_repeat.json")
    ]
    semantic_keys = ("class", "start", "count", "first_span", "last_span", "status_counts",
                     "sum_nodes", "max_nodes", "sum_constraints", "constraint_count_range",
                     "result_digest_fnv1a64", "digest_record_format")
    signatures = [{key: item[key] for key in semantic_keys} for item in repeats]
    require(signatures[0] == signatures[1] == signatures[2], "CUDA deterministic/cross-build replay mismatch")
    selftest = load_json(HERE / "raw_cuda" / "cuda_dpll_selftest.json")
    require(selftest["tests"] == 256 and selftest["errors"] == 0 and selftest["truth_table_mismatches"] == 0,
            "CUDA truth-table self-test failed")
    return summaries


def verify_residuals() -> dict:
    directory = HERE / "structural"
    metadata_files = sorted(directory.glob("residual_*.json"))
    require(len(metadata_files) == 4, "expected four residual metadata files")
    expected_spans = {(1, 2, 3), (2, 5, 7), (256, 512, 768), (256, 3072, 3328)}
    expected_stems = {
        "residual_1_2", "residual_2_5", "residual_256_512", "residual_256_3072"
    }
    exit_summary = load_json(directory / "survivor_solver_exit_codes.json")
    require(set(exit_summary) == expected_stems, "wrong residual exit-code summary entries")
    expected_exit_codes = {
        "cadical": 20,
        "cryptominisat5": 20,
        "drat_trim_original": 0,
        "drat_trim_core": 0,
    }
    spans: set[tuple[int, int, int]] = set()
    compressed_proof_bytes = 0
    for metadata_path in metadata_files:
        data = load_json(metadata_path)
        stem = metadata_path.stem
        require(stem in expected_stems, f"unexpected residual instance {stem}")
        require(exit_summary[stem] == expected_exit_codes, f"bad solver/proof exit codes for {stem}")
        span = tuple(data["span"])
        spans.add(span)
        require(data["class"] == 0 and data["constraints"] == 105,
                f"bad residual metadata in {metadata_path.name}")
        cnf_path = metadata_path.with_suffix(".cnf")
        require(cnf_path.is_file(), f"missing {cnf_path.name}")
        require(sha256(cnf_path) == data["cnf_sha256"], f"CNF hash mismatch for {cnf_path.name}")
        for solver in ("cadical", "cryptominisat"):
            log = directory / f"{stem}.{solver}.log"
            require(log.is_file(), f"missing solver log {log.name}")
            text = log.read_text(encoding="utf-8", errors="replace")
            require("UNSATISFIABLE" in text, f"no terminal UNSAT marker in {log.name}")
        for suffix in ("drat-trim.log", "core-drat-trim.log", "core-check.log"):
            log = directory / f"{stem}.{suffix}"
            require(log.is_file(), f"missing proof-check log {log.name}")
            text = log.read_text(encoding="utf-8", errors="replace")
            require("s VERIFIED" in text, f"no terminal VERIFIED marker in {log.name}")
        proof = directory / f"{stem}.core.drat.gz"
        require(proof.is_file(), f"missing compressed proof {proof.name}")
        decompressed = 0
        with gzip.open(proof, "rb") as stream:
            for block in iter(lambda: stream.read(1 << 20), b""):
                decompressed += len(block)
        require(decompressed > 0, f"empty compressed proof {proof.name}")
        compressed_proof_bytes += proof.stat().st_size
    require(spans == expected_spans, "wrong four residual kernel spans")
    return {
        "instances": 4,
        "constraints_each": 105,
        "compressed_proof_bytes": compressed_proof_bytes,
        "spans": sorted(map(list, spans)),
    }


def verify_manifest() -> dict:
    require(MANIFEST.is_file(), "audit manifest is missing")
    manifest = load_json(MANIFEST)
    require(manifest["schema"] == "etzion-silberstein-e6-f2-release-v2", "unexpected manifest schema")
    files = manifest["files"]
    for relative, record in files.items():
        path = REPO / relative
        require(path.is_file(), f"manifest file missing: {relative}")
        require(path.stat().st_size == record["bytes"], f"size mismatch: {relative}")
        require(sha256(path) == record["sha256"], f"SHA-256 mismatch: {relative}")
    actual = {
        path.relative_to(REPO).as_posix()
        for path in HERE.rglob("*")
        if path.is_file()
        and path != MANIFEST
        and "__pycache__" not in path.parts
        and "build" not in path.parts
    }
    require(actual == set(files), "manifest file set differs from release file set")
    return {"files": len(files)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-manifest", action="store_true", help="verify mathematical records before sealing hashes")
    args = parser.parse_args()
    result = {
        "dimension_11": verify_dimension_11(),
        "small_examples": verify_small_examples(),
        "row_cone_family": verify_row_cone_family(),
        "classification": verify_classification(),
        "primary": verify_primary(),
        "independent": verify_independent(),
        "raw_cuda": verify_raw_cuda(),
        "residuals": verify_residuals(),
    }
    if not args.skip_manifest:
        result["manifest"] = verify_manifest()
    result["status"] = "VERIFIED"
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
