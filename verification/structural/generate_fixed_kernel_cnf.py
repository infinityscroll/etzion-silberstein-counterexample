#!/usr/bin/env python3
"""Generate a compact independent CNF for one fixed kernel/lift instance."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from generate_cnf import CNF, encode_parity


D = {
    0: [33825, 17048, 10692, 40034, 35397, 42200, 19754, 53988, 36015, 51832, 42828, 29802],
    2: [32771, 16398, 8196, 4104, 2058, 1029, 524, 260, 140, 72, 47, 22],
    3: [32781, 16396, 8206, 4106, 2054, 1029, 516, 268, 140, 70, 41, 20],
}


def parity(x: int) -> int:
    return x.bit_count() & 1


def vector_basis(vectors: list[int]) -> list[int]:
    pivots: dict[int, int] = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot in pivots:
                vector ^= pivots[pivot]
            else:
                pivots[pivot] = vector
                break
    return [pivots[pivot] for pivot in sorted(pivots, reverse=True)]


def matrix_of(message: int, basis: list[int]) -> int:
    matrix = 0
    for bit, generator in enumerate(basis):
        if (message >> bit) & 1:
            matrix ^= generator
    return matrix


def rows_of(matrix: int) -> list[int]:
    return [(matrix >> (4 * row)) & 15 for row in range(4)]


def rank4(matrix: int) -> int:
    rows = rows_of(matrix)
    rank = 0
    for bit in range(4):
        pivot = next((i for i in range(rank, 4) if (rows[i] >> bit) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(rank + 1, 4):
            if (rows[i] >> bit) & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def null_basis(matrix: int) -> list[int]:
    null = [v for v in range(1, 16) if all(parity(row & v) == 0 for row in rows_of(matrix))]
    result = vector_basis(null)
    assert len(result) == 2
    return result


def constraints(code_class: int, l1: int, l2: int) -> list[tuple[int, int]]:
    ambient = D[code_class]
    words = [x for x in range(4096) if parity(l1 & x) == parity(l2 & x) == 0]
    basis = vector_basis(words)
    assert len(words) == 1024 and len(basis) == 10
    coordinates = {}
    for u in range(1024):
        x = 0
        for i, generator in enumerate(basis):
            if (u >> i) & 1:
                x ^= generator
        coordinates[x] = u
    result = set()
    for x in words:
        matrix = matrix_of(x, ambient)
        if x == 0 or rank4(matrix) != 2:
            continue
        masks = []
        for null in null_basis(matrix):
            mask = 0
            for component in range(4):
                if (null >> component) & 1:
                    mask |= coordinates[x] << (10 * component)
            masks.append(mask)
        result.add(tuple(sorted(masks)))
    return sorted(result)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("code_class", type=int, choices=(0, 2, 3))
    parser.add_argument("l1", type=int)
    parser.add_argument("l2", type=int)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    clauses = constraints(args.code_class, args.l1, args.l2)
    cnf = CNF(args.output)
    variables = [cnf.var() for _ in range(40)]
    for left, right in clauses:
        left_var = encode_parity(cnf, [variables[i] for i in range(40) if (left >> i) & 1])
        right_var = encode_parity(cnf, [variables[i] for i in range(40) if (right >> i) & 1])
        cnf.clause(left_var, right_var)
    variable_count, clause_count = cnf.variables, cnf.clauses
    cnf.finish()
    metadata = {
        "class": args.code_class, "span": [args.l1, args.l2, args.l1 ^ args.l2],
        "constraints": len(clauses), "variables": variable_count, "clauses": clause_count,
        "cnf_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
    }
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata, sort_keys=True))


if __name__ == "__main__":
    main()
