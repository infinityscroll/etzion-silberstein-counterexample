#!/usr/bin/env python3
"""Independent Python checker for a fixed dual-spread kernel/lift instance."""

from __future__ import annotations

import argparse
import json
import time

from xor_brancher_core import XorBrancher, form_mask, parity, vector_basis, verify_assignment

D2 = [32771, 16398, 8196, 4104, 2058, 1029, 524, 260, 140, 72, 47, 22]
D3 = [32781, 16396, 8206, 4106, 2054, 1029, 516, 268, 140, 70, 41, 20]


def matrix_of(message: int, basis: list[int]) -> int:
    matrix = 0
    for i, generator in enumerate(basis):
        if (message >> i) & 1:
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


def solve(code_class: int, l1: int, l2: int) -> dict[str, object]:
    ambient = D2 if code_class == 2 else D3
    words = [x for x in range(4096) if parity(l1 & x) == parity(l2 & x) == 0]
    kernel_basis = vector_basis(words)
    assert len(words) == 1024 and len(kernel_basis) == 10
    coordinates = {}
    for u in range(1024):
        x = 0
        for i, b in enumerate(kernel_basis):
            if (u >> i) & 1:
                x ^= b
        coordinates[x] = u
    constraints = set()
    for x in words:
        matrix = matrix_of(x, ambient)
        if x == 0 or rank4(matrix) != 2:
            continue
        n1, n2 = null_basis(matrix)
        a, b = form_mask(coordinates[x], n1), form_mask(coordinates[x], n2)
        constraints.add(tuple(sorted((a, b))))
    started = time.time()
    solver = XorBrancher(sorted(constraints))
    result = solver.solve()
    record: dict[str, object] = {
        "class": code_class, "l1": l1, "l2": l2, "constraints": len(constraints),
        "nodes": solver.nodes, "conflicts": solver.conflicts, "elapsed_s": time.time() - started,
        "status": "SAT" if result is not None else "UNSAT",
    }
    if result is not None:
        assignment, _ = result
        assert verify_assignment(assignment, sorted(constraints))
        record["lift_assignment"] = assignment
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("code_class", type=int, choices=(2, 3))
    parser.add_argument("l1", type=int)
    parser.add_argument("l2", type=int)
    args = parser.parse_args()
    print(json.dumps(solve(args.code_class, args.l1, args.l2), sort_keys=True))


if __name__ == "__main__":
    main()
