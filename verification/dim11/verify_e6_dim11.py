#!/usr/bin/env python3
"""Verify the explicit binary [E_6,11,3] Ferrers-diagram code."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CERTIFICATE = ROOT / "e6_dim11_certificate.json"


def rank_columns(columns: list[int], row_count: int) -> int:
    """Exact column rank over F_2 for bit-packed columns."""
    work = columns[:]
    rank = 0
    for bit in range(row_count):
        pivot = next(
            (index for index in range(rank, len(work)) if (work[index] >> bit) & 1),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(len(work)):
            if index != rank and ((work[index] >> bit) & 1):
                work[index] ^= work[rank]
        rank += 1
    return rank


def main() -> None:
    certificate = json.loads(CERTIFICATE.read_text())
    basis = certificate["basis_column_masks"]
    heights = certificate["support_column_heights"]
    dimension = certificate["dimension"]
    row_count = max(heights)

    assert certificate["field"] == "F2"
    assert dimension == len(basis) == 11
    assert heights == [5, 5, 5, 5, 1, 1]

    for generator in basis:
        assert len(generator) == len(heights)
        for column, height in zip(generator, heights):
            assert 0 <= column < (1 << height)

    words: set[tuple[int, ...]] = set()
    distribution: Counter[int] = Counter()
    for coefficients in range(1 << dimension):
        word = [0] * len(heights)
        for index, generator in enumerate(basis):
            if (coefficients >> index) & 1:
                word = [left ^ right for left, right in zip(word, generator)]
        words.add(tuple(word))
        rank = rank_columns(word, row_count)
        distribution[rank] += 1
        if coefficients:
            assert rank >= certificate["minimum_nonzero_rank"]

    assert len(words) == 1 << dimension
    expected = {int(rank): count for rank, count in certificate["rank_distribution"].items()}
    assert dict(sorted(distribution.items())) == expected
    assert sum(expected.values()) == 1 << dimension
    assert min(rank for rank in distribution if rank) == 3

    print(json.dumps({
        "status": "VERIFIED",
        "dimension": dimension,
        "codewords": len(words),
        "minimum_nonzero_rank": 3,
        "rank_distribution": dict(sorted(distribution.items())),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
