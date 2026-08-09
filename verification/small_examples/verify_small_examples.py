#!/usr/bin/env python3
"""Verify the explicit binary codes for E_4 and E_5 exactly."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def rank_columns(columns: list[int]) -> int:
    """Return the exact column rank over F_2 for bit-packed columns."""
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


def singleton_values(heights: list[int], distance: int) -> tuple[int, ...]:
    return tuple(
        sum(max(0, height - distance + cut + 1) for height in heights[cut:])
        for cut in range(distance)
    )


def verify(filename: str, expected_distribution: Counter[int]) -> None:
    certificate = json.loads((ROOT / filename).read_text())
    basis = certificate["basis_column_masks"]
    heights = certificate["support_column_heights"]
    dimension = certificate["dimension"]
    distance = certificate["minimum_nonzero_rank"]

    assert certificate["field"] == "F2"
    assert len(basis) == dimension
    assert min(singleton_values(heights, distance)) == dimension

    for generator in basis:
        assert len(generator) == len(heights)
        assert all(0 <= value < (1 << height)
                   for value, height in zip(generator, heights))

    words: set[tuple[int, ...]] = set()
    distribution: Counter[int] = Counter()
    for coefficients in range(1 << dimension):
        word = [0] * len(heights)
        for index, generator in enumerate(basis):
            if (coefficients >> index) & 1:
                word = [left ^ right for left, right in zip(word, generator)]
        words.add(tuple(word))
        distribution[rank_columns(word)] += 1

    assert len(words) == 1 << dimension
    assert distribution == expected_distribution
    assert min(rank for rank in distribution if rank) == distance
    print(json.dumps({
        "certificate": filename,
        "dimension": dimension,
        "minimum_nonzero_rank": distance,
        "rank_distribution": dict(sorted(distribution.items())),
        "singleton_values": singleton_values(heights, distance),
        "status": "VERIFIED",
    }, sort_keys=True))


def main() -> None:
    verify("e4_dim2_certificate.json", Counter({0: 1, 3: 3}))
    verify("e5_dim6_certificate.json", Counter({0: 1, 3: 57, 4: 6}))


if __name__ == "__main__":
    main()
