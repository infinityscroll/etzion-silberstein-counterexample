#!/usr/bin/env python3
"""Exact finite checks for the recursively coned E6 counterexample family."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


CERT = Path(__file__).resolve().parents[1] / "dim11" / "e6_dim11_certificate.json"


def rank_columns(columns: list[int]) -> int:
    pivots: dict[int, int] = {}
    for value in columns:
        x = value
        while x:
            p = x.bit_length() - 1
            if p not in pivots:
                pivots[p] = x
                break
            x ^= pivots[p]
    return len(pivots)


def nus(heights: tuple[int, ...], d: int) -> tuple[int, ...]:
    return tuple(
        sum(max(0, h - (d - 1 - j)) for h in heights[j:])
        for j in range(d)
    )


def cone_heights(heights: tuple[int, ...], pad: int = 12) -> tuple[int, ...]:
    return tuple(h + 1 for h in heights) + (1,) * pad


def cone_basis(basis: list[list[int]], pad: int = 12) -> list[list[int]]:
    assert len(basis) <= pad
    ans = []
    for i, row in enumerate(basis):
        appended = [0] * pad
        appended[i] = 1
        ans.append([x << 1 for x in row] + appended)
    return ans


def distribution(basis: list[list[int]]) -> Counter[int]:
    out: Counter[int] = Counter()
    for coeffs in range(1 << len(basis)):
        columns = [0] * len(basis[0])
        for i, gen in enumerate(basis):
            if (coeffs >> i) & 1:
                columns = [x ^ y for x, y in zip(columns, gen)]
        out[rank_columns(columns)] += 1
    return out


def main() -> None:
    raw = json.loads(CERT.read_text())
    basis = raw["basis_column_masks"]
    heights = tuple(raw["support_column_heights"])

    for t in range(11):
        d = 3 + t
        got = distribution(basis)
        expected = Counter({0: 1, 3 + t: 605, 4 + t: 1098, 5 + t: 344})
        assert got == expected
        assert min(nus(heights, d)) == 12
        assert max((x.bit_length() for row in basis for x in row), default=0) <= max(
            heights
        )
        for row in basis:
            assert all(mask < (1 << h) for mask, h in zip(row, heights))
        print(
            json.dumps(
                {
                    "t": t,
                    "rows": max(heights),
                    "columns": len(heights),
                    "distance": d,
                    "singleton_values": nus(heights, d),
                    "rank_distribution": dict(sorted(got.items())),
                    "status": "verified",
                },
                sort_keys=True,
            )
        )
        basis = cone_basis(basis)
        heights = cone_heights(heights)


if __name__ == "__main__":
    main()
