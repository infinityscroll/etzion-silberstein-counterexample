#!/usr/bin/env python3
"""Exhaustive F2 sanity check of the row-cone theorem on tiny diagrams.

Every subspace is enumerated once via its unique reduced-row-echelon basis.
The calculation is intentionally independent of the E6 verification code.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations


def rank_rows(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for value in rows:
        x = value
        while x:
            p = x.bit_length() - 1
            if p not in pivots:
                pivots[p] = x
                break
            x ^= pivots[p]
    return len(pivots)


def word_ranks(heights: tuple[int, ...]) -> list[int]:
    cells = [(r, c) for c, h in enumerate(heights) for r in range(h)]
    out = []
    for word in range(1 << len(cells)):
        rows = [0] * max(heights)
        for bit, (r, c) in enumerate(cells):
            if (word >> bit) & 1:
                rows[r] |= 1 << c
        out.append(rank_rows(rows))
    return out


def rref_bases(n: int, k: int):
    if k == 0:
        yield ()
        return
    for pivots in combinations(range(n), k):
        free = [
            (i, j)
            for i, p in enumerate(pivots)
            for j in range(p + 1, n)
            if j not in pivots
        ]
        for mask in range(1 << len(free)):
            rows = [1 << p for p in pivots]
            for bit, (i, j) in enumerate(free):
                if (mask >> bit) & 1:
                    rows[i] |= 1 << j
            yield tuple(rows)


def gaussian_binomial(n: int, k: int) -> int:
    numerator = 1
    denominator = 1
    for i in range(k):
        numerator *= (1 << (n - i)) - 1
        denominator *= (1 << (k - i)) - 1
    return numerator // denominator


def valid_span(basis: tuple[int, ...], ranks: list[int], d: int) -> bool:
    span = [0]
    for v in basis:
        span += [x ^ v for x in span]
    return all(ranks[x] >= d for x in span[1:])


@lru_cache(maxsize=None)
def kappa(heights: tuple[int, ...], d: int) -> int:
    nbits = sum(heights)
    ranks = word_ranks(heights)
    for k in range(nbits, -1, -1):
        if any(valid_span(basis, ranks, d) for basis in rref_bases(nbits, k)):
            return k
    raise AssertionError("zero subspace was not found")


def nus(heights: tuple[int, ...], d: int) -> tuple[int, ...]:
    return tuple(
        sum(max(0, h - (d - 1 - j)) for h in heights[j:])
        for j in range(d)
    )


def partitions(max_height: int, max_columns: int, max_cells: int):
    def rec(prefix: tuple[int, ...], ceiling: int):
        if prefix:
            yield prefix
        if len(prefix) == max_columns:
            return
        for h in range(ceiling, 0, -1):
            if sum(prefix) + h <= max_cells:
                yield from rec(prefix + (h,), h)

    yield from rec((), max_height)


def main() -> None:
    for n in range(9):
        for k in range(n + 1):
            assert sum(1 for _ in rref_bases(n, k)) == gaussian_binomial(n, k)
    print("RREF enumeration counts verified through ambient dimension 8")

    tested = 0
    for heights in partitions(max_height=4, max_columns=4, max_cells=6):
        for d in range(1, min(max(heights), len(heights)) + 1):
            b = min(nus(heights, d))
            if b <= 0:
                continue
            cone = tuple(h + 1 for h in heights) + (1,) * b
            # Exhaustively enumerate both ambient subspace lattices only when
            # the coned matrix space is at most eight-dimensional.
            if sum(cone) > 8:
                continue
            a = kappa(heights, d)
            cone_a = kappa(cone, d + 1)
            cone_b = min(nus(cone, d + 1))
            assert cone_b == b
            assert cone_a == a
            tested += 1
            print(
                f"D={heights} d={d} cells={sum(heights)} "
                f"bound={b} kappa={a} R={cone} "
                f"Rbound={cone_b} Rkappa={cone_a}"
            )
    print(f"VERIFIED {tested} tiny diagram-distance pairs exhaustively")


if __name__ == "__main__":
    main()
