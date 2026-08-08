"""Small exact GF(2) brancher used by ``audit_dual_spread.py``."""

from __future__ import annotations


def parity(value: int) -> int:
    return value.bit_count() & 1


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


def form_mask(coordinates: int, null_vector: int) -> int:
    mask = 0
    for component in range(4):
        if (null_vector >> component) & 1:
            mask |= coordinates << (10 * component)
    return mask


class XorBrancher:
    def __init__(self, constraints: list[tuple[int, int]]) -> None:
        self.constraints = constraints
        self.nodes = 0
        self.conflicts = 0

    @staticmethod
    def reduce(
        mask: int, rhs: int, equations: dict[int, tuple[int, int]]
    ) -> tuple[int, int]:
        while mask:
            pivot = mask.bit_length() - 1
            equation = equations.get(pivot)
            if equation is None:
                break
            mask ^= equation[0]
            rhs ^= equation[1]
        return mask, rhs

    @classmethod
    def add(
        cls,
        equations: dict[int, tuple[int, int]],
        mask: int,
        rhs: int,
    ) -> dict[int, tuple[int, int]] | None:
        mask, rhs = cls.reduce(mask, rhs, equations)
        if mask == 0:
            return equations if rhs == 0 else None
        result = equations.copy()
        result[mask.bit_length() - 1] = (mask, rhs)
        return result

    def solve(
        self, equations: dict[int, tuple[int, int]] | None = None
    ) -> tuple[int, dict[int, tuple[int, int]]] | None:
        if equations is None:
            equations = {}
        self.nodes += 1
        chosen: tuple[int, int, int, int] | None = None
        chosen_key: tuple[int, int] | None = None
        for left, right in self.constraints:
            reduced_left, value_left = self.reduce(left, 0, equations)
            reduced_right, value_right = self.reduce(right, 0, equations)
            if (reduced_left == 0 and value_left == 1) or (
                reduced_right == 0 and value_right == 1
            ):
                continue
            if (
                reduced_left == 0
                and value_left == 0
                and reduced_right == 0
                and value_right == 0
            ):
                self.conflicts += 1
                return None
            key = (
                int(reduced_left != 0) + int(reduced_right != 0),
                reduced_left.bit_count() + reduced_right.bit_count(),
            )
            if chosen_key is None or key < chosen_key:
                chosen_key = key
                chosen = (reduced_left, value_left, reduced_right, value_right)

        if chosen is None:
            assignment = 0
            for pivot in sorted(equations):
                mask, rhs = equations[pivot]
                value = rhs ^ parity(assignment & (mask & ((1 << pivot) - 1)))
                if value:
                    assignment |= 1 << pivot
            return assignment, equations

        left, value_left, right, value_right = chosen
        if left == 0:
            next_equations = self.add(equations, right, 1 ^ value_right)
            return None if next_equations is None else self.solve(next_equations)
        if right == 0:
            next_equations = self.add(equations, left, 1 ^ value_left)
            return None if next_equations is None else self.solve(next_equations)

        left_branch = self.add(equations, left, 1 ^ value_left)
        if left_branch is not None:
            result = self.solve(left_branch)
            if result is not None:
                return result
        right_branch = self.add(equations, left, value_left)
        if right_branch is not None:
            right_branch = self.add(right_branch, right, 1 ^ value_right)
            if right_branch is not None:
                result = self.solve(right_branch)
                if result is not None:
                    return result
        return None


def verify_assignment(assignment: int, constraints: list[tuple[int, int]]) -> bool:
    return all(
        parity(assignment & left) or parity(assignment & right)
        for left, right in constraints
    )
