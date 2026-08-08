#!/usr/bin/env python3
"""Compare original NETLIB MPS files with their p_* presolved counterparts."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def parse_mps(filename: Path) -> dict:
    section = ""
    objective_row = None
    row_type: dict[str, str] = {}
    columns: dict[str, dict[str, float]] = defaultdict(dict)
    objective: dict[str, float] = defaultdict(float)
    rhs: dict[str, float] = defaultdict(float)
    bounds: dict[str, list[float | bool]] = {}
    for raw in filename.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("*"):
            continue
        fields = line.split()
        keyword = fields[0].upper()
        if keyword in {"NAME", "ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS", "ENDATA"}:
            section = keyword
            continue
        if section == "ROWS":
            kind, name = fields[0], fields[1]
            if kind == "N":
                objective_row = name
            else:
                row_type[name] = kind
        elif section == "COLUMNS":
            if len(fields) >= 3 and fields[1].strip("'") == "MARKER":
                continue
            column = fields[0]
            for position in range(1, len(fields) - 1, 2):
                row = fields[position]
                value = float(fields[position + 1])
                if row == objective_row:
                    objective[column] += value
                else:
                    columns[column][row] = columns[column].get(row, 0.0) + value
        elif section == "RHS":
            for position in range(1, len(fields) - 1, 2):
                rhs[fields[position]] = float(fields[position + 1])
        elif section == "BOUNDS":
            kind = fields[0].upper()
            column = fields[2]
            value = float(fields[3]) if len(fields) > 3 else 0.0
            lower, upper, has_lower, has_upper = bounds.get(
                column, [0.0, math.inf, True, False]
            )
            if kind == "FR":
                lower, upper, has_lower, has_upper = -math.inf, math.inf, False, False
            elif kind in {"LO", "LI"}:
                lower, has_lower = value, True
            elif kind in {"UP", "UI"}:
                upper, has_upper = value, True
            elif kind == "FX":
                lower, upper, has_lower, has_upper = value, value, True, True
            elif kind == "MI":
                lower, has_lower = -math.inf, False
            elif kind == "PL":
                upper, has_upper = math.inf, False
            elif kind == "BV":
                lower, upper, has_lower, has_upper = 0.0, 1.0, True, True
            bounds[column] = [lower, upper, has_lower, has_upper]
    all_columns = set(columns) | set(objective) | set(bounds)
    for column in all_columns:
        columns.setdefault(column, {})
        bounds.setdefault(column, [0.0, math.inf, True, False])
    rows: dict[str, dict[str, float]] = {name: {} for name in row_type}
    for column, entries in columns.items():
        for row, value in entries.items():
            if value != 0.0 and row in rows:
                rows[row][column] = value
    return {
        "rows": rows,
        "row_type": row_type,
        "columns": dict(columns),
        "objective": dict(objective),
        "rhs": dict(rhs),
        "bounds": bounds,
    }


def different(left: float, right: float) -> bool:
    return abs(left - right) > 1e-10 * (1.0 + abs(left) + abs(right))


def forcing_rows(model: dict) -> tuple[set[str], set[str]]:
    """Return statically forcing rows and the columns fixed by their extreme."""
    rows: set[str] = set()
    columns: set[str] = set()
    for row_name, entries in model["rows"].items():
        if not entries:
            continue
        minimum = maximum = 0.0
        minimum_finite = maximum_finite = True
        for column, coefficient in entries.items():
            lower, upper, has_lower, has_upper = model["bounds"][column]
            if coefficient > 0.0:
                if has_lower:
                    minimum += coefficient * lower
                else:
                    minimum_finite = False
                if has_upper:
                    maximum += coefficient * upper
                else:
                    maximum_finite = False
            else:
                if has_upper:
                    minimum += coefficient * upper
                else:
                    minimum_finite = False
                if has_lower:
                    maximum += coefficient * lower
                else:
                    maximum_finite = False
        rhs = model["rhs"].get(row_name, 0.0)
        kind = model["row_type"][row_name]
        scale = 1.0 + abs(rhs)
        if minimum_finite:
            scale += abs(minimum)
        if maximum_finite:
            scale += abs(maximum)
        force_minimum = (
            minimum_finite
            and kind in {"L", "E"}
            and abs(minimum - rhs) <= 1e-9 * scale
        )
        force_maximum = (
            maximum_finite
            and kind in {"G", "E"}
            and abs(maximum - rhs) <= 1e-9 * scale
        )
        if force_minimum or force_maximum:
            rows.add(row_name)
            columns.update(entries)
    return rows, columns


def zero_cost_singleton_inequality_columns(model: dict) -> set[str]:
    """Find singleton inequality columns fixable at a relaxing finite bound."""
    eligible: set[str] = set()
    for column, entries in model["columns"].items():
        if len(entries) != 1 or model["objective"].get(column, 0.0) != 0.0:
            continue
        row, coefficient = next(iter(entries.items()))
        kind = model["row_type"].get(row)
        lower, upper, has_lower, has_upper = model["bounds"][column]
        del lower, upper
        if kind == "L":
            finite = has_lower if coefficient > 0.0 else has_upper
        elif kind == "G":
            finite = has_upper if coefficient > 0.0 else has_lower
        else:
            finite = False
        if finite:
            eligible.add(column)
    return eligible


def singleton_equality_classes(model: dict) -> tuple[set[str], set[str], set[str]]:
    """Classify equality singleton columns by the number of surviving bounds."""
    removable: set[str] = set()
    projectable: set[str] = set()
    ranged: set[str] = set()
    for column, entries in model["columns"].items():
        if len(entries) != 1:
            continue
        row_name, pivot = next(iter(entries.items()))
        if model["row_type"].get(row_name) != "E" or pivot == 0.0:
            continue
        constant = model["rhs"].get(row_name, 0.0) / pivot
        minimum = maximum = constant
        minimum_finite = maximum_finite = True
        for other, coefficient in model["rows"][row_name].items():
            if other == column:
                continue
            multiplier = -coefficient / pivot
            lower, upper, has_lower, has_upper = model["bounds"][other]
            if multiplier > 0.0:
                if has_lower:
                    minimum += multiplier * lower
                else:
                    minimum_finite = False
                if has_upper:
                    maximum += multiplier * upper
                else:
                    maximum_finite = False
            else:
                if has_upper:
                    minimum += multiplier * upper
                else:
                    minimum_finite = False
                if has_lower:
                    maximum += multiplier * lower
                else:
                    maximum_finite = False
        lower, upper, has_lower, has_upper = model["bounds"][column]
        scale = 1.0 + abs(constant)
        if minimum_finite:
            scale += abs(minimum)
        if maximum_finite:
            scale += abs(maximum)
        lower_implied = not has_lower or (
            minimum_finite and minimum >= lower - 1e-9 * scale
        )
        upper_implied = not has_upper or (
            maximum_finite and maximum <= upper + 1e-9 * scale
        )
        if lower_implied and upper_implied:
            removable.add(column)
        elif lower_implied or upper_implied:
            projectable.add(column)
        else:
            ranged.add(column)
    return removable, projectable, ranged


def fingerprint(original: dict, reduced: dict, name: str) -> dict[str, int | str]:
    original_rows = set(original["rows"])
    reduced_rows = set(reduced["rows"])
    original_columns = set(original["columns"])
    reduced_columns = set(reduced["columns"])
    removed_rows = original_rows - reduced_rows
    removed_columns = original_columns - reduced_columns
    common_rows = original_rows & reduced_rows
    common_columns = original_columns & reduced_columns

    forcing, forcing_columns = forcing_rows(original)
    zero_cost_singletons = zero_cost_singleton_inequality_columns(original)
    singleton_removable, singleton_projectable, singleton_ranged = (
        singleton_equality_classes(original)
    )

    fixed_removed = 0
    singleton_removed = 0
    free_singleton_removed = 0
    zero_cost_removed = 0
    for column in removed_columns:
        lower, upper, has_lower, has_upper = original["bounds"][column]
        degree = len(original["columns"][column])
        fixed_removed += int(has_lower and has_upper and lower == upper)
        singleton_removed += int(degree == 1)
        free_singleton_removed += int(degree == 1 and not has_lower and not has_upper)
        zero_cost_removed += int(original["objective"].get(column, 0.0) == 0.0)

    removed_row_singletons = sum(len(original["rows"][row]) == 1 for row in removed_rows)
    removed_row_doubletons = sum(len(original["rows"][row]) == 2 for row in removed_rows)
    changed_rhs = sum(
        different(original["rhs"].get(row, 0.0), reduced["rhs"].get(row, 0.0))
        for row in common_rows
    )
    changed_objective = sum(
        different(
            original["objective"].get(column, 0.0),
            reduced["objective"].get(column, 0.0),
        )
        for column in common_columns
    )
    changed_bounds = sum(
        any(
            left != right
            if isinstance(left, bool)
            else different(left, right)
            for left, right in zip(
                original["bounds"][column], reduced["bounds"][column]
            )
        )
        for column in common_columns
    )
    changed_coefficients = 0
    added_coefficients = 0
    deleted_common_coefficients = 0
    changed_rows = 0
    for row in common_rows:
        left = original["rows"][row]
        right = reduced["rows"][row]
        row_changed = False
        for column in common_columns & (set(left) | set(right)):
            a = left.get(column, 0.0)
            b = right.get(column, 0.0)
            if a == 0.0 and b != 0.0:
                added_coefficients += 1
                row_changed = True
            elif a != 0.0 and b == 0.0:
                deleted_common_coefficients += 1
                row_changed = True
            elif different(a, b):
                changed_coefficients += 1
                row_changed = True
        changed_rows += int(row_changed)

    original_nnz = sum(len(row) for row in original["rows"].values())
    reduced_nnz = sum(len(row) for row in reduced["rows"].values())
    return {
        "model": name,
        "original_rows": len(original_rows),
        "presolved_rows": len(reduced_rows),
        "removed_rows": len(removed_rows),
        "original_columns": len(original_columns),
        "presolved_columns": len(reduced_columns),
        "removed_columns": len(removed_columns),
        "original_nnz": original_nnz,
        "presolved_nnz": reduced_nnz,
        "fixed_removed_columns": fixed_removed,
        "singleton_removed_columns": singleton_removed,
        "free_singleton_removed_columns": free_singleton_removed,
        "zero_cost_removed_columns": zero_cost_removed,
        "removed_singleton_rows": removed_row_singletons,
        "removed_doubleton_rows": removed_row_doubletons,
        "forcing_rows": len(forcing),
        "removed_forcing_rows": len(forcing & removed_rows),
        "forcing_columns": len(forcing_columns),
        "removed_forcing_columns": len(forcing_columns & removed_columns),
        "zero_cost_singleton_inequality_columns": len(zero_cost_singletons),
        "removed_zero_cost_singleton_inequality_columns": len(
            zero_cost_singletons & removed_columns
        ),
        "singleton_equality_removable_columns": len(singleton_removable),
        "removed_singleton_equality_removable_columns": len(
            singleton_removable & removed_columns
        ),
        "singleton_equality_projectable_columns": len(singleton_projectable),
        "removed_singleton_equality_projectable_columns": len(
            singleton_projectable & removed_columns
        ),
        "singleton_equality_ranged_columns": len(singleton_ranged),
        "removed_singleton_equality_ranged_columns": len(
            singleton_ranged & removed_columns
        ),
        "changed_rhs_rows": changed_rhs,
        "changed_objective_columns": changed_objective,
        "changed_bound_columns": changed_bounds,
        "changed_common_rows": changed_rows,
        "changed_coefficients": changed_coefficients,
        "added_coefficients": added_coefficients,
        "deleted_common_coefficients": deleted_common_coefficients,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("original_dir", type=Path)
    parser.add_argument("presolved_dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = []
    for presolved_file in sorted(args.presolved_dir.glob("p_*.mps")):
        name = presolved_file.stem[2:]
        original_file = args.original_dir / f"{name}.mps"
        if not original_file.exists():
            continue
        records.append(
            fingerprint(parse_mps(original_file), parse_mps(presolved_file), name)
        )
    if not records:
        raise SystemExit("no original/presolved pairs found")
    fieldnames = list(records[0])
    if args.output:
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
    writer = csv.DictWriter(__import__("sys").stdout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(records)


if __name__ == "__main__":
    main()
