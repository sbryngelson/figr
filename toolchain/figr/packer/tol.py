"""Tolerance-based comparison of two Pack objects."""

import math
from typing import Optional, Tuple

from .pack import Pack


class Tolerance:
    """Absolute and relative tolerance for comparing packed values."""

    def __init__(self, abs_tol: float, rel_tol: float):
        self.abs_tol = abs_tol
        self.rel_tol = rel_tol


def compare(result: Pack, golden: Pack, tol: Tolerance) -> Tuple[int, Optional[str]]:
    """Compare *result* against *golden* entry-by-entry.

    Returns (0, None) on success, or (1, message) on the first mismatch.
    """
    for filepath, g_entry in golden.entries.items():
        r_entry = result.find(filepath)
        if r_entry is None:
            return 1, f"Missing file in result: {filepath}"

        if len(r_entry.doubles) != len(g_entry.doubles):
            return 1, (
                f"{filepath}: length mismatch "
                f"(result={len(r_entry.doubles)}, golden={len(g_entry.doubles)})"
            )

        for i, (r, g) in enumerate(zip(r_entry.doubles, g_entry.doubles)):
            if math.isnan(r) or math.isnan(g):
                if math.isnan(r) and math.isnan(g):
                    continue
                return 1, f"{filepath}[{i}]: NaN mismatch (result={r}, golden={g})"

            threshold = tol.abs_tol + tol.rel_tol * max(abs(r), abs(g))
            if abs(r - g) > threshold:
                return 1, (
                    f"{filepath}[{i}]: |{r} - {g}| = {abs(r - g):.6e} "
                    f"> tol {threshold:.6e}"
                )

    # Check for extra files in result not in golden
    for filepath in result.entries:
        if golden.find(filepath) is None:
            return 1, f"Extra file in result not in golden: {filepath}"

    return 0, None
