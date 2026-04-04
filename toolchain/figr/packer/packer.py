"""Pack simulation output (.dat files) into a comparable in-memory structure."""

import math
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class PackEntry:
    """One data file's numeric values."""

    def __init__(self, filepath: str, doubles: List[float]):
        self.filepath = filepath
        self.doubles = doubles

    def __repr__(self) -> str:
        return f"{self.filepath} {' '.join([str(d) for d in self.doubles])}"


class Pack:
    """Collection of PackEntry objects keyed by filepath."""

    def __init__(self):
        self.entries: Dict[str, PackEntry] = {}

    def has_bad_values(self) -> bool:
        for entry in self.entries.values():
            for v in entry.doubles:
                if math.isnan(v) or math.isinf(v):
                    return True
        return False

    def save(self, filepath: str) -> None:
        with open(filepath, "w") as f:
            for entry in sorted(self.entries.values(), key=lambda x: x.filepath):
                f.write(f"{entry}\n")

    def find(self, filepath: str) -> Optional[PackEntry]:
        return self.entries.get(filepath)

    def set(self, entry: PackEntry) -> None:
        self.entries[entry.filepath] = entry

    def remove(self, entry: PackEntry) -> None:
        self.entries.pop(entry.filepath, None)


def _extract_doubles(s: str) -> List[float]:
    """Parse a string of whitespace-separated numbers into a list of floats."""
    return [float(e) for e in re.sub(r"[\n\t\s]+", " ", s).strip().split(" ")]


def pack(casepath: str) -> Tuple[Pack, Optional[str]]:
    """Scan casepath/D/ for text .dat files and extract values.

    Each .dat file has lines of the form: x [y [z]] value
    We discard coordinate columns and keep only the value column.
    """
    case_dir = os.path.dirname(casepath) if os.path.isfile(casepath) else casepath
    d_dir = os.path.join(case_dir, "D")

    if not os.path.isdir(d_dir):
        return Pack(), f"Directory {d_dir} does not exist."

    result = Pack()

    for filepath in sorted(Path(d_dir).rglob("*.dat")):
        short_filepath = str(filepath).replace(f"{case_dir}", "")[1:].replace("\\", "/")

        try:
            with open(filepath) as f:
                content = f.read()
        except Exception as e:
            return Pack(), f"Failed to read {short_filepath}: {e}"

        try:
            # Each line is: x [y [z]] value
            # Determine dimensionality from the first line
            ndims = len(_extract_doubles(content.split("\n", 1)[0])) - 1
            # Extract only the value column (every ndims+1'th element, starting at ndims)
            doubles = _extract_doubles(content)[ndims::ndims + 1]
        except (ValueError, IndexError) as e:
            return Pack(), f"Failed to parse {short_filepath}: {e}"

        result.set(PackEntry(short_filepath, doubles))

    if not result.entries:
        return result, f"No .dat files found in {d_dir}"

    return result, None


def load(filepath: str) -> Pack:
    """Load a golden/pack .txt file into a Pack."""
    result = Pack()
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(" ")
            path = parts[0]
            doubles = [float(v) for v in parts[1:]]
            result.set(PackEntry(path, doubles))
    return result


# Alias for backward compatibility
compile = pack
