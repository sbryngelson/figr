"""Pack simulation output (.dat files) into a comparable in-memory structure."""

import math
import os
import struct
from typing import Dict, List, Optional, Tuple


class PackEntry:
    """One data file's numeric values."""

    def __init__(self, filepath: str, doubles: List[float]):
        self.filepath = filepath
        self.doubles = doubles


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
            for path in sorted(self.entries):
                entry = self.entries[path]
                vals = " ".join(repr(v) for v in entry.doubles)
                f.write(f"{entry.filepath} {vals}\n")

    def find(self, filepath: str) -> Optional[PackEntry]:
        return self.entries.get(filepath)

    def set(self, entry: PackEntry) -> None:
        self.entries[entry.filepath] = entry

    def remove(self, entry: PackEntry) -> None:
        self.entries.pop(entry.filepath, None)


def _read_dat(filepath: str) -> List[float]:
    """Read a Fortran unformatted sequential binary file and return its values as doubles."""
    with open(filepath, "rb") as f:
        data = f.read()

    if len(data) == 0:
        return []

    # Try Fortran unformatted sequential: 4-byte record marker, data, 4-byte record marker
    if len(data) >= 8:
        (header,) = struct.unpack("i", data[:4])
        if 0 < header <= len(data) - 8:
            (footer,) = struct.unpack("i", data[4 + header : 4 + header + 4])
            if header == footer:
                payload = data[4 : 4 + header]
                # Detect precision from record length
                if header % 8 == 0:
                    count = header // 8
                    return list(struct.unpack(f"{count}d", payload))
                elif header % 4 == 0:
                    count = header // 4
                    return [float(v) for v in struct.unpack(f"{count}f", payload)]

    # Fallback: raw binary (MPI-IO parallel output, no record markers)
    if len(data) % 8 == 0:
        count = len(data) // 8
        return list(struct.unpack(f"{count}d", data))
    elif len(data) % 4 == 0:
        count = len(data) // 4
        return [float(v) for v in struct.unpack(f"{count}f", data)]

    return []


def pack(dirpath: str) -> Tuple[Pack, Optional[str]]:
    """Scan dirpath for .dat files in D/ subdirectory and load them."""
    d_dir = os.path.join(dirpath, "D")
    if not os.path.isdir(d_dir):
        return Pack(), f"Directory {d_dir} does not exist."

    result = Pack()
    for root, _, files in os.walk(d_dir):
        for fname in sorted(files):
            if not fname.endswith(".dat"):
                continue
            full = os.path.join(root, fname)
            rel = os.path.relpath(full, dirpath)
            try:
                doubles = _read_dat(full)
            except Exception as e:
                return Pack(), f"Failed to read {rel}: {e}"
            result.set(PackEntry(rel, doubles))

    if not result.entries:
        return result, f"No .dat files found in {d_dir}"

    return result, None


def load(filepath: str) -> Pack:
    """Load a golden.txt file into a Pack."""
    result = Pack()
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            path = parts[0]
            doubles = [float(v) for v in parts[1:]]
            result.set(PackEntry(path, doubles))
    return result


# Alias used by check_case_optimization_output.py
compile = pack
