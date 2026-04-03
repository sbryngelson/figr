"""
Parse Fortran namelist definitions to extract valid parameters for each target.

Reads the actual Fortran source files (src/*/m_start_up.fpp) so the Python
toolchain stays in sync with what Fortran accepts.
"""

import re
from pathlib import Path
from typing import Dict, Set


def parse_namelist_from_file(filepath: Path) -> Set[str]:
    content = filepath.read_text()

    namelist_match = re.search(
        r"namelist\s+/user_inputs/\s*(.+?)(?=\n\s*\n|\n\s*!(?!\s*&)|\n\s*[a-zA-Z_]+\s*=)",
        content, re.DOTALL | re.IGNORECASE
    )
    if not namelist_match:
        raise ValueError(f"Could not find namelist /user_inputs/ in {filepath}")

    namelist_text = namelist_match.group(1)
    namelist_text = re.sub(r"&\s*\n\s*", " ", namelist_text)
    namelist_text = re.sub(r"#:.*", "", namelist_text)
    namelist_text = re.sub(r"!.*", "", namelist_text)

    found_params = set()
    for match in re.finditer(r"\b([a-zA-Z_][a-zA-Z0-9_]*)\b", namelist_text):
        name = match.group(1)
        if name.lower() not in {"namelist", "user_inputs", "if", "endif", "not"}:
            found_params.add(name)
    return found_params


def get_figr_root() -> Path:
    return Path(__file__).resolve().parent.parent.parent.parent


_TARGET_PARAMS_CACHE: Dict[str, Set[str]] = {}


def get_target_params() -> Dict[str, Set[str]]:
    if not _TARGET_PARAMS_CACHE:
        root = get_figr_root()
        targets = {
            "pre_process": root / "src" / "pre_process" / "m_start_up.fpp",
            "simulation":  root / "src" / "simulation"  / "m_start_up.fpp",
            "post_process": root / "src" / "post_process" / "m_start_up.fpp",
        }
        for name, path in targets.items():
            _TARGET_PARAMS_CACHE[name] = parse_namelist_from_file(path)
    return _TARGET_PARAMS_CACHE
