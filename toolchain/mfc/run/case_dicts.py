"""
Case parameter dictionaries — uses Fortran namelist parsing as source of truth.
"""

import re
from ..params.namelist_parser import get_target_params
from ..state import ARG

ALL_PARAMS = get_target_params()

# Union of all target params — used to check if a param is known at all
_ALL_KNOWN = set().union(*ALL_PARAMS.values())

IGNORE = []
CASE_OPTIMIZATION = []  # GPU case-opt params; not used in IGR-only build

_BASE_RE = re.compile(r"^([a-zA-Z_][a-zA-Z0-9_]*)")


class _AllParamMapping:
    """Checks if a param name is known by any target."""
    def __contains__(self, name):
        m = _BASE_RE.match(name)
        base = m.group(1) if m else name
        return base in _ALL_KNOWN

    def keys(self):
        return _ALL_KNOWN


class _TargetKeySet:
    """Set-like object that checks if a param is valid for a target."""
    def __init__(self, target_name):
        self._valid = ALL_PARAMS.get(target_name, set())

    def __contains__(self, name):
        m = _BASE_RE.match(name)
        base = m.group(1) if m else name
        return base in self._valid


ALL = _AllParamMapping()


def get_validator():
    return None  # Schema validation removed


def get_input_dict_keys(target_name):
    return _TargetKeySet(target_name)
