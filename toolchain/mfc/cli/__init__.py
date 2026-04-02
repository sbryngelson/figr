"""
CLI Schema and Auto-Generation Module.

This package provides the single source of truth for MFC's CLI definition,
with generators for argparse, shell completions, and documentation.

Usage:
    from mfc.cli.commands import MFC_CLI_SCHEMA
    from mfc.cli.argparse_gen import generate_parser
    # completion_gen removed for mini-app
"""

from .schema import (
    ArgAction,
    Argument,
    CLISchema,
    Command,
    CommonArgumentSet,
    Completion,
    CompletionType,
    Example,
    MutuallyExclusiveGroup,
    Positional,
)

__all__ = [
    "CLISchema",
    "Command",
    "Argument",
    "Positional",
    "Example",
    "CommonArgumentSet",
    "MutuallyExclusiveGroup",
    "ArgAction",
    "CompletionType",
    "Completion",
]
