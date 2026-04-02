"""
figr Clean Command - Remove build artifacts.
"""

import os
import shutil

from .common import FIGR_BUILD_DIR
from .printer import cons


def clean():
    """Remove the build directory and all build artifacts."""
    if os.path.isdir(FIGR_BUILD_DIR):
        cons.print(f"Removing [bold magenta]{FIGR_BUILD_DIR}[/bold magenta]...")
        try:
            shutil.rmtree(FIGR_BUILD_DIR)
            cons.print("[bold green]Build directory cleaned successfully.[/bold green]")
        except OSError as e:
            cons.print(f"[bold red]Error cleaning build directory:[/bold red] {e}")
    else:
        cons.print("[yellow]Build directory does not exist, nothing to clean.[/yellow]")
