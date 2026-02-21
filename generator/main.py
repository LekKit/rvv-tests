from __future__ import annotations

"""Main driver for the RVV test generator.

Invokes every gen module's ``generate()`` function, collects the list of
produced ``.S`` files, and writes a manifest (``tests/MANIFEST``) so the
build system and test runner know exactly which tests exist.
"""

import shutil
import sys
import time
from pathlib import Path

from .gen import (
    config,
    load_store,
    indexed_load_store,
    segment_load_store,
    int_arith,
    int_cmp,
    int_widening,
    int_macc,
    int_extension,
    int_adc,
    fixed_point,
    float_arith,
    float_widening,
    reduction,
    mask,
    permutation,
    edge_cases,
)

# Ordered so that foundational tests (config, load/store) come first.
_GENERATORS: list[tuple[str, object]] = [
    ("config",          config),
    ("load_store",      load_store),
    ("indexed_ls",      indexed_load_store),
    ("segment_ls",      segment_load_store),
    ("int_arith",       int_arith),
    ("int_cmp",         int_cmp),
    ("int_widening",    int_widening),
    ("int_macc",        int_macc),
    ("int_extension",   int_extension),
    ("int_adc",         int_adc),
    ("fixed_point",     fixed_point),
    ("float_arith",     float_arith),
    ("float_widening",  float_widening),
    ("reduction",       reduction),
    ("mask",            mask),
    ("permutation",     permutation),
    ("edge_cases",      edge_cases),
]


def _clean_tests_dir(tests_dir: Path, *, verbose: bool = False) -> None:
    """Remove the ``tests/`` directory entirely before regenerating.

    Safety: only removes if the directory is actually named ``tests``.
    Everything inside is generated and will be recreated immediately after.
    """
    if not tests_dir.exists():
        return
    if tests_dir.name != "tests":
        raise ValueError(f"Refusing to remove {tests_dir}: not named 'tests'")

    shutil.rmtree(tests_dir)
    if verbose:
        print("  cleaned tests/ directory", file=sys.stderr)


def run(base_dir: Path, *, verbose: bool = False) -> list[str]:
    """Run all generators and return the list of produced files.

    Parameters
    ----------
    base_dir:
        Project root (the directory that contains ``tests/``).
    verbose:
        Print progress to stderr.

    Returns
    -------
    list[str]
        Sorted list of paths (relative to *base_dir*) that were written.
    """
    # Remove tests/ completely before regenerating (safe: only known file types).
    _clean_tests_dir(base_dir / "tests", verbose=verbose)

    all_files: list[str] = []
    total_t0 = time.monotonic()

    for name, mod in _GENERATORS:
        t0 = time.monotonic()
        try:
            files = mod.generate(base_dir)  # type: ignore[attr-defined]
        except Exception as exc:
            print(f"ERROR generating {name}: {exc}", file=sys.stderr)
            raise
        elapsed = time.monotonic() - t0
        if verbose:
            print(
                f"  {name:20s}  {len(files):4d} files  ({elapsed:.2f}s)",
                file=sys.stderr,
            )
        all_files.extend(files)

    total_elapsed = time.monotonic() - total_t0

    # Make paths relative to base_dir for the manifest.
    rel_files: list[str] = []
    for f in all_files:
        try:
            rel = str(Path(f).relative_to(base_dir))
        except ValueError:
            rel = f
        rel_files.append(rel)
    rel_files.sort()

    # Write manifest.
    manifest = base_dir / "tests" / "MANIFEST"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text("\n".join(rel_files) + "\n")

    if verbose:
        print(
            f"\nTotal: {len(rel_files)} test files generated in {total_elapsed:.2f}s",
            file=sys.stderr,
        )
        print(f"Manifest written to {manifest}", file=sys.stderr)

    return rel_files
