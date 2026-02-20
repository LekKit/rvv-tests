#!/usr/bin/env python3
from __future__ import annotations

"""Entry point for generating all RVV 1.0 test assembly files.

Usage:
    python generate_tests.py [-v | --verbose] [BASE_DIR]

BASE_DIR defaults to the directory containing this script.
"""

import argparse
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate RVV 1.0 test assembly files."
    )
    parser.add_argument(
        "base_dir",
        nargs="?",
        default=str(Path(__file__).resolve().parent),
        help="Project root directory (default: directory of this script)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print progress information",
    )
    args = parser.parse_args()
    base_dir = Path(args.base_dir).resolve()

    # Import the generator package.
    from generator.main import run

    files = run(base_dir, verbose=args.verbose)
    if not args.verbose:
        print(f"{len(files)} test files generated.")


if __name__ == "__main__":
    main()
