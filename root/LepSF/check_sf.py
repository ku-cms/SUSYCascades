#!/usr/bin/env python3
"""Validate a correctionlib JSON file and report which corrections fail."""

import argparse
import json
import sys
from pathlib import Path

from correctionlib.schemav2 import Correction, CorrectionSet


def main():
    parser = argparse.ArgumentParser(
        description="Validate a correctionlib JSON file and locate broken corrections."
    )
    parser.add_argument(
        "json_file",
        type=Path,
        help="Path to the correctionlib JSON file to validate.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print full pydantic error messages for each failing correction.",
    )
    parser.add_argument(
        "--inspect",
        type=int,
        metavar="IDX",
        help="Pretty-print the raw JSON of correction at index IDX and exit.",
    )
    args = parser.parse_args()

    if not args.json_file.exists():
        print(f"ERROR: file not found: {args.json_file}", file=sys.stderr)
        sys.exit(1)

    with open(args.json_file) as f:
        cset = json.load(f)

    corrections = cset.get("corrections", [])
    print(f"Loaded {args.json_file} with {len(corrections)} corrections.\n")

    # Inspect mode: dump one correction and exit
    if args.inspect is not None:
        if not 0 <= args.inspect < len(corrections):
            print(f"ERROR: index {args.inspect} out of range "
                  f"(0..{len(corrections) - 1})", file=sys.stderr)
            sys.exit(1)
        print(json.dumps(corrections[args.inspect], indent=2))
        return

    # Try the whole file first
    try:
        CorrectionSet.parse_obj(cset)
        print("✓ Full CorrectionSet validates cleanly.")
        return
    except Exception as e:
        print("✗ Full CorrectionSet validation failed. "
              "Checking each correction individually...\n")

    # Per-correction validation
    failures = []
    for i, corr in enumerate(corrections):
        name = corr.get("name", "<unnamed>")
        try:
            Correction.parse_obj(corr)
        except Exception as e:
            failures.append((i, name, e))
            if args.verbose:
                print(f"--- correction {i} ({name}) ---")
                print(e)
                print()
            else:
                print(f"  [{i:3d}] {name}: {type(e).__name__}")

    print()
    if failures:
        print(f"Summary: {len(failures)} / {len(corrections)} corrections failed validation.")
        print(f"Failing indices: {[i for i, _, _ in failures]}")
        print("\nRe-run with --verbose for full error messages, "
              "or --inspect IDX to dump the raw JSON of one correction.")
        sys.exit(1)
    else:
        print("All individual corrections validated, but the full set did not. "
              "Check top-level fields (schema_version, description, corrections list shape).")
        sys.exit(1)


if __name__ == "__main__":
    main()