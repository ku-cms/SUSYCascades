#!/usr/bin/env python3
"""
Convert integers to floats in 'edges' (pt, eta) and 'content' arrays
within all multibinning nodes of a correctionlib JSON file.
"""

import json
import sys
from pathlib import Path
from typing import Optional


def convert_multibinning_node(node: dict) -> None:
    """
    Recursively walk the node tree.
    When a multibinning node is found, convert all int values in
    'edges' (all sub-lists) and 'content' to floats, in-place.
    """
    if not isinstance(node, dict):
        return

    if node.get("nodetype") == "multibinning":
        # Convert edges: list of lists
        if "edges" in node:
            node["edges"] = [
                [float(v) for v in edge_list]
                for edge_list in node["edges"]
            ]

        # Convert content: flat list of values
        if "content" in node:
            node["content"] = [float(v) for v in node["content"]]

    # Recurse into all dict values and list items
    for value in node.values():
        if isinstance(value, dict):
            convert_multibinning_node(value)
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    convert_multibinning_node(item)


def process_file(input_path: str, output_path: Optional[str] = None) -> None:
    input_file = Path(input_path)
    if not input_file.exists():
        print(f"Error: file not found: {input_path}")
        sys.exit(1)

    with input_file.open("r") as f:
        data = json.load(f)

    corrections = data.get("corrections", [])
    if not corrections:
        print("No corrections found in file.")
        return

    for correction in corrections:
        if "data" in correction:
            convert_multibinning_node(correction["data"])

    # Write output
    out_file = Path(output_path) if output_path else input_file.with_suffix(".converted.json")
    with out_file.open("w") as f:
        json.dump(data, f, indent=2)

    print(f"Done. Written to: {out_file}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python convert_json_ints.py <input.json> [output.json]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    process_file(input_path, output_path)