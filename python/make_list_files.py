#!/usr/bin/env python3

import os
import glob

base_dir = "samples/NANO"
lists_dir = os.path.join(base_dir, "Lists")

# Create Lists directory if it does not exist
os.makedirs(lists_dir, exist_ok=True)

# Loop over subdirectories in samples/NANO
for entry in sorted(os.listdir(base_dir)):

    subdir_path = os.path.join(base_dir, entry)

    # Skip non-directories and the Lists directory itself
    if not os.path.isdir(subdir_path):
        continue
    if entry == "Lists":
        continue

    # Find all .txt files in this subdirectory
    txt_files = sorted(glob.glob(os.path.join(subdir_path, "*.txt")))

    # Output list filename
    output_file = os.path.join(lists_dir, f"{entry}.list")

    # Write file list
    with open(output_file, "w") as f:
        for txt in txt_files:
            f.write(txt + "\n")

    print(f"Wrote {len(txt_files)} entries to {output_file}")
