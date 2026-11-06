#!/usr/bin/env python3
"""
Extract 20-character barcodes from all count matrices in a directory.

The script extracts barcodes following the pattern:
<20_char_barcode>.<rest>
Example: TGGACGTGACCCGTGCAACC.DP800007365BR_L01_501
         -> TGGACGTGACCCGTGCAACC
"""

import pandas as pd
import re
import sys
from pathlib import Path
import os
from collections import defaultdict


def extract_barcode(column_name):
    """
    Helper function. Extract 20-character barcode from column name.

    Handles multiple patterns:
    - <barcode>.<flowcell ID><Lane>  (e.g., TGGACGTGACCCGTGCAACC.DP800007365BR_L01_501)
    - <barcode>  (e.g., TGGACGTGACCCGTGCAACC)
    - <barcode>_[A-Z]  (e.g., TGGACGTGACCCGTGCAACC_A)
    """
    match = re.search(r"^([ACGT]{20})(?:\.|_|$)", column_name)
    if match:
        return match.group(1)
    return None


def extract_cell_line_name(filename):
    """
    Helper function. Extract cell line name from filename.

    Pattern: <CellLine>.count_mtx.tsv.txt
    Example: Hs578T.count_mtx.tsv.txt -> hs578t
    """
    # Remove extension and get the cell line name
    name = filename.replace(".count_mtx.tsv.txt", "").lower()
    return name


def detect_separator(input_file):
    """
    Detect the separator used in the file (tab, comma, space, etc.).

    Reads the first line and tries common separators.
    """
    try:
        with open(input_file, "r") as f:
            first_line = f.readline()

        # Try common separators
        for sep in ["\t", ",", ";", " ", "|"]:
            if sep in first_line:
                return sep

        # Default to tab
        return "\t"
    except Exception:
        return "\t"


def process_count_matrix(input_file, output_dir):
    """
    Process a single count matrix file and extract barcodes from column names.

    Column format: <barcode>.<flowcell ID><Lane>
    """
    try:
        # Generate output filename based on input filename
        cell_line_name = extract_cell_line_name(os.path.basename(input_file))

        print(f"\nProcessing file: {input_file}")
        print(f"Detected cell line: {cell_line_name}")

        # Read only the header (column names)
        df = pd.read_csv(input_file, sep=detect_separator(input_file), nrows=0)

        # Get column names, skip first column (usually gene names)
        columns = df.columns.tolist()

        if columns[0] == "ID":
            print(
                f"Warning: First column '{columns[0]}' does not look like a barcode. Skipping it."
            )
            columns = columns[1:]  # Skip the first column if it's 'ID'

        # Extract barcodes from column names
        barcodes = [extract_barcode(col) for col in columns]

        # Filter out failed extractions
        extracted_barcodes = [bc for bc in barcodes if bc is not None]

        # Generate output filename based on input filename
        cell_line_name = extract_cell_line_name(os.path.basename(input_file))

        output_file = os.path.join(output_dir, f"{cell_line_name}_barcodes.tsv")

        # Save barcodes to file, one per line
        with open(output_file, "w") as f:
            for barcode in extracted_barcodes:
                f.write(barcode + "\n")

        print(
            f"Extracted {len(extracted_barcodes)} barcodes from {len(columns)} total columns"
        )
        print(f"Results saved to: {output_file}")

        return True, len(extracted_barcodes)
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        return False, 0


def main():
    # Directory containing the count matrices
    input_dir = "/Users/nikoegorov/Documents/internship/data/CSE0000121_processed"
    output_dir = "output/barcodes"

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Find all count matrix files
    count_matrix_files = sorted(Path(input_dir).glob("*.count_mtx.tsv.txt"))

    if not count_matrix_files:
        print(f"No count matrix files found in {input_dir}")
        return

    print(f"Found {len(count_matrix_files)} count matrix files")
    print("=" * 60)

    total_barcodes = 0
    successful = 0

    # Process each file
    for input_file in count_matrix_files:
        success, barcode_count = process_count_matrix(str(input_file), output_dir)
        if success:
            successful += 1
            total_barcodes += barcode_count

    print("\n" + "=" * 60)
    print(
        f"Processing complete: {successful}/{len(count_matrix_files)} files processed successfully"
    )
    print(f"Total barcodes extracted: {total_barcodes}")

def find_non_unique_barcodes(input_path="output/barcodes/CL200136255_L01/"):
    # ---- Phase 1: load barcodes from each file ----
    barcode_files = sorted(Path(input_path).glob("*.tsv"))
    # We'll build a parallel list: barcodes_by_cell_line[i] holds the set of barcodes
    # found in barcode_files[i].
    barcodes_by_cell_line = []
    for cell_line in barcode_files:
        with open(cell_line, "r") as f:
            barcodes = f.read().splitlines()
            barcodes = set(barcodes)  # Remove duplicates
            barcodes_by_cell_line.append(barcodes)
    # ---- Phase 2: invert the mapping ----
    # Goal: go from {file -> {barcodes}} to {barcode -> {files it appears in}}.
    # The inverted form is what makes the collision check trivial: a collision is
    # just a key whose value-set has size >= 2.
    barcode_to_lines = defaultdict(set)
    for path, bcs in zip(barcode_files, barcodes_by_cell_line):
        # strips the directory and the final extension
        # PosixPath("output/.../hs578t_barcodes.tsv").stem -> "hs578t_barcodes"
        cell_line = path.stem 
        for bc in bcs:
            barcode_to_lines[bc].add(cell_line)
    # ---- Phase 3: filter to collisions and report ----
    # Dict comprehension: keep only barcodes seen in >=2 distinct cell lines.
    # Anything with len(lines) == 1 is unique to its cell line and uninteresting.
    collisions = {bc: lines for bc, lines in barcode_to_lines.items() if len(lines) >= 2}
    print(f"Found {len(collisions)} barcodes shared across cell lines:")
    # Sort collisions by severity
    for bc, lines in sorted(collisions.items(), key=lambda x: -len(x[1])):
        print(f"  {bc}: {len(lines)} cell lines -> {sorted(lines)}")

if __name__ == "__main__":
    find_non_unique_barcodes()
