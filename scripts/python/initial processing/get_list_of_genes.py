"""Find genes detected in every count matrix."""

import sys
import os
import csv
import pandas as pd
from glob import glob

def load_count_matrix(mfile):
    """Load a count matrix from a file, automatically detecting the separator."""
    # Peek at the first line to detect separator:
    with open(mfile) as f:
        first_line = f.readline()
    
    # Try common separators
    for sep in ["\t", ",", ";", " ", "|"]:
        if sep in first_line:
            return pd.read_csv(mfile, sep=sep, index_col=0)

    # Default to tab
    return pd.read_csv(mfile, sep="\t", index_col=0)

def extract_genes(mfile):
    """Extract genes detected from a count matrix."""
    # The index of the DataFrame contains gene names:
    cell_line = os.path.basename(mfile).replace(".count_mtx.tsv.txt", "")
    df = load_count_matrix(mfile)
    # Count cells with nonzero expression per gene:
    has_expression = df > 0
    cells_per_gene = has_expression.sum(axis=1)
    detected = set(cells_per_gene[cells_per_gene >= MIN_CELLS].index)
    if not detected:
        print(f"Warning: No genes detected in {cell_line} with the current threshold.")
        print("Check empty_cell_lines.txt for details.")
        with open(f"output/gene list/empty_cell_lines.txt", "w") as f:
            f.write(cell_line + "\n")
        print("-" * 40)
    else:
        print(f"{len(detected):,} genes detected in ≥{MIN_CELLS} cells")
        print("Writing detected genes to file...")
        with open(f"output/gene list/{cell_line}_genes.txt", "w") as f:
            for g in sorted(detected):
                f.write(g + "\n")
        print("-" * 40)
    return detected

def extract_most_expressed_genes(file, top_n=5):
    cell_line = os.path.basename(file).replace(".count_mtx.tsv.txt", "")
    df = load_count_matrix(file)
    total_expression = df.sum(axis=1)
    top_genes = total_expression.nlargest(top_n)
    print(f"Top {top_n} most expressed genes in {cell_line}:")
    os.makedirs("output/gene list/top genes", exist_ok=True)
    with open(f"output/gene list/top genes/{cell_line}_top{top_n}_genes.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene", "Total_Expression"])
        for gene, expression in top_genes.items():
            print(f"  {gene}: {expression:,.0f}")
            writer.writerow([gene, expression])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python get_list_of_genes.py /path/to/count_mtx.tsv.txt")
        sys.exit(1)

    count_matrix = sys.argv[1]
    if not os.path.isfile(count_matrix):
        print(f"Error: file not found: {count_matrix}")
        sys.exit(1)

    df = load_count_matrix(count_matrix)
    total_expression = df.sum(axis=1)
    sorted_genes = total_expression.sort_values(ascending=False)

    cell_line = os.path.basename(count_matrix).replace(".count_mtx.tsv.txt", "")
    os.makedirs("output/cell line genes", exist_ok=True)
    output_file = f"output/cell line genes/{cell_line}_genes_sorted_by_expression.txt"

    with open(output_file, "w") as f:
        for gene, expression in sorted_genes.items():
            f.write(f"{gene}\t{expression}\n")

    print(f"Wrote {len(sorted_genes):,} genes to {output_file}")