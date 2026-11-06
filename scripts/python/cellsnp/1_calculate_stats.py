# This script computes statistics for cellSNP output, 
# including counts and percentages of positions with sufficient depth (DP) and allele depth (AD) for each cell line. 
# It reads the ground truth VCF files, DP and AD matrices, and generates summary statistics and top variant lists.

import os
import scipy.io
import pandas as pd

# Create output directory if it doesn't exist
OUTPUT_DIR = "output/cellsnp"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Read cell line information
INPUT = "input/cellsnp"
CELL_LINES = pd.read_csv("input/cell_lines_modelids.tsv", sep="\t")
LINES = CELL_LINES["Cell line"]
MODEL_IDS = CELL_LINES["Model_ID"]

# Loop through each cell line and compute statistics
for line, model_id in zip(LINES, MODEL_IDS):
    # Read ground truth VCF for the cell line
    # This VCF contains the known somatic mutations for the cell line
    GROUND_TRUTH = pd.read_csv(
                            f"output/SomaticMutationsFiltered/{line}/{line}_SomaticMutations_expressed.vcf", 
                            sep="\t", 
                            comment="#", 
                            header=None,
                            names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
                        )
    # Get the total number of known mutations for the cell line
    OUT_LINE, N_TRUTH = f"{INPUT}/{line}", len(GROUND_TRUTH)
    # Read DP and AD matrices for the cell line
    # row col value triplet per non-zero entry, 1-indexed. 
    # For example, 3 120 5 means SNP 3 has 5 ALT-supporting reads in cell 120.
    DP = scipy.io.mmread(f"{OUT_LINE}/cellSNP.tag.DP.mtx").tocsr()
    AD = scipy.io.mmread(f"{OUT_LINE}/cellSNP.tag.AD.mtx").tocsr()
    # Compute the number of positions with DP and AD above thresholds
    dp_pos = DP.sum(axis=1).A1                                  
    ad_pos = AD.sum(axis=1).A1                                  
    # Create statistics dataframe
    stats_data = []
    # Loop through DP thresholds and compute counts and percentages
    for k in (1, 5, 10, 20):
        n_dp = int((dp_pos >= k).sum())
        n_adk = int((ad_pos >= k).sum())
        n_ad1 = int((ad_pos[dp_pos >= k] >= 1).sum())

        pct_dp        = 100 * n_dp  / N_TRUTH
        pct_ad1       = 100 * n_ad1 / N_TRUTH
        pct_adk       = 100 * n_adk / N_TRUTH
        pct_ad1_cond  = 100 * n_ad1 / n_dp if n_dp else float("nan")
        # print(f"DP>={k:<2}: {n_dp:3d}/{N_TRUTH}  ({pct:.1f}%)")
        stats_data.append(
            {
                "DP_threshold": k,
                "count DP": n_dp,
                "count AD>=k": n_adk,
                "count AD>=1": n_ad1,
                "total": N_TRUTH,
                "percentage DP": pct_dp,
                "percentage AD>=k": pct_adk,
                "percentage AD>=1": pct_ad1,
                "percentage AD>=1|DP>=k": pct_ad1_cond
             }
            )
    # Create a DataFrame from the statistics data
    stats_df = pd.DataFrame(stats_data)
    # Save statistics to CSV
    os.makedirs(f"output/cellsnp/stats/{line}", exist_ok=True)
    stats_df.to_csv(f"output/cellsnp/stats/{line}/{line}_stats.csv", index=False)
    # Read the base VCF and add DP and AD columns
    base = pd.read_csv(f"{OUT_LINE}/cellSNP.base.vcf.gz", sep="\t", comment="#", header=None,
                    usecols=[0,1,3,4], names=["CHROM","POS","REF","ALT"])
    base["DP"] = dp_pos
    top_variants = base.sort_values("DP", ascending=False).head(20)
    # Save output to CSV
    base.to_csv(f"output/cellsnp/stats/{line}/{line}_variants.csv", index=False)
    top_variants.to_csv(f"output/cellsnp/stats/{line}/{line}_top_variants.csv", index=False)