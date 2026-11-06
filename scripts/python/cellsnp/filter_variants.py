import os
import sys
import pandas as pd


def filter_variants(depmap_profile, query_df):
    """
    Filter DepMap mutation profile for a respective cell line
    Inputs:
        depmap_profile: DataFrame containing DepMap mutation profile for cell lines of interest
        query_df: DataFrame containing genes and their expression levels for a respective cell line
    Output: 
        cl_expressed_vars: DataFrame that contains only expressed variants
        cl_not_expressed_vars
    """
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    expressed_query = query_df[query_df["Expression"] > 0]
    not_expressed_query = query_df[query_df["Expression"] == 0]
    # Expressed
    variants_expressed = depmap_profile.merge(expressed_query, on="HugoSymbol")
    variants_expressed["Chrom"] = pd.Categorical(variants_expressed["Chrom"], categories=chrom_order, ordered=True)
    variants_expressed = variants_expressed.sort_values(["Chrom", "Pos"])
    # Not expressed
    variants_not_expressed = depmap_profile.merge(not_expressed_query, on="HugoSymbol")
    variants_not_expressed["Chrom"] = pd.Categorical(variants_not_expressed["Chrom"], categories=chrom_order, ordered=True)
    variants_not_expressed = variants_not_expressed.sort_values(["Chrom", "Pos"])
    # Not in depmap
    not_in_depmap = depmap_profile[~depmap_profile["HugoSymbol"].isin(query_df["HugoSymbol"])]
    not_in_depmap = not_in_depmap.sort_values("HugoSymbol").reset_index(drop=True)
    variants_not_expressed["Chrom"] = pd.Categorical(variants_not_expressed["Chrom"], categories=chrom_order, ordered=True)
    variants_not_expressed = variants_not_expressed.sort_values(["Chrom", "Pos"])

    return variants_expressed, variants_not_expressed, not_in_depmap

if __name__ == "__main__":
    # Load input files
    MUTATIONS_DIR = "input/OmicsSomaticMutations.csv"
    SomaticMutations = pd.read_csv(MUTATIONS_DIR, 
                                   sep=",",
                                   usecols=["ModelID", "IsDefaultEntryForModel", "Chrom", "Pos", "Ref", "Alt", "VariantType", "HugoSymbol", "EnsemblGeneID"]
                                   )
    # Clean the mutation profile, keep only default entries and SNVs
    SomaticMutations = SomaticMutations[(SomaticMutations["IsDefaultEntryForModel"] == "Yes") & 
                                        (SomaticMutations["VariantType"] == "SNV")
                                        ]
    SomaticMutations = SomaticMutations.drop(["IsDefaultEntryForModel", "VariantType"], axis=1)

    LINES_DIR = "input/cell_lines_modelids.tsv"
    LINES_TSV = pd.read_csv(LINES_DIR,sep="\t")
    cell_lines = LINES_TSV["Cell line"]
    model_ids = LINES_TSV["Model_ID"]

    # load each cell line expression profile
    EXPR_DIR  = "input/ExpressionData"
    for (cl, model_id) in zip(cell_lines, model_ids):
        expression_data = pd.read_csv(f"{EXPR_DIR}/{cl}/{cl}_genes_sorted_by_expression.txt",
                                      sep="\t", header=None, names=["HugoSymbol", "Expression"]
                                      )
        CL_SomaticMutations = SomaticMutations[SomaticMutations["ModelID"] == model_id]
        cl_expressed_vars, cl_not_expressed_vars, not_in_depmap = filter_variants(CL_SomaticMutations, expression_data)

        # save output 
        #outdir = f"output/FilteredSomaticMutations/{cl}"
        #os.makedirs(f"{outdir}", exist_ok=True)
        #CL_SomaticMutations.to_csv(f"{outdir}/all_variants.csv", sep="\t", index=False)
        #cl_expressed_vars.to_csv(f"{outdir}/expressed_variants.csv", sep="\t", index=False)
        #cl_not_expressed_vars.to_csv(f"{outdir}/not_expressed_variants.csv", sep="\t", index=False)
        #not_in_depmap.to_csv(f"{outdir}/not_in_depmap_variants.csv", sep="\t", index=False)

        
    