import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def parse_info(vcf):
    vcf = vcf.copy()
    for key in ("AD", "DP"):
        vcf[key] = vcf["info"].str.extract(rf"{key}=(\d+)", expand=False).astype(int)
    return vcf.drop(columns="info")


def plot_ecdf(df_all, title, save_path):
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.ecdfplot(data=df_all, x="dist", hue="group", log_scale=True, ax=ax, legend=False)
    ax.set_xlabel("Distance to 3′ end (transcript bp)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    # Load input data
    mtx_dir = "input/cellsnp"
    dist_dir = "output/distances for all exons"
    cell_lines_tsv = pd.read_csv("input/cell_lines_modelids.tsv", sep="\t")
    # Plot dir
    plot_dir = "output/distance plots for all exons (no legends)"
    os.makedirs(plot_dir, exist_ok=True)

    cols_to_use = ["chr", "pos", "dist"]
    # for each cell line
    for line in cell_lines_tsv["Cell line"].tolist():
        # Load distance matrixes and cellsnp callability stats
        expr_dist = pd.read_csv(f"{dist_dir}/{line}/{line}_SomaticMutations_distance.csv", sep="\t", header=0, names=cols_to_use)
        unexpr_dist = pd.read_csv(f"{dist_dir}/{line}/{line}_SomaticMutations_distance_unexpressed.csv", sep="\t", header=0, names=cols_to_use)
        vcf = pd.read_csv(f"{mtx_dir}/{line}/cellSNP.base.vcf.gz", 
                                    compression='gzip', 
                                    sep="\t", 
                                    comment="#", 
                                    header=None, 
                                    names=["chr", "pos", "id", "ref", "alt", "qual", "filter", "info"],
                                    usecols=["chr", "pos", "info"]
                            )
        # Reformat the last column of cellsnp matrix
        vcf = parse_info(vcf)

        # Divide distances of expressed variants by thresholds k = [1,5] and k >= 16
        low = vcf[vcf["DP"].between(1, 5)].merge(expr_dist, on=["chr", "pos"])
        high = vcf[vcf["DP"] >= 16].merge(expr_dist, on=["chr", "pos"])
        df_all = pd.concat(
        [
            unexpr_dist.assign(group="not expressed"),
            high.assign(group="expressed, k>=16"),
            low.assign(group="expressed, 1<=k<=5"),
        ],
        ignore_index=True,
        )
        # Plot ECDF 
        title = f"ECDF for {line}"
        plot_ecdf(df_all,title, f"{plot_dir}/{line}_ecdf.png")