import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
 
# ---- config ----------------------------------------------------------------
LINES_TSV = "input/cell_lines_modelids.tsv"
STATS_DIR = "output/cellsnp/stats"
OUT_DIR   = "output/cellsnp/plots_nolegends"
os.makedirs(OUT_DIR, exist_ok=True)
 
GRAY = "#898781"   # denominator / context
 
# input column -> short internal name
COLS = {
    "count DP":    "dp",       # DP >= k
    "count AD>=1": "dp_ad1",   # DP >= k & AD >= 1
    "count AD>=k": "ad_k",     # AD >= k (implies DP >= k)
    "total":       "total",
}
# bars per threshold: (internal name, legend label, color)
SERIES = [
    ("dp",     "DP $\\geq$ k",                "#2a78d6"),
    ("dp_ad1", "DP $\\geq$ k & AD $\\geq$ 1", "#e8702a"),
    ("ad_k",   "AD $\\geq$ k",                "#7a4fa3"),
]
COUNT_COLS = [c for c, _, _ in SERIES]
 
 
def sanitize(name):
    """Return a safe filename for a cell line name."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(name))
 
 
def add_pct(df):
    """
    Add percentage columns to a stats DataFrame.
    Input: df with columns in COUNT_COLS plus "total"
    """
    df = df.copy()
    for c in COUNT_COLS:
        df[f"{c}_pct"] = 100 * df[c] / df["total"]
    # share of DP-callable SNVs where the variant is actually seen
    df["ad1_given_dp_pct"] = 100 * df["dp_ad1"] / df["dp"].where(df["dp"] > 0)
    return df
 
 
# ---- load ------------------------------------------------------------------
cell_lines = pd.read_csv(LINES_TSV, sep="\t")
LINES = cell_lines["Cell line"].tolist()
REQUIRED = {"DP_threshold", *COLS}
 
frames, missing = [], []
for line in LINES:
    f = f"{STATS_DIR}/{line}/{line}_stats.csv"
    if not os.path.exists(f):
        missing.append(line)
        continue
    df = pd.read_csv(f)
    absent = REQUIRED - set(df.columns)
    if absent:
        raise SystemExit(f"[error] {f} is missing columns: {sorted(absent)}")
    df = df.rename(columns=COLS)[["DP_threshold", *COLS.values()]]
    df["cell_line"] = line
    frames.append(df)
 
if missing:
    print(f"[warn] no stats file for: {', '.join(missing)}")
if not frames:
    raise SystemExit(f"[error] no callability stats found under {STATS_DIR}/")


# Combine all per-line stats into one DataFrame, sorted by line and DP threshold
combined = pd.concat(frames, ignore_index=True).sort_values(["cell_line", "DP_threshold"])
present = combined["cell_line"].unique().tolist()
combined = add_pct(combined)

 
# ---- per-line grouped bar plot ---------------------------------------------
def plot_line(df, title, fname):
    df = df.sort_values("DP_threshold")
    total = int(df["total"].iloc[0])
    thr = df["DP_threshold"].tolist()
 
    fig, ax = plt.subplots(figsize=(8, 4.5))
 
    ax.bar([0], [total], color=GRAY, width=0.5, zorder=3)
    ax.text(0, total * 1.02, f"{total}", ha="center", va="bottom", fontsize=9)
 
    x = np.arange(1, len(thr) + 1)
    n = len(SERIES)
    w = 0.8 / n
    for i, (col, label, color) in enumerate(SERIES):
        counts = df[col].tolist()
        off = (i - (n - 1) / 2) * w
        bars = ax.bar(x + off, counts, width=w, color=color, zorder=3, label=label)
        for bar, c in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + total * 0.02,
                    f"{c}\n({100 * c / total:.0f}%)",
                    ha="center", va="bottom", fontsize=7)
 
    ax.set_xticks([0] + list(x))
    ax.set_xticklabels(["Expressed\ntruth SNVs"] + [f"k = {t}" for t in thr])
    ax.set_ylabel("Truth-set SNVs (count)")
    ax.set_title(title, fontsize=11)
    ax.set_ylim(0, total * 1.25)
    ax.grid(axis="y", alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    #ax.legend(fontsize=8, frameon=False, loc="upper right")
 
    fig.tight_layout()
    fig.savefig(fname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {fname}")
 
 
for line in present:
    plot_line(combined[combined["cell_line"] == line],
              title=f"{line} \u2014 somatic SNV callability (expressed genes)",
              fname=f"{OUT_DIR}/{sanitize(line)}_callability.png")
 
 
# ---- pooled across lines ---------------------------------------------------
pooled = combined.groupby("DP_threshold")[[*COUNT_COLS, "total"]].sum().reset_index()
pooled = add_pct(pooled)
 
# unweighted mean of per-line percentages: every line counts equally
pct_cols = [f"{c}_pct" for c in COUNT_COLS] + ["ad1_given_dp_pct"]
line_means = (combined.groupby("DP_threshold")[pct_cols].mean()
              .add_prefix("mean_line_").reset_index())
pooled = pooled.merge(line_means, on="DP_threshold")
 
pooled.to_csv(f"{OUT_DIR}/callability_pooled.csv", index=False)
 
 
# ---- decay curves across lines ---------------------------------------------
stats = combined.groupby("DP_threshold")[pct_cols].agg(["min","max"])
x = stats.index.to_numpy()

fig, ax = plt.subplots(figsize=(5.5, 3.5))          # ~ final printed size
for col, label, color in SERIES:
    # min-max range across cell lines
    ax.fill_between(x, stats[(f"{col}_pct", "min")], stats[(f"{col}_pct", "max")],
                    color=color, alpha=0.15, lw=0)
    # equal-weight mean across cell lines
    ax.plot(pooled["DP_threshold"], pooled[f"mean_line_{col}_pct"],
            marker="o", lw=1.8, color=color, label=label)

ax.set_xticks(x)
ax.set_xlabel("Depth threshold $k$ (UMIs)", fontsize=10)
ax.set_ylabel("% of SNVs in expressed genes", fontsize=10)
ax.set_ylim(0, stats.xs("max", axis=1, level=1).max().max() * 1.1)
ax.tick_params(labelsize=9)
ax.legend(fontsize=9, frameon=False)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/callability_overlay.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"saved: {OUT_DIR}/callability_overlay.png")
 
 
# ---- summary table ---------------------------------------------------------
print()
print(combined[["cell_line", "DP_threshold", "total",
                "dp", "dp_pct", "dp_ad1", "dp_ad1_pct", "ad_k", "ad_k_pct",
                "ad1_given_dp_pct"]].round(1).to_string(index=False))
print("\ndone.")