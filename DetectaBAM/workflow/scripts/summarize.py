import pandas as pd

inputs = snakemake.input
thr    = snakemake.params.thr
outpath = snakemake.output[0]

rows = []
for f in inputs:
    line = Path(f).stem.replace("_depth", "")
    d = pd.read_csv(f, sep="\t", header=None, names=["chrom","pos","depth"])
    n = len(d); row = {"cell_line": line, "n_truth_positions": n}
    for t in thr:
        c = int((d.depth >= t).sum())
        row[f"callable_ge{t}"] = c
        row[f"rate_ge{t}"] = round(c/n, 3) if n else 0.0
    rows.append(row)

pd.DataFrame(rows).to_csv(outpath, sep="\t", index=False)