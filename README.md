# How Much of the Somatic Mutation Landscape Is Callable from 3′ scRNA-seq?

Analysis code for a feasibility study: **what fraction of known somatic driver mutations can actually be detected in 3′-end single-cell RNA-seq (scRNA-seq) data?**

Calling variants from scRNA-seq is a popular route to study genetic heterogeneity in cancer. But 3′ protocols only sequence a short stretch at the end of each transcript, so most of a gene body is never read. Before building any variant-based heterogeneity analysis on this data, you need to know how much of the mutational signal is observable at all. This repository answers that question and checks one proposed explanation: that callability depends on how far a mutation sits from the 3′ end of its transcript.


---

## Contents

- [Study design](#study-design)
- [Repository layout](#repository-layout)
- [Data](#data)
- [Setup](#setup)
- [Workflow](#workflow)
- [Key definitions and metrics](#key-definitions-and-metrics)
- [Results](#results)
- [Tests](#tests)
- [Known issues and caveats](#known-issues-and-caveats)
- [Data sources and references](#data-sources-and-references)

---

## Study design

1. **Truth set.** For each of 8 cancer cell lines, take the known somatic SNVs from DepMap.
2. **Observation.** Count alleles at those positions in the scRNA-seq BAMs with `cellsnp-lite`, using UMI-collapsed depth.
3. **Callability.** Report what fraction of truth-set SNVs reach a given depth (DP ≥ k) and show alt-allele support (AD ≥ 1).
4. **Explanation.** Compute each variant's distance to the 3′ end of its transcript and compare callable vs. non-callable variants.

| Cell line | Cancer type | DepMap ModelID |
|---|---|---|
| A-253 | Head and neck | ACH-000740 |
| FaDu | Head and neck | ACH-000846 |
| SCC-4 | Head and neck | ACH-000238 |
| MCF-7 | Breast | ACH-000019 |
| MDA-MB-231 | Breast | ACH-000768 |
| BT-474 | Breast | ACH-000927 |
| HCT-116 | Colorectal | ACH-000971 |
| HCT-15 | Colorectal | ACH-000997 |

The cell line / ModelID mapping lives in `input/cell_lines_modelids.tsv` and drives every per-line loop.

---

## Repository layout

```
.
├── config/
│   └── config.R                      # Shared constants for the R coverage plots
├── scripts/
│   ├── python/
│   │   ├── initial processing/       # Barcodes, gene lists, sample sheets
│   │   ├── cellsnp/                  # Truth-set filtering, callability stats and plots
│   │   └── distance processing/      # 3′-end distance calculation, ECDF plots, unit tests
│   └── R/
│       └── plot_track.R              # Gviz coverage track with the mutation highlighted
├── input/                            # Not tracked
└── output/                           # Not tracked
```

### Scripts

| Script | Purpose |
|---|---|
| `initial processing/extract_barcodes_from_cellid.py` | Extracts the 20-nt cell barcode from each count-matrix column name (`<barcode>.<flowcell><lane>`), writes one barcode list per cell line, and reports barcodes shared between cell lines. |
| `initial processing/get_list_of_genes.py` | Ranks genes by total UMI count for one count matrix and writes `<line>_genes_sorted_by_expression.txt`. |
| `initial processing/get_cell_lines.py` | Lists the flowcell IDs found for each cell line. |
| `initial processing/create_csv.py` | Builds a cell line sample sheet (`resources/cell_lines.tsv`) with barcode and gene-list paths for the upstream pipeline. |
| `cellsnp/filter_variants.py` | Filters DepMap to default-entry SNVs per line and splits them into expressed / not expressed / not in the expression table. |
| `cellsnp/1_calculate_stats.py` | Loads `cellsnp-lite` DP/AD matrices and computes callability at DP thresholds 1, 5, 10, 20. |
| `cellsnp/2_plot.py` | Per-line bar plots, the pooled summary table, and the decay-curve overlay. |
| `distance processing/calculate_distance.py` | Core function: distance from a variant to the transcript 3′ end, measured in transcript (spliced) coordinates. |
| `distance processing/main.py` | Runs the distance calculation for every variant in every cell line. |
| `distance processing/process_gtf.ipynb` | Parses the GENCODE GTF into an exon table (all exons and `basic`-tagged transcripts only). |
| `distance processing/plot_results.py` | ECDF of 3′ distance for not-expressed variants vs. expressed variants with low depth (1–5) and high depth (≥16). |
| `R/plot_track.R` | Gviz plot of a coverage bigWig over a gene, with the mutation highlighted in red. |

---

## Data

| Input | Expected location | Source |
|---|---|---|
| Cell line / ModelID table | `input/cell_lines_modelids.tsv` | Tab-separated, columns `Cell line`, `Model_ID` |
| DepMap somatic mutations | `input/OmicsSomaticMutations.csv` | [DepMap portal](https://depmap.org/portal/) |
| Per-cell-line expression ranking | `input/ExpressionData/<line>/<line>_genes_sorted_by_expression.txt` | Output of `get_list_of_genes.py` |
| `cellsnp-lite` output | `cellSNP.base.vcf.gz`, `cellSNP.tag.DP.mtx`, `cellSNP.tag.AD.mtx`, `cellSNP.tag.OTH.mtx`, `cellSNP.samples.tsv` | Upstream pipeline |
| Gene annotation | `input/gtf/genes.gtf` | GENCODE v32 (Ensembl 98), GRCh38 |
| Coverage tracks | `input/bw/<line>.bw` | `deepTools bamCoverage` |
| Count matrices | `<line>.count_mtx.tsv.txt` |  |

- **Reference genome:** hg38, with `chr`-prefixed contigs.

---

## Setup

Python ≥ 3.12 and R ≥ 4.x.

```bash
python -m venv .venv
source .venv/bin/activate
pip install pandas numpy scipy matplotlib seaborn pysam pytest
```

R packages for `plot_track.R`:

```r
install.packages("data.table")
BiocManager::install(c("Gviz", "BSgenome.Hsapiens.UCSC.hg38",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))
```

**Run every script from the repository root.** Paths such as `input/...` and `output/...` are relative.

---

## Workflow

### 1. Prepare per-line inputs

```bash
# Rank genes by total UMIs for a count matrix
python "scripts/python/initial processing/get_list_of_genes.py" path/to/cell_line.count_mtx.tsv.txt

# Barcode lists per cell line (edit input_dir in main() first)
python "scripts/python/initial processing/extract_barcodes_from_cellid.py"
```

### 2. Build the truth sets

```bash
python scripts/python/cellsnp/filter_variants.py
```

Filters DepMap to `IsDefaultEntryForModel == "Yes"` and `VariantType == "SNV"`, then splits each line's SNVs by whether the gene is detected in that line's expression table. The CSV exports are currently commented out in the script; uncomment the `to_csv` block to write `output/FilteredSomaticMutations/<line>/*.csv`.

### 3. Compute callability

```bash
python scripts/python/cellsnp/1_calculate_stats.py
python scripts/python/cellsnp/2_plot.py
```

- `1_calculate_stats.py` writes `output/cellsnp/stats/<line>/<line>_stats.csv`, `_variants.csv` and `_top_variants.csv`.
- `2_plot.py` writes per-line bar plots, `callability_overlay.png` and `callability_pooled.csv`.
- Step 3 reads truth-set VCFs from `output/SomaticMutationsFiltered/<line>/<line>_SomaticMutations_expressed.vcf`.

### 4. Distance to the 3′ end

1. Run `scripts/python/distance processing/process_gtf.ipynb` to produce annotation file that includes only exons (and the `basic`-transcript subset).
2. Compute distances:
   ```bash
   python "scripts/python/distance processing/main.py"
   ```
   `main.py` calls `main("unexpressed")` at the bottom. Edit it to `main("expressed")` and run again for the other group.
3. Plot ECDFs:
   ```bash
   python "scripts/python/distance processing/plot_results.py"
   ```

**How the distance is defined.** For each variant, all transcripts whose genomic span contains the position are looked up. If the position falls inside an exon, the distance is the number of exonic bases between the variant and the transcript's 3′ end (introns are excluded; strand is respected). The final value is the mean over all overlapping transcripts. Variants that fall in no exon get `None`.

### 5. Coverage track for a single locus

Edit the gene name, coordinates and mutation position at the top of `scripts/R/plot_track.R`, point `BAM_DIR` in `config/config.R` at a bigWig, then:

```bash
Rscript scripts/R/plot_track.R
```

---

## Key definitions and metrics

- **Truth set.** DepMap SNVs flagged `IsDefaultEntryForModel == "Yes"`. DepMap stores one row per variant per sequencing run and runs disagree, so filtering to the default entry gives one coherent experiment per line instead of a union across runs.
- **Two denominators.**
  - `all`: every SNV in the truth set.
  - `expressed`: SNVs in genes detected in that line's expression table. This is the denominator used in the plots from `2_plot.py` ("SNVs in expressed genes").
- **DP (depth).** UMI-collapsed molecule count at the position, summed over all cells. This is the operative callability metric. Raw read depth (bigWig/`bamCoverage`) is for visualisation only, because PCR duplicates inflate it at highly expressed loci.
- **AD.** UMI-collapsed alt-allele count, summed over cells.
- **Reported thresholds.** For k ∈ {1, 5, 10, 20}:
  - `DP ≥ k`: the position is covered.
  - `DP ≥ k & AD ≥ 1`: covered *and* the alternate allele is seen.
  - `AD ≥ k`: the alt allele alone reaches k.
  - `AD ≥ 1 | DP ≥ k`: of covered positions, the share where the alt allele is observed.

---

## Data sources and references

- Zhu et al. (2023), *Nature Communications*. scRNA-seq of pooled cancer cell lines on the DNBelab C4 platform (3′-end tag sequencing).
- DepMap: `OmicsSomaticMutations.csv`.
- GENCODE v32 / Ensembl 98 annotation, GRCh38.
- Tools: [`cellsnp-lite`](https://github.com/single-cell-genetics/cellsnp-lite), [PISA](https://github.com/shiquan/PISA), [deepTools](https://deeptools.readthedocs.io/), [Gviz](https://bioconductor.org/packages/Gviz/), pandas, SciPy, seaborn.
