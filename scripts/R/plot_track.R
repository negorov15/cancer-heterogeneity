# This script plots coverage tracks for a given gene using Gviz.
# See this blog post for an example: https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/

library(data.table)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
source("config/config.R")

# Specify the range to plot
names <- c("TUFM")
chromosomes <- c("chr16")
positions <- list(
    c(28840411, 28848348)
)
mutations_pos <- list(
    28844030
)

name <- names[1]
thechr <- chromosomes[1]
st <- positions[[1]][1]
en <- positions[[1]][2]
mut_pos <- mutations_pos[[1]]

covtrack <- DataTrack(
    range      = BAM_DIR,
    genome     = GENOME,
    type       = "a",            # or "a" / "histogram" for the coverage profile
    chromosome = thechr,
    name       = "Coverage (CPM)",
    fontcolor.title="black", 
    col.axis="black"
    )

itrack <- IdeogramTrack(
    genome = "hg38", chromosome = thechr
)

grtrack <- GeneRegionTrack(
    txdb,
    chromosome = thechr, start = st, end = en,
    showId = TRUE,
    collapseTranscripts = "longest",
    transcriptAnnotation = "symbol",
    background.title = "transparent"
)

# knownGene stores Entrez IDs in the 'gene' column; map them to symbols
entrez <- gene(grtrack)
uniq   <- unique(entrez[!is.na(entrez)])
lut    <- mapIds(
    org.Hs.eg.db,
    keys      = uniq,
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
)
syms <- lut[entrez]                       # vectorized lookup, preserves length
syms[is.na(syms)] <- entrez[is.na(syms)]  # fall back to the ID where no symbol
symbol(grtrack) <- unname(syms)

ht <- HighlightTrack(
    trackList  = list(itrack, covtrack, grtrack), #itrack, gtrack,
    start      = mut_pos - 1,   # widen slightly so it's visible (see note)
    end        = mut_pos + 1,
    chromosome = thechr,
    col        = "red",          # band border
    fill       = "#FF000022"     # translucent red fill
)

out_file <- file.path(RESULTS_DIR, paste0(name, "_coverage", ".png"))

png(
    filename = out_file,
    width = 1000, height = 500, res = 200      # res=200 keeps text crisp
)

plotTracks(
    ht, 
    from = st, to = en,
    sizes = c(0.2, 0.8, 0.2),
    fontsize = 8        # base size everything scales from
    # sizes = c(0.2, 1, 0.8, 0.2)
)

dev.off()