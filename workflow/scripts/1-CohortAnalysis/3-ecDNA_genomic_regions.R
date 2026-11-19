# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
})

PROCDATA_DIR <- "data/procdata/"
RESULTS_DIR <- "data/results/"

source("workflow/scripts/utils/palettes.R")
source("workflow/scripts/utils/genomic_regions.R")

###########################################################
# Load in data
###########################################################

# load in AA results
AA_results <- readRDS(paste0(PROCDATA_DIR, "AA_processed/AA_results.rds"))
all_samples <- readRDS(paste0(PROCDATA_DIR, "AA_processed/sample_df.rds"))

# ------------------------- removing BCCA for now!!!!
AA_results <- AA_results[AA_results$Centre == "PM2C",]
all_samples <- all_samples[all_samples$Centre == "PM2C",]

# keep ecDNA
ecDNA <- AA_results[AA_results$Classification == "ecDNA",]

###########################################################
# Get ecDNA segment positions
###########################################################

segments <- parse_ecDNA_segments(ecDNA)

###########################################################
# Evaluate bin enrichment
###########################################################

RESULTS_DATA_DIR <- paste0(RESULTS_DIR, "data/1-CohortAnalysis/3-ecDNA_genomic_regions/")
enrichment_file <- paste0(RESULTS_DATA_DIR, "enrichment.rds")
RERUN_PERM <- FALSE

# only reruns permutations if file does not exist, manual override to rerun
if (RERUN_PERM || !file.exists(enrichment_file)) {
    res <- bin_enrichment(segments, label = "ecDNA", n_perm = 1000)
    dir.create(RESULTS_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(res, file = enrichment_file)
} else {
    res <- readRDS(enrichment_file)
}

###########################################################
# Bin frequency and hotspot regions
###########################################################

n_samples <- length(unique(all_samples$sample))
freq <- bin_frequency(res, n_samples)
regions <- get_hotspot_regions(freq, segments, ecDNA, n_samples)
write.csv(regions, paste0(RESULTS_DATA_DIR, "hotspot_regions.csv"), row.names = FALSE)

###########################################################
# Plot bin frequency
###########################################################

p <- plot_frequency(freq, regions)
filename <- "data/results/figures/1-CohortAnalysis/ecDNA_genomic_frequency.png"
ggsave(filename, p, width = 10, height = 2)
