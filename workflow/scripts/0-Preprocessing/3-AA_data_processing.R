# clean AA_output for easier load in after

PROCDATA_DIR <- "data/procdata/"
RESULTS_DIR <- "data/results/"

source("workflow/scripts/0-Preprocessing/utils.R")
source("workflow/scripts/utils/palettes.R")
source("workflow/scripts/utils/cohorts.R")

###########################################################
# Load in data
###########################################################

# load in AA results
AA_results <- readRDS(paste0(PROCDATA_DIR, "AA_output/AA_results.rds")) |> as.data.frame()

# load in all samples
all_samples <- readRDS(paste0(PROCDATA_DIR, "AA_output/sample_df.rds"))

# process cohort names
AA_results <- process_cohort_names(AA_results)
all_samples <- process_cohort_names(all_samples)

###########################################################
# Preprocessing
###########################################################

# remove duplicates
all_samples <- all_samples[-which(all_samples$duplicated == "duplicated_removed"),]

# remove missing amplicons
all_samples <- all_samples[-which(all_samples$result_table == "missing"),]

# remove NA amplicons
AA_results <- AA_results[!is.na(AA_results$'AA amplicon number'),]
all_samples <- all_samples[all_samples$sample %in% AA_results$'Sample name',]

###########################################################
# Save for easier access
###########################################################

saveRDS(AA_results, file = paste0(PROCDATA_DIR, "AA_processed/AA_results.rds"))
saveRDS(all_samples, file = paste0(PROCDATA_DIR, "AA_processed/sample_df.rds"))
