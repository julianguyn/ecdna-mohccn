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
# Removed samples
###########################################################
# see notes in oncotree for reason

# load in clinical data
oncotree <- read.csv("data/rawdata/clinical/mohda2026003_oncotree.csv")
to_rm <- oncotree$samples[oncotree$oncotree_root == ""]

AA_results <- AA_results[-which(AA_results$'Sample name' %in% to_rm),]
all_samples <- all_samples[-which(all_samples$sample %in% to_rm),]

# 8399 AA_results, 1212 samples

###########################################################
# Catch the one BTC without cancer type
###########################################################

PM_results <- AA_results[AA_results$Centre == "PM2C",]

# catch the one BTC with no matches
tobind <- data.frame(
    submitter_donor_id = "",
    samples = PM_results$'Sample name'[-which(PM_results$'Sample name' %in% oncotree$samples)],
    cohort = "BTC",
    cancer_type_code_oncotree = "CHOL",
    oncotree_root = "BILIARY_TRACT",
    notes = ""
)
oncotree <- rbind(oncotree, tobind)

###########################################################
# Save for easier access
###########################################################

saveRDS(AA_results, file = paste0(PROCDATA_DIR, "AA_processed/AA_results.rds"))
saveRDS(all_samples, file = paste0(PROCDATA_DIR, "AA_processed/sample_df.rds"))
write.csv(oncotree, file = paste0(PROCDATA_DIR, "clinical/oncotree.csv"))
