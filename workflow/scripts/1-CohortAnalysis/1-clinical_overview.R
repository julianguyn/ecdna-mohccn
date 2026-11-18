# only have PM2C cancer types at the moment

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
})

PROCDATA_DIR <- "data/procdata/"
RESULTS_DIR <- "data/results/"

source("workflow/scripts/utils/palettes.R")

###########################################################
# Load in data
###########################################################

# load in AA results
AA_results <- readRDS(paste0(PROCDATA_DIR, "AA_processed/AA_results.rds"))
all_samples <- readRDS(paste0(PROCDATA_DIR, "AA_processed/sample_df.rds"))

# load in ID mappping file
map <- read.table("data/rawdata/ID_mapping_20260925/MOHCCN_PM2C_ID_map.tsv", header = TRUE)
map$match_id <- sub("^(([^-]+-){2}[^-]+)-.*", "\\1", map$PM2C_SAMPLE_ID)

# load in cancer type
oncotree <- read.csv(paste0(PROCDATA_DIR, "clinical/oncotree.csv"))

# load in age
age <- read.csv("data/rawdata/clinical/20261111_year_frommatthew/csv/20261111_julia_ages.csv")

###########################################################
# Check data matches with cancer type
###########################################################

PM_results <- AA_results[AA_results$Centre == "PM2C",]

table(PM_results$'Sample name' %in% oncotree$samples)
table(oncotree$samples %in% PM_results$'Sample name')

###########################################################
# Number of samples per cancer types
###########################################################

toPlot <- table(oncotree$oncotree_root) |> as.data.frame()
toPlot <- toPlot[order(toPlot$Freq, decreasing = TRUE),]
toPlot$Var1 <- factor(toPlot$Var1, levels = rev(toPlot$Var1))

p <- ggplot(toPlot, aes(y = Var1, x = Freq)) +
    geom_bar(stat = "identity", color = "black", fill = olive) +
    geom_text(
            aes(label = Freq),
            hjust = -0.3,
            color = "black",
            size = 3
        ) +
    theme_classic() +
    labs(y = "Cancer Type", x = "Number of Participants")
filename <- "data/results/figures/1-CohortAnalysis/cancertype_count.png"
ggsave(filename, p, width = 7.5, height = 7)

###########################################################
# Bin distribution of age
###########################################################

# TODO: map samples first
table(age$submitter_donor_id %in% map$match_id)
# missing mapping: 138

toPlot <- data.frame(
    Bin = c("1-20", "21-40", "41-60", "61-80", "81-100"),
    Count = c(
        nrow(age[age$age_at_diagnosis_min <= 20,]),
        nrow(age[age$age_at_diagnosis_min > 20 & age$age_at_diagnosis_min <= 40,]),
        nrow(age[age$age_at_diagnosis_min > 40 & age$age_at_diagnosis_min <= 60,]),
        nrow(age[age$age_at_diagnosis_min > 60 & age$age_at_diagnosis_min <= 80,]),
        nrow(age[age$age_at_diagnosis_min > 80,])
    )
)

p <- ggplot(toPlot, aes(x = Bin, y = Count)) +
    geom_bar(stat = "identity", color = "black", fill = olive) +
    geom_text(
            aes(label = Count),
            vjust = -0.3,
            color = "black",
            size = 3
        ) +
    theme_classic() +
    labs(x = "Age Range", y = "Number of Participants")
filename <- "data/results/figures/1-CohortAnalysis/age_count.png"
ggsave(filename, p, width = 3, height = 4)