# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(scales)
    library(tidyverse)
    library(patchwork)
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

# ------------------------- removing BCCA for now!!!!
AA_results <- AA_results[AA_results$Centre == "PM2C",]
all_samples <- all_samples[all_samples$Centre == "PM2C",]

# load in ID mappping file
map <- read.table("data/rawdata/ID_mapping_20260925/MOHCCN_PM2C_ID_map.tsv", header = TRUE)
map$match_id <- sub("^(([^-]+-){2}[^-]+)-.*", "\\1", map$PM2C_SAMPLE_ID)

# load in cancer type
oncotree <- read.csv(paste0(PROCDATA_DIR, "clinical/oncotree.csv"))
all_samples$cancer <- oncotree$oncotree_root[match(all_samples$sample, oncotree$samples)]
all_samples$cancer_code <- oncotree$cancer_type_code_oncotree[match(all_samples$sample, oncotree$samples)]

# load in age
age <- read.csv("data/rawdata/clinical/20261111_year_frommatthew/csv/20261111_julia_ages.csv")

###########################################################
# Number of samples per cancer types
###########################################################

toPlot <- table(all_samples$cancer) |> as.data.frame()
toPlot <- toPlot[order(toPlot$Freq, decreasing = TRUE),]
toPlot$Var1 <- factor(toPlot$Var1, levels = rev(toPlot$Var1))

# keep order of cancer types
order_cancer <- toPlot$Var1

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

# -------- SEND TO MATTHEW
table(all_samples$sample %in% map$OICR_SAMPLE_FILE) #43 samples missing the PM2C_SAMPLE_ID
t1 <- all_samples$sample[-which(all_samples$sample %in% map$OICR_SAMPLE_FILE)]
table(age$submitter_donor_id %in% map$match_id) #138 samples missing OICR_SAMPLE_ID
t2 <- age$submitter_donor_id[-which(age$submitter_donor_id %in% map$match_id)]

write.csv(t1, file = "data/temp/missing_PM2C_SAMPLE_ID.csv", quote = F, row.names = F)
write.csv(t2, file = "data/temp/missing_OICR_SAMPLE_ID.csv", quote = F, row.names = F)


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

###########################################################
# Keep only ecDNA
###########################################################

ecDNA <- AA_results[AA_results$Classification == "ecDNA",]
ecDNA$cancer <- oncotree$oncotree_root[match(ecDNA$'Sample name', oncotree$samples)]

###########################################################
# Proportion with ecDNA per cancer type
###########################################################

# format toPlot
toPlot <- all_samples[,c("sample", "cancer", "cancer_code")]
toPlot$ecDNA <- ifelse(toPlot$sample %in% ecDNA$'Sample name', "ecDNA+", "ecDNA-")
toPlot$cancer <- factor(toPlot$cancer, levels = order_cancer)

toPlot <- toPlot %>%
  group_by(cancer_code, cancer) %>%
  summarise(
    n_samples = n(),
    n_ecDNA_pos = sum(ecDNA == "ecDNA+"),
    prop_ecDNA = n_ecDNA_pos / n_samples * 100,
    plot_label = paste0(first(cancer_code), " (n=", n_samples, ")"),
    .groups = "drop"
  ) %>%
  arrange(cancer, cancer_code) %>%
  mutate(
    plot_label = as.character(plot_label)
  ) %>% as.data.frame()

# helper function (split into two rows)
plot_prop <- function(p) {
    toPlot <- toPlot[toPlot$cancer %in% p,]
    p <- ggplot(toPlot, aes(x = plot_label, y = prop_ecDNA, fill = cancer)) +
    geom_hline(yintercept = c(25, 50, 75), color = "gray", alpha = 0.5) +
    geom_col() +
    facet_grid(~ cancer, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = palette_20[2:20]) +
    scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        expand = c(0, 0)
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 5),
        panel.spacing.x = unit(0.4, "lines"),
        axis.title.x = element_blank(),
        legend.position = "none"
    ) +
    labs(y = "ecDNA+ (%)")
    return(p)
}

p1 <- plot_prop(order_cancer[1:10])
p2 <- plot_prop(order_cancer[11:length(order_cancer)])

p <- p1 / p2

filename <- "data/results/figures/1-CohortAnalysis/ecDNA_cancercode_prop.png"
ggsave(filename, p, width = 8, height = 3.5)

###########################################################
# Count ecDNA per cancer type
###########################################################

toPlot <- table(ecDNA$'Sample name') |> as.data.frame()
toPlot$cancer <- oncotree$oncotree_root[match(toPlot$Var1, oncotree$samples)]
toPlot <- toPlot[order(toPlot$Freq, decreasing = TRUE),]
toPlot$Var1 <- factor(toPlot$Var1, levels = rev(toPlot$Var1))

p <- ggplot(toPlot, aes(x = cancer, y = Freq)) +
    geom_boxplot(color = "black", fill = olive) +
    geom_jitter(alpha = 0.5, width = 0.2) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Cancer Type", y = "Number of ecDNA Molecules")
filename <- "data/results/figures/1-CohortAnalysis/ecDNA_cancertype_count.png"
ggsave(filename, p, width = 5, height = 5)
