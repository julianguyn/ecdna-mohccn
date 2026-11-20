# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(scales)
    library(tidyverse)
    library(patchwork)
    library(tidytext)
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
all_samples$cohort <- droplevels(all_samples$cohort)

# load in ID mapping file
map <- read.csv("data/procdata/clinical/ID_mapping.csv")

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

p <- ggplot(toPlot, aes(y = Var1, x = Freq, fill = Var1)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    geom_text(
            aes(label = Freq),
            hjust = -0.3,
            color = "black",
            size = 3
        ) +
    scale_fill_manual(values = cancer_pal) +
    scale_x_continuous(expand = expansion(add = c(0, 20))) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(y = "Cancer Type", x = "Number of Participants")
filename <- "data/results/figures/1-CohortAnalysis/cancertype_count.png"
ggsave(filename, p, width = 7, height = 6.5)

###########################################################
# Number of samples per cohort
###########################################################

toPlot <- table(all_samples$cohort, all_samples$cancer) |> as.data.frame()
toPlot <- toPlot[order(toPlot$Freq, decreasing = TRUE),]

to_order <- table(all_samples$cohort) |> as.data.frame()
to_order <- to_order[order(to_order$Freq, decreasing = TRUE),]

toPlot$Var1 <- factor(toPlot$Var1, levels = rev(to_order$Var1))
toPlot$Var2 <- factor(toPlot$Var2, levels = rev(names(cancer_pal)))
toPlot$label <- to_order$Freq[match(toPlot$Var1, to_order$Var1)]

p <- ggplot(toPlot, aes(y = Var1, x = Freq, fill = Var2)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    geom_text(
        data = toPlot %>%
            group_by(Var1) %>%
            summarise(
                Freq = sum(Freq),
                label = first(label),
                .groups = "drop"
            ),
        aes(y = Var1, x = Freq, label = label),
        inherit.aes = FALSE, hjust = -0.3, color = "black", size = 3
    ) +
    scale_fill_manual("Cancer Type", values = cancer_pal, breaks = rev(levels(toPlot$Var2))) +
    scale_x_continuous(expand = expansion(add = c(0, 20))) +
    theme_classic() +
    theme(
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0)
    ) +
    labs(y = "Cohort", x = "Number of Participants")
filename <- "data/results/figures/1-CohortAnalysis/cohort_count.png"
ggsave(filename, p, width = 6, height = 6)

###########################################################
# Bin distribution of age
###########################################################

table(all_samples$sample %in% map$OICR_SAMPLE_FILE)

to_match <- map$match_id[map$OICR_SAMPLE_FILE %in% all_samples$sample]
age <- age[age$submitter_donor_id %in% to_match,]

missing_age <- nrow(all_samples) - nrow(age)
# missing age: 90

toPlot <- data.frame(
    Bin = c("1-20", "21-40", "41-60", "61-80", "81-100", "NA"),
    Count = c(
        nrow(age[age$age_at_diagnosis_min <= 20,]),
        nrow(age[age$age_at_diagnosis_min > 20 & age$age_at_diagnosis_min <= 40,]),
        nrow(age[age$age_at_diagnosis_min > 40 & age$age_at_diagnosis_min <= 60,]),
        nrow(age[age$age_at_diagnosis_min > 60 & age$age_at_diagnosis_min <= 80,]),
        nrow(age[age$age_at_diagnosis_min > 80,]),
        missing_age
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
ggsave(filename, p, width = 4, height = 4)

###########################################################
# Keep only ecDNA
###########################################################

ecDNA <- AA_results[AA_results$Classification == "ecDNA",]
ecDNA$cancer <- oncotree$oncotree_root[match(ecDNA$'Sample name', oncotree$samples)]

# format toPlot
toPlot <- all_samples[,c("sample", "cancer", "cancer_code")]
toPlot$ecDNA <- ifelse(toPlot$sample %in% ecDNA$'Sample name', "ecDNA+", "ecDNA-")
toPlot$cancer <- factor(toPlot$cancer, levels = names(cancer_pal))

###########################################################
# Proportion with ecDNA per cancer type
###########################################################

plot_df <- toPlot %>%
  group_by(cancer) %>%
  summarise(
    n = n(),
    n_ecDNA_pos = sum(ecDNA == "ecDNA+"),
    prop_ecDNA_pos = n_ecDNA_pos / n,
    .groups = "drop"
  )
p <- ggplot(plot_df, aes(x = cancer, y = prop_ecDNA_pos, fill = cancer)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(
      y = prop_ecDNA_pos / 2,
      label = ifelse(prop_ecDNA_pos == 0, "", paste0(round(prop_ecDNA_pos * 100, 0), "%")),
      color = ifelse(cancer %in% c("PANCREAS", "UTERUS"), "white", "black")
    ),
    size = 3
  ) +
  scale_color_identity() +
  scale_fill_manual(values = cancer_pal) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(y = "ecDNA+ (%)")
filename <- "data/results/figures/1-CohortAnalysis/ecDNA_cancertype_prop.png"
ggsave(filename, p, width = 8, height = 3.5)

###########################################################
# Proportion with ecDNA per cancer code
###########################################################

toPlot <- toPlot %>%
  group_by(cancer_code, cancer) %>%
  summarise(
    n_samples = n(),
    n_ecDNA_pos = sum(ecDNA == "ecDNA+"),
    prop_ecDNA = n_ecDNA_pos / n_samples * 100,
    plot_label = paste0(first(cancer_code), " (n=", n_samples, ")"),
    .groups = "drop"
  ) %>%
  group_by(cancer) %>%
  mutate(
    plot_label = reorder(plot_label, n_samples)
  ) %>%
  ungroup() %>%
  arrange(cancer, desc(n_samples)) %>%
  mutate(
    plot_label = as.character(plot_label)
  ) %>%
  as.data.frame()

toPlot <- toPlot %>%
  mutate(plot_label = reorder_within(plot_label, -n_samples, cancer))

# helper function (split into two rows)
plot_prop <- function(p) {
    toPlot <- toPlot[toPlot$cancer %in% p,]
    p <- ggplot(toPlot, aes(x = plot_label, y = prop_ecDNA, fill = cancer)) +
    geom_hline(yintercept = c(25, 50, 75), color = "gray", alpha = 0.5) +
    geom_col() +
    facet_grid(~ cancer, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cancer_pal) +
    scale_x_reordered() +
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

p1 <- plot_prop(names(cancer_pal)[1:10])
p2 <- plot_prop(names(cancer_pal)[11:length(names(cancer_pal))])

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
toPlot$cancer <- factor(toPlot$cancer, levels = names(cancer_pal))

p <- ggplot(toPlot, aes(x = cancer, y = Freq, fill = cancer)) +
    geom_jitter(alpha = 0.5, width = 0.2) +
    geom_boxplot(color = "black", outlier.shape = NA) +
    scale_fill_manual(values = cancer_pal) +
    theme_classic() + 
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none"
    ) +
    labs(y = "Number of ecDNA Molecules")
filename <- "data/results/figures/1-CohortAnalysis/ecDNA_cancertype_count.png"
ggsave(filename, p, width = 5, height = 5)

###########################################################
# Count ecDNA per cancer type
###########################################################

ecDNA$cancer <- factor(ecDNA$cancer, levels = names(cancer_pal))

p <- ggplot(ecDNA, aes(x = cancer, y = log2(.data[["Feature median copy number"]]), fill = cancer)) +
    geom_jitter(alpha = 0.5, width = 0.2) +
    geom_boxplot(color = "black", outlier.shape = NA) +
    scale_fill_manual(values = cancer_pal) +
    theme_classic() + 
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none"
    ) +
    labs(y = "log2 Median Copy Number")
filename <- "data/results/figures/1-CohortAnalysis/ecDNA_cancertype_cn.png"
ggsave(filename, p, width = 5, height = 5)
