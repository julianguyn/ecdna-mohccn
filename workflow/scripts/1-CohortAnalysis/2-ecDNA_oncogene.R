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

# ------------------------- removing BCCA for now!!!!
AA_results <- AA_results[AA_results$Centre == "PM2C",]
all_samples <- all_samples[all_samples$Centre == "PM2C",]

# load in cancer type
oncotree <- read.csv(paste0(PROCDATA_DIR, "clinical/oncotree.csv"))

# keep ecDNA
ecDNA <- AA_results[AA_results$Classification == "ecDNA",]
ecDNA$cancer <- oncotree$oncotree_root[match(ecDNA$'Sample name', oncotree$samples)]


###########################################################
# Get oncogenes
###########################################################

# get oncogenes
amplified_onco <- ecDNA[ecDNA$Oncogene != "[]",] # 2139 amplicons with oncogene
oncogenes <- unlist(
  strsplit(gsub("\\[|\\]|'", "", amplified_onco$Oncogenes), ",\\s*")
) |> unique() # 813 unique amplified oncogenes

# create new dataframe
oncogene_df <- data.frame(matrix(nrow=0, ncol=4))
for (oncogene in oncogenes) {
    samples <- amplified_onco[grep(oncogene, amplified_onco$Oncogene),c("cohort", "Sample name", "cancer")]
    samples$oncogene <- oncogene
    oncogene_df <- rbind(oncogene_df, samples)
}

###########################################################
# Plot oncogenes
###########################################################

# get top 25
tt <- table(oncogene_df$oncogene) |> as.data.frame()
tt <- tt[order(tt$Freq, decreasing = T),]
top_genes <- tt$Var1[1:25]

toPlot <- oncogene_df[oncogene_df$oncogene %in% top_genes,]
toPlot$oncogene <- factor(toPlot$oncogene, levels = rev(top_genes))
toPlot$cancer <- factor(toPlot$cancer, levels = rev(names(cancer_pal)))

p <- ggplot(toPlot, aes(x = oncogene, fill = cancer)) +
  geom_bar(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual("Cancer Type", values = cancer_pal, breaks = rev(levels(toPlot$cancer))) +
  theme_classic() +
  theme(
    legend.key.size = unit(4, "mm"),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0)
  ) +
  coord_flip() +
  labs(y = "Number of ecDNA occurrences", x = "Oncogene") +
  scale_y_continuous(position = "right", expand = c(0, 0))
filename <- "data/results/figures/1-CohortAnalysis/top_oncogenes_ecdna.png"
ggsave(filename, p, width = 4, height = 5)

