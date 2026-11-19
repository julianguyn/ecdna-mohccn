# preprocessing of metadata

###########################################################
# Sample Mapping file
###########################################################

map <- read.table("data/rawdata/ID_mapping_20260925/MOHCCN_PM2C_ID_map.tsv", header = TRUE)
missing_map <- read.csv("data/rawdata/ID_mapping_20260925/missing_PM2C_SAMPLE_ID_matched_11182026.csv")


map <- data.frame(
    PM2C_SAMPLE_ID = c(map$PM2C_SAMPLE_ID, missing_map$submitter_specimen_id),
    OICR_sID = c(map$OICR_sID, missing_map$oicr_sid),
    OICR_SAMPLE_FILE = c(map$OICR_SAMPLE_FILE, missing_map$julia_id)
)
map$match_id <- sub("^(([^-]+-){2}[^-]+)-.*", "\\1", map$PM2C_SAMPLE_ID)

write.csv(map, file = "data/procdata/clinical/ID_mapping.csv", quote = F, row.names = F)
