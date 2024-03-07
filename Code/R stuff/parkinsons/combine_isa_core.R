#!/usr/bin/env Rscript
library(tidyverse)
library(dplyr)
library(tibble)

## Combining the ASVs from the core microbiome and indicator species analysis
load("core_farm_ASV.RData")
load("isa_frame.RData")
view(isa_frame)
view(farm_ASV_core)


# change the taxonomy table from phyloseq to a data frame
farm_ASV_core_df <- as.data.frame(farm_ASV_core)
view(farm_ASV_core_df)

# name the first column in the taxonomy table "ASV" to match the name in the isa_frame
farm_ASV_core_df <- rownames_to_column(farm_ASV_core_df, var = "row_names")
names(farm_ASV_core_df)[1] <- "ASV"

#combine the 
isa_core_farm <- inner_join(isa_frame, farm_ASV_core_df, by = "ASV")
view(isa_core_farm)

write.csv(isa_core_farm, "isa_core_farm.csv", row.names = FALSE)

# 3 species were identified
