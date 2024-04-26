#!/usr/bin/env Rscript
library(tidyverse)
library(dplyr)
library(tibble)

#### Combining the ASVs from the core microbiome farm subset and indicator species analysis ####
load("farm_ASV_unique_df_FIXED.RData")
load("isa_frame_species_FIXED.RData")
view(isa_frame_all)
view(farm_ASV_unique_df)


# change the taxonomy table from phyloseq to a data frame
farm_ASV_unique_df <- as.data.frame(farm_ASV_unique_df)
view(farm_ASV_unique_df)

#combine the tables using the inner_join function
isa_core_farm <- inner_join(isa_frame_all, farm_ASV_unique_df, by = "ASV")
view(isa_core_farm)

write.csv(isa_core_farm, "isa_core_farm_FIXED.csv", row.names = FALSE)

# generated 4 species, all p < 0.05
