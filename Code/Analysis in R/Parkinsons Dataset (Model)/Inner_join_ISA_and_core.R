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


##### IGNORE ######
# name the first column in the taxonomy table "ASV" to match the name in the isa_frame
farm_ASV_core_df <- rownames_to_column(farm_ASV_core_df, var = "row_names")
names(farm_ASV_core_df)[1] <- "ASV"

########

#combine the tables using the inner_join function
isa_core_farm <- inner_join(isa_frame_all, farm_ASV_unique_df, by = "ASV")
view(isa_core_farm)
# 3 species were identified

write.csv(isa_core_farm, "isa_core_farm_FIXED.csv", row.names = FALSE)

# generated 4 species, all p < 0.05 and IV > 0.39

#### IGNORE (done above) ###########
#Combine core microbiome species unique to the farm group and indicator species

load("farm_ASV_unique_df.RData")
isa_core_farm_unique <- inner_join(isa_frame, farm_ASV_unique_df, by = "ASV")
view(isa_core_farm_unique)
# 2 species identified!