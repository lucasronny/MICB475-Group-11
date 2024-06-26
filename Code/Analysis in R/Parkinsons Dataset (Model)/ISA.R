#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
# again, non-rarefied

healthy <- subset_samples(phylobj, Disease == "Control")
healthy_metadata <- sample_data(healthy)
view(healthy_metadata)

#Indicator Species/Taxa Analysis

##### Glom by Genus #####
# glom to Genus (grouping based on genus level - the species were not well defined)
# convert to RA by mathematical function
parkinsons_genus <- tax_glom(phylobj, "Genus", NArm = FALSE)
parkinsons_genus_RA <- transform_sample_counts(parkinsons_genus, fun=function(x) x/sum(x))

#ISA
# multipatt function requires that the otu table is oriented a specific way, use t() to put sample names in rows, ASV IDs in columns
isa_parkinsons <- multipatt(t(otu_table(parkinsons_genus_RA)), cluster = sample_data(parkinsons_genus_RA)$`Lived_on_farm`)
summary(isa_parkinsons)

# create a results table by combining with the tax table
isa_taxtable <- tax_table(phylobj) %>% as.data.frame() %>% rownames_to_column(var="ASV")
view(isa_taxtable)


# consider that your table is only going to be resolved up to the species level
isa_frame <- isa_parkinsons$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(isa_taxtable) %>%
  filter(p.value<0.05)

save(isa_frame, file = "isa_frame_all.RData")
save(isa_frame_high, file = "isa_frame_high.RData")
write.csv(isa_frame, "isa_frame.csv", row.names = FALSE)

# too few taxa are identified as indicators, change parameters (see below) 


##### Glom by species
# too few indicators generated when gloming by genus, attempt to glom by species (make criteria for analysis less stringent)
parkinsons_species <- tax_glom(healthy, "Species", NArm = FALSE)
parkinsons_species_RA <- transform_sample_counts(parkinsons_species, fun=function(x) x/sum(x))

#ISA
# multipatt function requires that the otu table is oriented a specific way, use t() to put sample names in rows, ASV IDs in columns
# cluster --> define predictor
# multipatt --> calculate all IVs and create a comprehensive table
isa_parkinsons <- multipatt(t(otu_table(parkinsons_species_RA)), cluster = sample_data(parkinsons_species_RA)$`Lived_on_farm`)

# create a nice table by combining with the tax table
isa_taxtable <- tax_table(phylobj) %>% as.data.frame() %>% rownames_to_column(var="ASV")


# generate a table containing all ASVs with IV>0.2 and p<0.05
isa_frame_high <- isa_parkinsons$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(isa_taxtable) %>%
  filter(p.value<0.05) %>%
  filter(stat>0.2)

save(isa_frame_high, file = "isa_frame_high_FIXED.RData")
write.csv(isa_frame_high, "isa_frame_high_FIXED.csv", row.names = FALSE)
#this file is for reference, not used in downstream analysis

isa_frame_all <- isa_parkinsons$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(isa_taxtable) %>%
  filter(p.value<0.05)


save(isa_frame_all, file = "isa_frame_species_FIXED.RData")
write.csv(isa_frame_all, "isa_frame_species_FIXED.csv", row.names = FALSE)
#these are the files used in downstream analysis

