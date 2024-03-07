#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
# again, non-rarefied
load("phyloseqobject_parkinsons.RData")

#### Indicator Species/Taxa Analysis ####

### Glom by Genus ###
# glom to Genus (grouping based on genus level - the species were not well defined)
# convert to RA by function
parkinsons_genus <- tax_glom(phylobj, "Genus", NArm = FALSE)
parkinsons_genus_RA <- transform_sample_counts(parkinsons_genus, fun=function(x) x/sum(x))

#ISA
# multipatt function requires that the otu table is oriented a specific way, so t() is used here (sample names in rows, ASVs in columns)
# cluster --> define predictor
# multipatt --> calculate all IVs and create a comprehensive table
# p-values --> from random sampling, are the IVs actually characteristic of the different groups within the predictor or are they random
# if indicative, p < 0.05
isa_parkinsons <- multipatt(t(otu_table(parkinsons_genus_RA)), cluster = sample_data(parkinsons_genus_RA)$`Lived_on_farm`)
view(isa_parkinsons)
# create a nice table by combining with the tax table
# now, unlike in the phyloseq, we want the rows to be indexed as per normal and and the row names to be ASVs
isa_taxtable <- tax_table(phylobj) %>% as.data.frame() %>% rownames_to_column(var="ASV")
view(isa_taxtable)


# consider that your table is only going to be resolved up to the species level
isa_frame <- isa_parkinsons$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)

save(isa_frame_high, file = "isa_frame_high.RData")
  
write.csv(isa_frame_high, "isa_frame_high.csv", row.names = FALSE)




#### Glom by species ####
# too few indicators generated when gloming by genus, attempt to glom by genus (make criteria for analysis less stringent)
parkinsons_species <- tax_glom(phylobj, "Species", NArm = FALSE)
parkinsons_species_RA <- transform_sample_counts(parkinsons_species, fun=function(x) x/sum(x))

#ISA
isa_parkinsons <- multipatt(t(otu_table(parkinsons_species_RA)), cluster = sample_data(parkinsons_species_RA)$`Lived_on_farm`)
view(isa_parkinsons)

# create a nice table by combining with the tax table
isa_taxtable <- tax_table(phylobj) %>% as.data.frame() %>% rownames_to_column(var="ASV")
view(isa_taxtable)

# generate a table containing all ASVs with IV>0.2 and p<0.05
isa_frame_high <- isa_parkinsons$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(isa_taxtable) %>%
  filter(p.value<0.05) %>%
  filter(stat>0.2)

save(isa_frame_high, file = "isa_frame_high.RData")

write.csv(isa_frame_high, "isa_frame_high.csv", row.names = FALSE)


