#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
# again, non-rarefied
load("phyloseq.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus (grouping based on genus level - the species were not well defined)
# convert to RA by function
mpt_genus <- tax_glom(mpt_final, "Genus", NArm = FALSE)
mpt_genus_RA <- transform_sample_counts(mpt_genus, fun=function(x) x/sum(x))

#ISA
# multipatt function requires that the otu table is oriented a specific way, so t() is used here (sample names in rows, ASVs in columns)
# cluster --> define predictor
# multipatt --> calculate all IVs and create a comprehensive table
# p-values --> from random sampling, are the IVs actually characteristic of the different groups within the predictor or are they random
# if indicative, p < 0.05
isa_mpt <- multipatt(t(otu_table(mpt_genus_RA)), cluster = sample_data(mpt_genus_RA)$`reported.antibiotic.usage`)
summary(isa_mpt)
# create a nice table by combining with the tax table
# now, unlike in the phyloseq, we want the rows to be indexed as per normal and and the row names to be ASVs
taxtable <- tax_table(mpt_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")
view(taxtable)

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
isa_mpt$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()