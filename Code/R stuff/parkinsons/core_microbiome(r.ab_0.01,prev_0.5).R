#!/usr/bin/env RScript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("phyloseqobject_parkinsons.RData")
# unrarefied that has been filtered - nonrarefied because the analyses below take into consideration that the different datasets come with different sequencing depths

#### "core" microbiome ####

# Convert otu table to relative abundance
# this is done because the analyses work with the relative, not absolute, abundance
parkinsons_RA <- transform_sample_counts(phylobj, fun=function(x) x/sum(x))

# Filter dataset by agricultural exposure
# create vector of ASV ids with only exposed samples
# create vector with non-exposed samples
phylobj_farm <- subset_samples(parkinsons_RA, `Lived_on_farm`=="Yes")
phylobj_nofarm <- subset_samples(parkinsons_RA, `Lived_on_farm`=="No")

# What ASVs are found in more than 70% of samples in each antibiotic usage category?
# trying changing the prevalence to see what happens
# detection = include anything that's more than 1 percent abundant
# prevalence = 0.7 is HIGH - include only if present in over 50% samples
farm_ASVs <- core_members(phylobj_farm, detection=0.01, prevalence = 0.5)
nofarm_ASVs <- core_members(phylobj_nofarm, detection=0.01, prevalence = 0.5)

# only 1 ASV for the farm subset, 3 for the no_farm subset

# What are these ASVs? you can code it in two different ways to see the same things
farm_ASVs_table <- prune_taxa(farm_ASVs,phylobj) %>%
  tax_table()
nofarm_ASVs_table <- tax_table(prune_taxa(nofarm_ASVs,phylobj))
view(farm_ASVs_table)
view(nofarm_ASVs_table)

# can plot those ASVs' relative abundance
farm_ASVs_barplot <- prune_taxa(farm_ASVs,parkinsons_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Lived_on_farm`, scales ="free") +
  labs(x="Sample, Agricultural exposure", y="Relative Abundance")

ggsave("core_ASVs_barplot(r.ab_0.01,prev_0.5).png", farm_ASVs_barplot, width = 20, height = 8)

# Notice that in this dataset, there are very few CORE microbiome members. This is common

### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?

# combine the vectors above into a list for ease of coding later on
farm_list_full <- list(Farm = farm_ASVs, No_Farm = nofarm_ASVs)

# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
core_venn <- ggVennDiagram(x = farm_list_full) + coord_flip()

ggsave("venn_core(r.ab_0.01,prev_0.5).png", core_venn)

### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
farm_list <- core_members(phylobj_farm, detection=0.001, prevalence = 0.10)
nofarm_list <- core_members(phylobj_nofarm, detection=0.001, prevalence = 0.10)

# combine the vectors above into a list for ease of coding later on
farm_list_full <- list(Farm = farm_list, No_Farm = nofarm_list)

# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
all_venn <- ggVennDiagram(x = farm_list_full) + coord_flip()

ggsave("venn_antibiotic.png", all_venn)
