#!/usr/bin/env Rscript
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


# Filter dataset by exposure
# create vector of ASV ids with only exposed samples for both males and females
# create vector with unexposed samples for both males and females
parkinsons_farm <- subset_samples(parkinsons_RA, `Lived_on_farm`=="Yes")
parkinsons_male <- subset_samples(parkinsons_RA, `Sex`=="Male")
parkinsons_female <- subset_samples(parkinsons_RA, `Sex`=="Female")
parkinsons_nofarm <- subset_samples(parkinsons_RA, `Lived_on_farm`=="No")

# What ASVs are found in more than 70% of samples in each antibiotic usage category?
# trying changing the prevalence to see what happens
# detection = 0 means that no matter how abundant it is within the group, include it
# prevalence = 0.7 is HIGH - include only if present in over 70% samples

farm_ASVs <- core_members(parkinsons_farm, detection=0.01, prevalence = 0.5)
nofarm_ASVs <- core_members(parkinsons_nofarm, detection=0.01, prevalence = 0.5)
female_ASVs <- core_members(parkinsons_female, detection=0.01, prevalence = 0.5)
male_ASVs <- core_members(parkinsons_male, detection=0.01, prevalence = 0.5)


# What are these ASVs? you can code it in two different ways to see the same things
prune_taxa(farm_ASVs,phylobj) %>%
  tax_table()

tax_table(prune_taxa(nofarm_ASVs,phylobj))

prune_taxa(female_ASVs,phylobj) %>%
  tax_table()

tax_table(prune_taxa(male_ASVs,phylobj))

# can plot those ASVs' relative abundance
gender_barplot <- prune_taxa(nofarm_ASVs,phylobj) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Sex`, scales ="free") +
  labs (x = "Samples, Gender", y = "Relative Abundance")

ggsave("core_ASVs_barplot_gender(r.ab_0.01,prev_0.5).png", gender_barplot, width = 20, height = 8)

# Notice that in this dataset, there are very few CORE microbiome members. This is common


### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?

# combine the vectors above into a list for ease of coding later on
exposure_list_core <- list(Farm = farm_ASVs, No_Farm = nofarm_ASVs, Female = female_ASVs, Male = male_ASVs)

# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
core_venn <- ggVennDiagram(x = exposure_list_core)

ggsave("venn_core_gender(r.ab_0.01,prev_0.5).png", core_venn)


# Create a Venn diagram of all species
farm_list <- core_members(parkinsons_farm, detection=0.001, prevalence = 0.10)
nofarm_list <- core_members(parkinsons_nofarm, detection=0.001, prevalence = 0.10)

# combine the vectors above into a list for ease of coding later on
exposure_list_full <- list(Farm1 = farm_list, No_Farm1 = nofarm_list)

# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
all_venn <- ggVennDiagram(x = exposure_list_full) + coord_flip()

ggsave("venn_antibiotic.png", first_venn)
