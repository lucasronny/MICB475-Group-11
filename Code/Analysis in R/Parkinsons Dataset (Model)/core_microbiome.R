#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("phyloseqobject_parkinsons.RData")
# unrarefied that has been filtered - nonrarefied because the analyses below take into consideration that the different datasets come with different sequencing depths

load("phyloseqobject_parkinsons.RData")
metadata <- sample_data(phylobj)
view(metadata)

# 04.03.2024 DATASET HAS NOT BEEN FILTERED! MUST RE-DO CORE, ISA, AND ALL DOWNSTREAM STEPS

healthy <- subset_samples(phylobj, Disease == "Control")
view(sample_data(healthy))

#### "core" microbiome ####

# Convert otu table to relative abundance
# this is done because the analyses work with the relative, not absolute, abundance
parkinsons_RA <- transform_sample_counts(healthy, fun=function(x) x/sum(x))


# Filter dataset by exposure
# create vector of ASV ids with only exposed samples
# create vector unexposed samples
parkinsons_farm <- subset_samples(parkinsons_RA, `Lived_on_farm`=="Yes")
parkinsons_nofarm <- subset_samples(parkinsons_RA, `Lived_on_farm`=="No")

# What ASVs are found in more than 70% of samples in each exposure category?
# trying changing the prevalence to see what happens
# detection = 0 means that no matter how abundant it is within the group, include it
# prevalence = 0.7 is HIGH - include only if present in over 70% samples

farm_ASVs <- core_members(parkinsons_farm, detection=0, prevalence = 0.1)
nofarm_ASVs <- core_members(parkinsons_nofarm, detection=0, prevalence = 0.1)


# What are these ASVs? you can code it in two different ways to see the same things
farm_ASV_core <- prune_taxa(farm_ASVs,phylobj) %>%
  tax_table()
view(farm_ASV_core)
save(farm_ASV_core, file = "core_farm_ASV.RData")

tax_table(prune_taxa(nofarm_ASVs,phylobj))

# can plot those ASVs' relative abundance
barplot <- prune_taxa(farm_ASVs,phylobj) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Lived_on_farm`, scales ="free") +
  labs (x = "Samples (grouped by agricultural exposure)", y = "Relative Abundance") + 
  theme(legend.key.width = unit(1, "cm")) +
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 1, ncol = 3))
barplot

ggsave("core_ASVs_barplot(r.ab_0,prev_0.1).png", barplot, width = 20, height = 8)

# combine the vectors above into a list for ease of coding later on
exposure_list_core <- list(Farm = farm_ASVs, No_Farm = nofarm_ASVs)

# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
core_venn <- ggVennDiagram(x = exposure_list_core) + coord_flip() +
  labs(title = "Total Count of Shared and Unique ASVs Per Group",  # Modify the title
  x = "Group",                                     # Modify the x-axis label
  y = "Count",                                     # Modify the y-axis label
  fill = "Count") +                                # Modify the legend label
  theme(plot.title = element_text(hjust = 0)) +  # Center the title
  theme(plot.title = element_text(size = 30),      # Increase title font size
        text = element_text(size = 25))            # Increase label font size

ggsave("venn_core(r.ab_0,prev_0.1)_FIXED.png", core_venn, height = 7, width = 10)



#### Subset all species that are unique to just the farm dataset ####

#conver the taxonomy table to a dataframe and transfer ASVs from the indexing column into a legitimate column and name it ASV
farm_ASV_core_df <- as.data.frame(farm_ASV_core)
view(farm_ASV_core_df)
farm_ASV_core_df <- rownames_to_column(farm_ASV_core_df, var = "row_names")
names(farm_ASV_core_df)[1] <- "ASV"

nofarm_ASVs_remove <- c(nofarm_ASVs)

for (r in 1:nrow(farm_ASV_core_df)) {
  farm_ASV_unique_df <- farm_ASV_core_df %>%
    filter(!ASV %in% nofarm_ASVs_remove)
}
view(farm_ASV_unique_df)

save(farm_ASV_unique_df, file = "farm_ASV_unique_df_FIXED.RData")


######### IGNORE ###############
# Create a Venn diagram of all species
farm_list <- core_members(parkinsons_farm, detection=0.001, prevalence = 0.10)
nofarm_list <- core_members(parkinsons_nofarm, detection=0.001, prevalence = 0.10)

# combine the vectors above into a list for ease of coding later on
exposure_list_full <- list(Farm1 = farm_list, No_Farm1 = nofarm_list)

# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
all_venn <- ggVennDiagram(x = exposure_list_full) + coord_flip()

ggsave("venn_antibiotic.png", first_venn)
