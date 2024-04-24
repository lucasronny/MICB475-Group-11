#Load in packages
library(phyloseq)
library(tidyverse)
library(picante)
library(ggplot2)
library(gridExtra)
library(vegan)



#Loading phyloseq object for healthy individuals
load("phyloseqobject_parkinsons_rarefied.RData") #loaded as "healthy"

#subset the healthy individuals
healthy <- subset_samples(phylobj_raref, Disease == "Control")

###### Filtering phyloseq object #########
#Filtering for the healthy individuals who have lived on a farm
healthyfarm <- subset_samples(healthy, Lived_on_farm == "Yes")
healthyfarm #24 samples in healhtyfarm after rarefaction

#Filtering for the healthy individuals who have not lived on a farm
healthynonfarm <- subset_samples(healthy, Lived_on_farm == "No")
healthynonfarm #64 samples in healthynonfarm after rarefaction



####### Alpha diversity analysis #########
#Generating richness plots for the healthy individuals
gg_richness <- plot_richness(healthy, x = "Lived_on_farm" ,measures = c("Observed")) +
  xlab("Group") +
  ylab("Richness")+
  #labs(title= "Observed Richness")+
  theme(plot.title = element_text(size = 10),
        axis.text.x = element_text(angle= 0, vjust = 1, hjust=0.5))+
  geom_boxplot()+
  scale_x_discrete(labels = c("Yes" = "Farm", "No" = "No-farm"))


gg_richness1 <- plot_richness(healthy, x = "Lived_on_farm" ,measures = c("Shannon")) +
  xlab("") +
  ylab("Shannon's Index")+
  #labs(title= "Observed Richness")+
  theme(plot.title = element_text(size = 10),
        axis.text.x = element_text(angle= 0, vjust = 1, hjust=0.5))+
  geom_boxplot()+
  scale_x_discrete(labels = c("Yes" = "Farm", "No" = "No-farm"))

# Faith's Phylogenetic distance analysis
phylo_dist <- pd(t(otu_table(healthy)), phy_tree(healthy),
                 include.root=F) #unrooted tree

# add faith's PD to metadata table to make graphs with ggplot
sample_data(healthy)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(healthy), aes(Lived_on_farm, PD)) + 
  geom_point()+
  geom_boxplot()+
  scale_x_discrete(labels = c("Yes" = "Farm", "No" = "No-farm"))+
  xlab("") +
  ylab("Phylogenetic Diversity") +
  #labs(title= "Observed Richness")+
  theme(plot.title = element_text(size = 10),
        axis.text.x = element_text(angle= 0, vjust = 1, hjust=0.5))
  

# Combine alpha diversity plots
alpha_plots <- list(gg_richness1, gg_richness, plot.pd)
alpha_figure <- grid.arrange(grobs = alpha_plots, nrow = 1)

#Save alhpa diversity figure
ggsave("Alpha_Diversity_Plots.png", alpha_figure, height = 7, width = 7)




########## Statistical analysis ############


#Observed richness
# Estimate richness for healthy individuals who have lived on a farm
richness_farm <- estimate_richness(healthyfarm, measures = "Observed")$Observed

# Estimate richness for healthy individuals who have not lived on a farm
richness_nonfarm <- estimate_richness(healthynonfarm, measures = "Observed")$Observed

# Perform Wilcoxon/Mann-Whitney U test for Observed richness
wilcox_observed <- wilcox.test(richness_farm, richness_nonfarm)

# Print the test result
print(wilcox_observed)



#Shannon diversity
# Estimate shannon diversity for healthy individuals who have lived on a farm
shannon_farm <- estimate_richness(healthyfarm, measures = "Shannon")$Shannon

# Estimate shannon diversity for healthy individuals who have not lived on a farm
shannon_nonfarm <- estimate_richness(healthynonfarm, measures = "Shannon")$Shannon

# Perform Wilcoxon/Mann-Whitney U test for Shannon Diversity
wilcox_shannon <- wilcox.test(shannon_farm, shannon_nonfarm)

# Print the test result
print(wilcox_shannon)


#Faith's PD
#Estimate faith's PD for the farm group
faiths_farm <- subset_samples(healthy, Lived_on_farm == "Yes")
faiths_farm_pd <- sample_data(faiths_farm)$PD

#Estimate faith's PD for the non-farm group
faiths_nonfarm <- subset_samples(healthy, Lived_on_farm == "No")
faiths_nonfarm_pd <- sample_data(faiths_nonfarm)$PD

#Perform Wilcoxon/Mann-Whitney U test for Faith's PD
wilcox_pd <- wilcox.test(faiths_farm_pd, faiths_nonfarm_pd)

#Print the test results
print(wilcox_pd)





####### Beta diversity plots ##########
#Bray Curtis Beta diversity
#create a distance metric for bray curtis

colors <- c("Yes" = "blue", "No"="red")


braycurt <- distance(healthy, method="bray")

#generating a PCoA matrix using the bray curtis distances
pcoa_braycurt <- ordinate(healthy, method="PCoA", distance=braycurt)

#Generating a PCoA plot from the bray curtis distances
gg_bpcoa <- plot_ordination(healthy, pcoa_braycurt, color = "Lived_on_farm") +
  labs(col="")+
  ggtitle("Bray Curtis")+
  theme_bw()+
  scale_color_manual(values = colors, 
                     labels = c("Yes" = "Farm", "No" = "No-farm")) +
  labs(color = "Group")+
  theme(plot.title = element_text(hjust = 0 ),
        legend.key = element_rect(fill = "lightgray", color = "lightgray"))

  

#Unifrac Beta diversity
#create a distance metric for unifrac distances
unif <- distance(healthy, method="unifrac")

#generating a PCoA matrix using the unifrac distances
pcoa_unifrac <- ordinate(healthy, method="PCoA", distance=unif)

#Generating a PCoA plot from the unifrac distances
gg_upcoa <- plot_ordination(healthy, pcoa_unifrac, color = "Lived_on_farm") +
  labs(col="")+
  ggtitle("Unifrac")+
  scale_color_manual(values = colors, 
                     labels = c("Yes" = "Farm", "No" = "No-farm")) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() +
  labs(color = "Group")

#Weighted Unifrac Beta diversity
#create a distance metric for Weighted Unifrac distances
wunif <- distance(healthy, method="wunifrac")

#generating a PCoA matrix using the Weighted Unifrac distances
pcoa_wunifrac <- ordinate(healthy, method="PCoA", distance=wunif)

#Generating a PCoA plot from the Weighted Unifrac distances
gg_wupcoa <- plot_ordination(healthy, pcoa_wunifrac, color = "Lived_on_farm") +
  labs(col="")+
  ggtitle("Weighted Unifrac")+
  scale_color_manual(values = colors, 
                     labels = c("Yes" = "Farm", "No" = "No-farm")) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  labs(color = "Group")


#Jaccard Beta diversity
#create a distance metric for Jaccard distances
jacc <- distance(healthy, method="jaccard")

#generating a PCoA matrix using the Jaccard distances
pcoa_jacc <- ordinate(healthy, method="PCoA", distance=jacc)


#Generating a PCoA plot from the Jaccard distances
gg_japcoa <- plot_ordination(healthy, pcoa_jacc, color = "Lived_on_farm") +
  labs(col="")+
  ggtitle("Jaccard")+
  scale_color_manual(values = colors, 
                     labels = c("Yes" = "Farm", "No" = "No-farm")) +
  theme(plot.title.position = "plot") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  labs(color = "Group")



# Combine beta diversity plots
beta_plots <- list(gg_bpcoa, gg_upcoa, gg_wupcoa, gg_japcoa)
beta_figure <- grid.arrange(grobs = beta_plots, nrow = 2)


#Formatting of the plots into a multi-plot
# Remove legends from individual plots
gg_bpcoa1 <- gg_bpcoa + theme(legend.position = "none")
gg_upcoa1 <- gg_upcoa + theme(legend.position = "none")
gg_wupcoa1 <- gg_wupcoa + theme(legend.position = "none")
gg_japcoa1 <- gg_japcoa + theme(legend.position = "none")

# Combine beta diversity plots without legends
combined_plots <- grid.arrange(gg_bpcoa1, gg_upcoa1, gg_wupcoa1, gg_japcoa1, nrow = 2)

# Create a custom legend for the overall plot
custom_legend <- cowplot::get_legend(gg_bpcoa)
  

# Add a single legend to the combined plot
final_plot <- grid.arrange(combined_plots, custom_legend, ncol = 2, widths = c(4, 1))

# Save combined beta diversity figure with single legend
ggsave("Combined_Beta_Diversity_Plots_with_Legend.png", final_plot, height = 7, width = 9)


######### Statistical Analysis ###########

#Making dataframe with the metadata
dat <- data.frame(sample_data(healthy))

#Bray Curtis distance matrix
dm_braycurtis <- vegdist(t(otu_table(healthy)), method="bray") 
#running permanova for Bray Curtis
adonis2(dm_braycurtis ~ Lived_on_farm, data=dat)


#Unifrac distance matrix
dm_unifrac <- UniFrac(healthy, weighted=FALSE)
#Running permanova
adonis2(dm_unifrac ~ Lived_on_farm, data=dat)


#Weighted Unifrac distance matrix
dm_wunifrac <- UniFrac(healthy, weighted=TRUE)
#Running permanova
adonis2(dm_wunifrac ~ Lived_on_farm, data=dat)


#Jaccard distance matrix
dm_jaccard <- vegdist(t(otu_table(healthy)), method="jaccard") # Bray-curtis
#Running permanova
adonis2(dm_jaccard ~ Lived_on_farm, data=dat)


