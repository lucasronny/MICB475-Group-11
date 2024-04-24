#Load in packages
library(phyloseq)
library(tidyverse)
library(ape)
library(picante)
library(ggpubr)
library(ggplot2)
library(gridExtra)



#Loading phyloseq object for healthy individuals
load("phyloseqobject_parkinsons_rarefied.RData") #loaded as "phylobj_raref"

#subset the healthy individuals
healthy <- subset_samples(phylobj_raref, Disease == "Control")

###### Filtering phyloseq object #########
#Filtering for the healthy individuals who have lived on a farm
healthyfarm <- subset_samples(healthy, Lived_on_farm == "Yes")
healthyfarm #24 samples in healhtyfarm after rarefaction

#Filtering for the healthy individuals who have not lived on a farm
healthynonfarm <- subset_samples(healthy, Lived_on_farm == "No")
healthynonfarm #64 samples in healthy nonfarm after rarefaction



####### Alpha diversity analysis #########
#Generating richness plots for the healthy individuals
gg_richness <- plot_richness(healthy, x = "Lived_on_farm" ,measures = c("Observed")) +
  xlab("Lived on a farm") +
  ylab("Richness")+
  #labs(title= "Observed Richness")+
  theme(plot.title = element_text(size = 10))+
  geom_boxplot()


gg_richness1 <- plot_richness(healthy, x = "Lived_on_farm" ,measures = c("Shannon")) +
  xlab("Lived on a farm") +
  ylab("Shannon's Index")+
  #labs(title= "Shannon's Diversity Index")+
  theme(plot.title = element_text(size = 10))+
  geom_boxplot()

# Faith's Phylogenetic distance analysis
phylo_dist <- pd(t(otu_table(healthy)), phy_tree(healthy),
                 include.root=F) #unrooted tree

# add faith's PD to metadata table to make graphs with ggplot
sample_data(healthy)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(healthy), aes(Lived_on_farm, PD)) + 
  geom_point()+
  geom_boxplot() +
  xlab("Lived on a farm") +
  ylab("Phylogenetic Diversity") +
  #labs(title= "Faith's Phylogenetic Diversity")+
  theme(plot.title = element_text(size = 10))
  

# Combine alpha diversity plots
alpha_plots <- list(gg_richness1, gg_richness, plot.pd)
alpha_figure <- grid.arrange(grobs = alpha_plots, nrow = 1)

#Save alhpa diversity figure
ggsave("Alpha_Diversity_Plots.png", alpha_figure, height = 7, width = 7)






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
  scale_color_manual(values = colors)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()
  

#Unifrac Beta diversity
#create a distance metric for unifrac distances
unif <- distance(healthy, method="unifrac")

#generating a PCoA matrix using the unifrac distances
pcoa_unifrac <- ordinate(healthy, method="PCoA", distance=unif)

#Generating a PCoA plot from the unifrac distances
gg_upcoa <- plot_ordination(healthy, pcoa_unifrac, color = "Lived_on_farm") +
  labs(col="")+
  ggtitle("Unifrac")+
  scale_color_manual(values = colors)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()

#Weighted Unifrac Beta diversity
#create a distance metric for Weighted Unifrac distances
wunif <- distance(healthy, method="wunifrac")

#generating a PCoA matrix using the Weighted Unifrac distances
pcoa_wunifrac <- ordinate(healthy, method="PCoA", distance=wunif)

#Generating a PCoA plot from the Weighted Unifrac distances
gg_wupcoa <- plot_ordination(healthy, pcoa_wunifrac, color = "Lived_on_farm") +
  labs(col="")+
  ggtitle("Weighted Unifrac")+
  scale_color_manual(values = colors)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()

#Jaccard Beta diversity
#create a distance metric for Jaccard distances
jacc <- distance(healthy, method="jaccard")

#generating a PCoA matrix using the Jaccard distances
pcoa_jacc <- ordinate(healthy, method="PCoA", distance=jacc)


#Generating a PCoA plot from the Jaccard distances
gg_japcoa <- plot_ordination(healthy, pcoa_jacc, color = "Lived_on_farm") +
  labs(col="")+
  ggtitle("Jaccard")+
  scale_color_manual(values = colors)+
  theme(plot.title.position = "plot") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()


# Combine beta diversity plots
beta_plots <- list(gg_bpcoa, gg_upcoa, gg_wupcoa, gg_japcoa)
beta_figure <- grid.arrange(grobs = beta_plots, nrow = 2)

# Save beta diversity figure
ggsave("Beta_Diversity_Plots.png", beta_figure, height = 7, width = 7)
