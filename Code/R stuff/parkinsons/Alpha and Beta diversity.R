#Load in packages
library(phyloseq)
library(tidyverse)
library(ape)
library(picante)
library(ggpubr)


#Loading phyloseq object for healthy individuals
load("phyloseqobject_parkinsons_rarefied.RData") #loaded as "healthy"

###### Filtering phyloseq object #########
#Filtering for the healthy individuals who have lived on a farm
healthyfarm <- subset_samples(healthy, Lived_on_farm == "Yes")
healthyfarm #24 samples in healhtyfarm after rarefaction

#Filtering for the healthy individuals who have not lived on a farm
healthynonfarm <- subset_samples(healthy, Lived_on_farm == "No")
healthynonfarm #64 samples in healthynonfarm after rarefaction


####### Alpha diversity analysis #########

#Generating richness plots for the healthy individuals
gg_richness <- plot_richness(healthy, x = "Lived_on_farm", measures = c("Shannon","Observed")) +
  xlab("Lived_on_farm") +
  geom_boxplot() +
  geom_point() +
  geom_jitter()
gg_richness

#save plot
ggsave(filename = "plot_richness.png"
       , gg_richness
       , height=4, width=6)




# Faith's Phylogenetic distance analysis
phylo_dist <- pd(t(otu_table(healthy)), phy_tree(healthy),
                 include.root=F) #unrooted tree

# add faith's PD to metadata table to make graphs with ggplot
sample_data(healthy)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(healthy), aes(Lived_on_farm, PD)) + 
  geom_boxplot() +
  xlab("Lived on a farm") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd

ggsave(filename = "Faith's_PD_yes_and_no_farm.png"
       , plot.pd
       , height =5, width =5)

# Need to extract information:
# alphadiv: generate a table of all diversity measures in terms of richness
# samp_dat: extract the metadta out of the rarefied dataset
# samp_dat_wdiv: combine the above into 1 dataframe
alphadiv <- estimate_richness(healthy)
healthy_meta <- sample_data(healthy)
healthy_meta_wdiv <- data.frame(healthy_meta, alphadiv)
print(colnames(healthy_meta_wdiv))
# testing data normality and equal variances to know which t-test to run...
shapiro.test(healthy_meta_wdiv$Observed)
leveneTest(Observed ~ Lived_on_farm, data=healthy_meta_wdiv)
# running 2-tailed t-tests on observed and shannon's div metrices
t.test(healthy_meta_wdiv$Shannon ~ healthy_meta_wdiv$Lived_on_farm, alternative = "two.sided", var.equal = TRUE)
t.test(healthy_meta_wdiv$Observed ~ healthy_meta_wdiv$Lived_on_farm, alternative = "two.sided", var.equal = TRUE)


######### Beta diversity plots ##########
# Check methods for distance function
distanceMethodList


#Bray Curtis Beta diversity
#create a distance metric for bray curtis
braycurt <- distance(healthy, method="bray")

#generating a PCoA matrix using the bray curtis distances
pcoa_braycurt <- ordinate(healthy, method="PCoA", distance=braycurt)

#Generating a PCoA plot from the bray curtis distances
gg_bpcoa <- plot_ordination(healthy, pcoa_braycurt, color = "Lived_on_farm") +
  labs(col="Lived On a Farm")+
  ggtitle("Bray Curtis")
gg_bpcoa

#Saving the plot
ggsave("Bray_curtis_pcoa.png"
       , gg_bpcoa
       , height=4, width=5)



#Unifrac Beta diversity
#create a distance metric for unifrac distances
unif <- distance(healthy, method="unifrac")

#generating a PCoA matrix using the unifrac distances
pcoa_unifrac <- ordinate(healthy, method="PCoA", distance=unif)

#Generating a PCoA plot from the unifrac distances
gg_upcoa <- plot_ordination(healthy, pcoa_unifrac, color = "Lived_on_farm") +
  labs(col="Lived On a Farm")+
  ggtitle("Unifrac")
gg_upcoa

#Saving the plot
ggsave("unifrac_pcoa.png"
       , gg_upcoa
       , height=4, width=5)


#Weighted Unifrac Beta diversity
#create a distance metric for Weighted Unifrac distances
wunif <- distance(healthy, method="wunifrac")


#generating a PCoA matrix using the Weighted Unifrac distances
pcoa_wunifrac <- ordinate(healthy, method="PCoA", distance=wunif)

#Generating a PCoA plot from the Weighted Unifrac distances
gg_wupcoa <- plot_ordination(healthy, pcoa_wunifrac, color = "Lived_on_farm") +
  labs(col="Lived On a Farm")+
  ggtitle("Weighted Unifrac")
gg_wupcoa

#Saving the plot
ggsave("weighted_unifrac_pcoa.png"
       , gg_wupcoa
       , height=4, width=5)



#Jaccard Beta diversity
#create a distance metric for Jaccard distances
jacc <- distance(healthy, method="jaccard")

#generating a PCoA matrix using the Jaccard distances
pcoa_jacc <- ordinate(healthy, method="PCoA", distance=jacc)

#Generating a PCoA plot from the Jaccard distances
gg_japcoa <- plot_ordination(healthy, pcoa_jacc, color = "Lived_on_farm") +
  labs(col="Lived On a Farm")+
  ggtitle("Jaccard")
gg_japcoa

#Saving the plot
ggsave("jaccard_pcoa.png"
       , gg_japcoa
       , height=4, width=5)


colnames(tax_table(healthy))


####### Taxonomic bar plots ########


#Generating a taxonomic bar plot for Phyla
# Convert healthy phyloseq object to relative abundances
healthy_RA <- transform_sample_counts(healthy, function(x) x/sum(x))


# To remove black bars, "glom" by phylum first
taxaplot <- tax_glom(healthy_RA, taxrank = "Phylum", NArm=TRUE)

gg_taxa <- plot_bar(taxaplot, fill="Phylum") + 
  facet_wrap(.~Lived_on_farm, scales = "free_x")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", linewidth=0.005)

gg_taxa

ggsave("plot_taxonomy_Phylum.png"
       , gg_taxa
       , height=8, width =12)

#Generating a taxonomic bar plot for Classes
# To remove black bars, "glom" by class first
taxaplot_class <- tax_glom(healthy_RA, taxrank = "Class", NArm=TRUE)

gg_taxa_class <- plot_bar(taxaplot_class, fill="Class") + 
  facet_wrap(.~Lived_on_farm, scales = "free_x")+
  geom_bar(aes(color=Class, fill=Class), stat="identity", linewidth=0.005)

gg_taxa_class

ggsave("plot_taxonomy_Class.png"
       , gg_taxa_class
       , height=8, width =12)
