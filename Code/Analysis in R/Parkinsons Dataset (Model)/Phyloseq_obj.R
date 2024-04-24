#!/usr/bin/env RScript

#Load in packages 'phyloseq', 'tidyverse', 'vegan', and 'ape'
library(phyloseq)
library(tidyverse)
library(ape)
library(vegan)



########## Importing Files ##############
meta <- read_delim("parkinsons_metadata.txt", delim = "\t")

otu <- read_delim("feature-table.txt", delim = "\t", skip=1)

tax <- read_delim("taxonomy.tsv", delim= "\t")

tree <- read.tree("tree.nwk")

#generating a new otu table with only the samples from non-PD individuals who lived on a farm
otu_hf <- otu %>% 
  select(all_of(c(meta$`#SampleID`))) %>%
  mutate(`#OTU ID` = otu$`#OTU ID`) %>%  #Adding the OTU column back
  select(`#OTU ID`, everything()) #moving the OTU column to the first column




####### Generating Phyloseq Object ##########

#Formating OTU table to make OTU object
otu_matrix <- as.matrix(otu_hf[,-1]) 
rownames(otu_matrix) <- otu_hf$`#OTU ID` #Formating the matrix to display the "#OTU ID" as the rownames
OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE) #make OTU table for phyloseq object
class(OTU) #checking that it's a phyloseq object

#metadata file formatting
meta_table <- as.data.frame(meta[,-1])
rownames(meta_table) <- meta$`#SampleID` #Formating the matrix to display the "#SampleID" as the rownames
META <- sample_data(meta_table) #Creating phyloseq object
class(META) #Checking it's a phyloseq object

#formatting taxonomy table
tax_matrix <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep=";", into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))%>% #use c() because the vector means you separate the different objects into each column
  as.matrix()

tax_matrix <- as.matrix(tax_matrix[,-1])
rownames(tax_matrix) <- tax$`Feature ID`
TAX <- tax_table(tax_matrix)#Make phyloseq object
class(TAX)

#Make final Phyloseq object
phylobj <- phyloseq(OTU, META, TAX, tree)

#Checking the object to make sure it's correct
otu_table(phylobj)
sample_data(phylobj)
tax_table(phylobj)
phy_tree(phylobj)

######## Processing of Phyloseq object ############
#Generating rarefaction curve for all samples
rarecurve(t(as.data.frame(otu_table(phylobj))), cex=0.0001)

rm(list=ls())

#Generating rarefaction curve for healthy individuals (non-PD)
healthy <- subset_samples(phylobj, Disease == "Control")
rarecurve(t(as.data.frame(otu_table(healthy))), cex=0.0005)



#generating new phyloseq object which has been rarefied, cutoff decided from rarefaction curve above
phylobj_raref <- rarefy_even_depth(phylobj, rngseed = 1, sample.size = 8478)

?rarecurve()


#Filtering for the healthy individuals
healthy <- subset_samples(phylobj_raref, Disease == "Control")

healthy

#This should remove any NA's

#Saving data
save(phylobj, file="phyloseqobject_parkinsons.RData")
save(phylobj_raref, file="phyloseqobject_parkinsons_rarefied.RData")
save(healthy, file="phyloseqobject_healthy.RData")

