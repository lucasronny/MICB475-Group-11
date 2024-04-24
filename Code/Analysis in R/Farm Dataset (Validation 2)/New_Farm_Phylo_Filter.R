#loading in the necessary libraries
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#Loading metadata, feature table, taxonomy, and phylogenetic tree
meta <- read_delim("farm_metadata.tsv", delim = "\t")

otu <- read_delim("table_export/feature-table.txt",
                  delim = "\t", skip = 1)

tax <- read_delim("taxonomy_export/taxonomy.tsv", delim = "\t")

phylotree <- read.tree("root_tree_export/tree.nwk")

#adjusting metadata to phyloseq format and saving it into SAMP
samp_df <- as.data.frame(meta[,-1])    #save everything except sample-id as a data frame
rownames(samp_df) <- meta$'sample-id'    #make sample-id the rownames
SAMP <- sample_data(samp_df)
class(SAMP) 

#adjusting feature table to phyloseq format and saving it into OTU
otu_matrix <- as.matrix(otu[,-1])    #save everything except #OTU ID as a matrix
rownames(otu_matrix) <- otu$'#OTU ID'    #make #OTU ID the row names
OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)
class(OTU)

#adjusting taxonomy table to phyloseq format and saving it into TAX
tax_matrix <- tax |>    #convert taxon into separate taxa rank columns
  select(-Confidence) |>
  separate(col = Taxon, sep = "; ",
           into = c("Domain", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species")) |>
  as.matrix()

tax_matrix <- tax_matrix[,-1]   #save everything except feature ID
rownames(tax_matrix) <- tax$'Feature ID'    
TAX <- tax_table(tax_matrix)
class(TAX)

#merge OTU, SAMP, TAX, and phylotree into a phyloseq object 
#file named newfarm_phyloseq
newfarm_phyloseq <- phyloseq(OTU, SAMP, TAX, phylotree)

#Filter to find species of interest
#Must be repeated per species of interest
#Example for Erysipelotrichaceae holdemanella filtering below

#List Species of Interest as c("o__Order", "f__Family", "g__Genus")
interestSpecies <- c("o__Erysipelotrichales", "f__Erysipelotrichaceae", "g__Holdemanella")

#Run code below in order to print out filtered species of interest in phyloseq object
spcOrder <- subset_taxa(nepal_phyloseq, Order==interestSpecies[1]) #Filter phyloseq object for specified Order
spcFamily <- subset_taxa(spcOrder, Family==interestSpecies[2]) #Filter phyloseq object again for specified Family
spcGenus <- subset_taxa(spcFamily, Genus==interestSpecies[3]) #Filter phyloseq object again for specified Genus
tax_table(spcGenus)  # Print out the tax table
#Save object as CSV for ease of BLAST searching
fileName <- "Erysi.csv"
write.csv(tax_table(spcGenus), file = fileName)
