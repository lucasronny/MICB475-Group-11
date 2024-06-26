#loading in the necessary libraries
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#reading in the metadata, feature table, taxonomy, and phylogenetic tree
meta <- read_delim("colombia_metadata.txt", delim = "\t")
otu <- read_delim("colombia_table_export/feature-table.txt",
                  delim = "\t", skip = 1)
tax <- read_delim("colombia_taxonomy_export/taxonomy.tsv", delim = "\t")
phylotree <- read.tree("colombia_tree_export/tree.nwk")

#adjusting metadata to phyloseq format and saving it into SAMP
samp_df <- as.data.frame(meta[,-1])    #save everything except #SampleID as a data frame
rownames(samp_df) <- meta$'#SampleID'    #make #SampleID the rownames
SAMP <- sample_data(samp_df)
class(SAMP)    #just to check if it was done properly, should return "phyloseq"

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

tax_matrix <- tax_matrix[,-1]    #save everything except feature ID
rownames(tax_matrix) <- tax$'Feature ID'    #make sample ID the row names
TAX <- tax_table(tax_matrix)
class(TAX)

#merge OTU, SAMP, TAX, and phylotree into a phyloseq object 
#called colombia_phyloseq
colombia_phyloseq <- phyloseq(OTU, SAMP, TAX, phylotree)

#making 5 phyloseq objects, one for each city
bogota <- subset_samples(colombia_phyloseq, city == "Bogota")
medellin <- subset_samples(colombia_phyloseq, city == "Medellin")
cali <- subset_samples(colombia_phyloseq, city == "Cali")
barranquilla <- subset_samples(colombia_phyloseq, city == "Barranquilla")
bucaramanga <- subset_samples(colombia_phyloseq, city == "Bucaramanga")
