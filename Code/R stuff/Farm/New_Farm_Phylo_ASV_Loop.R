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

#Keegan -> ASV finding for loop

#Create list of target ASVs from analyses (Random ASVs used below)
target_asvs <- c("7696986a6aaff60166f5aa47d67f02d8", 
                 "cbd66fa79a3e7fa1c33d29ed908de048", 
                 "3b092a59cd0e2d57483420cc994ec49c",
                 "9c15e7e160a34bddcf638ffa932425fd",
                 "3f38509af22157cdcb4006457529d25e")

#Can be used to double check if certain ASVs are within the phyloseq object
phylo_rownames <- rownames(otu_table(newfarm_phyloseq))
grep("3f38509af22157cdcb4006457529d25e", phylo_rownames)


# Loop through each ASVs to find target ASVs
for (asv in target_asvs) {     
  if (asv %in% rownames(otu_table(newfarm_phyloseq))) { # Check if the ASV exists in the phyloseq object
    asv_count <- sum(rownames(otu_table(newfarm_phyloseq)) == asv)# If the ASV exists, count occurrences
    cat("ASV", asv, "occurs", asv_count, "times.\n")    # Print the ASV name and count
  } else {
    message(paste("ASV", asv, "not found in the dataset."))
  }
}

# Code for comparing ASVs from Parkinsons and New Farm data
not_found_count <- 0
found_count <- 0
for (asv in rownames(otu_table(phylobj))) {     
  if (asv %in% rownames(otu_table(newfarm_phyloseq))) { # Check if the ASV exists in the phyloseq object
    asv_count <- sum(rownames(otu_table(newfarm_phyloseq)) == asv)# If the ASV exists, count occurrences
    cat("ASV", asv, "occurs", asv_count, "times.\n")    # Print the ASV name and count
    found_count <- found_count + 1
  } else {
    not_found_count <- not_found_count + 1
    message(paste("ASV", asv, "not found in the dataset."))
  }
}
cat("Total ASVs not found:", not_found_count, "\n")
cat("Total ASVs found:", found_count, "\n")




