#Load in packages
library(phyloseq)
library(tidyverse)
library(ape)
library(picante)
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(readr)

########## Importing Files ##############
meta <- read_delim("farm_metadata.tsv", delim = "\t")

otu <- read_delim("feature-table.txt", delim = "\t", skip=1)

tax <- read_delim("taxonomy_farm.tsv", delim= "\t")

tree <- read.tree("tree_farm.nwk")


####### Generating Phyloseq Object ##########

#Formating OTU table to make OTU object
otu_matrix <- as.matrix(otu[,-1])
rownames(otu_matrix) <- otu$`#OTU ID` #Formating the matrix to display the "#OTU ID" as the rownames
OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE) #make OTU table for phyloseq object
class(OTU) #checking that it's a phyloseq object


#metadata file formatting
meta_table <- as.data.frame(meta[,-1])
rownames(meta_table) <- meta$`sample-id` #Formating the matrix to display the "#SampleID" as the rownames
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
newfarm_phyloseq <- phyloseq(OTU, META, TAX, tree)

#Checking the object to make sure it's correct
ooooo <- otu_table(newfarm_phyloseq)
sample_data(newfarm_phyloseq)
tax_table(newfarm_phyloseq)
phy_tree(newfarm_phyloseq)







########## New Farm ASV counts ################

#Import dataset with ASVs and annotated species
annotated_asvs <- read_delim("target_asvs_annotated_species.csv", delim = ";")

#Can be used to double check if certain ASVs are within the phyloseq object
phylo_rownames <- rownames(otu_table(newfarm_phyloseq))
grep("b55ddaf95d9d2217ea087863a7afbbad", phylo_rownames)


#Create list of target ASVs from analyses and remove NAs (Random ASVs used below)
target_asvs <- c(annotated_asvs$ASV)
target_asvs <- na.omit(target_asvs)


#Count the number of occurences per ASV from the whole dataset
asv_counts <- rowSums(otu_matrix)


#Counts and prints the nr of occurences per asv from target_asvs
for (asv in target_asvs) {
  if (asv %in% names(asv_counts)) {
    count_per_asv <- asv_counts[asv]
    cat("ASV", asv, "occurs", count_per_asv, "times.\n")
  } else {
    message(paste("ASV", asv, "not found in the dataset."))
  }
}


####Create a matrix "asv_matrix" with the ASVs and the counts for each ASV from target_asvs, and sort by species
asv_matrix <- matrix(nrow = length(target_asvs), ncol = 3)
colnames(asv_matrix) <- c("ASV_ID", "Count", "Species")
asv_matrix[,3] <- annotated_asvs$Species

for (i in seq_along(target_asvs)) {
  if (asv %in% target_asvs) {
    asv <- target_asvs[i]
    count_per_asv <- asv_counts[asv]
    asv_matrix[i, 1] <- asv
    asv_matrix[i, 2] <- count_per_asv
  } else {
    asv_matrix[i, 1] <- asv
    asv_matrix[i, 2] <- NA
    print("ASV", asv, "not found in datasets")
  }
  
}




####### Creating bubble plot for abundance of each ASV #######
total_counts_df <- as.data.frame(asv_matrix) %>%
  group_by(Species) %>%
  summarize(Count = sum(as.numeric(Count)))

# Create the bubble plot
bub_plot <- ggplot(total_counts_df, 
       aes(x = 0, y = Species, size = as.numeric(Count))) +
  geom_point(shape = 21, fill = "gray", color = "black") +
  scale_size_continuous(range = c(3, 40), breaks = c(1000, 3000, 6000)) +  # Adjust the range of bubble sizes
  labs(title = "Relative abundance of indicator ASVs", x = "", y = "Species", size = "ASV counts per species")+  # Add labels and title
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())+
  theme_bw()# Remove x-axis ticks # Add labels and title


ggsave("Bubble_plot_counts.png", bub_plot, width=8, height= 8, units="in", dpi=1000 )




####### Creating bubble plot for number of samples with each indicator ASV #######

#creating matrix with the number of samples with each indicator ASV 
samp_w_asv <- matrix(nrow = length(target_asvs), ncol = 3)
colnames(samp_w_asv) <- c("ASV_ID", "samples_with_asv", "Species")
samp_w_asv[,3] <- annotated_asvs$Species

for (i in seq_along(target_asvs)) {
  if (asv %in% target_asvs) {
    asv <- target_asvs[i]
    sampcounts <- sum(otu_matrix[asv, ] !=0)
    samp_w_asv[i, 1] <- asv
    samp_w_asv[i, 2] <- sampcounts
  }else {
    samp_w_asv[i, 1] <- NA
    samp_w_asv [i, 2] <- NA
    print("ASV", asv, "not found in datasets")
  }
}

  



total_sampcounts_df <- as.data.frame(samp_w_asv) %>%
  group_by(Species) %>%
  summarize(Count = sum(as.numeric(samples_with_asv)))

# Create the bubble plot
bub_samp_plot <- ggplot(total_sampcounts_df, 
                   aes(x = 0, y = Species, size = Count, fill=total_counts_df$Count)) +
  geom_point(shape = 21, color = "black") +
  geom_text(aes(label = total_counts_df$Count), size=total_sampcounts_df$Count*0.15) +
  scale_fill_gradient(low = "#F1F1F1", high = "#929191", name = "Total ASV abundance") +  # Add legend for fill color and change color
  scale_size_continuous(range = c(min(total_sampcounts_df$Count), max(total_sampcounts_df$Count)), breaks=c(min(total_sampcounts_df$Count), 40, max(total_sampcounts_df$Count))) +  # Adjust the range of bubble sizes
  labs(title = "#Samples with Indicator ASV & total ASV abundance in the Farm Dataset", x = "", y = "Species", size = "#Samples with ASV", fill="Total ASV counts in dataset")+
  theme_bw()+# Add labels and title
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank())
  #theme(axis.text.y = element_text(angle = 45, hjust = 0.5, face = "italic"))




ggsave("Bubble_plot_sample_counts.png", bub_samp_plot, width=8, height= 8, units="in", dpi=1000 )

