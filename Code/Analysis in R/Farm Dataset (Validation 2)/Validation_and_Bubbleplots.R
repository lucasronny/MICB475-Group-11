#Load in packages
library(phyloseq)
library(tidyverse)
library(ape)
library(ggtext)
library(glue)



########## Importing Files ##############
meta <- read_delim("farm_metadata.tsv", delim = "\t")

otu <- read_delim("feature-table-farm.txt", delim = "\t", skip=1)

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
newfarm_phyloseq1 <- phyloseq(OTU, META, TAX, tree)

#Checking the object to make sure it's correct
otu_table(newfarm_phyloseq1)
sample_data(newfarm_phyloseq1)
tax_table(newfarm_phyloseq1)
phy_tree(newfarm_phyloseq1)


# removing samples from American non-farmers
newfarm_phyloseq <- subset_samples(newfarm_phyloseq1, Country == "Nepal")


########## New Farm ASV counts ################

#Import dataset with target ASVs and annotated species
annotated_asvs <- read_delim("target_asvs_annotated_species.csv", delim = ";")


#Create list of target ASVs from analyses and remove NAs (Random ASVs used below)
target_asvs <- c(annotated_asvs$ASV)
target_asvs <- na.omit(target_asvs)


#Count the number of occurences per ASV from the whole dataset
asv_counts <- rowSums(otu_table(newfarm_phyloseq))

#Counts and prints the nr of occurences per asv from target_asvs
for (asv in target_asvs) {
  if (asv %in% names(asv_counts) & asv_counts[asv] > 0) {
    count_per_asv <- asv_counts[asv]
    cat("ASV", asv, "occurs", count_per_asv, "times.\n")
  } else {
  }
}


####Create a matrix "asv_matrix" with the ASVs and the counts for each ASV from target_asvs, and sort by species
asv_matrix <- matrix(nrow = length(target_asvs), ncol = 3)
colnames(asv_matrix) <- c("ASV_ID", "Count", "Species")
asv_matrix[,3] <- annotated_asvs$Species

#adds each ASV corresponding to an indicator taxa to a table, including species name and total ASV count
for (i in seq_along(target_asvs)) {
  asv <- target_asvs[i]
  if (asv %in% names(asv_counts) & asv_counts[asv] > 0) {
    count_per_asv <- asv_counts[asv]
    asv_matrix[i, 1] <- asv
    asv_matrix[i, 2] <- count_per_asv
  } else {
    asv_matrix[i, 1] <- asv
    asv_matrix[i, 2] <- NA
  }
}


# Filter out rows with NA values in the second column
asv_matrix <- asv_matrix[!is.na(asv_matrix[, 2]), ]


#Summarizing the data from the previous for loop and group by species.
total_counts_df <- as.data.frame(asv_matrix) %>%
  group_by(Species) %>%
  summarize(Count = sum(as.numeric(Count)))



#creating matrix with the number of samples with each indicator ASV 
samp_w_asv <- matrix(nrow = length(target_asvs), ncol = 3)
colnames(samp_w_asv) <- c("ASV_ID", "samples_with_asv", "Species")
samp_w_asv[,3] <- annotated_asvs$Species

for (i in seq_along(target_asvs)) {
  asv <- target_asvs[i]
  if (asv %in% target_asvs) {
    sampcounts <- sum(otu_table(newfarm_phyloseq)[asv, ] !=0)
    if (sampcounts > 0) {
      samp_w_asv[i, 1] <- asv
      samp_w_asv[i, 2] <- sampcounts
    } else {
      samp_w_asv[i, 1] <- NA
      samp_w_asv[i, 2] <- NA
      print(paste("ASV", asv, "not found in datasets"))
    }
  } else {
    samp_w_asv[i, 1] <- NA
    samp_w_asv[i, 2] <- NA
    print(paste("ASV", asv, "not found in datasets"))
  }
}

# Filter out rows with NA values in the second column
samp_w_asv <- samp_w_asv[!is.na(samp_w_asv[, 2]), ]


#total nr of samples in dataset
a <- ncol((otu_table(newfarm_phyloseq)))


#Summarizing the prevalence of ASVs in the dataset
total_sampcounts_df <- as.data.frame(samp_w_asv) %>%
  group_by(Species) %>%
  summarize(Count = sum(as.numeric(samples_with_asv)))

#Removing any ASVs with only 1 count
total_sampcounts_df <- total_sampcounts_df[(total_sampcounts_df[, 2]) > 1 , ]

#Calculating the percentage prevalence of each ASV
total_sampcounts_df <- mutate(total_sampcounts_df, Count = (Count/a))

#Making the abundance counts fit the same ASVs as the prevalence counts
total_counts_df <- filter(total_counts_df, Species %in% total_sampcounts_df$Species)

#subsetting genus and species to format itallics (Redundant code but has to stay for rest of code to work)
total_sampcounts_df$Genus <- gsub("^(\\S+)\\s.*", "\\1", total_sampcounts_df$Species)
total_sampcounts_df$Species_name <- gsub("^\\S+\\s(.*)", "\\1", total_sampcounts_df$Species)

#formatting itallics for all lines except for the Flavobacteriaceae which is a family
total_sampcounts_df <- total_sampcounts_df %>% mutate(
  name = glue("{Genus} {Species_name}"),
  name = if_else(row_number() == 4, glue("{Genus} {Species_name}"), name)
)


# Create the bubble plot
bub_samp_plot <- ggplot(total_sampcounts_df, 
                        aes(x = 0, y = name, size = Count, fill=total_counts_df$Count)) +
  geom_point(shape = 21, color = "black") +
  geom_text(data = total_sampcounts_df, 
            aes(label = total_counts_df$Count), 
            size = total_sampcounts_df$Count * 10.15)+
  xlim(-0.5, 0.5)+
  #annotate("text", label=total_counts_df$Count[4], x=0.2, y="Uncultured Flavobacteriaceae", size=5)+
  scale_fill_gradient(low = "#F1F1F1", high = "#929191", name = "Total ASV abundance") +  # Add legend for fill color and change color
  scale_size_continuous(range = c(70*min(total_sampcounts_df$Count), 70*max(total_sampcounts_df$Count)), breaks=c(min(total_sampcounts_df$Count), 0.5, max(total_sampcounts_df$Count)), 
                        labels = scales::percent_format()) +  # Adjust the range of bubble sizes
  labs(x = "", y = "Taxa", size = "% of Samples with ASV", fill="Total ASV counts in dataset")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "italic"))+
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank(),
        axis.text.y = element_markdown())



ggsave("Bubble_plot_sample_counts_no_1s.png", bub_samp_plot, width=8, height= 8, units="in", dpi=1000 )

