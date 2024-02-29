# MICB475-Group-11: Project 2
A repository for documentation, coding scripts for our project.<br>
For meeting minutes and agendas, see [Meeting notes](https://github.com/lucasronny/MICB475-Group-11/tree/3647d9e5b9e84cc3018353200aefe5036ffafb39/notes).

## Team
Ekaterina Galysheva <br> Keegan McDonald <br> Alix Najera Mazariegos <br> Lucas Rönn <br> Leonardo Wu

## Table of Contents
[Summary](#summary) <br>
[Project Aims](#project-Aims)<br>
[Lab Notebook](#lab-Notebook)<br>

## Summary
The project aims to analyze the microbiome of healthy individuals living on a farm to create a predictive model for determining whether a test individual of unknown background lives on a farm or not. First, we will determine whether the microbiome diversity metrics are different between people living in an agricultral region compared to those that live in urban centers. Then we will select bugs that are differentially abundant in the individuals coming from agricultral communities to include in our predictive model. We will test the accuracy of our model by applying it to test subjects with known background. If successful, our model should predict, based on the richness and abudance of the microbiome, whether a test subject lives on a farm or not.

## Project Aims
1. Initial Processing of the Microbiome Datasets in QIIME2: process the parkinsons and columbia datasets in QIIME2 to prepare data for downstream analysis. Quality control of sample reads and proper sampling depth for each dataset will be determined here. After denoising and clustering, export the taxonomy, table, metadata, and the phylogenetic tree into a local drive to perform diversity metrics analysis in R.
2. Generate microbiome diversity metrics in R studio: for the parkinson's data, filter out all healthy patients and remove n/a values in the Lived_on_farm column, and generate different &alpha; and &beta; diversity metrics, comparing the samples from donors that lived on a farm to those that did not.
3. Core microbiome analysis
4. Indicator taxon analysis
5. Model validation

## Lab Notebook

### Processing of the Parkinson's dataset in QIIME2 (INSERT DATE HERE):

Results:
- 25 samples with YES-farm, NO-PD

To inquire in meeting: 
- Should I have filtered the mitochondria + chlorplast?
- Was 251 a good cutoff for filtering (i.e the last base pair, ~25 quality score)? 
- All the "NO DISEASE & FARM=YES" have 0 in constipation?
- Help with rarefaction curve - should we cut 1 of the samples off?
- Further steps?
- organize Github (dates, progress made)
- NEW: I recall in the last meeting that we might be pooling the microbiomes in the Colombia dataset by city (for when we test our model). Each city has a different number of samples, so should we rarefy the Colombia dataset too?

### Core Microbiome Analysis of the Farm dataset (healhty idividuals from the Parkinson's dataset) | February 28, 2024:

Procedure:
- load the phyloseq object containing the Parkinson's dataset filtered for only healhty patients
- convert the frequencies in the OTU table
- create vectors for all ASVs found in the samples that had been exposed to an agricultural setting and in those that had not been exposed
- core microbiome analysis of either vector
	- relative abundance threshold: 1% (only take the ASV if its relative abundance is higher than 1%; this filters out all the rare ASVs and retains all abundant and non-abundant ASVs)
	- prevalence: 50% (only accept the ASV if it is found in over 50% of the samples)
- subset a new table for both exposed and unexposed groups by isolating the ASVs found in the core microbiome analysis
	- create 1 table per group
- generate a barplot of ASVs identified in either group
- generate a Venn diagram representing the unique and shared ASVs between groups

Results:
- only 3 core ASVs have been identified in the unexposed group, while 4 were identified in the exposed group
- exposed group ASV taxonomic classification:

|ASV code                        |Domain  |Phylum      |Class      |Order          |Family           |Genus             |Species               |
|--------------------------------|--------|------------|-----------|---------------|-----------------|------------------|----------------------|
|d71fe67e73a45e78b36ab2586bb5dcc8|Bacteria|Firmicutes  |Clostridia |Oscillospirales|_Ruminococcaceae_|_Faecalibacterium_|_NA_                  |
|96c7cedbae877361d087a8b39206d408|Bacteria|Firmicutes  |Clostridia |Oscillospirales|_Ruminococcaceae_|_Faecalibacterium_|_NA_                  |
|eb1cb381c314ecbbf46265b2e2351ee8|Bacteria|Bacteroidota|Bacteroidia|Bacteroidales  |_Bacteroidaceae_ |_Bacteroides_     |_NA_                  |
|76d6be90e19a660c9a84b35fa31d801c|Bacteria|Bacteroidota|Bacteroidia|Bacteroidales  |_Bacteroidaceae_ |_Bacteroides_     |_Bacteroides_vulgatus_|

- unexposed group ASV taxonomic classification:

|ASV code                        |Domain  |Phylum      |Class      |Order          |Family           |Genus             |Species               |
|--------------------------------|--------|------------|-----------|---------------|-----------------|------------------|----------------------|
|d71fe67e73a45e78b36ab2586bb5dcc8|Bacteria|Firmicutes  |Clostridia |Oscillospirales|_Ruminococcaceae_|_Faecalibacterium_|_NA_                  |
|96c7cedbae877361d087a8b39206d408|Bacteria|Firmicutes  |Clostridia |Oscillospirales|_Ruminococcaceae_|_Faecalibacterium_|_NA_                  |
|76d6be90e19a660c9a84b35fa31d801c|Bacteria|Bacteroidota|Bacteroidia|Bacteroidales  |_Bacteroidaceae_ |_Bacteroides_     |_Bacteroides_vulgatus_|

- as a result of the taxonomic analysis, one ASV unique to the exposed group was identified. This ASV belongs to the genus _Bacteroides_. The taxonomic analysis could not identify the species to which this species belongs to.
- _Bacteroides_ genus is commonly found in colon and has both positive and negative effects on overall health depending on their interactions with other microbes in the gut microbiome [Zafar, 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7872030/).

- the bar plot of the distribution of the identified ASVS
	- overall, the 2 ASVs show a higher relative abundance in the unexposed group
	- there is higher variability of the _Bacteroides_ genus in both the exposed and unexposed groups compared to the _Faecalibacterium_
<br>![Core Microbiome Analysis](https:/Code/R_stuff/parkinsons/core_ASVs_barplot(r.ab_0.01,prev_0.5).png)

- Venn diagram generated for the core microbiome analysis
<br>![Core Microbiome Analysis Venn Diagram](https:/Code/R_stuff/parkinsons/venn_core(r.ab_0.01,prev_0.5).png)

