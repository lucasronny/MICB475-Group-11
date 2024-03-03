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


### Core Microbiome Analysis of the Farm dataset (healhty idividuals from the Parkinson's dataset) | February 28, 2024, EG:

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
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Code/R%20stuff/parkinsons/core_ASVs_barplot(r.ab_0.01%2Cprev_0.5).png">

- Venn diagram generated for the core microbiome analysis
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Code/R%20stuff/parkinsons/venn_core(r.ab_0.01%2Cprev_0.5).png" width = "675" height = "450">

CONT [BELOW](#core-Microbiome-Analysis,-continued:-changing-parameters-to-increase-the-number-of-identified-core-microbiome-species-|-March-3,-2024,-EG:)


### Indicator Species Analysis of the Farm dataset (healhty idividuals from the Parkinson's dataset) | February 29, 2024, EG:

Protocol:

- packages required: tidyverse, phyloseq, indicspecies
- load the non-rarefied phyloseq object for the Farm dataset
- glom the tax table to Genus, make sure to include the NA rows
- conver the frequences for the glommed Genus by transform_sample_counts()
- use multipatt() to calculate indicator values for the dataset and cluster by the `Lived_on_farm` column
- combine the selected ASVs with the data in the tax table
	- assign the name "ASV" to the first column in the tax table
	- filter for the rows with a p-value < 0.05

Results:
The following list of ASVs was generated in the ISA:

|ASV code                        |Domain  |Phylum           |Class          |Order             |Family                     |Genus                   |Species               |Lived on farm?|Stat |p-value|
|--------------------------------|--------|-----------------|---------------|------------------|---------------------------|------------------------|----------------------|--------------|-----|-------|
|2fbf02f0d8043728d6e93f0bfa432d85|Bacteria|Actinobacteriota |Coriobacteriia |Coriobacteriales  |_Eggerthellaceae_          |_NA_                    |_NA_                  |Yes           |0.209|0.035  |
|23667231b1722eb51bc7ccfaa5e5c38a|Bacteria|Bacteroidota     |Bacteroidia    |Bacteroidales     |_Prevotellaceae_           |_Prevotellaceae_UCG-001_|_uncultured bacterium_|Yes           |0.233|0.010  |
|898f2833fcd964a4cfae006e66bf10b6|Bacteria|Bacteroidota     |Bacteroidia    |Bacteroidales     |_Barnesiellaceae_          |_NA_                    |_NA_                  |Yes           |0.184|0.040  |
|455cd7d5df48c77aedb33f14590d1e78|Bacteria|Firmicutes       |Bacilli        |Erysipelotrichales|_NA_                       |_NA_                    |_NA_                  |Yes           |0.228|0.035  |
|cbd66fa79a3e7fa1c33d29ed908de048|Bacteria|Firmicutes       |Bacilli        |Erysipelotrichales|_Erysipelatoclostridiaceae_|_Catenibacterium_       |_NA_                  |Yes           |0.329|0.045  |
|b7fb7b8de1b4e013e7255cec428893a7|Bacteria|Verrucomicrobiota|Lentisphaeria  |Victivallales     |_Victivallaceae_           |_Victivallaceae_        |_uncultured rumen_    |Yes           |0.325|0.015  |

- According to the literature, these ASVs are present abudantly in healthy human microbiome
- None of the ASVs were identified to the Species level, but 3 Genera were identified: _Prevotllaceae_, Cotenibacterium_, and _Victivallaceae_
- the Stat values are only considered meaningful when > 0.7. The highest Stat value in this analysis is 0.329 for the ASV in the _Catenibacterium_ genus

CONT BElOW


### Core Microbiome Analysis, continued: changing parameters to increase the number of identified core microbiome species | March 3, 2024, EG:

BEGIN [ABOVE](#core-Microbiome-Analysis-of-the-Farm-dataset-(healhty-idividuals-from-the-Parkinson's-dataset)-|-February-28,-2024,-EG:)

- to increase the number of species in the core mircobiome, make the detection and prevalence criteria less stringent
	- set detection = 0, prevalence = 0.1

Results:
- there are 311 species in the exposed group's core microbiome, and 280 species in the unexposed group's core microbiome 
- bar graph plot generated:

### Indicator Species Analysis of the Farm dataset (healhty idividuals from the Parkinson's dataset) | February 29, 2024, EG:

- only 6 ASVs were identified as indictor species; these are too few for downstream analysis
- addtionally, the Stat values for these ASVs were below, which is considered a ________________

Protocol:
- perform the same analysis, but glom by species
- any ASVs that could not be identified up to the species level, use the table.qzv file generated in QIIME2 analysis of the dataset to blast these ASVs manually and fill in the table

Results:
- 19 indicator species were identified
