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

### Processing of the Parkinson's dataset in QIIME2 and Rstudio:
January 29, 2024, LR

Purpose: To denoise, filter and rarefy the Parkinson's dataset to retrieve only healthy individuals for the comparison between individuals who have lived on a farm and individuals who have not. Then to export the processed OTU table, phylogenetic tree, metadata and the taxonomy table to create a phyloseq object in R for further analysis. 

Procedure:
- Combine manifest file with the sequence reads .qza file to demultiplex samples
- Using the Quality scores per base pair plot (from below), trim and denoise the samples to a read depth of n=251 using DADA2
- Filter out the non-PD individuals from the .qza file to leave only the healthy individuals
- Generate the 4 files needed to create a phyloseq object; taxonomy.tsv, feature-table.txt, metadata.txt and tree.nwk; Then export
- Import the files into Rstudio and generate a phyloseq object
- Generate alpha rarefaction curve in Rstudio and in QIIME2 to determine rarefaction depth
- Rarefaction depth of 8478 was chosen as the optimal rarefaction depth, retaining the highest number of samples while retaining most features.

Quality score for each base in the sequence reads.
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/images%20and%20files/parkinsons/filter_depth.png">


Alpha rarefaction curve generated in R, for only healthy individuals
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/images%20and%20files/parkinsons/rarefaction//rarefactioncurve.jpg">



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

### Alpha and Beta Diversity of the Parkinson's dataset in R:
February 10, 2024, LR

### Core Microbiome Analysis of the Farm dataset (healhty idividuals from the Parkinson's dataset)
February 28, 2024, EG:

Purpose: <br> Having selected all the healthy individuals in the Parkinson's dataset, we can now proceed to create a list of species characteristic of the mcirobiota of individuals that have been exposed to an agricultural environment (in our case, a farm). One of the analysis we decided to perform is the core microbiome analysis, which allows to determine the shared and unique species between selected groups of samples. The analysis will be performed in R using the microbiome package. <br> A high number of shared species will indicate that the individuals in the exposed and unexposed groups have a similar gut microbiome composition, whereas a diffirence between the 2 groups would be indicated by a low number of shared species. The group with a higer \u03B1-diversity will have a higher number of unique species compared to the other group. This will enable us to qualitatively assess the extent of similarity and divergence between the 2 groups. <br> Downstream, we will perform indicator species analysis, determine which ASVs are common between the indicator species and the species present uniquely in exposed individuals. This will comprise the basis of the predictive taxonomic model.

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
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/images%20and%20files/parkinsons/core/core_ASVs_barplot(r.ab_0.01%2Cprev_0.5).png">

- Venn diagram generated for the core microbiome analysis
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/images%20and%20files/parkinsons/core/venn_core(r.ab_0.01%2Cprev_0.5).png">

continued below: [Core Microbiome Analysis, continued](#Core-Microbiome-Analysis,-continued)


### Indicator Species Analysis of the Farm dataset (healhty idividuals from the Parkinson's dataset)
February 29, 2024, EG:

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

|ASV code                        |Domain  |Phylum           |Class          |Order             |Family                     |Genus                   |Species               |Lived on farm?|Indicator value (IV)|_p_    |
|--------------------------------|--------|-----------------|---------------|------------------|---------------------------|------------------------|----------------------|--------------|--------------------|-------|
|2fbf02f0d8043728d6e93f0bfa432d85|Bacteria|Actinobacteriota |Coriobacteriia |Coriobacteriales  |_Eggerthellaceae_          |_NA_                    |_NA_                  |Yes           |0.209               |0.035  |
|23667231b1722eb51bc7ccfaa5e5c38a|Bacteria|Bacteroidota     |Bacteroidia    |Bacteroidales     |_Prevotellaceae_           |_Prevotellaceae_UCG-001_|_uncultured bacterium_|Yes           |0.233               |0.010  |
|898f2833fcd964a4cfae006e66bf10b6|Bacteria|Bacteroidota     |Bacteroidia    |Bacteroidales     |_Barnesiellaceae_          |_NA_                    |_NA_                  |Yes           |0.184               |0.040  |
|455cd7d5df48c77aedb33f14590d1e78|Bacteria|Firmicutes       |Bacilli        |Erysipelotrichales|_NA_                       |_NA_                    |_NA_                  |Yes           |0.228               |0.035  |
|cbd66fa79a3e7fa1c33d29ed908de048|Bacteria|Firmicutes       |Bacilli        |Erysipelotrichales|_Erysipelatoclostridiaceae_|_Catenibacterium_       |_NA_                  |Yes           |0.329               |0.045  |
|b7fb7b8de1b4e013e7255cec428893a7|Bacteria|Verrucomicrobiota|Lentisphaeria  |Victivallales     |_Victivallaceae_           |_Victivallaceae_        |_uncultured rumen_    |Yes           |0.325               |0.015  |

- According to the literature, these ASVs are present abudantly in healthy human microbiome
- None of the ASVs were identified to the Species level, but 3 Genera were identified: _Prevotllaceae_, Cotenibacterium_, and _Victivallaceae_
- Indicator values indicate how strongly an ASV is associated with a group (here, with the farm group)
		- IV = 100 x relative abundance x relative frequency
				- relative abundance: how many individuals in the group have the ASV
				- relative frequency: **proportion of sites in the group that have the individual**
		- values are only considered meaningful when > 0.7. The highest IV in this analysis is 0.329 for the ASV in the _Catenibacterium_ genus

continued below: [Indicator Species Analysis, continued](#Indicator-species-analysis.+continued)

### Core Microbiome Analysis, continued
Changing parameters to increase the number of identified core microbiome species
March 3, 2024, EG:

- to increase the number of species in the core mircobiome, make the detection and prevalence criteria less stringent
	- set detection = 0, prevalence = 0.1

Results:
- there are 311 species in the exposed group's core microbiome, and 280 species in the unexposed group's core microbiome 
- bar graph plot generated:
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/images%20and%20files/parkinsons/core/core_ASVs_barplot(r.ab_0%2Cprev_0.1).png">
- there are more ASVs included in the current analysis because the parameters are less stringent

- and the Venn diagram
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/images%20and%20files/parkinsons/core/venn_core(r.ab_0%2Cprev_0.1).png">
- overall, there are more ASVs included into the analysis with the less stringent parameters, consistent with the bar plot generated
- the farm group has more unique species compared to the non-farm group, meaning that 

### Indicator Species Analysis, continued
February 29, EG

- only 6 ASVs were identified as indictor species; these are too few for downstream analysis
- addtionally, the IVs for these ASVs were below 0.7, which indicates that 

Protocol:
- perform the same analysis, but glom by species
- subset ASVs with IV>0.2 and p<0.05
- blast the species selected using the NIH BLAST&reg; tool

Results:
- 19 indicator species were identified
- see this [file](https://github.com/lucasronny/MICB475-Group-11/blob/main/images%20and%20files/parkinsons/isa/isa_frame_high1.csv) for tables listing the taxonomic classification of ASVs performed using the representative sequences in QIIME2 and NIH BLAST &reg; tool classification

### Core Microbiome and Indicator Species Analysis: Combining Results to Identify Species for the Predictive Model
March 4, EG:

Purpose: <br>
The core microbiome indicator species analyses have been performed to the necessary level of stringency. To generate the predictive taxonomic model, we will combine the species identified 

Procedure:
- tables containing the core microbiome species and indicator species for the farm-exposed group were joined using the inner_join() function, basing the join on the ASV column in either table (see code here)
- the same operations were performed for the indicator species and the core microbiome members unique to the exposed group

Results:
- there were 3 species in common betwen the total ASVs belonging to the exposed group and the indicator species:

|ASV code                        |Domain  |Phylum           |Class           |Order             |Family            |Genus                         |Species                    |Lived on farm?|Indicator value (IV)|_p_    |
|--------------------------------|--------|-----------------|----------------|------------------|------------------|------------------------------|---------------------------|--------------|--------------------|-------|
|bfab02a86c5187ed451db10d9b81b8d5|Bacteria|Firmicutes       |Clostridia      |Lachnospirales    |_Lachnospiraceae_ |_Eubacterium ventriosum group_|_uncultured bacterium_     |Yes           |0.420               |0.005  |
|**cbd66fa79a3e7fa1c33d29ed908de048**|Bacteria|Bacillota        |Erysipelotrichia|Erysipelotrichales|_Coprobacillaceae_|_Catenibacterium_             |_Catenibacterium misuokai_ |Yes           |0.327               |0.015  |
|**3f9fcc47fb363e0fdfa75dc56bb107f5**|Bacteria|Bacteroidota     |Bacteroidia     |Bacteroidales     |_Tannerellaceae_  |_Parabacteroides_             |_Parabacteroides johnsonii_|Yes           |0.303               |0.030  |

- all the sequences were blasted using the NIH BLAST&reg; tool
- all ASVs were identified down to the species level, except for the species belonging to the _Lachnospiraceae_ family
- additionally the uncultured bacterium ASV and _Parabacteroides johnsonii_ were in common between the unique core microbiome species in the exposed group and the indicator species
- when comparing the species from core microbiome analysis unique to and the indicator species for the farm individuals, only 2 ASVs are identified: _Catenibacterium misuokai_ and _Parabacteroides johnsonii_ (the 2 bolded in the table) 

Conclusion: 
- the common species are very few
- additionally, the IV are low (IV<0.42), indicating that the species identified are not unique to the farm group
	- note that all the ubiquotous species (i.e., common to both groups) are excluded by the _p_-value threshold
- only the species for the 2 bolded ASVs will be used as the species for the model
	- this makes for a rather weak model, considering there are only 2 species and the indicator values are quite low

 ### Ranking Cities by Agricultural Exposure to Predict Which Will Be Identified as "Farming Populations" by The Predictive Model
March 13, ANM, LW
<br>Purpose: To build a hypothesis to use in testing our model. Before looking for the model species in the Colombian population (grouped by city), we propose which cities are more or less likely to be identified as agricultural communities by our model. We will calculate **(fill in)** to assign each selected city a ranking. The highest rankings will be given to the cities which we believe to be likely agricultural communities, the lower ranking will be given to the species that are the least likely.

Procedure:
- Data taken from the Censo Nacional Agropecuario 2014
- Cities in the Colombia dataset were classified based on their associated Department.
- Using the Census we looked and the data for Land Cover (the amount of land in hectares that is used for agricultural activities in a specific Department). Then we looked and the tables covering the amount of agricultural workers for each Department.
- The total amount of hectares in land cover and the total amount of workers respectively, were used to calculate the proportion of Agricultural land and workers for each Department.
- City ranking made by looking at the proportion results.


Results:
| City        | Associated Department | Land Cover | Proportion of Agricultural Land | Total Workers | Proportion of workers |
|-------------|------------------------|------------|----------------------------------|---------------|------------------------|
| Bogota      | DC and Cundinamarca   | 1,482,641  | 63.07%                           | 285,275       | 2.70%                  |
| Barranquilla| Atlantico             | 192,664    | 65.70%                           | 18,905        | 0.75%                  |
| Medellin    | Antioquia             | 3,426,582  | 54.75%                           | 247,384       | 4.14%                  |
| Cali        | Valle del Cauca       | 1,100,046  | 54.11%                           | 88,327        | 2.33%                  |
| Bucaramanga | Santander             | 1,946,487  | 64.70%                           | 170,329       | 8.48%                  |

Bogota data collected by adding census results from both Bogota D.C. and Cundinamarca since officially the city of Bogota is considered to be part of Cundinamarca but was treated separately in the census. 

### Rankings (most to least agricultural):

1. Bucaramanga
2. Medellin
3. Bogota
4. Cali
5. Barranquilla


### Applying The Model to The Colombia Dataset
March 13, KM
<br>Purpose: To look for the model species in each city in the Columbia dataset to determine whether 

Procedure:

Results:


