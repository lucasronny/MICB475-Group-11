# MICB475-Group-11: Project 2
A repository for documentation and coding scripts for our project.<br>
For meeting minutes and agendas, see [Meeting notes](https://github.com/lucasronny/MICB475-Group-11/tree/3647d9e5b9e84cc3018353200aefe5036ffafb39/notes).

## Team
Ekaterina Galysheva <br> Alix Najera Mazariegos <br> Keegan McDonald <br> Lucas Rönn <br> Leonardo Wu

## Table of Contents
[Summary](#summary) <br>
[Project Aims](#project-Aims)<br>
[Lab Notebook](#lab-Notebook)<br>

## Summary
The project aims to analyze the microbiome of healthy individuals living on a farm to create a predictive model for determining whether a test individual of unknown background lives on a farm or not. First, we will determine whether the microbiome diversity metrics are different between people living in an agricultral region compared to those that live in urban centers. Then we will select bugs that are differentially abundant in the individuals coming from agricultral communities to include in our predictive model. We will test the accuracy of our model by applying it to test subjects with known background. If successful, our model should predict, based on the richness and abudance of the microbiome, whether a test subject lives on a farm or not. <br>

Upon completing core microbiome analysis, we decided to alter the overall research aim of our project to designing a methodology for the development and training of taxonomic predictive models based on microbiome composition that would allow to identify taxa characteristic of different environments.

## Project Aims
1. Initial Processing of the Microbiome Datasets in QIIME2: process the Parkinsons and Columbia datasets in QIIME2 to prepare data for downstream analysis. Quality control of sample reads and proper sampling depth for each dataset will be determined here. After denoising and clustering, we will export the taxonomy, table, metadata, and the phylogenetic tree into a local drive to perform diversity metrics analysis in R.
2. Generate microbiome diversity metrics in R studio: for the Parkinson's data, filter out all healthy patients and remove n/a values in the Lived_on_farm column, and generate different &alpha; and &beta; diversity metrics, comparing the samples from donors that lived on a farm (Farm group) to those that did not (No-farm group).
3. Core microbiome analysis: determine the unique and shared taxa between the two groups of interest.
4. Indicator taxon analysis: determine the indicator taxa for the Farm group.
5. Reconciliation: keep only the taxa common to the core microbiome subset "Unique to Farm group" and "Farm Indicator Species". These will be the model taxa.
7. Model validation: calculatig the abundance and prevalence for the model taxa in a validation dataset or subset which contains only individuals who live in agricultural settings.

## Lab Notebook

### Processing of the Parkinson's dataset in QIIME2 and Rstudio:
January 29, 2024, LR

__Purpose:__ To denoise, filter and rarefy the Parkinson's dataset to retrieve only healthy individuals for the comparison between individuals who have lived on a farm and individuals who have not. Then to export the processed OTU table, phylogenetic tree, metadata and the taxonomy table to create a phyloseq object in R for further analysis. 

__Procedure:__
- Combine manifest file with the sequence reads .qza file to demultiplex samples
- Using the Quality scores per base pair plot (see below), trim and denoise the samples to a read depth of n=251 using DADA2 (effectively no trimming was done because all Q values were above the set threshold of 32).
- Filter out the non-PD individuals from the .qza file to leave only the healthy individuals
- Generate the 4 files needed to create a phyloseq object: taxonomy.tsv, feature-table.txt, metadata.txt and tree.nwk. Export thses files as a folder to a local drive.
- Import the files into Rstudio and generate a phyloseq object
- Generate alpha rarefaction curve in Rstudio and in QIIME2 to determine rarefaction depth.
- Rarefaction depth of 8478 was chosen as the optimal rarefaction depth, retaining the highest number of samples while retaining most features. This was decided using the alpha rarefaction curve below, and using the QIIME2 online viewer.

Quality score for each base in the sequence reads:
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)/filter_depth.png">

Alpha rarefaction curve generated in R, for only healthy individuals:
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)/rarefaction/rarefactioncurve.jpg">

__Results__ after denoising and rarefaction:
- 24 samples with YES-farm, NO-PD retained
- 64 samples with NO-farm, NO-PD retained


### Alpha and Beta Diversity of the Parkinson's dataset in Rstudio:
February 10, 2024, LR

__Purpose:__ To compare the diversity metrics between Farm and No-farm individuals to see if there are any noticable differences between the two metadata categories. In order to be able to create a model, the two sets of data have to be significantly different, and diversity metrics are a good way to estimate this.

__Procedure:__
- Separate the Farm and No-farm individuals into separate phyloseq objects
- Create alpha diversity plots;
	- Use the plot_richness() function with the "Observed" and "Shannon" measures.
	- Also create a Faith's PD alpha diversity plot using the pd() function, and an unrooted tree
	
- Create beta diversity plots;
	- Use the distance() and ordinate() functions to generate beta diversity PCoA plots of the following metrics:
 		- Bray Curtis, Unifrac, Weighted Unifrac and Jaccard metrics
- R packages used: tidyverse, phyloseq, vegan, picante
- Calculate significance using the mann-whitney U test for Alpha diversity, and the PERMANOVA test for the Beta diversity

__Results:__
The alpha and beta diversity plots can be seen below.
None of the diversity metrics show significant difference in diversity when comparing the Farm and No-farm individuals

Alpha Diversity Plots for the Farm vs No-farm individuals:
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)/Diversity_plots/Alpha_Diversity_Plots.jpg">

Beta Diversity Plots for the Farm vs No-farm individuals:
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)/Diversity_plots/Combined_Beta_Diversity_Plots_with_Legend.png">


### Core Microbiome Analysis of the Farm dataset (healhty idividuals from the Parkinson's dataset)
February 28, 2024, EG:

__Purpose:__ Having selected all the healthy individuals in the Parkinson's dataset, we can now proceed to create a list of species characteristic of the mcirobiota of individuals that have been exposed to an agricultural environment (in our case, a farm). One of the analysis we decided to perform is the core microbiome analysis, which allows to determine the shared and unique species between selected groups of samples. The analysis will be performed in R using the microbiome package. <br> 
A high number of shared species will indicate that the individuals in the exposed and unexposed groups have a similar gut microbiome composition, whereas a diffirence between the 2 groups would be indicated by a low number of shared species. The group with a higer \u03B1-diversity will have a higher number of unique species compared to the other group. This will enable us to qualitatively assess the extent of similarity and divergence between the 2 groups. <br> Downstream, we will perform indicator species analysis, determine which ASVs are common between the indicator species and the species present uniquely in exposed individuals. This will comprise the basis of the predictive taxonomic model.

__Procedure:__
- load the phyloseq object containing the Parkinson's dataset filtered for only healhty patients
- convert the frequencies in the OTU table
- create vectors for all ASVs found in the samples in the Farm group and the No-farm group
- core microbiome analysis of either vector
	- relative abundance threshold: 1% (only take the ASV if its relative abundance is higher than 1%; this filters out all the rare ASVs and retains all abundant and non-abundant ASVs)
	- prevalence: 50% (only accept the ASV if it is found in over 50% of the samples)
- subset a new table for both exposed and unexposed groups by isolating the ASVs found in the core microbiome analysis
	- create 1 table per group
- generate a barplot of ASVs identified in either group
- generate a Venn diagram representing the unique and shared ASVs between groups

__Results:__
- 3 core ASVs have been identified in the unexposed group, while 4 were identified in the exposed group
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
- _Bacteroides_ genus is commonly found in colon and has both positive and negative effects on overall health depending on their interactions with other microbes in the gut microbiome ([Zafar, 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7872030/)).

- the bar plot of the distribution of the identified ASVS
	- overall, the 2 ASVs show a higher relative abundance in the unexposed group
	- there is higher variability of the _Bacteroides_ genus in both the exposed and unexposed groups compared to the _Faecalibacterium_
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)//core/core_ASVs_barplot(r.ab_0.01%2Cprev_0.5).png">

- Venn diagram generated for the core microbiome analysis
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)/core/venn_core(r.ab_0.01%2Cprev_0.5).png">

continued below: [Core Microbiome Analysis, continued](#Core-Microbiome-Analysis,-continued)


### Indicator Species Analysis of the Farm dataset (healhty idividuals from the Parkinson's dataset)
February 29, 2024, EG <br>

__Purpose:__ The indicator species analysis will allow to identify any taxa that are indicators of the Farm and No-farm group. This is the analysis that is performed regularly in the field to identify taxa that are commonly found in a specific environment. In our method, the results for the group of interst are checked against the core microbiome-identified taxa unique to the Farm group. 

__Protocol:__
- packages required: tidyverse, phyloseq, indicspecies
- load the non-rarefied phyloseq object for the Farm dataset
- glom the tax table to Genus, make sure to include the NA rows
- conver the frequences for the glommed Genus by transform_sample_counts()
- use multipatt() to calculate indicator values for the dataset and cluster by the `Lived_on_farm` column
- combine the selected ASVs with the data in the tax table
	- assign the name "ASV" to the first column in the tax table
	- filter for the rows with a p-value < 0.05
 - used the phyloseq, tidyverse, and indicspecies packages in R

__Results:__
The following list of ASVs was generated in the ISA:

|ASV code                        |Domain  |Phylum           |Class          |Order             |Family                     |Genus                   |Species               |Lived on farm?|Indicator value (IV)|_p-value_|
|--------------------------------|--------|-----------------|---------------|------------------|---------------------------|------------------------|----------------------|--------------|--------------------|---------|
|2fbf02f0d8043728d6e93f0bfa432d85|Bacteria|Actinobacteriota |Coriobacteriia |Coriobacteriales  |_Eggerthellaceae_          |_NA_                    |_NA_                  |Yes           |0.209               |0.035    |
|23667231b1722eb51bc7ccfaa5e5c38a|Bacteria|Bacteroidota     |Bacteroidia    |Bacteroidales     |_Prevotellaceae_           |_Prevotellaceae_UCG-001_|_uncultured bacterium_|Yes           |0.233               |0.010    |
|898f2833fcd964a4cfae006e66bf10b6|Bacteria|Bacteroidota     |Bacteroidia    |Bacteroidales     |_Barnesiellaceae_          |_NA_                    |_NA_                  |Yes           |0.184               |0.040    |
|455cd7d5df48c77aedb33f14590d1e78|Bacteria|Firmicutes       |Bacilli        |Erysipelotrichales|_NA_                       |_NA_                    |_NA_                  |Yes           |0.228               |0.035    |
|cbd66fa79a3e7fa1c33d29ed908de048|Bacteria|Firmicutes       |Bacilli        |Erysipelotrichales|_Erysipelatoclostridiaceae_|_Catenibacterium_       |_NA_                  |Yes           |0.329               |0.045    |
|b7fb7b8de1b4e013e7255cec428893a7|Bacteria|Verrucomicrobiota|Lentisphaeria  |Victivallales     |_Victivallaceae_           |_Victivallaceae_        |_uncultured rumen_    |Yes           |0.325               |0.015    |

- According to the literature, these ASVs are present abudantly in healthy human microbiome
- None of the ASVs were identified to the Species level, but 3 Genera were identified: _Prevotllaceae_, _Catenibacterium_, and _Victivallaceae_
- Indicator values indicate how strongly an ASV is associated with a group (here, with the farm group)
		- Calculation of IV = 100 x relative abundance x relative frequency
				- relative abundance: how many individuals in the group have the ASV
				- relative frequency: **proportion of sites in the group that have the individual**
		- values are only considered meaningful when > 0.7. The highest IV in this analysis is 0.329 for the ASV in the _Catenibacterium_ genus

continued below: [Indicator Species Analysis, continued](#Indicator-species-analysis.+continued)


### Core Microbiome Analysis, continued
March 3, 2024, EG:

__Purpose:__ Changing parameters to increase the number of identified core microbiome species. The previous analysis (see results above) did not yield enough species for to proceed as only 2 species were unique to the Farm group. Two taxa are not enough for a strong predictive model. To increase the number of species in the core mircobiome, the detection and prevalence criteria were made less stringent.

__Procedure:__ Perform the same analysis as above, but change the detection parameter to 0, and the prevalence to 0.1.

__Results:__
- there are 311 species in the Farm group's core microbiome, and 280 species in the unexposed group's core microbiome 
- bar graph plot generated:
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)/core/core_ASVs_barplot(r.ab_0%2Cprev_0.1).png">
- there are more ASVs included in the current analysis because the parameters are less stringent

- and the Venn diagram
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)/core/venn_core(r.ab_0%2Cprev_0.1).png">
- overall, there are more ASVs included into the analysis with the less stringent parameters, consistent with the bar plot generated
- the farm group has more unique species compared to the non-farm group, meaning that 


### Indicator Species Analysis, continued
February 29, EG

__Purpose:__ In the previous attempt at the ISA, only 6 species had ben identified as indicators of the Farm group; these are too few for downstream analysis. Additionally, the IVs are not very high. While making the parameters less stringent will not increase the IV for the identified indicators, but at least the number of indicator taxa will be increased. This will improve our chances of finding taxa that are unique to the Farm group (identified through the core microbiome analysis above) and are indicators of the Farm group.

__Protocol:__
- perform the same analysis, but glom by species
- subset ASVs with IV>0.2 and p<0.05
- blast the species selected using the NIH BLAST&reg; tool

__Results:__
- 19 indicator species were identified
- 12 ASVs were subsetted given the parameters above
- see this [file](https://github.com/lucasronny/MICB475-Group-11/blob/main/images%20and%20files/parkinsons/isa/isa_frame_high1.csv) for tables listing the taxonomic classification of ASVs performed using the representative sequences in QIIME2 and NIH BLAST &reg; tool classification


### Core Microbiome and Indicator Species Analysis: Combining Results to Identify Species for the Predictive Model
March 4, EG

__Purpose:__ The core microbiome and indicator species analyses have been performed to the necessary level of stringency. To generate the predictive taxonomic model, we will combine the species identified. This will mitigate the reduced stringency of either analysis to produce a stronger model.

__Procedure:__
- tables containing the core microbiome species and indicator species for the Farm group were joined using the inner_join() function, basing the join on the ASV column in either table (see code here)
- the same operations were performed for the indicator species and the core microbiome members unique to the exposed group

__Results:__
- there were 3 species in common betwen the total ASVs belonging to the exposed group and the indicator species:

|ASV code                            |Domain  |Phylum           |Class           |Order             |Family            |Genus                         |Species                    |Lived on farm?|Indicator value (IV)|_p-value_|
|------------------------------------|--------|-----------------|----------------|------------------|------------------|------------------------------|---------------------------|--------------|--------------------|---------|
|bfab02a86c5187ed451db10d9b81b8d5    |Bacteria|Firmicutes       |Clostridia      |Lachnospirales    |_Lachnospiraceae_ |_Eubacterium ventriosum group_|_uncultured bacterium_     |Yes           |0.420               |0.005    |
|**cbd66fa79a3e7fa1c33d29ed908de048**|Bacteria|Bacillota        |Erysipelotrichia|Erysipelotrichales|_Coprobacillaceae_|_Catenibacterium_             |_Catenibacterium misuokai_ |Yes           |0.327               |0.015    |
|**3f9fcc47fb363e0fdfa75dc56bb107f5**|Bacteria|Bacteroidota     |Bacteroidia     |Bacteroidales     |_Tannerellaceae_  |_Parabacteroides_             |_Parabacteroides johnsonii_|Yes           |0.303               |0.030    |

- all the sequences were blasted using the NIH BLAST&reg; tool
- all ASVs were identified down to the species level, except for the species belonging to the _Lachnospiraceae_ family
- additionally the uncultured bacterium ASV and _Parabacteroides johnsonii_ were in common between the unique core microbiome species in the exposed group and the indicator species
- when comparing the species from core microbiome analysis unique to and the indicator species for the farm individuals, only 2 ASVs are identified: _Catenibacterium misuokai_ and _Parabacteroides johnsonii_ (**the 2 bolded in the table**) 

__Conclusion:__
- the common species are very few
- additionally, the IV are low (IV<0.42), indicating that the species identified are not unique to the farm group
	- note that all the ubiquotous species (i.e., common to both groups) are excluded by the _p-value_ threshold
- only the species for the __2 bolded ASVs will be used as the species for the model__
	- this makes for a rather weak model, considering there are only 2 species and the indicator values are quite low


### Ranking Cities by Agricultural Exposure to Predict Which Will Be Identified as "Farming Populations" by The Predictive Model
March 13, ANM, LW

<br> __Purpose:__ To build a hypothesis to use in testing our model. Before looking for the model species in the Colombian population (grouped by city), we propose which cities are more or less likely to be identified as agricultural communities by our model. We will calculate the proportion of all workers tht are involved in agricultural activities. The highest rankings will be given to the cities which we believe to be likely agricultural communities, the lower ranking will be given to the species that are the least likely.

__Procedure:__
- Data taken from the [Censo Nacional Agropecuario 2014](https://www.dane.gov.co/files/images/foros/foro-de-entrega-de-resultados-y-cierre-3-censo-nacional-agropecuario/CNATomo3-Mapas.pdf)
- Cities in the Colombia dataset were classified based on their associated Department.
- Using the Census we looked and the data for Land Cover (the amount of land in hectares that is used for agricultural activities in a specific Department). Then we looked and the tables covering the amount of agricultural workers for each Department.
- The total amount of hectares in land cover and the total amount of workers respectively, were used to calculate the proportion of Agricultural land and workers for each Department.
- City ranking made by looking at the proportion results.

__Results:__
| City        | Associated Department| Land Cover| Proportion of Agricultural Land| Total Workers| Proportion of workers|
|-------------|----------------------|-----------|--------------------------------|--------------|----------------------|
| Bogota      | DC and Cundinamarca  | 1,482,641 | 63.07%                         | 285,275      | 2.70%                |
| Barranquilla| Atlantico            | 192,664   | 65.70%                         | 18,905       | 0.75%                |
| Medellin    | Antioquia            | 3,426,582 | 54.75%                         | 247,384      | 4.14%                |
| Cali        | Valle del Cauca      | 1,100,046 | 54.11%                         | 88,327       | 2.33%                |
| Bucaramanga | Santander            | 1,946,487 | 64.70%                         | 170,329      | 8.48%                |

Bogota data collected by adding census results from both Bogota D.C. and Cundinamarca since officially the city of Bogota is considered to be part of Cundinamarca but was treated separately in the census. 

__Rankings (most to least agricultural):__

1. Bucaramanga
2. Medellin
3. Bogota
4. Cali
5. Barranquilla


### Applying The Model to The Colombia Dataset
March 13, KM

<br> __Purpose:__ To look for the model species in each city in the Colombia dataset to determine whether or if these species exist in another dataset.

__Procedure:__
- Loaded required phyloseq packages in R (tidyverse, phyloseq, ape)
- Used previous code to build phyloseq object using the Colombia dataset
- Created a loop that would go through a list of target ASVs to determine if ASVs were found in the colombia dataset
- Loop determined if ASV was present then added to a counter to print out how many times the ASV was found
- If ASV was not present, then it would print "ASV not found in dataset"

__Results:__
- No ASVs of our defined species were found in the Colombia dataset
- New loop went through all exisiting ASVs in each dataset and compared to find no in common ASVs between the two datasets
- Misunderstanding that ASVs are unique to datasets made all results above invalid


### Filtering Phyloseq Objects to Find Specific Species
March 28, KM

<br> __Purpose:__ To re-work the initial purpose of the for loop to determine how many species of interests exist in our validation/test dataset

__Procedure:__
- Loaded required phyloseq packages in R
- Used previous code to build phyloseq object of our test dataset
- Wrote a vector including the specific names/filter depth of the species of interest 
	- eg: c("o__Flavobacteriales", "f__Flavobacteriaceae", "g__uncultured")
- Created a new filtered phyloseq object with subset_taxa of the name/filter of species of interest until the genus level
- Print and save new taxa table with filtered species as CSV
- R packages needed: tidyverse, phyloseq

__Results:__
- Code was able to successfully filter and find species of interests selected in the vector
- [R Script here](https://github.com/lucasronny/MICB475-Group-11/blob/main/Code/Analysis%20in%20R/Farm%20Dataset%20(Validation%202)/New_Farm_Phylo_Filter.R)


### New Validation Dataset Processing in QIIME2 - Farm Dataset:
March 23, EG <br>

__Purpose:__ As the original validation dataset did not prove to be useful (see Applying The Model to The Colombia Dataset), we decided to attempt validating the model using another dataset collected to study the differences between the inhabitants of agricultural and urban communities. The new dataset is one that compares stool samples between American industrialized and Himalayan agricultural populations. 

__Procedure:__
- Demultiplex the samples using the manifest file and single-end sequence reads (.qza file)
- Determine trimming parameters using the quality scores per base pair plot (see below),
- Trim and denoise the samples to a read depth of n=150 using the DADA2 algorithm (used Q = 32, all scores above the threshold)
- Generate the table.qzv and rep-seqs.qzv files and transfer to the local server for visualization
- Generate taxonomy.qza file for taxonomic classification of ASVs clustered in previous steps
- Filter out mitochondrial and chloroplast DNA
- Generate the unrooted tree using the rep-seqs.qza file
- Export all the necessary files for a phyloseq object
- all the files have been uploaded to farm_export folder

__Results:__
- the files generated are stored in the individual server and on a local computer
- the quality score for every base pair in the sequence is as follows:
<br> < img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images%20and%20files/Farm%20Dataset%20(Validation%202)/demux_quality_scores.png">

### New Validation Dataset Phyloseq Object
March 24th, AN  <br>
__Purpose:__  To take the QIIME2 data processing files for the new validation dataset and convert them into a phyloseq object to find our indicator species.  <br>
__Protocol:__ 
- Loaded necessary packages: phyloseq, ape, tidyverse and vegan.
- Loaded metadata, feature_table, taxonomy and phylogenetic tree produced from Qiime processing of the New Farm dataset.
- Adjusted metadata, table and taxonomy to phyloseq format.

__Results__ 
- Code succesfully generated phyloseq object.
- The file generated was exported to GitHub to continue with the validation of the model.

### EDIT: Core microbiome analysis, ISA, and reconciliation
April 4th, EG <br>

__Purpose:__ While the analyses above were performed correctly (see core microbiome, ISA, and reconciliation sections above), the team missed the fact that the phyloseq object used for those analysis was not propoerly filtered. The object included samples from Parkinson's disease patients, whereas we aimed to include samples only from healhty individuals who lived on farms. The effect of Parkinson's disease would be a third variable which we do not want to account for. Consequently, we need to control for the overall health state of the included samples by including only healthy subjects.

__Protocol:__
- in the code for ISA and core microbiome analysis, include a function to filter the phyloseq object to make it contain only samples with "Control" output in the Disease metadata column
- repeat the analysis using previously written code (see sections on ISA and core microbiome analysis)
- packages needed: tidyverse and phyloseq (for both core and ISA), indicspecies (ISA), microbiome (core)

__Results:__
- the new outputs generated in the analyses
	Core Microbiome Analysis done on the newly filtered phyloseq object:
	<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Parkinsons Dataset (Model)/core/venn_core(r.ab_0%2Cprev_0.1)_FIXED.png">
	File path in the repository: images and files/parkinsons/core/venn_core(r.ab_0,prev_0.1)_FIXED.png_

	Table of Indicator Taxa(generated by the ISA on the newly filtered phyloseq object), not annotated using NCBI BLAST:
|ASV code                        |Domain  |Phylum          |Class           |Order             |Family                                 |Genus                                  |Species                            |Lived on farm?|Indicator value (IV)|_p-value_|
|--------------------------------|--------|----------------|----------------|------------------|---------------------------------------|---------------------------------------|-----------------------------------|--------------|--------------------|---------|
|54da3e844b42d39369127a368dcced1e|Bacteria|Firmicutes      |Bacilli         |RF39              |_RF39_                                 |_RF39_                                 |_human gut_                        |Yes           |0.380               |0.040    |
|fdbdbe81dd5e0e2b71693fe87db9691d|Bacteria|Firmicutes      |Bacilli         |Erysipelotrichales|_Erysipelotrichaceae_                  |_Holdemanella_                         |_uncultured bacterium_             |Yes           |0.409               |0.045    |
|cbd66fa79a3e7fa1c33d29ed908de048|Bacteria|Firmicutes      |Bacilli         |Erysipelotrichales|_ Erysipelatoclostridiaceae_           |_Catenibacterium_                      |_NA_                               |Yes           |0.398               |0.025    |
|80dbd6cdd5798902a9038eb121f7f640|Bacteria|Firmicutes      |Clostridia      |Clostridia UCG-014|_ Clostridia UCG-014_                  |_Clostridia UCG-014_                   |_uncultured bacterium_             |No            |0.501               |0.030    |
|ce7e172e3d64bc42dddb2f4aa26e3a03|Bacteria|Desulfobacterota|Desulfovibrionia|Desulfovibrionales|_Desulfovibrionaceae_                  |_Desulfovibrio_                        |_NA_                               |Yes           |0.553               |0.020    |
|9b695df7c620cdc44eb078b53cef23de|Bacteria|Firmicutes      |Clostridia      |Oscillospirales   |_[Eubacterium] coprostanoligenes group_|_[Eubacterium] coprostanoligenes group_|_Eubacterium coprostanoligenes_    |Yes           |0.277               |0.040    |
|8b7dd3e4e3f153a95c528f452b48a2bd|Bacteria|Bacteroidota    |Bacteroidia     |Bacteroidales     |_Rikenellaceae_                        |_Rikenellaceae RC9 gut group_          |_uncultured bacterium_             |Yes           |0.387               |0.045    |
|dc04cb67b3295173275e20a8a7a3d16a|Bacteria|Bacteroidota    |Bacteroidia     |Flavobacteriales  |_Flavobacteriaceae_                    |_uncultured_                           |_gut metagenome_                   |Yes           |0.408               |0.040    |
	File path to the .csv file containing this table: images and files/parkinsons/isa/isa_frame_high_FIXED.csv

	Table of Model Taxa generated by merging the core microbiome and ISA performed on the filtered dataset (annotated using NCBI BLAST and ordered by decreasing IV):
|ASV code                        |Domain  |Phylum          |Class           |Order             |Family                                 |Genus                                  |Species                            |Lived on farm?|Indicator value (IV)|_p-value_|
|--------------------------------|--------|----------------|----------------|------------------|---------------------------------------|---------------------------------------|-----------------------------------|--------------|--------------------|---------|
|ce7e172e3d64bc42dddb2f4aa26e3a03|Bacteria|Desulfobacterota|Desulfovibrionia|Desulfovibrionales|_Desulfovibrionaceae_                  |_Desulfovibrio_                        |_NA_                               |Yes           |0.553               |0.020    |
|fdbdbe81dd5e0e2b71693fe87db9691d|Bacteria|Firmicutes      |Bacilli         |Erysipelotrichales|_Erysipelotrichaceae_                  |_Holdemanella_                         |_uncultured bacterium_             |Yes           |0.409               |0.045    |
|dc04cb67b3295173275e20a8a7a3d16a|Bacteria|Bacteroidota    |Bacteroidia     |Flavobacteriales  |_Flavobacteriaceae_                    |_uncultured_                           |_gut metagenome_                   |Yes           |0.408               |0.040    |
|cbd66fa79a3e7fa1c33d29ed908de048|Bacteria|Firmicutes      |Bacilli         |Erysipelotrichales|_ Erysipelatoclostridiaceae_           |_Catenibacterium_                      |_NA_                               |Yes           |0.398               |0.025    |
	Path to the .csv file containing this table: images and files/parkinsons/core_isa/isa_core_farm_FIXED_annotated.csv

__Conclusion:__ <br>
Filtering the phyloseq object to contain only the healthy individuals:
- identified a higher number of species unique to the Farm group (core microbiome analysis)
- identified a higher number of indicator taxa (ISA)
- allowed to include more species with higher IVs as model taxa

Compared to previous results, these improvements better demonstrate that our method can generate model taxa that are good predictors of exposure to a certain environment.


### Results from model validation of indicator ASVs on the Farm dataset:
April 4, LR <br>

__Purpose:__ The list of indicator ASVs was applied to the new Farm dataset to validate if the indicator ASVs are actually indicators for agricultural exposure. This step is thus done to assess the abundance and prevalence of each indicator ASV in the Farm dataset to see how valid the indicator ASVs are.

__Procedure:__
- Multiple ASVs were corresponding to the same taxa, thus ASVs for the same taxa were kept separate for preventing overrepresentation of each taxa
- Remove non-farmers (Americans) from the dataset, only keeping Nepalese farmers
- Scan through Farm dataset for the presence and abundance of each indicator ASV
- A bubble plot was generated, displaying total ASV abundance as the color of each bubble, and the percentage prevalence is shown as the bubble size.
- the packages needed for this analysis: tidyverse, phyloseq, ape

__Results:__
- Some ASVs were only present in one sample, and were thus removed from the figure. This can be seen in comparing the two figures below.
- _Catenibacterium mitsuokai_ ASV 1 and _Holdemanellea biformis_ ASV 2 were present in >70% of samples.
- _Holdemanellea biformis_ ASV 1 was only present in 35% of samples.
- _Desulfovibrio piger_ was present in 48% of samples.

Bubble plot including ASVs which appeared only once:
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Farm Dataset (Validation 2)/Bubble_plot_sample_counts.png">

Bubble plot excluding ASVs which appeared only once:
<br><img src = "https://github.com/lucasronny/MICB475-Group-11/blob/main/Images and files/Farm Dataset (Validation 2)/Bubble_plot_sample_counts_no_1ss.png">
