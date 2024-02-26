# MICB475-Group-11: Project 2
A repository for documentation, coding scripts for our project.<br>

### Upcoming Meeting Agenda (Feb/29/2024)
- Item 1
- Item 2
- Item 3

For meeting minutes, see [Meeting notes](https://github.com/lucasronny/MICB475-Group-11/tree/3647d9e5b9e84cc3018353200aefe5036ffafb39/notes).

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
- I recall in the last meeting that we might be pooling the microbiomes in the Colombia dataset by city (for when we test our model). Each city has a different number of samples, so should we rarefy the Colombia dataset too?
