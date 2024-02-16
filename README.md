# MICB475-Group-11: Project 2

### Meeting notes (Feb 8):
- decided on parkinson & columbia datasets to make a model for predicting the relationship between farm land and microbiome indicators.
- for Thursday, February 15
-   See if there's enough samples with NO PARKINSON & YES FARM
-   Paired end vs single end?

### Processing of the Parkinson's dataset in QIIME2 (INSERT DATE HERE):

Results:
- 25 samples with YES-farm, NO-PD

To inquire in meeting: 
- Should I have filtered the mitochondria + chlorplast?
- Was 251 a good cutoff for filtering (i.e the last base pair, ~25 quality score)? 
- All the "NO DISEASE & FARM=YES" have 0 in constipation?
- Help with rarefaction curve - should we cut 1 of the samples off?
- Further steps?

### Meeting (Thursday, February 15):
- some things that we might want to do as intermediate steps in our project (unless other papers already did it in which case we can just cite their data): 
- check whether the beta diversity between farm and non-farm people is significantly different
- possibly check alpha diversity of the farm people
- indicator species analuysis for the farm people

To do for the next meeting (Thursday, Febraury 29):
- SEND TABLE.QZV FILES (before rarefaction)
- make phyloseq object for Colombia (don't need to rarefy)

- alpha (Faith's, Shannon's) and beta (Bray Curtis, Jaccard) diversity and include figures (can post figure on Github)

- organize Github (dates, progress made)
