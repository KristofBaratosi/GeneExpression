# GeneExpression
## Deseq_workflow: 
Differential Gene Expression  and Principal Component Analysis were performed on a dataset with primary human airway smooth muscle cell lines. <br>
There are 4 cell lines, every cell line has an untreated and a treated(with dexamethasone) sample. <br>
The following plots were made:
- MA Plot
- PCA Plot (PC1 and PC2)
- Volcano Plot
- Heatmap
### Data availability: 
The data was obtained from: https://github.com/kpatel427/YouTubeTutorials/tree/main/files
## Machine_learning/bowel_disease: 
### Data availability: 
The data was obtained from: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11349 <br>
Experimental factors: normal, ulcerative colitis, Crohn's disease <br>
The labels for the machine learning algorithms were encoded in the following way: <br>
- normal was encoded with 0 
-  ulcerative colitis and Crohn's disease were encoded with 1.

### Files:
```data_prep.R```: This code was used to prepare the sample information and the counts matrix for machine learning. There is no need to run this code, as its outputs are already saved: (counts_matrix.csv and sample_information.csv) <br>
```machine_learning.R```: Genes that don't have at least 75 reads in total are removed, the raw counts are normalized first with the DESeq2 package, then quantile normalization is performed. The data is split into a train(80%) and test(20%) set. PCA is done on the training set, then the test set is transformed to PC's. Then, the  first 35 principal components are selected for machine learning. <br>
Four machine learning algorithms are tested:
- Random Forest Classsifier
- Logistic Regression
- Support Vector Machine
- Naive Bayes classifier
 <br>
Then, a roc plot and confusion matrices are created to assess the performance of the models.
