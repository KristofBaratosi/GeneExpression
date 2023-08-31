library(tidyverse)
library(DESeq2)
library(ggplot2)
#reading in counts data 
counts_data <- read.csv("counts_data.csv")
#head(counts_data)
dim(counts_data)
#reading in sample information 
colData <- read.csv("sample_info.csv")
#checking whether the sample info rows and counts data columns match 
all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data) ==rownames(colData))
#Deseq dataset object 
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~dexamethasone)
dds
#filtering genes that dont't have at least 10 reads in total 
keep <- rowSums(counts(dds))>=10
keep
dds <- dds[keep,]
dim(dds)
#setting factor level, the reference will be the untreated and we will compare it against the treated samples
dds$dexamethasone<- relevel(dds$dexamethasone,ref="untreated")
dds$dexamethasone
#running deseq
dds <- DESeq(dds)
#exploring results
res <- results(dds)
res
#summary of res 
summary(res)
#changing p-value for the results 
res.001= results(dds, alpha = 0.01)
res.001
summary(res.001)
#MAPlot
plotMA(res.001)
library(EnhancedVolcano)
# Creating a volcano plot using EnhancedVolcano
evolcano <- EnhancedVolcano(res,
                            lab = rownames(res.001),
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            xlab = bquote(~Log[2]~ 'fold change'),
                            pCutoff = 10e-32,
                            FCcutoff = 2.0,
                            pointSize = 4.0,
                            labSize = 6.0,
                            colAlpha = 1,
                            legendPosition = 'right',
                            legendLabSize = 12,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 0.75,
                            subtitle = "")

# Printing volcano plot
print(evolcano)
# creating a dataframe for the normalized counts
df <- counts(dds, normalized=T)
head(df)
#PCA
pca <- prcomp(t(df),scale. = TRUE)
#scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
#Plot of PC1 and PC2
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])

ggplot(data = pca.data, aes(x = X, y = Y, color = colData$dexamethasone)) +
  geom_point() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA Plot",subtitle = "PCA1(x-axis) and PCA2(y-axis)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  labs(color = "Treatment")
#HEATMAP
#for converting gene ids
library(org.Hs.eg.db)
#extracting only significant genes from the res df object
sigs <- na.omit(res.001)
sigs<- sigs[sigs$padj < 0.01, ]
sigs
#converting sigs to df, extracting genes with basemean >100 and asbolute value of log2foldchange >2.3
sigs.df <- as.data.frame(sigs)
sigs.df <- sigs.df[(sigs.df$baseMean >100)&(abs(sigs.df$log2FoldChange)>2.3),]
#converting ensemble gene id's to symbols 
sigs.df$symbol <- mapIds(org.Hs.eg.db,keys = rownames(sigs.df), keytype = "ENSEMBL", column = "SYMBOL")
mat <- counts(dds, normalized=T)[rownames(sigs.df),]
dim(mat)
mat.z <- t(apply(mat,1,scale))
colnames(mat.z) <- colnames(df)
#creating the actual heatmap, (you might need to make the plot window bigger)
library(ComplexHeatmap)
Heatmap(mat.z, cluster_rows = T,cluster_columns = T, column_labels = colnames(mat.z), name = "Z-score",row_labels = sigs.df[rownames(mat.z),]$symbol)
