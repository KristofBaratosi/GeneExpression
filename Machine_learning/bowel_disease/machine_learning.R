library(tidyverse)
library(dplyr)
library(caret)
library(pROC)
library(pROC)
library(e1071)
library(randomForest)
library(DESeq2)
library(preprocessCore)
#reading in the raw counts and the sample information files
colData <-  read.csv("sample_information.csv",row.names = 1)
counts_data <-  read.csv("counts_matrix.csv",row.names = 1)
#converting the condition and outcome columns to factors 
colData[c( 'condition', 'outcome')] <- lapply(colData[c( 'condition', 'outcome')], factor)
#checking whether the order of the samples in colData is identical to counts_data
all(colnames(counts_data) %in% colData$sampleID)
all(colnames(counts_data) ==colData$sampleID)
#Deseq dataset object 
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~condition)
#Removing genes that don't have at least 75 reads total 
keep <- rowSums(counts(dds))>=75
dds <- dds[keep,]
dim(dds)
# Extracting normalized counts
dds <- estimateSizeFactors(dds)
df <- counts(dds, normalized = TRUE)
#quantile normalization
df_norm <- normalize.quantiles(df)
#restoring row and column names
rownames(df_norm) <- rownames(df)
colnames(df_norm) <- colnames(df)
#Transposing the df so we have the samples as rows, and genes as columns
ml_df <-  as.data.frame(t(df_norm))
dim(ml_df)
#assigning the target variable to the dataframe
ml_df$outcome <- colData$outcome
dim(ml_df)
#MACHINE LEARNING
#setting seed for reproducibility
set.seed(5) 
#Creating training and test sets- train: 80%, test: 20%
split_index <- sample(1:nrow(ml_df), 0.8 * nrow(ml_df))
train_data <- ml_df[split_index, ]
test_data <- ml_df[-split_index, ]

train_features <- train_data[, -ncol(train_data)]  
train_target <- train_data[, ncol(train_data)]   

test_features <- test_data[, -ncol(test_data)]
test_target <- test_data[, ncol(test_data)]
#PCA on train features, transforming  test_features to PC's
#(Using first 35 PC's for the models)
pca <- prcomp(train_features,scale. = TRUE)
train_features <- as.data.frame(pca$x[,1:35])
test_features <- predict(pca, newdata = test_features)
test_features <- test_features[,1:35]
test_features <- as.data.frame(test_features)
dim(train_features)
dim(test_features)
#Creating dataframe to plot the first and second principal components
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
#Extracting variation explained by PC's
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
#Scree plot
barplot(pca.var.per[1:15], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
#Plotting PC1 and PC2 
ggplot(data = pca.data, aes(x = X, y = Y, color = train_target)) +
  geom_point() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA Plot")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color = "Disease")
# RANDOM FOREST CLASSIFIER
rf_model <- randomForest(train_features, train_target, ntree = 100)
predictions <- predict(rf_model, test_features)
prediction_probs <- predict(rf_model, test_features,type = "prob")
rfc_accuracy <- sum(predictions == test_target) / length(test_target)
print(paste("Accuracy of Random Forest Classifier on  the test set:", round(rfc_accuracy * 100, 2), "%"))
#LOGISTIC REGRESSION
logit_model <- glm(train_target ~ ., data = train_features, family = "binomial")
test_predictions <- predict(logit_model, newdata = test_features, type = "response")
test_binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
logreg_accuracy <- sum(test_binary_predictions == test_target) / length(test_target)
print(paste("Accuracy of Logistic Regression on the test set:", round(logreg_accuracy * 100, 2), "%"))
#NAIVE BAYES
nb_model <- naiveBayes(train_features, train_target)
nb_predictions <- predict(nb_model, test_features)
nb_accuracy <- sum(nb_predictions == test_target) / length(test_target)
cat("Accuracy of Naive Bayes on the test set:", nb_accuracy, "\n")
#Support Vector Machine 
svm_model <- svm(train_target ~ ., data = train_features)
svm_predictions <- predict(svm_model, newdata = test_features)
svm_accuracy <- sum(svm_predictions == test_target) / length(test_target)
cat("Accuracy of SVM on the test set:", svm_accuracy, "\n")
svm_model <- svm(train_target ~ ., data = train_features, probability = TRUE)
svm_probabilities <- predict(svm_model, newdata = test_features, probability = TRUE)
# Extracting probabilities for the positive class 
svm_positive_probs <- attr(svm_probabilities, "probabilities")[, "1"]
# Creating a data frame that contains the accuracy of the models
accuracy_df <- data.frame(model = c("Random Forest Classifier", "Logistic Regression","Naive Bayes","Support Vector Machine"),
                          accuracy = c(rfc_accuracy,logreg_accuracy,nb_accuracy,svm_accuracy) )
#Creating a bar plot for model accuracy 
ggplot(accuracy_df, aes(x = reorder(model,-accuracy), y = accuracy, fill = model)) +
  geom_bar(stat = "identity", color = "black",width = 0.4) +
  labs(title = "Accuracy of the Models", x = "Model", y = "Accuracy") +
  ylim(0, 1)+
  theme_minimal()
#ROC PLOT 
#This prevents the axes from going beyond 1 
par(pty="s")
#Calculating roc and auc for logistic regression, random forest classifier and support vector machine
roc(test_target, test_predictions,plot= TRUE, legacy.axes = TRUE,
                  col= "#377eb8",lwd=4,print.auc=TRUE, percent=TRUE,xlab="False Positive Percentage", ylab="True Postive Percentage")
plot.roc(test_target,prediction_probs[,2],col="#4daf4a", lwd=4, percent=TRUE,print.auc=TRUE, add=TRUE, print.auc.y=42)
plot.roc(test_target,svm_positive_probs,col="#de2d26", lwd=4, percent=TRUE,print.auc=TRUE, add=TRUE, print.auc.y=35)
#Creating legend 
legend("bottomright", legend=c("Logisitic Regression", "Random Forest", "Support Vector Machine"), col=c("#377eb8", "#4daf4a","#de2d26"), lwd=3,cex = 0.7)
#Setting  graphics back to normal 
par(pty = "m")
#Creating confusion  matrices
svm_cm <- confusionMatrix(svm_predictions,test_target)
rfc_cm <- confusionMatrix(predictions,test_target)
logreg_cm <- confusionMatrix(as.factor(test_binary_predictions),test_target)
nb_cm <- confusionMatrix(nb_predictions,test_target)
#Viewing the matrices
svm_cm
rfc_cm
logreg_cm
nb_cm