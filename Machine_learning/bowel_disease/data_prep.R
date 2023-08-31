#DATA PREPERATION
#Reading in the files
library(tidyverse)
library(dplyr)
raw_counts <- read.table(file="ArrayExpress-raw.csv",
                             sep=",",
                             header=T,
                             fill=T,
                             check.names=F)
dim(raw_counts)
sample_info <- read.table(file="E-MTAB-11349.sdrf.txt", sep="\t", header=T, fill=T, check.names=F)
dim(sample_info)

glimpse(sample_info)
#Selecting important columns
sample_info_selected <- select(sample_info,
                                   'Source Name',
                                   'Characteristics[age]',
                                   'Characteristics[sex]',
                                   'Characteristics[disease]'
)

sample_info_selected <- rename(sample_info_selected,
                                   'sampleID' = 'Source Name',
                                   'age' = 'Characteristics[age]',
                                   'sex' = 'Characteristics[sex]',
                                   'condition' = 'Characteristics[disease]'
)
#Renaming patients with chrons disease
sample_info_selected$condition[agrep("Crohns", sample_info_selected$condition)] <- "crohns_disease"
#encoding "normal" patients with 0 and every other patient with 1
sample_info_selected <- sample_info_selected %>%
  mutate(outcome = if_else(condition == 'normal', 0, 1))
#Changing spaces to _
sample_info_selected <- sample_info_selected %>%
  mutate(sampleID = gsub(" ", "_", sampleID),
                condition = gsub(" ", "_", condition))
sample_info_selected[c('sex', 'condition', 'outcome')] <- lapply(sample_info_selected[c('sex', 'condition', 'outcome')], factor)
glimpse(sample_info_selected)

#Cleaning up the counts matrix 
counts_matrix <- raw_counts[-1,-1]
rownames(counts_matrix) <- NULL
counts_matrix <-  counts_matrix %>% column_to_rownames('read')
colnames(counts_matrix) <- gsub(" ", "_", colnames(counts_matrix))
#all(sample_info_selected$sampleID==colnames(counts_matrix))

#Saving them as csv files 
write.csv(counts_matrix, file="/Users/baratosikristof/Desktop/Main/bioinfo/chrons_raw_files/counts_matrix.csv")
write.csv(sample_info_selected, file="/Users/baratosikristof/Desktop/Main/bioinfo/chrons_raw_files/sample_information.csv")

