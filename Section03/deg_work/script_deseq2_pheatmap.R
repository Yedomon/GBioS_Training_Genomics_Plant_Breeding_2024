### Set the working directory


setwd("C:/Users/ange_/Downloads/GBIOS/Section03")



### Librairies

library(DESeq2)
library(pheatmap)


### Data importation

data_count = read.table("C:/Users/ange_/Downloads/GBIOS/Section03/sesame_count.txt", header = TRUE, sep = "\t")



data_count = read.csv("C:/Users/ange_/Downloads/GBIOS/Section03/sesame_count.csv", header = TRUE, sep = ",", row.names = 1)


head(data_count)



### Rename the columns
colnames(data_count) <- c("ho_1", "ho_2", "ho_3", "lo_1", "lo_2", "lo_3")


### Verify the change
head(data_count)


### Load the sample information


sample_info = read.csv("coldata_table.csv", h=T, sep=",", row.names = 1)


# Making sure the row names in coldata matches the column names in data_count 

all(colnames(data_count)%in% rownames(sample_info)) 


# Double check if they are in the same order

all(colnames(data_count) == rownames(sample_info)) 


### Construct DESEQDataSet Object

dds <- DESeqDataSetFromMatrix(countData=data_count, 
                              colData=sample_info, 
                              design=~replication)
