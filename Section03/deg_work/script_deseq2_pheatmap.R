### Set the working directory


setwd("C:/Users/ange_/Downloads/GBIOS/Section03")



### Librairies

library(DESeq2)
library(pheatmap)


### Data importation

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
                              design=~oil_content)



# pre-filtering: removing rows with low gene counts by keeping rows that have at least 10 reads in total

keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]



## Now we’re ready to run DESEQ function

dds <- DESeq(dds)


# Take a look at the results table
res <- results(dds)

head(results(dds)) #let's look at the results table

# Explore results

summary(res)

res0.05 = results(dds, alpha = 0.05)
summary(res0.05)

# Contrast

resultsNames(dds)


results(dds, contrast = c ("oil_content", "low", "high"))

### Save your results

## Order by adjusted p-value

res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

## Write results
write.csv(resdata, file="diffexpr-results.csv",quote = FALSE,row.names = F)


### PCA and Heatmap

# Regularized log transformation for clustering/heatmaps, etc

rld <- rlogTransformation(dds)


rld


# Principal Components Analysis
plotPCA(rld, intgroup = "oil_content")


# Package

library(pheatmap) 


# Data: We will top 10 with highest log fold change for illustration

data_set = read.csv("diffexpr-results_for_heatmap.csv", h = T, sep = ",", row.names = 1)



# Make a matrix

data_matrix = as.matrix(data_set)


# Scale the data


data_matrix_scaled = scale(data_matrix)


# Render the heatmap

pheatmap(data_matrix,  
         scale = "row")


##### BONUS: Want to play with color ?

# Load required libraries
library(pheatmap)
library(RColorBrewer)

# Define your custom color palette
# For example, using a blue-to-red gradient with 50 colors
my_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Generate the heatmap with custom colors
pheatmap(data_matrix, 
         scale = "row", 
         color = my_colors)



