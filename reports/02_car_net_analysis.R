# car net analysis

### load libraries ################
library(tidyverse)
library(readr)

### load data #####################
asvs <- read.delim("cleaned-data/asvtable.txt",sep=" ")
data <- read.delim("cleaned-data/data.txt",sep=" ")
taxonomy <- read.delim("cleaned-data/taxonomytable.txt",sep=" ")

### aligning data for analysis ###################
names(asvs) = gsub(pattern = "X*", replacement = "", x = names(asvs)) # remove the Xs that have emerged in the column headers that are numerical
asvs <- asvs[, grepl("_QS$", names(asvs))] # only keep qiasymphony (automatic) purified samples
names(asvs) <-  gsub(pattern = "_QS*", replacement = "", x = names(asvs)) # remove purification method from the data
names(asvs)

keep <- data$PCRID
asvs <- asvs[, (names(asvs) %in% keep)] # subset asvtable to only contain samples with metadata

# 35 samples were in the sequencing run. Data needs to be updated
keep <- names(asvs)
data <- data %>% filter(PCRID %in% keep)

# do all samples have seuqnces?
colSums(asvs)
asvs <- asvs[,colSums(asvs) > 0]

# remove empty asvs
rowSums(asvs)
asvs <- asvs[rowSums(asvs) > 0, ]

# remove the asvs from the taxonomy
keep <- rownames(asvs)
taxonomy <- taxonomy[rownames(taxonomy) %in% keep, ]

# remove the two samples with no sequences
keep <- names(asvs)
data <- data %>% filter(PCRID %in% keep)

# remove the negative with low reads
table(data$Sample)
data <- data %>% filter(!PCRID == "IM17_47")
keep <- data$PCRID
asvs <- asvs[, (names(asvs) %in% keep)]

### preliminary data exploration #####
table(data$Date)
