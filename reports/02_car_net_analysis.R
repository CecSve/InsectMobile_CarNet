# car net analysis

### load libraries ################
library(tidyverse)
library(readr)

### load data #####################
asvs <- read.delim("cleaned-data/asvtable.txt",sep=" ")
data <- read.delim("cleaned-data/data.txt",sep=" ")
taxonomy <- read.delim("cleaned-data/taxonomytable.txt",sep=" ")

### subsetting the asv table ###################
names(asvs) = gsub(pattern = "X*", replacement = "", x = names(asvs)) # remove the Xs that have emerged in the column headers that are numerical
asvs <- asvs[, grepl("_QS$", names(asvs))] # only keep qiasymphony (automatic) purified samples
names(asvs) <-  gsub(pattern = "_QS*", replacement = "", x = names(asvs)) # remove purification method from the data
names(asvs)

keep <- data$PCRID
asvs <- asvs[, (names(asvs) %in% keep)] # subset asvtable to only contain samples with metadata
