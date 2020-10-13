# combining fwh data from 2018

### load libraries ######################
library(tidyverse)
library(data.table)
library(seqinr)

### fwh primer ##########################

# sequence data for libraries 10, 13, 16, 19, 22, 25, 28, 34, 36, 38, 40, 42, 44, 46 (library 4 not included)
lulified_nochim <-
  readRDS("H:/Documents/Insektmobilen/Analysis/DADA2/fwh/all_2018/data/lulified_nochim.RDS") # read in the lulufied RDS file for library 10, 13 and 16
otus <- lulified_nochim[["curated_table"]] # extract the otutable
otus <- rownames_to_column(otus, var = "otuid")
asvtable <- otus[,!grepl("IM17",colnames(otus))] # remove samples from 2017 - German samples still in there, but will be removed when merging with metadata eventually

# the fasta file has been subset to only contain sequences with a length over 200 bp, so the asv table should only have asvs with the correct read length
fastas <- read.fasta("H:/Documents/Insektmobilen/Analysis/DADA2/fwh/all_2018/data/lulufied_otus_lengthcorrected.fasta")
keep <- names(fastas)
test <- asvtable %>% column_to_rownames(var = "otuid")
asvtable <- test[(rownames(test) %in% keep), ]

# save only 2018 data - here blanks will be removed and German data as well
#test <- asvtable %>% column_to_rownames(var = "otuid")
asvs <- asvtable[, grepl("IM18_*", names(asvtable))] # only keep Danish samples from 2018 - be aware this removes blanks and negatives as well
names(asvs)

# get metadata
data <- read.delim("cleaned-data/DK_metadata_2018.txt",sep=" ")

data <- data %>% filter(PCRID %in% colnames(asvs)) # only retain metadata for the samples that are in the sequence table
keep <- data$PCRID
asvs <- asvs[, (names(asvs) %in% keep)] # subset metadata to match sequence data

# is there taxa in all samples?
colSums(asvs)
min(colSums(asvs))

# is some taxa not present in some samples?
rowSums(asvs)
min(rowSums(asvs)) # some taxa are not the retained samples and should be removed

asvs <- asvs[rowSums(asvs) > 0, ]

colSums(asvs)
min(colSums(asvs))

# taxonomy data 
taxonomy_1 <- read.delim("H:/Documents/Insektmobilen/Analysis/DADA2/fwh/all_2018/data/blastresult_1.csv", sep = ",")
taxonomy_2 <- read.delim("H:/Documents/Insektmobilen/Analysis/DADA2/fwh/all_2018/data/blastresult_2.csv", sep = ",")
taxonomy_3 <- read.delim("H:/Documents/Insektmobilen/Analysis/DADA2/fwh/all_2018/data/blastresult_3.csv", sep = ",")

# merge taxonomy data
taxonomy <- rbind(taxonomy_1, taxonomy_2)
taxonomy <- rbind(taxonomy, taxonomy_3)

# match asv ids to taxonomy
keep <- rownames(asvs)
taxonomy <- taxonomy[taxonomy$occurrenceId %in% keep, ]

taxa <- taxonomy
taxonomy <- taxonomy %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package

# choose only OTUs that have 99 % or higher ID match and choose only the OTUs that are assigned to class Insecta. To make sure the otuids are not deleted, we nedd to make the rownames into a column and then revert back to rownames
taxonomy_class <-
  taxonomy %>% filter(class == "Insecta") %>% column_to_rownames('occurrenceId') 

taxonomy_99 <-
  taxonomy %>% filter(class == "Insecta") %>% filter(identity >= 99) %>% column_to_rownames('occurrenceId')

# save outputs
write.table(asvs,file="cleaned-data/DK_asvtable_2018_data.txt",sep="\t")
write.table(data,file="cleaned-data/DK_metadata_2018_sequenced.txt",sep="\t")
write.table(taxonomy_class,file="cleaned-data/DK_taxonomy_Insecta.txt",sep="\t")
write.table(taxonomy_99,file="cleaned-data/DK_taxonomy_Insecta_99.txt",sep="\t")

# examining taxonomy data
#table(taxonomy$order)
#lepidoptera <- taxonomy %>% filter(order == "Lepidoptera")
#table(lepidoptera$family)


