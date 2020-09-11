# preparing data for analysis

### load libraries ################
library(tidyverse)
library(readr)
library(lubridate)

### load & subset metadata #################
pilot_metadata <-
  read_delim(
    "raw-data/Metadata_Pilotstudy2017.txt",
    "\t",
    escape_double = FALSE,
    locale = locale(encoding = "WINDOWS-1252"),
    na = "0",
    trim_ws = TRUE
  )

metadata_validation <-
  read_delim(
    "raw-data/Metadata_Validation2018.txt",
    "\t",
    escape_double = FALSE,
    locale = locale(encoding = "WINDOWS-1252"),
    trim_ws = TRUE
  )

metadata_validation <- metadata_validation %>% filter(Sample == "Car net") # only keep car net samples

### load sequence data ###########
lulified_fwh <- readRDS("raw-data/lulified_fwh_nochim.RDS")
asvs <- lulified_fwh[["curated_table"]]

### load & subset lab data ###############
IM2017_labdata <-
  read_delim(
    "raw-data/IM2017_labdata.txt",
    "\t",
    escape_double = FALSE,
    locale = locale(encoding = "WINDOWS-1252"),
    trim_ws = TRUE
  )

IMvalidation_labdata <- read_delim("raw-data/IMvalidation_labdata.txt", "\t", escape_double = FALSE, locale = locale(encoding = "WINDOWS-1252"), trim_ws = TRUE)

### load taxonomy ###############
taxonomy <- read_csv("raw-data/blastresult_GBIF_lulufied_fwh_novaseq.csv")

### filter taxonomy & ASVs ###################
# filter taxonomy based on 99% identity match or higher and only keep insects
taxonomy <- taxonomy %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package

taxonomy <- 
  taxonomy %>% filter(identity >= 99) %>% filter(class == "Insecta") %>% column_to_rownames('occurrenceId') # choose only OTUs that have 99 % or higher ID match and choose only the OTUs that are assigned to class Insecta. To make sure the otuids are not deleted, we nedd to make the rownames into a column and then revert back to rownames
write.table(taxonomy, file = "cleaned-data/taxonomytable.txt")

# subset asv table to fit taxonomy
keep <- rownames(taxonomy)
asvs <- asvs %>% rownames_to_column(var = "otuids")
asvs <- asvs %>% filter(otuids %in% keep) %>% column_to_rownames(var = "otuids")
write.table(asvs, file = "cleaned-data/asvtable.txt")

### merge metadata & lab data ##########
# lab data
names(IM2017_labdata)
labdata_pilot <- IM2017_labdata %>% select(PCRID, `Sample ID`, `Size fraction (mm)`, Sample, Notes) %>% drop_na(PCRID)

names(IMvalidation_labdata)
labdata_validation <- IMvalidation_labdata %>% select(PCRID, SampleID, `Size fraction (mm)`, Sample, Notes) %>% filter(Sample == "Car net")

labdata_pilot <- labdata_pilot %>% rename(SampleID = `Sample ID`)
labdata <- rbind(labdata_pilot, labdata_validation)

# metadata
names(metadata_validation)
metadata_validation <- metadata_validation %>% select(PCRID, SampleID, Locality, Route, Lat, Long, LandUseType, Date, TimeStart, TimeEnd, `Wind (m/s)`, `Temperature (celsius)`, SamplingNotes)

names(pilot_metadata)
morph_pilot <- pilot_metadata[, c(1, 15:140)] # save morphology in a separate table
pilot_metadata <- pilot_metadata[, c(1:13)]
pilot_metadata <- pilot_metadata %>% drop_na(Dato)
names(pilot_metadata) <- c("SampleID", "Locality", "Route", "Lat", "Long", "Quadrant", "LandUseType", "Date", "TimeStart", "TimeEnd", "Wind (m/s)", "Temperature (celsius)", "SamplingNotes")
pilot_metadata$PCRID <- NA

pilot_metadata <- pilot_metadata %>% select(PCRID, SampleID, Locality, Route, Lat, Long, LandUseType, Date, TimeStart, TimeEnd, `Wind (m/s)`, `Temperature (celsius)`, SamplingNotes)

pilot_metadata$Date <- as.Date(pilot_metadata$Date)
metadata_validation$Date <- as.Date(metadata_validation$Date, "%d-%m-%Y")
str(metadata_validation)
str(pilot_metadata)

metadata <- merge(metadata_validation, pilot_metadata, all = TRUE)
metadata <- metadata %>% drop_na(Lat)
table(metadata$LandUseType)

test <- left_join(labdata, metadata, by = "SampleID")
names(test)
data <- test[, c(1:5, 7:17)]
names(data)
data <- data %>% rename(PCRID = PCRID.x)
data <- data[!duplicated(data[,c('PCRID')]),] # only keep unique PCRIDs

write.table(data, file = "cleaned-data/data.txt")
table(data$LandUseType)
