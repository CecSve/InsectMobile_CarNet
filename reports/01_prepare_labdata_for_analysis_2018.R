# metadata preparation

### load libraries ######################
library(tidyverse)
library(data.table)
library(plyr)
library(lubridate)
library(reshape2)

### lab data ############################
labdata_2018 <- read_delim("raw-data/labdata_2018.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
imagerecognition_labdata_2018 <- read_delim("raw-data/imagerecognition_labdata_2018.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

str(labdata_2018)
str(imagerecognition_labdata_2018)

# select only the relevant columns, in this case we only need the sample names and no extra variables will be used
data1 <- labdata_2018 %>% dplyr::select(PID, SampleID_size, SampleID, PCRID3)
data2 <- imagerecognition_labdata_2018 %>% dplyr::select(PID, SampleID, PilotTripID, PCRID3)

# change column headers 
data2 <-
  data2 %>% dplyr::rename(
    SampleID_size = SampleID,
    SampleID = PilotTripID,
  )


labdata <- rbind(data1, data2)
labdata <- droplevels(labdata)

# trying to merge data with duplicate sampleids due to size sorting
str(labdata)

# change column headers 
labdata <-
  labdata %>% dplyr::rename(
    PCRID = PCRID3
  )

### sampling event data ############################
SamplingEvent <- read.csv("H:/Documents/Insektmobilen/Analysis/InsectMobile_Diversity/data/samplingevent/SamplingEvent.csv", sep=";")
routeID <- read.delim("H:/Documents/Insektmobilen/Analysis/InsectMobile_Diversity/data/samplingevent/DK_routeID.txt")
metadata <- merge(SamplingEvent, routeID, by = "SampleID")
setdiff(labdata$SampleID, metadata$SampleID) # the difference is extraction blanks, ethanol extraction tests, car net ethanol extraction tests (so see if disinfecting the net worked), and tests of different ethanol volumes on DNA concentration - they are not necessary to retain
metadata <- merge(metadata, labdata, by = c("SampleID", "PID"))
str(metadata) 
metadata <-
  metadata %>% dplyr::select(
    PCRID,
    SampleID,
    SampleID_size,
    PID,
    TripID,
    RouteID,
    LandUSeType,
    subLandUseType,
    Date,
    StartTime,
    EndTime,
    Wind,
    Temperature,
    FullySampled
  ) # only select relevant columns for analysis

table(metadata$LandUSeType)

# Rename the values from Danish to English
metadata$LandUSeType <- mapvalues(
  metadata$LandUSeType,
  from = c(
    'mark',
    'skov',
    'skov, tør',
    'tør',
    'tør, mark',
    'tør, våd',
    'urban',
    'våd',
    'våd, tør'
  ),
  to = c(
    'Farmland',
    'Forest',
    'Forest_dry',
    'Dryland',
    'Dry_agriculture',
    'Dry_wet',
    'Urban',
    'Wetland',
    'Dry_wet'
  )
)

table(metadata$FullySampled) # all are fully sampled

data<-metadata[!(metadata$PCRID=="#N/A"),] 

# get summaries of how many samples there is for each variable and their levels
length(unique(data[["RouteID"]])) # count how many routes were sampled - but notice some have received new routes
length(unique(data[["PID"]])) # count how many pilots that carried out the sampling
length(unique(data[["SampleID"]])) # count how many samples
data.frame(table(data$Wind)) # how often were the different wind categories registered
data.frame(table(data$Temperature)) # how many samples were collected at different temperature intervals
data.frame(table(data$Date)) # how many samples per day


write.table(data, file = "cleaned-data/DK_metadata_2018.txt")
