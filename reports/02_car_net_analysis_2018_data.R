
# all DK 2018 data with fwh primer
# car net analysis
# data input from Insect Diversity Manuscript 

### load libraries ################
# importing and data wrangling
library(tidyverse)
library(readr)
library(data.table)

# map construction
library(sp)
library(raster)
library(rgdal)
library(maps)
library(mapproj)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rnaturalearthhires)
library(ggsn)
library(ggsci)
library(reshape2)
library(geonames)
library(ISOcodes)
library(purrr)
library(places)

# stats
library(vegan)
library(ade4)
library(iNEXT)
library(phyloseq)
library(fossil)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ampvis2)

# visualisation
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
library(grid)
library(gridExtra)
library(jpeg)
library(gtable)
library(ggimage)
library(wesanderson)
library(sjPlot)

# GBIF comparison
library(rgbif)
library(dplyr)
library(purrr)
library(readr)  
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid_

### colour scheme ##################
iwanthue_19 <- palette(c(rgb(255,75,203, maxColorValue=255),rgb(35,124,0, maxColorValue=255), rgb(233,69,247, maxColorValue=255), rgb(77,108,0, maxColorValue=255),rgb(52,102,255, maxColorValue=255), rgb(216,141,0, maxColorValue=255), rgb(0,48,187, maxColorValue=255), rgb(236,86,0, maxColorValue=255), rgb(1,82,175, maxColorValue=255), rgb(255,70,95, maxColorValue=255), rgb(1,130,77, maxColorValue=255), rgb(237,129,250, maxColorValue=255), rgb(172,175,111, maxColorValue=255), rgb(187,154,251, maxColorValue=255), rgb(121,65,0, maxColorValue=255), rgb(138,170,243, maxColorValue=255), rgb(134,0,80, maxColorValue=255),rgb(214,140,143, maxColorValue=255), rgb(255,103,168, maxColorValue=255)))

custcol <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey")

custcol2 <- c("#8F1542", "#DBCB21", "#21B7DB", "#DB0B57", "#00738F", "darkgrey")

### load data #####################
asvs <- read.delim("cleaned-data/DK_asvtable_2018_data_totalsamples.txt", sep="\t")
data <- read.delim("cleaned-data/DK_metadata_2018_sequenced_totalsamples.txt",sep="\t")
#taxonomy_filtered <- read.delim("cleaned-data/DK_taxonomy.txt",sep="\t")
taxonomy_insects <- read.delim("cleaned-data/DK_taxonomy_Insecta.txt",sep="\t")
taxonomy_insects99 <- read.delim("cleaned-data/DK_taxonomy_Insecta_99.txt",sep="\t")
environ <- read.delim("cleaned-data/allInsects_environdata.txt",sep="\t")

# load centroid coordinates for each route
coords <- read.delim("cleaned-data/DK_ruter2018_pkt_koordinater_withregions.txt")
coords <- coords %>% dplyr::select(routeID, utm_x, utm_y, regionNavn)

### standardize data input #########

# format date
data$yDay <- yday(data$Date)
str(data)
#centering
data$cyDay <- data$yDay - median(data$yDay)

# format time
#sort time data to standard each around the time band
data$numberTime <- as.numeric(lubridate::hms(data$StartTime))#Denmark
#transform to minutes
data$numberTime <- data$numberTime/60 

# adding utm coordinates for route centroids
mergedData <-
  merge(data, coords, by.x = "RouteID", by.y = "routeID")

# first replace commas with points for decimals
mergedData$utm_x <- sapply(mergedData$utm_x, gsub, pattern = ",", replacement= ".")
mergedData$utm_y <- sapply(mergedData$utm_y, gsub, pattern = ",", replacement= ".")
str(mergedData)

# change from character to numeric
mergedData$utm_x <- as.numeric(mergedData$utm_x)
mergedData$utm_y <- as.numeric(mergedData$utm_y)
str(mergedData)

data <- mergedData

### sampling map #######

# extract Latitude and Longitude variables and put them into simple data frame called lat.long.df.
lat.long.df <- data.frame(data$utm_x, data$utm_y) 
str(lat.long.df)

coordinates(lat.long.df) <-  ~data.utm_x + data.utm_y #Function coordinates() creates a spatial object
str(lat.long.df) # at this point, this dataset doesn’t have CRS. Spatial data can’t exist without CRS(Coordinates Reference System). Spatial data requires, at least, coordinates and CRS.

proj4string(lat.long.df) <- CRS("+init=epsg:25832")
#this is +proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 

head(lat.long.df)

# Now, dataset has coordinates and CRS. Next thing is to convert this to Longitude-Latitude data.
dist.location <- spTransform(lat.long.df, CRS("+init=epsg:4326"))
dist.location

### adding regions to data based on coordinates ####
regiondata <- 
  data.frame(region = data$regionNavn,
             lat = dist.location$data.utm_x,
             lng = dist.location$data.utm_y)

# use this option for administrative regions
#options(geonamesUsername="yourusername")

#countrysub <- map2_dfr(regiondata[["lng"]], regiondata[["lat"]], .f = GNcountrySubdivision, radius = "10")

#data <- cbind(data, countrysub)

region.map <- 
  data.frame(Region = data$regionNavn,
             lat = dist.location$data.utm_x,
             long = dist.location$data.utm_y)

worldmap <- ne_countries(scale = 'large', type = 'map_units',
                         returnclass = 'sf')

denmark <- worldmap[worldmap$name == 'Denmark',]

# map the data and colour by administrative regions - note that this region subsetting does not correspond to the way taxonomists usually divide regions in Denmark when comparing species distributions, e.g. combining the Island of Bornholm with the capital region (blue) does not make much sense, since the capital region is much more urbanised compared to Bornholm
mapplot <-
  denmark %>% ggplot() + geom_sf(data = denmark,
                                 fill = "white",
                                 colour = "black") + coord_sf() + geom_point(
                                   data = region.map,
                                   aes(x = lat, y = long, colour = Region),
                                   size = 4,
                                   show.legend = T
                                 ) + theme_minimal() + scale_colour_manual(
                                   values = wes_palette(name = "Moonrise3"),
                                   labels = c("Bornholm", "Funen", "West Jutland", "East Jutland", "Sealand")
                                 ) + scalebar(
                                   denmark,
                                   dist = 25,
                                   dist_unit = "km",
                                   transform = T,
                                   model = "WGS84",
                                   st.size = 4
                                 ) + north(denmark, symbol = 4, scale = 0.07) + theme(
                                   plot.subtitle = element_text(face = "bold", size = 20),
                                   plot.margin = margin(0, 0, 0, 0, "cm"),
                                   legend.position = "bottom",
                                   legend.text = element_text(face = "bold"),
                                   legend.title = element_text(face = "bold")
                                 ) + panel_border() + labs(y = "Longitude\n", x = "\nLatitude", colour = "Biogeographic region") + guides(colour = guide_legend(nrow = 2))

ggsave("plots/Sampling_map_DK_1streview.jpg", height = 10, width = 12, dpi = 600) # remember to increase DPI for publication

### analysis  #####
#how many routes
length(unique(data$RouteID))

# how many samples
length(unique(data$SampleID)) # total/complete samples
#length(unique(data$SampleID_size)) # size sorted sample IDs

# how many pilots?
length(unique(data$PID))

# how many regions
length(unique(data$regionNavn))

# How many unique in each taxonomic level?
taxonomy <- taxonomy_insects
taxonomy %>% summarise_all(n_distinct) # for order
taxonomy %>% drop_na(family) %>% summarise_all(n_distinct) # family
taxonomy %>% drop_na(genus) %>% summarise_all(n_distinct) # genus
taxonomy %>% drop_na(species) %>% summarise_all(n_distinct) # species
table(taxonomy$order)

### insect order - Denmark comparison ##############
# load species data on insects from Denmark (downloaded manually, not updated since May 27th 2020)

allearter <- read_csv("raw-data/allearter_19092020.csv")

dkarter <-  allearter %>% group_by(Orden) %>% dplyr::summarise(value = n()) # count how many species are known for each order

# remove non-flying insect orders
remove <- c("Phthiraptera", "Siphonaptera", "Zygentoma", "Microcoryphia")
dkarter <-  dkarter %>% dplyr::filter(!Orden %in% remove)

allearter %>% dplyr::filter(!Orden %in% remove) %>% group_by(Rødlistestatus) %>% dplyr::summarise(value = n())

### get summaries for the taxonomy assigned with different filters are rename Psocodea to Psocoptera ########
allinsects <- taxonomy %>% group_by(order) %>% dplyr::summarise(value = n()) %>% 
  mutate(order = as.character(order)) %>% 
  mutate(order = replace(order, order == 'Psocodea', 'Psocoptera')) # count unique species in insect order for the data with filter: class = insecta
insects99 <- taxonomy_insects99 %>% drop_na(order) %>% group_by(order) %>% dplyr::summarise(value = n()) %>% 
  mutate(order = as.character(order)) %>% 
  mutate(order = replace(order, order == 'Psocodea', 'Psocoptera'))

# how many uniquely named species in the dataset?
unique_species <- taxonomy_insects99 %>% distinct(species, .keep_all = T) %>% drop_na(species)
uniquenames <- unique_species %>% drop_na(species) %>% group_by(order) %>% dplyr::summarise(value = n()) %>% 
  mutate(order = as.character(order)) %>% 
  mutate(order = replace(order, order == 'Psocodea', 'Psocoptera'))

# apply filter column
allinsects$filter <- "classInsecta"
insects99$filter <- "classInsecta99"
uniquenames$filter <- "uniquenames"

### Red List ###########

redlist <- allearter %>% dplyr::select(`Videnskabeligt navn`, Rødlistestatus, Fredningsstatus, `Bonn-konventionen`, `Bern-konventionen`, CITES, `NOBANIS-arter`, `NOBANIS (herkomst)`, `NOBANIS (etableringsstatus)`, `NOBANIS (invasiv optræden)`, `Øvrige forvaltningskategorier`)

redlistcar <- merge(unique_species, redlist, by.x = "species", by.y = "Videnskabeligt navn")
table(redlistcar$Rødlistestatus)
#write.table(redlistcar, file = "cleaned-data/redlistspecies.txt")

# match columns prior to merge
dkarter <- dkarter %>% dplyr::rename("order" = Orden) # first rename column header
dkarter$filter <- "Denmark"

# merge data for plotting
ordercomp <- bind_rows(allinsects, insects99, uniquenames, dkarter)

# all small insect orders will be renamed to 'other'
table(ordercomp$order)

ordercomp <-
  ordercomp %>% dplyr::mutate(
    final_order = dplyr::recode(
      order,
      "Dermaptera" = "Other",
      "Dictyoptera" = "Other",
      "Mecoptera" = "Other",
      "Megaloptera" = "Other",
      "Neuroptera"= "Other",
      "Odonata"= "Other",
      "Orthoptera"= "Other",
      "Plecoptera"= "Other",
      "Psocoptera"= "Other",
      "Raphidioptera"= "Other",
      "Strepsiptera"= "Other",
      "Ephemeroptera"= "Other",
      "Thysanoptera"= "Other",
      "Trichoptera"= "Other"
    )
  )

# which insect orders were not detected with the car net
setdiff(dkarter$order, allinsects$order)

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

### Fig 1: bar plot of insect orders ###########

plot_stack <-
  ordercomp %>% mutate(filter = fct_relevel(
    filter,
    "Denmark",
    "classInsecta",
    "classInsecta99",
    "uniquenames"
  )) %>% mutate(
    final_order = fct_relevel(
      final_order,
      "Diptera",
      "Hymenoptera",
      "Coleoptera",
      "Lepidoptera",
      "Hemiptera",
      "Other"
    )
  ) %>% ggplot(aes(
    x = reorder(filter, -value),
    y = value,
    fill = final_order
  )) + geom_bar(stat = "identity",
                width = 0.7,
                position = "dodge")  + scale_fill_grey() + labs(x = "", y = "Species richness\n", fill = "Insect order") + theme_cleveland() + theme(
                  title = element_text(face = "bold"),
                  legend.position = "bottom",
                  plot.margin = margin(0, 0, 0, 0, "cm"),
                  axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.y = element_text(size = 14, face = "bold"),
                  panel.background = element_rect(
                    fill = "white",
                    colour = "white",
                    size = 0.5,
                    linetype = "solid"
                  )
                ) + guides(fill = guide_legend(nrow = 1)) + scale_x_discrete(
                  labels = c(
                    "Denmark" = "Known species \nin Denmark",
                    "classInsecta" = "Insect ASVs",
                    "classInsecta99" = "Insect ASVs with a 99% or higher \nmatch to reference database",
                    "uniquenames" = "Unique species names \nin reference database"
                  )
                ) 

ggsave("plots/barplots_review.jpg", height = 10, width = 12, dpi = 600) 

### richness per sample, region ####
tasvs <- t(asvs)
patasvs <- decostand(tasvs, method = "pa")

patasvs <- as.data.frame(patasvs)
patasvs$richness <- rowSums(patasvs)
samplerich <- patasvs %>% rownames_to_column(var = "SampleID") %>% dplyr::select(SampleID, richness)
data <- merge(data, samplerich, by = "SampleID")

# richness per sample
sample_rich_plot <-
  ggplot(data = data, mapping = aes(x = reorder(SampleID,-richness), y = richness)) + geom_bar(
    stat = "identity",
    width = 0.5,
    position = "dodge",
    show.legend = F
  ) + theme_cleveland() +  theme(
    title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0, "cm"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.background = element_rect(
      fill = "white",
      colour = "white",
      size = 0.5,
      linetype = "solid"
    )
  ) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(limits = c(0, 100, 200, 300, 400),expand = c(0, 0)) + labs(x = "Samples (ordered by descending richness)",  y = "Species richness\n") + geom_hline( yintercept = mean(data$richness),linetype = 2, color = "black" ) + geom_hline(yintercept = median(data$richness), color = "black") +
  annotate(
    "text",
    x = "P42.2A",
    y = 133,
    label = "median = 126 \nmean = 131.2",
    vjust = -0.5
  )

### richness per administrative region ####
test <- taxonomy_insects %>% rownames_to_column(var = "otuid") %>% dplyr::select(otuid, order)
test2 <- asvs 
test2 <- decostand(test2, method = "pa") %>% rownames_to_column(var = "otuid")
test3 <- merge(test, test2, by = "otuid")
test3 <- test3 %>% dplyr::mutate(
    order_final = dplyr::recode(
      order,
      "Dermaptera"= "Other",
      "Dictyoptera" = "Other",
      "Mecoptera" = "Other",
      "Megaloptera" = "Other",
      "Neuroptera"= "Other",
      "Odonata"= "Other",
      "Orthoptera"= "Other",
      "Plecoptera"= "Other",
      "Psocodea"= "Other",
      "Raphidioptera"= "Other",
      "Strepsiptera"= "Other",
      "Ephemeroptera"= "Other",
      "Thysanoptera" = "Other",
      "Trichoptera"= "Other"
    )
  )

test3 <- dplyr::select(test3, -order)

test4 <- test3 %>% column_to_rownames(var = "otuid")

test5 <- test4 %>% group_by(order_final) %>%
  summarise_all(.funs = sum) %>% 
  pivot_longer(-order_final, names_to = "SampleID", values_to = "order_richness")

test6 <- data %>% dplyr::select(SampleID, regionNavn)
test7 <- merge(test6, test5, by = "SampleID")

# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
tgc <- summarySE(test7, measurevar="order_richness", groupvars=c("regionNavn","order_final"))
tgc

# Error bars represent standard error of the mean
tgc %>% ggplot(mapping = aes(x=reorder(regionNavn,-order_richness), y=order_richness, fill = order_final)) + geom_bar(stat="identity", position = position_dodge(), show.legend = F)  +
  geom_errorbar(aes(ymin=order_richness-se, ymax=order_richness+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + theme_cleveland() +
  theme(axis.text.x=element_text(),
        axis.ticks.x=element_blank(), axis.title.x = element_text(size = 12, face = "bold"),  axis.title.y = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid")) + scale_x_discrete(expand = c(0, 0), labels = c("Funen", "West Jutland", "East Jutland", "Sealand", "Bornholm")) + scale_y_discrete(limits = c(0, 25, 50, 75, 100), expand = c(0, 0)) + scale_fill_manual(values = wes_palette(name = "IsleofDogs1")) + labs(
          x = "",
          y = "Mean flying insect richness\n",
          fill = "Biogeographic region \nin Denmark")

plot_region <- data %>% ggplot(mapping = aes(x=reorder(regionNavn,-richness), y=richness, fill = regionNavn)) + geom_bar(stat="identity", width=0.5, position = "dodge", show.legend = F) + theme_cleveland() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.x = element_text(size = 12, face = "bold"),  axis.title.y = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid")) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(limits = c(0, 100, 200, 300, 400), expand = c(0, 0)) + scale_fill_manual(values = wes_palette(name = "Moonrise3")) + labs(
          x = "",
          y = "Flying insect richness\n",
          fill = "Biogeographic region \nin Denmark")

data <- merge(data, environ, by = "SampleID")
data[ ,c(22:133)] <- data[,c(22:133)]*100 # scale by 100 to get percentage cover

# Create new column that picks out samples that are at the extreme in representing a single land cover type (>60% or >80% of a single land cover type) at the buffer zone with the largest effect size 

data$hab50 = 'Mix' # samples that do not have more than 60% of one specific land type
data$hab50[data$Agriculture_1000>=50]<-'Farmland50'
data$hab50[data$Forest_1000>=50]<-'Forest50'
data$hab50[data$Urban_1000>=50]<-'Urban50'
table(data$hab50) # notice the variation in sample size

### richness land cover ####
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
statsdata <- summarySE(data, measurevar="richness", groupvars=c("regionNavn", "hab50"))
statsdata
statsdata <- statsdata %>% filter(!regionNavn == "Bornholm")

# Error bars represent standard error of the mean
plot_cover <- statsdata %>% mutate(hab50 = fct_relevel(hab50, "Farmland50", "Forest50", "Mix", "Urban50")) %>% ggplot(mapping = aes(x=regionNavn, y= richness, fill = hab50)) + geom_bar(stat="identity", position = position_dodge(), show.legend = F) + geom_errorbar(aes(ymin=richness-se, ymax=richness+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + theme_cleveland() +
  theme(axis.text.x=element_text(),
        axis.ticks.x=element_blank(), axis.title.x = element_text(size = 12, face = "bold"),  axis.title.y = element_text(size = 14, face = "bold"), panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),plot.margin = margin(0, 0, 0, 0, "cm")) + scale_x_discrete(expand = c(0, 0), labels = c("Funen", "West Jutland", "East Jutland", "Sealand")) + scale_y_discrete(limits = c(0, 50, 100, 150, 200), expand = c(0, 0)) + scale_fill_manual(values =custcol[c(2,5,6,1)]) + labs(
          x = "",
          y = "Mean flying insect richness\n",
          fill = ">50% land cover")

#library(indicspecies)
#abund <- t(pa) # based on presence absence
#land <- allInsects_totsample$hab50

# example with syrphids
tax_syr <- taxonomy_insects %>% rownames_to_column(var = "otuid") %>% dplyr::filter(family == "Syrphidae")

syrphids <- tax_syr$otuid

otu_syr <- patasvs[, (names(patasvs) %in% syrphids)] 
rowSums(otu_syr)
otu_syr <- otu_syr[rowSums(otu_syr) > 0, ]
syr_samples <- rownames(otu_syr)
data_syr <- data %>% filter(SampleID %in% syr_samples)
data_syr$syr_richness <- rowSums(otu_syr)

### mean richness plots (not used) ####
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
tgc <- summarySE(data_syr, measurevar="syr_richness", groupvars=c("regionNavn","hab50"))
tgc
tgc <- tgc %>% dplyr::filter(N > 1)
tgc <- tgc %>% filter(!regionNavn == "Bornholm")

# Error bars represent standard error of the mean
tgc %>% mutate(hab50 = fct_relevel(hab50, "Farmland50", "Forest50", "Mix", "Urban50")) %>% ggplot(mapping = aes(x=reorder(regionNavn,-syr_richness), y=syr_richness, fill = hab50)) + geom_bar(stat="identity", position = position_dodge(), show.legend = F)  +
  geom_errorbar(aes(ymin=syr_richness-se, ymax=syr_richness+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + theme_cleveland() +
  theme(axis.text.x=element_text(),
        axis.ticks.x=element_blank(), axis.title.x = element_text(size = 12, face = "bold"),  axis.title.y = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid")) + scale_x_discrete(expand = c(0, 0), labels = c("Funen", "West Jutland", "East Jutland", "Sealand")) + scale_y_discrete(limits = c(0, 1, 2, 3, 4, 5), expand = c(0, 0)) + scale_fill_manual(values = custcol[c(2,5,6,1)]) + labs(
          x = "",
          y = "Mean flying hoverfly richness\n",
          fill = "Biogeographic region \nin Denmark")

hoverfly <- png::readPNG("H:/Documents/Insektmobilen/Billeder og videoer/Billeder af insekt silhouetter/svirreflue_noback.png")

### total hoverfly richness by region #####
syr_region <-
  data_syr %>% ggplot(mapping = aes(
    x = reorder(regionNavn, -syr_richness),
    y = syr_richness,
    fill = regionNavn
  )) + geom_bar(
    stat = "identity",
    width = 0.5,
    position = "dodge",
    show.legend = F
  ) + theme_cleveland() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    panel.background = element_rect(
      fill = "white",
      colour = "white",
      size = 0.5,
      linetype = "solid"
    )
  ) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(limits = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                                            expand = c(0, 0)) + scale_fill_manual(values = wes_palette(name = "Moonrise3")) + labs(x = "", y = "Syrphidae species richness\n", fill = "Biogeograpgic region \nin Denmark")

syr_region_img <- ggdraw() + draw_plot(syr_region) +
  draw_image(hoverfly,  x = 0.3, y = 0.30, scale = .4) 
syr_region_img

#bottom_row <- plot_grid(plot_region, syr_region_img, labels = c('B', 'C'), label_size = 12, ncol = 1)
#fig1_plot <- plot_grid(mapplot, bottom_row, labels = c('A', ''), label_size = 12, ncol = 2)
#save_plot("plots/fig1_1streview.jpg", fig1_plot, base_height = 10, base_width = 14)
#ggsave('plots/figure1_review.jpg', plot = fig1_plot, width=450,height=550, units = "mm", dpi=600)

plots <- align_plots(mapplot, plot_region, align = 'h', axis = 'l')
# then build the bottom row
bottom_row <- plot_grid(plots[[2]], syr_region_img, labels = c('B', 'C'), label_size = 12)

# then combine with the top row for final plot
figure1 <- plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 2)
save_plot("plots/fig1_1streview.jpg", figure1, base_height = 8, base_width = 16)

data_syr %>% group_by(regionNavn) %>% dplyr::summarise(value = n())

### hoverfly richness across land cover types ####
syr_hab <-
  data_syr %>% mutate(hab50 = fct_relevel(hab50, "Farmland50", "Forest50", "Mix", "Urban50")) %>% ggplot(aes(
    x = hab50,
    y = syr_richness,
    fill = hab50,
    group = 1
  )) + geom_bar(
    stat = "identity",
    width = 0.5,
    position = "dodge",
    show.legend = T
  ) + theme_cleveland() + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.background = element_rect(
      fill = "white",
      colour = "white",
      size = 0.5,
      linetype = "solid"
    ),
    legend.position = "bottom"
  ) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(limits = c(0, 2, 4, 6, 8, 10), expand = c(0, 0)) + scale_fill_manual(
    values = custcol[c(2, 5, 6, 1)],
    labels = c(">50% Farmland", ">50% Forest", "Mix", ">50% Urban")
  ) + labs(x = "", y = "Syrphidae species richness\n", fill = "Land cover")

syr_hab_img <- ggdraw() + draw_plot(syr_hab) +
  draw_image(hoverfly,  x = 0.4, y = 0.30, scale = .35) 
syr_hab_img

### not used ####
# Calculates mean, median, sd, se and IC
my_sum <- data_syr %>%
  group_by(hab50) %>%
  summarise( 
    n=n(),
    mean=mean(syr_richness),
    median=median(syr_richness),
    sd=sd(syr_richness)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

#limits <- aes(ymax = resp + se, ymin = resp - se)

# Median & Confidence Interval
my_sum %>% mutate(hab50 = fct_relevel(hab50, "Farmland50", "Forest50", "Mix", "Urban50")) %>% ggplot(aes(x=hab50, y= median, fill = hab50, group = 1)) + geom_bar(stat="identity", width=0.5, position = "dodge", show.legend = T) + geom_errorbar(aes(x= as.factor(hab50), ymin=median, ymax=median+se), width=0.4, colour="darkgrey", alpha=0.9, size=1.5) + geom_errorbar(aes(x= as.factor(hab50), ymin=median, ymax=median-se), width=0.4, colour="white", alpha=0.9, size=1.5) + theme_cleveland() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.x = element_text(size = 12, face = "bold"),  axis.title.y = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid")) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(limits = c(0, 1, 2, 3, 4, 5), expand = c(0, 0)) + scale_fill_manual(values =custcol[c(2,5,6,1)], labels = c("Farmland", "Forest", "Mix", "Urban")) + labs(
          x = ">50% land cover",
          y = "Median Syrphidae species richness\n",
          fill = "Land cover")

### Fig 2: proportional data & land cover ########

prow <- plot_grid(
  plot_cover, # + theme(legend.position="none"),
  syr_hab_img, #+ theme(legend.position="none"),
  align = 'h',
  labels = c('C', 'D'), 
  hjust = -1,
  label_size = 12,
  ncol = 1
)
prow

zrow <- plot_grid(
  plot_stack, #+ theme(legend.position="none"),
  sample_rich_plot, #+ theme(legend.position="none"),
  align = 'h',
  labels = c('A', 'B'), 
  hjust = -1,
  label_size = 12,
  ncol = 1
)
zrow

fig2 <- plot_grid(zrow, prow, align = "hv", rel_widths = c(1.7,1), scale = 0.9)

save_plot("plots/fig2_1streview.jpg", fig2, base_height = 10, base_width = 20)

### richness analysis #####
hist(data$richness)
hist(log(data$richness))
hist(sqrt(data$richness))

# all insects
options(na.action = "na.fail")
fm1 <- lmer(sqrt(richness) ~ 
              regionNavn + 
              hab50 +
              Wind +
              Temperature +
              (1|RouteID) + (1|PID), data=data)

dd <- dredge(fm1)
subset(dd, delta < 4)
#'Best' model
summary(get.models(dd, 1)[[1]])

lme1000 <- lmer(sqrt(richness) ~ 
                  regionNavn + 
                  hab50 +
                  #Wind +
                  #Temperature +
                  (1|RouteID) + (1|PID), data=data)
# data=subset(allInsects, Open.uncultivated.land_1000 < 0.2)
summary(lme1000)
qqnorm(resid(lme1000))
AICc(lme1000)

tab_model(lme1000, pred.labels = c("Intercept (>50% farmland cover)", "Funen", "West Jutland", "East Jutland", "Sealand", ">50% forest cover", "Mixed cover", ">50% urban cover"))

# emmeans
library(emmeans)

ems <- emmeans(lme1000, specs = pairwise ~ regionNavn, type = "response")
ems$contrasts
ems$emmeans
plot(ems, comparisons = TRUE)

ems <- emmeans(lme1000, specs = pairwise ~ hab50, type = "response")
ems$contrasts
ems$emmeans
plot(ems, comparisons = TRUE)

# syrphidae
hist(data_syr$syr_richness)
hist(log(data_syr$syr_richness))
hist(sqrt(data_syr$syr_richness)) # not normally distributed by transformation

lme1000 <- lmer(sqrt(syr_richness) ~ 
                  regionNavn + 
                  hab50 +
                  #Wind +
                  #Temperature +
                  (1|RouteID) + (1|PID), data=data_syr)

summary(lme1000)
qqnorm(resid(lme1000))
AICc(lme1000)

tab_model(lme1000, pred.labels = c("Intercept (>50% farmland cover)", "Funen", "West Jutland", "East Jutland", "Sealand", ">50% forest cover", "Mixed cover", ">50% urban cover"))

ems <- emmeans(lme1000, specs = pairwise ~ regionNavn, type = "response")
ems$contrasts
ems$emmeans
plot(ems, comparisons = TRUE)

ems <- emmeans(lme1000, specs = pairwise ~ hab50, type = "response")
ems$contrasts
ems$emmeans
plot(ems, comparisons = TRUE)

### new species for DK? ###############
# how many species were detected by the net and is not present in the Danish species database? 
newspecies <- setdiff(taxonomy_insects99$species, allearter$`Videnskabeligt navn`) 

# comparison of car net and DK database
newspeciesDK <- taxonomy_insects99 %>% rownames_to_column(var = "ASVID") %>% filter(species %in% newspecies) 
#write.table(newspeciesDK, file = "cleaned-data/newspecies_DK_alldata.txt", sep = "\t", row.names = F)
#newspeciesDK <- newspeciesDK %>% column_to_rownames(var = "ASVID")
newdistinct_species <- newspeciesDK %>% distinct(species, .keep_all = T) %>% drop_na(species)
newspeciesDK <- newspeciesDK %>% group_by(order, species) %>% distinct(species) %>% drop_na(species)
#write.table(newspeciesDK, file = "cleaned-data/newspecies_DK.txt", sep = "\t", row.names = F)

### gbif comparison ###############

# run crossed out lines if you want to generate new comparison to GBIF 

# fill in your gbif.org credentials 
#user <- "" # your gbif.org username 
#pwd <- "" # your gbif.org password
#email <- "" # your email 

# match the names 
#gbif_taxon_keys <- 
#  newspeciesDK %>% 
#  pull("species") %>% # use fewer names if you want to just test 
#  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
#  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
#  bind_rows() %T>% # combine all data.frames into one
#  readr::write_tsv(path = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
#filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
#  filter(class == "Insecta") %>% # remove anything that might have matched to a non-insect
#  pull(usagekey) # get the gbif taxonkeys

# comparison to DK
# gbif_taxon_keys should be a long vector like this c(2977832,2977901,2977966,2977835,2977863)
# !!very important here to use pred_in!!
# use matched gbif_taxon_keys from above 
#occ_download(
#  pred_in("taxonKey", gbif_taxon_keys),
#  pred("country", "DK"),
#  format = "SIMPLE_CSV",
#  user=user,pwd=pwd,email=email
#) # this generates a file for your user where the matches are

# comparison to region (Germany, Norway and Sweden)
# use matched gbif_taxon_keys from above 
#occ_download(
#  pred_in("taxonKey", gbif_taxon_keys),
#  pred_in("country", c("SE", "NO", "DE")),
#  format = "SIMPLE_CSV",
#  user=user,pwd=pwd,email=email
#) # this generates a file for your user where the matches are

gbif_query <- read.delim("raw-data/0080813-200613084148143.csv") # GBIF.org (07 October 2020) GBIF Occurrence Download https://doi.org/10.15468/dl.2mker8 
gbif_query_scandi <- read.delim("raw-data/0080815-200613084148143.csv") # GBIF.org (07 October 2020) GBIF Occurrence Download https://doi.org/10.15468/dl.zvumvk 

#test <- gbif_query %>% distinct(species, .keep_all = T)

setdiff(newspeciesDK$species, gbif_query$species) # How many species are not found in DK
setdiff(newspeciesDK$species, gbif_query_scandi$species) # How many species are not found in the countries bordering DK 

nooccurrence_DK <- anti_join(newdistinct_species, gbif_query, by = "species")
nooccurrence_scandi <- anti_join(newdistinct_species, gbif_query_scandi, by = "species")

write.table(nooccurrence_DK, file = "cleaned-data/nooccurrence_DK.txt", sep = "\t", row.names = F)
write.table(nooccurrence_scandi, file = "cleaned-data/nooccurrence_scandi.txt", sep = "\t", row.names = F)

# which order does the new species belong to and what is the frequency 
newspeciesDK %>%
  group_by(order) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = (n / sum(n))*100) # get the proportion in percent

test <- allearter %>%
  group_by(Orden) %>%
  dplyr::summarise(n = n())

test2 <- taxonomy %>%
  group_by(order) %>%
  summarise(n = n())

test <- test %>% dplyr::rename("order" = Orden) %>% dplyr::rename("count_allearter" = n)
test2 <- test2  %>% dplyr::rename("count_carnet" = n) %>% mutate(order = replace(order, order == 'Psocodea', 'Psocoptera'))

test3 <- left_join(test2, test, by = "order")
test3$prop <- ((test3$count_carnet/test3$count_allearter)*100)

### species accumulation curves #########
otus <- decostand(asvs, method = "pa") # transform into presence absence

## Accumulation model

#Data has species as rows and sites as columns, but we need the opposite. Thus, transpose the original data
otusT <- t(otus)
#get richness estimators (for each sample, cumulative)
pool <- poolaccum(otusT, permutations = 1000)
plot(pool)

# using the fossil library to calculate richness estimates - including chao2
test <- spp.est(otus, rand = 100, abund = FALSE, counter = FALSE, max.est = 'all')
test2 <- as_tibble(test)

str(test2)

# the column names are bonkers so they need to be renamed
speciesest <- test2 %>% rename(sample = `# Samples`, Observed = S.obs, s.obs.upper = `S.obs(+95%)`, s.obs.lower = `S.obs(-95%)`,  chao2.upper = `Chao2(upper)`, chao2.lower = `Chao2(lower)`, ICE.upper = `ICE(lower)`, ICE.lower = `V10`, jack1.upper = `Jack1(lupper)`, jack1.lower = `Jack1(lower)`)

speciesest %>% tidyr::gather("id", "value", c(2, 5, 8, 11)) %>% 
  ggplot(., aes(sample, value))+
  geom_smooth(method = "loess", se=FALSE, color="black", size = 2)+
  facet_wrap(~id) + theme_minimal_grid() + labs(x = "Number of samples", y = "Number of insect ASVs") 

# rarefaction curve
#build the species accumulation curve & rarefaction curve (expected)
#otus.specaccum <- specaccum(otusT,method = "exact", permutations = 100) 
#plot the curve with some predefined settings
#plot(otus.specaccum,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#creating a dataframe for ggplot2
#data_specaccum <- data.frame(Sites=otus.specaccum$sites, Richness=otus.specaccum$richness, SD=otus.specaccum$sd)

#ggplot() + geom_point(data=data_specaccum, aes(x=Sites, y=Richness)) + geom_line(data=data_specaccum, aes(x=Sites, y=Richness)) + geom_ribbon(data=data_specaccum ,aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)

# plot more pretty
carnet <- png::readPNG("plots/car.png")

# test#
g <- rasterGrob(carnet, width=unit(80,"lines"), height=unit(60,"lines"), interpolate = T)


#build a expected curve (randomization for boxplot comparison)
#otus.specaccum.rand <- specaccum(otusT, "random")
#plot both curves ("observed" vs "randomized")
#plot(otus.specaccum,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
#boxplot(otus.specaccum.rand, col="yellow", add=TRUE, pch="+")

### iNEXT visualisation ####

otuspa <- iNEXT::as.incfreq(otus)# transform incidence raw data (a species by sites presence-absence matrix) to incidence frequencies data (iNEXT input format, a row-sum frequencies vector contains total number of sampling units)

# set a series of sample sizes (m) for R/E computation
t <- seq(1, 4000, by=50)

#apply `iNEXT` main function
otuspa.inext <- iNEXT(otuspa, q = c(0, 1, 2), datatype = "incidence_freq", size = t) # calculate for all three hill numbers
otuspa.inext$DataInfo # summarizing data information, returns basic data information including the reference sample size (n), observed species richness (S.obs), a sample coverage estimate (SC), and the first ten frequency counts (f1‐f10)
otuspa.inext$iNextEst # showing diversity estimates along with related statistics for a series of rarefied and extrapolated samples
otuspa.inext$AsyEst # showing asymptotic diversity estimates along with related statistics
ChaoRichness(otuspa, datatype = "incidence_freq", conf = 0.95)
#look at the data
otuspa.inext

out.inc <- iNEXT(otuspa, q=0, datatype="incidence_freq", size=t)
# Sample‐size‐based R/E curves
ggiNEXT(out.inc, type=1, color.var="site") + 
  theme_bw(base_size = 18) + 
  theme(legend.position="none")

# Sample completeness curves
ggiNEXT(out.inc, type=2, color.var="site") +
  ylim(c(0,1)) +
  theme_bw(base_size = 18) + 
  theme(legend.position="none")

# Coverage‐based R/E curves
ggiNEXT(out.inc, type=3, color.var ="site") + 
  xlim(c(0,1)) +
  theme_bw(base_size = 18) +
  theme(legend.position="bottom",
        legend.title=element_blank())

#plot the results
# Sample‐size‐based R/E curves - plots diversity estimates with confidence intervals (if se=TRUE) as a function of sample size up to double the reference sample size, by default, or a user‐specified endpoint.
plot.inext <-
  ggiNEXT(
    otuspa.inext,
    se = TRUE,
    type = 1,
    color.var = "order",
    grey = T
  ) + theme_pubclean() + theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    rect = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent"),
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16, face = "bold")
  ) + scale_fill_grey(start = 0, end = .4) +
  scale_colour_grey(start = .2, end = .2) + scale_shape_discrete(labels = c("ASV richness", "Shannon diversity", "Simpson diversity")) + guides(fill = FALSE, colour = FALSE) + labs(x = "\nNumber of samples", y = "Number of insect ASVs",subtitle = "C")  +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

#ggsave("plots/inext_accumulation.png", grid.draw(gList(rasterGrob(carnet, width = unit(1,"npc"), height = unit(1,"npc")), ggplotGrob(plot.inext))))

grid.draw(gList(rasterGrob(carnet, width=unit(40,"lines"), height=unit(30,"lines"), interpolate = T), 
                ggplotGrob(plot.inext)))

# Sample completeness curves - plots the sample coverage with respect to sample size for the same range described in sample-size-based R/E curves
ggiNEXT(otuspa.inext, se = T, type=2, color.var="order") + theme_cowplot()

# Coverage‐based R/E curves - plots the diversity estimates with confidence intervals (if se=TRUE) as a function of sample coverage up to the maximum coverage obtained from the maximum size described in sample-size-based R/E curves
ggiNEXT(otuspa.inext, se = TRUE, type=3, color.var ="order") + theme_cowplot()

# extract data for ggplot 
#accdata <- summary(pool, display = "chao")
#observed <- summary(pool, display = "S")
#data_plot <- as_tibble(accdata[["chao"]])
#observed_plot <- as_tibble(observed[["S"]])

#ggsave("plots/zoom_accumulation_chao1.png")

as_ggplot(gList(rasterGrob(carnet, width=unit(80,"lines"), height=unit(60,"lines"), interpolate = T), ggplotGrob(plot.inext)))

ggsave("plots/Fig1C_richness_DK.jpg", height = 10, width = 12, dpi = 600)

test <- as_ggplot(gList(rasterGrob(carnet, width=unit(52,"lines"), height=unit(42,"lines"), interpolate = T), ggplotGrob(plot.inext)))

g1 <- ggplotGrob(mapplot)
g2 <- ggplotGrob(stacked_plot) #small_leg_plot
g3 <- ggplotGrob(plot.inext)
#g3 <- ggplotGrob(acummulation_plot)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
grid.newpage()
grid.draw(g)

ggsave('plots/figure1.jpg', plot = g, width=450,height=550, units = "mm", dpi=600)

#figure1 <- ggarrange(mapplot, small_leg_plot, b, ncol = 1, align = "v")
#save_plot("plots/fig1_relabun_estimaterich.png", figure1, base_height = 10, base_width = 8)

### z-test  ##################
# how many unique observations in taxonomy
length(unique(taxonomy$order))
length(unique(taxonomy$family))
length(unique(taxonomy$genus))
length(unique(taxonomy$species))

# in the Danish species list
allearter %>% dplyr::filter(!Orden %in% remove) %>% group_by(Orden) %>% dplyr::summarise(value = n()) # count how many species are known for each order
dkarter_z <- allearter %>% dplyr::filter(!Orden %in% remove) # remove non-flying insect orders
length(unique(dkarter_z$Orden))
length(unique(dkarter_z$Familie))
length(unique(dkarter_z$Sl?gt))
length(unique(dkarter_z$`Videnskabeligt navn`))

# order
res <- prop.test(x = c(15, 19), n = c(19, 19))
# Printing the results
res

# family
res <- prop.test(x = c(241, 485), n = c(485, 485), alternative = "less")
# Printing the results
res

# genus
res <- prop.test(x = c(1274, 5467), n = c(5467, 5467), alternative = "less")
# Printing the results
res

# species
res <- prop.test(x = c(2115, 18791), n = c(18791, 18791), alternative = "less")
# Printing the results
res


test <- allearter %>% dplyr::filter(!Orden %in% remove) %>%
  group_by(Orden) %>%
  summarise(n = n())

taxonomy %>%
  group_by(order) %>%
  summarise(n = n())

# Is the proportion of insect orders car nets equal to the proportion of insect orders in Denmark
# Coleoptera
res <- prop.test(x = c(1014, 3866), n = c(8172, 18791))
# Printing the results
res 

# Dermaptera
res <- prop.test(x = c(1, 6), n = c(8172, 18791))
# Printing the results
res 

# Diptera
res <- prop.test(x = c(4828, 5086), n = c(8172, 18791))
# Printing the results
res 

# Ephemeroptera
res <- prop.test(x = c(12, 43), n = c(8172, 18791))
# Printing the results
res 

# Hemiptera  
res <- prop.test(x = c(656, 1495), n = c(8172, 18791))
# Printing the results
res 

# Hymenoptera     
res <- prop.test(x = c(1358, 5150), n = c(8172, 18791))
# Printing the results
res 

# Lepidoptera      
res <- prop.test(x = c(129, 2590), n = c(8172, 18791))
# Printing the results
res 

# Mecoptera         
res <- prop.test(x = c(3, 4), n = c(8172, 18791))
# Printing the results
res 

# Neuroptera        
res <- prop.test(x = c(2, 62), n = c(8172, 18791))
# Printing the results
res 

# Odonata           
res <- prop.test(x = c(9, 60), n = c(8172, 18791))
# Printing the results
res 
# Orthoptera        
res <- prop.test(x = c(10, 38), n = c(8172, 18791))
# Printing the results
res 
# Plecoptera        
res <- prop.test(x = c(3, 25), n = c(8172, 18791))
# Printing the results
res 
# Psocoptera (Psocodea)
res <- prop.test(x = c(49, 62), n = c(8172, 18791))
# Printing the results
res 
# Thysanoptera     
res <- prop.test(x = c(89, 113), n = c(8172, 18791))
# Printing the results
res 
# Trichoptera   
res <- prop.test(x = c(9, 172), n = c(8172, 18791))
# Printing the results
res 

### phyloseq #############
OTU = otu_table(asvs, taxa_are_rows = TRUE)
tax <- taxonomy_insects %>% rownames_to_column(var = "id") %>% dplyr::select(id, kingdom, phylum, class, order, family, genus, species) %>% column_to_rownames(var = "id")
tax <- as.matrix(tax, rownames.force = TRUE) # force the table into a matrix fromat so a phyloseq tax object can be created
TAX = tax_table(tax)
OTU
TAX

test <- data %>% dplyr::select(SampleID, Date, Wind, Temperature) %>% column_to_rownames(var = "SampleID")
metadata <- sample_data(test) # create the phyloseq sample data object, needs to be reordered to fit with PCRIDs as row names

# now we can create a phyloseq object
ps <- phyloseq(OTU, TAX, metadata)
ps

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(ps))

# mean, max and min of sample read counts
smin <- min(sample_sums(ps))
smean <- mean(sample_sums(ps))
smax <- max(sample_sums(ps))

### RAM analysis ######

library(RAM)

# To be able to merge the OTU table with the taxonomy table, they need to have a common column to call
otus <- rownames_to_column(asvs, var = "otuid")
taxonomy_insects <- read.delim("cleaned-data/DK_taxonomy_Insecta.txt",sep="\t")
taxonomy <- rownames_to_column(taxonomy_insects, var = "otuid") # remeber to remake the column into rownames for both datasets if you need to

# format taxonomy for RAM
taxonomy$kingdom <- paste0("k__", taxonomy$kingdom)
taxonomy$phylum <- paste0("; p__", taxonomy$phylum)
taxonomy$class <- paste0("; c__", taxonomy$class)
taxonomy$order <- paste0("; o__", taxonomy$order)
taxonomy$family <- paste0("; f__", taxonomy$family)
taxonomy$genus <- paste0("; g__", taxonomy$genus)
taxonomy$species <- paste0("; s__", taxonomy$species)

merged_taxonomy <- within(taxonomy, taxonomy <- paste(kingdom,phylum,class, order, family, genus, species))
merged_taxonomy <- merged_taxonomy %>% dplyr::select(otuid, taxonomy)

ramdata <- left_join(otus, merged_taxonomy, by = "otuid") # need to make otuid back to rownames
ramdata <- column_to_rownames(ramdata, var = "otuid")
ramdata <- as.data.frame(ramdata, stringsAsFactors = FALSE) # don't know if this is required

valid.taxonomy(data = list(data = ramdata)) # check if the format for taxonomy is correct
input.ranks <- c("kingdom", "phylum", "class", "order", 
                 "family", "genus", "species")
reform.data <- reformat.taxonomy(list(data=ramdata),
                                 input.ranks=input.ranks,
                                 sep="; ")[[1]]

# relative oTU abundance by taxonomic rank per sample
group.abundance(data = list(reform.data = reform.data), rank = "o", ggplot2 = TRUE)

# the order of the pcrid has to match the order of the sample columns in taxonomy
attach(data)
metadata <- data[order(PCRID),]
rownames(metadata) <- NULL
metadata <- metadata %>% column_to_rownames(var = "PCRID")
metadata$samplingmethod <- "carnet"
reform.data.test <- reform.data[,order(colnames(reform.data))]

### ampvis analysis - rarefaction curve ####
otus <- rownames_to_column(asvs, var = "otuid")
taxonomy <- rownames_to_column(taxonomy, var = "otuid")
taxonomy_amp <- taxonomy %>% dplyr::select(otuid, kingdom, phylum, class, order, family, genus, species)
otutable <- left_join(otus, taxonomy_amp, by = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames
otutable <-
  otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)
amp_abun <- amp_load(otutable = otutable, metadata = data)
amp_abun

# rarefaction curve
rarecurve <-
  amp_rarecurve(
    amp_abun,
    facet_by = "regionNavn",
    stepsize = 100,
    color_by = "regionNavn",
    facet_scales = "free"
  )

region_names <- list(
  'Bornholm'="Bornholm",
  'Fyn'="Funen",
  'Jylland_ATL'="West Jutland",
  'Jylland_KON'="East Jutland",
  'SjLolFal'="Sealand"
)

region_labeller <- function(variable,value){
  return(region_names[value])
}

rarfac <- rarecurve + scale_colour_manual(values = wes_palette(name = "Moonrise3")) + geom_line(size = 1, na.rm = TRUE) + scale_x_continuous(
  labels = scales::number,
  limits = c(0, 200000),
  breaks = seq(from = 0, to = 200000, by = 5000)
) + guides(colour=guide_legend(ncol=8)) + ylab("Number of observed ASVs") + xlab("Sequencing depth (reads)") + theme_minimal()+ theme(legend.position = "none", axis.text.x = element_text(angle = 90), legend.title = element_blank()) + facet_wrap(regionNavn ~ ., labeller=region_labeller, ncol = 2) 

cowplot::save_plot("plots/rarecurve_bioregion.png", rarfac, base_height = 10, base_width = 12)
