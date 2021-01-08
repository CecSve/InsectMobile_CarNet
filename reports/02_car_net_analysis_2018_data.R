
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

# stats
library(vegan)
library(ade4)
library(iNEXT)
library(phyloseq)
library(fossil)

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

### load data #####################
asvs <- read.delim("cleaned-data/DK_asvtable_2018_data_totalsamples.txt", sep="\t")
data <- read.delim("cleaned-data/DK_metadata_2018_sequenced_totalsamples.txt",sep="\t")
#taxonomy_filtered <- read.delim("cleaned-data/DK_taxonomy.txt",sep="\t")
taxonomy_insects <- read.delim("cleaned-data/DK_taxonomy_Insecta.txt",sep="\t")
taxonomy_insects99 <- read.delim("cleaned-data/DK_taxonomy_Insecta_99.txt",sep="\t")

# load centroid coordinates for each route
coords <- read.delim("cleaned-data/DK_ruter2018_pkt_koordinater.txt")

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
mergedData <-
  dplyr::select(mergedData,-OBJECTID) # remove objectid column since it is not needed  

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

landuse.map <- 
  data.frame(Landuse = data$LandUSeType,
             lat = dist.location$data.utm_x,
             long = dist.location$data.utm_y)

worldmap <- ne_countries(scale = 'large', type = 'map_units',
                         returnclass = 'sf')

denmark <- worldmap[worldmap$name == 'Denmark',]

mapplot <- denmark %>%
  ggplot() + 
  geom_sf(data = denmark, 
          fill="white", colour = "black") + 
  coord_sf() + 
  geom_point(data = landuse.map, 
             aes(x=lat, y = long, colour = "darkgrey"), size=4, show.legend = F) + theme_void() + scale_colour_manual(values = "darkgrey") + scalebar(denmark, dist = 25, dist_unit = "km", transform = T, model = "WGS84", st.size = 3) + labs(subtitle = "A") + north(denmark, symbol = 4, scale = 0.07) + theme(plot.subtitle = element_text(face = "bold", size = 20), plot.margin = margin(0, 0, 0, 0, "cm")) + panel_border()

ggsave("plots/Sampling_map_DK.jpg", height = 10, width = 12, dpi = 600) # remember to increase DPI for publication

### adding regions to coordinates ####
regiondata <- 
  data.frame(sampleID = data$SampleID,
             lat = dist.location$data.utm_x,
             long = dist.location$data.utm_y)

library(geonames)
library(ISOcodes)
library(purrr)

options(geonamesUsername="hvalli_2002")

my_ports_countrycode <- map2_dfr(
  regiondata[["lat"]], 
  regiondata[["long"]], 
  .f = GNcountryCode, 
  radius = 10
)

### analysis  #####
#how many routes
length(unique(data$RouteID))

# how many samples
length(unique(data$SampleID)) # total/complete samples
#length(unique(data$SampleID_size)) # size sorted sample IDs

# how many pilots?
length(unique(data$PID))

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

# match columns prior to merge
dkarter <- dkarter %>% dplyr::rename("order" = Orden) # first rename column header
dkarter$filter <- "Denmark"

# merge data for plotting
ordercomp <- bind_rows(allinsects, insects99, uniquenames, dkarter)

# which insect orders were not detected with the car net
setdiff(dkarter$order, allinsects$order)

### plot stacked bar plot of insect orders ###########

plot <- ordercomp %>% mutate(filter = fct_relevel(filter, 
                                                  "classInsecta", "classInsecta99", "uniquenames", "Denmark")) %>% ggplot(aes(filter, value))

stacked_plot <-
  plot + geom_bar(
    stat = "identity",
    aes(fill = order),
    show.legend = T,
    position = "fill"
  ) + labs(
    x = "",
    y = "Relative ASV/species richness\n",
    fill = "Insect order",
    subtitle = "A"
  ) + theme_classic2() + scale_fill_manual(values = iwanthue_19) + theme(
    title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) + guides(fill = guide_legend(nrow = 3)) + scale_x_discrete(labels=c("classInsecta" = "class = Insecta", "classInsecta99" = "class = Insecta & match = >99%", "uniquenames" = "Unique names"))# position = "fill"in geom_bar gives relative

ggsave("plots/barplots.jpg", height = 10, width = 12, dpi = 600) 

addSmallLegend <- function(myPlot, pointSize = 1, textSize = 7, spaceLegend = 0.7) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"), plot.margin = margin(0, 0, 0, 0, "cm"))
}

stacked_plot <-
  plot + geom_bar(
    stat = "identity",
    aes(fill = order),
    show.legend = T,
    position = "fill"
  ) + labs(
    x = "",
    y = "Relative ASV/species richness\n",
    fill = "Insect order",
    subtitle = "B"
  ) + theme_minimal() + scale_fill_manual(values = iwanthue_19) + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0, "cm"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16, face = "bold")
  ) + guides(fill = guide_legend(nrow = 3)) + scale_x_discrete(
    labels = c(
      "classInsecta" = "Class = \nInsecta",
      "classInsecta99" = "Class = Insecta & \nmatch = >99%",
      "uniquenames" = "Unique \nnames"
    )
  )# position = "fill"in geom_bar gives relative

ggsave("plots/Fig1B_Stacked_bars_DK.jpg", height = 10, width = 12, dpi = 600) # remember to increase DPI for publication

# Apply on original plot
small_leg_plot <- addSmallLegend(stacked_plot)

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

### DOES NOT WORK###########
# relative OTU abundance per sampling method
group.abundance.meta(data = list(reform.data = reform.data.test), top = 30, rank = "o", drop.unclassified = TRUE, meta = metadata, meta.factor = "LandUSeType")  # the sampling strategies does not capture similar taxa

# core taxa - to detect swarming insects on specific days
core.OTU(data = list(reform.data = reform.data.test), meta = metadata, meta.factor="LandUSeType", percent=0.5) # This function returns a list showing otus that present in a pre-defined percent of samples in eachlevel of a given metadata category

core.OTU.rank(data = list(reform.data = reform.data.test), rank="s", drop.unclassified=TRUE, meta = metadata, meta.factor="samplingmethod", percent=0.25) #This function returns a list showing otus that present in a pre-defined percent of samples in eachlevel of a given metadata category.

