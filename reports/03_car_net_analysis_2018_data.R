# all DK 2018 data with fwh primer
# car net analysis

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

# visualisation
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
library(grid)
library(gridExtra)

### colour scheme ##################
landuseCols <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73") # colour friendly, ordered by land cover 
landuseOrder <- c("Urban","Farmland","Bog","Forest")

### load data #####################
asvs <- read.delim("cleaned-data/DK_asvtable.txt", sep="\t")
data <- read.delim("cleaned-data/DK_rough_landuse_biomass.txt",sep=" ")
taxonomy <- read.delim("cleaned-data/DK_taxonomy.txt",sep="\t")

# load centroid coordinates for each route
coords <- read.delim("cleaned-data/DK_ruter2018_pkt_koordinater.txt")

### aligning data for analysis ###################
names(asvs) = gsub(pattern = "X*", replacement = "", x = names(asvs)) # remove the Xs that have emerged in the column headers that are numerical
asvs <- asvs[, grepl("IM18_*", names(asvs))] # only keep Dnaish samples from 2018 - be aware this removes blanks and negatives as well
names(asvs)

keep <- data$PCRID
asvs <- asvs[, (names(asvs) %in% keep)] # subset asvtable to only contain samples with metadata

keep <- names(asvs)
data <- data %>% filter(PCRID %in% keep)

# do all samples have seuqnces?
colSums(asvs) # yes
#asvs <- asvs[,colSums(asvs) > 0]

# remove empty asvs
rowSums(asvs)
asvs <- asvs[rowSums(asvs) > 0, ]

# remove the asvs from the taxonomy
keep <- rownames(asvs)
taxonomy <- taxonomy[rownames(taxonomy) %in% keep, ]

# rmeove asvs that does not fit with taxonomy of class = insecta, 99% identity score threshold etc
keep <- rownames(taxonomy)
asvs <- asvs[rownames(asvs) %in% keep, ]

# do all samples have seuqnces?
colSums(asvs) # no
asvs <- asvs[,colSums(asvs) > 0]

# subset data to the samples with reads
keep <- colnames(asvs)
data <- data %>% filter(PCRID %in% keep)

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

## haven't done this! rewrite to fit
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
  data.frame(Landuse = data$roughLand_use,
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
             aes(x=lat, y = long, colour = "darkgrey"), size=4, show.legend = F) + theme_void() + scale_colour_manual(values = "deepskyblue1") + scalebar(denmark, dist = 25, dist_unit = "km", transform = T, model = "WGS84", st.size = 3) + labs(subtitle = "A") + north(denmark, symbol = 4, scale = 0.07) + theme(plot.subtitle = element_text(face = "bold", size = 20), plot.margin = margin(0, 0, 0, 0, "cm")) + panel_border()

#ggsave("plots/Sampling_map_DK.png", height = 10) # remember to increase DPI for publication

### analysis  #####
#how many routes
length(unique(data$RouteID))

# how many samples
length(unique(data$SampleID))

### insect order ############
taxon <- dplyr::select(taxonomy, "order") # or choose other taxonomic levels (but remember morphology only has order level for all individuals)
otus <- decostand(asvs, "pa") # make asvtable into presence absence

# To be able to merge the OTU table with the taxonomy table, they need to have a common column to call
otus <- rownames_to_column(otus, var = "otuid")
taxon <- rownames_to_column(taxon, var = "otuid") # remeber to remake the column into rownames for both datasets if you need to

test <- otus %>% column_to_rownames(var = "otuid") 
min(colSums(test))

taxon_data <- left_join(otus, taxon, by = "otuid")
taxon_data <- column_to_rownames(taxon_data, var = "otuid")

# working with reshape2 to wrangle data from wide to long format
longdata <- melt(taxon_data) 
#longdata <- longdata[longdata$value!=0,]
longdata <- aggregate(. ~ order + variable, data = longdata, FUN = sum) # if values in the insect order column and the sample (variable) column are identical, then sum up how many unique otus there were in the sample from that insect order (richness)

longdata2 <- longdata %>% group_by(order) %>% 
  mutate(value_mean = mean(value)) %>% 
  arrange(value_mean) %>% ungroup() 

longdata2$order <- factor(longdata2$order, levels=unique(longdata2$order))
str(longdata2)

p1 <- ggplot(data = longdata2, aes(x = order , y = variable , fill = value)) 
order_raster_plot <- p1 + geom_raster() + scale_fill_gradient(low="midnightblue", high="darkgoldenrod2") + 
  labs(title = "A", fill = "Unique \nASVs") + xlab("Insect order") +
  theme_minimal_grid() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3, hjust = 1),
                               axis.text.y=element_blank(),
                               axis.ticks = element_blank(),
                               axis.title.y = element_blank(),
                               axis.title.x = element_text(face = "bold"),
                               plot.title=element_text(face = "bold")) + scale_x_discrete(limits = rev(levels(longdata2$order)))

# how frequent is the different orders detected
sum(longdata2$value)
aggregate(longdata2$value, by=list(order=longdata2$order), FUN=sum)
23280/35646 # diptera
4802/35646 # coleoptera
2969/35646 #hymenoptera
3111/35646 # hemiptera
1014/35646 # Thysanoptera

plot_data <- merge(data, longdata, by.x = "PCRID", by.y = "variable")
plot_data <- as_tibble(plot_data)
plot_data <- plot_data %>% dplyr::rename(richness = value)

plot_data %>% dplyr::group_by(Time_band, order) %>% dplyr::summarise(richness = sum(richness)) # richness by time band and insect order

# Psocodea needs to be renamed to Psocoptera
plot_data <- plot_data%>% 
  mutate(order = as.character(order)) %>% 
  mutate(order = replace(order, order == 'Psocodea', 'Psocoptera'))

custom_11 <- c("#771155", "#AA4488", "#EA6CC0", "#CC99BB", "#114477", "#4477AA","#1E78D2", "#77AADD", "#117777", "#44AAAA", "#3FE4E4")

custom_other = c("#77CCCC", "#117744","#44AA77", "#1ED278", "#88CCAA", "#771122", "#AA4455", "#D21E2C","#DD7788","#777711", "#AAAA44", "#D2D21E", "#DDDD77","#774411", "#AA7744", "#D2781E", "#DDAA77")

hist(log(plot_data$richness)) # non-normal and log and sqrt does not resolve this

str(plot_data)

# Stacked bar chart of ASV richness
l <- ggplot(plot_data, aes(Time_band, richness))
s <- ggplot(plot_data, aes(SampleID, richness))

p1 <- l + geom_bar(stat = "identity", aes(fill = order), show.legend = T, position = "fill") + labs(x = "Time band", y = "", fill = "Insect order", title = "B") + theme_classic2() + scale_fill_manual(values = c("red", viridis::viridis(15))) + theme(title = element_text(face = "bold")) + guides(fill = guide_legend(ncol = 6))

p2 <- s + geom_bar(stat = "identity", aes(fill = order), show.legend = F, position = "fill") + labs(x = "Sample", y = "", fill = "Insect order", title = "C") + theme_classic2() + scale_fill_manual(values = c("red", viridis::viridis(15))) + theme(axis.text.x = element_blank(), title = element_text(face = "bold"), axis.ticks.x = element_blank())

figure <- ggarrange(p1,p2, common.legend = T, legend = "bottom", nrow = 2, align = "hv")

# Annotate the figure by adding a common labels
fig1 <- annotate_figure(figure,
                        left = text_grob("Relative richness", rot = 90, face = "bold")
)

save_plot("plots/fig1_richness_landcover_timeband_routes.png", fig1, base_height = 12, base_width = 10)

### insect order - Denmark comparison ##############
# load species data on insects from Denmark (downloaded manually, not updated since May 27th 2020)

allearter <- read_csv("raw-data/allearter_19092020.csv")

dkarter <-  allearter %>% group_by(Orden) %>% dplyr::summarise(value = n()) # count how many species are known for each order

# remove non-flying insect orders
remove <- c("Phthiraptera", "Siphonaptera", "Zygentoma", "Microcoryphia")
dkarter <-  dkarter %>% dplyr::filter(!Orden %in% remove)

# prepare study data
taxon <- taxonomy %>% drop_na(species) # remove the obeservation were the ASV was not assigned species level taxonomy
taxon <- dplyr::select(taxonomy, "order") # or choose other taxonomic levels (but remember morphology only has order level for all individuals)
otus <- decostand(asvs, "pa") # make asvtable into presence absence

# To be able to merge the OTU table with the taxonomy table, they need to have a common column to call
otus <- rownames_to_column(otus, var = "otuid")
taxon <- rownames_to_column(taxon, var = "otuid") # remeber to remake the column into rownames for both datasets if you need to

taxon_data <- left_join(otus, taxon, by = "otuid")
taxon_data <- column_to_rownames(taxon_data, var = "otuid")

# working with reshape2 to wrangle data from wide to long format
longdata <- melt(taxon_data) 
longdata <- longdata[longdata$value!=0,] # remove zeros in this case
longdata <- aggregate(. ~ order + variable, data = longdata, FUN = sum) # if values in the insect order column and the sample (variable) column are identical, then sum up how many unique otus there were in the

# Psocodea needs to be renamed to Psocoptera
longdata <- longdata%>% 
  mutate(order = as.character(order)) %>% 
  mutate(order = replace(order, order == 'Psocodea', 'Psocoptera'))

# match columns prior to merge
dkarter <- dkarter %>% rename("order" = Orden) # first rename column header
dkarter$variable <- NA
dkarter$origin <- "Denmark"
longdata$origin <- "Car net"

str(longdata)
dkarter <- dkarter %>% dplyr::select(order, variable, value, origin)

ordercomp <- rbind(longdata, dkarter)

# which insect orders were not detected with the car net
setdiff(dkarter$order, longdata$order)

### plot stacked bar plot of insect orders ###########
c <- ggplot(ordercomp, aes(origin, value))
stacked_plot <- c + geom_bar(stat = "identity", aes(fill = order), show.legend = T, position = "fill") + labs(x = "", y = "Relative abundance", fill = "Insect order", subtitle = "A") + theme_classic2() + scale_fill_manual(values = c("red", viridis::viridis(20))) + theme(title = element_text(), legend.position = "bottom",axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) + guides(fill = guide_legend(nrow = 2))# position = "fill"in geom_bar gives relative

addSmallLegend <- function(myPlot, pointSize = 1, textSize = 7, spaceLegend = 0.7) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"), plot.margin = margin(0, 0, 0, 0, "cm"))
}

stacked_plot <- c + geom_bar(stat = "identity", aes(fill = order), show.legend = T, position = "fill") + labs(x = "", y = "Relative species richness", fill = "Insect order", subtitle = "B")  + theme_minimal() + scale_fill_manual(values = c("red", viridis::viridis(20))) + theme(plot.subtitle = element_text(size = 20, face = "bold"), legend.position = "bottom", plot.margin = margin(0, 0, 0, 0, "cm"), axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8)) + guides(fill = guide_legend(nrow = 3))# position = "fill"in geom_bar gives relative

# Apply on original plot
small_leg_plot <- addSmallLegend(stacked_plot)

#le1 <- cowplot::get_legend(small_leg_plot)

# + scale_fill_manual(values = c("red", viridis::viridis(20)))
# + scale_fill_viridis_d(option = "plasma") 

### new species for DK? ###############
# how many species were detected by the net and is not present in the Danish species database?
newspecies <- setdiff(taxonomy$species, allearter$`Videnskabeligt navn`) # the first observation is NA, so 338 species not in the Danish database

# comparison of car net and DK database
newspeciesDK <- taxonomy %>% rownames_to_column(var = "ASVID") %>% filter(species %in% newspecies) 
#write.table(newspeciesDK, file = "cleaned-data/newspecies_DK_alldata.txt", sep = "\t", row.names = F)
#newspeciesDK <- newspeciesDK %>% column_to_rownames(var = "ASVID")
newdistinct_species <- newspeciesDK %>% distinct(species, .keep_all = T)%>% drop_na(species)
newspeciesDK <- newspeciesDK %>% group_by(order, species) %>% distinct(species) %>% drop_na(species)
#write.table(newspeciesDK, file = "cleaned-data/newspecies_DK.txt", sep = "\t", row.names = F)

### gbif comparison ###############
library(rgbif)

# fill in your gbif.org credentials 
user <- "" # your gbif.org username 
pwd <- "" # your gbif.org password
email <- "" # your email 

library(dplyr)
library(purrr)
library(readr)  
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid_

# match the names 
gbif_taxon_keys <- 
  newspeciesDK %>% 
  pull("species") %>% # use fewer names if you want to just test 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_tsv(path = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
  #filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
  filter(class == "Insecta") %>% # remove anything that might have matched to a non-insect
  pull(usagekey) # get the gbif taxonkeys

# gbif_taxon_keys should be a long vector like this c(2977832,2977901,2977966,2977835,2977863)
# !!very important here to use pred_in!!
# use matched gbif_taxon_keys from above 
occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred("country", "DK"),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
) # this generates a file for your user where the matches are

# use matched gbif_taxon_keys from above 
occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in("country", c("SE", "NO", "DE")),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
) # this generates a file for your user where the matches are

gbif_query <- read.delim("raw-data/newspecies_DK.csv")
gbif_query_scandi <- read.delim("raw-data/newspecies_scandi.csv")

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
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n))*100) # get the proportion in percent

allearter %>%
  group_by(Orden) %>%
  summarise(n = n())

taxonomy %>%
  group_by(order) %>%
  summarise(n = n())

# Lepidoptera 
(2590/18882)*100 # almost 14% of insects in DK are Lepidoptera
butterfies <- taxonomy %>% filter(order == "Lepidoptera") %>% distinct(species)
dkbutterflies <- allearter %>% filter(Orden == "Lepidoptera") %>% distinct(`Videnskabeligt navn`)
setdiff(butterfies$species, dkbutterflies$`Videnskabeligt navn`) # all detected were already known in Denmark
setdiff(dkbutterflies$`Videnskabeligt navn`, butterfies$species)
(74/2590)*100 # we detected less than 3% of the known Lepidoptera in Denmark

(5086/18882)*100 # almost 27% of insects in DK are flies
flies <- taxonomy %>% filter(order == "Diptera") %>% distinct(species)
dkflies <- allearter %>% filter(Orden == "Diptera") %>% distinct(`Videnskabeligt navn`)
setdiff(flies$species, dkflies$`Videnskabeligt navn`) # 192 flies new to DK
setdiff(dkflies$`Videnskabeligt navn`, flies$species)
((2405-192)/5086)*100

### species accumulation curves #########
otus <- otus %>% column_to_rownames(var = "otuid")

## Accumulation model
pool <- poolaccum(otus, permutations = 1000)
plot(pool)

# extract data for ggplot 
accdata <- summary(pool, display = "chao")
data_plot <- as_tibble(accdata[["chao"]])

acummulation_plot <- data_plot %>% ggplot(aes(N, Chao)) +
  geom_line(aes(color = "black"), size = 2) + theme_minimal() + scale_colour_manual(values = "darkgrey") + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm")
  ) + geom_ribbon(
        aes(
          ymin = Chao-Std.Dev,
          ymax = Chao+Std.Dev
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        x = "Number of samples",
        y = "Estimated Richness (Chao)",
        subtitle = "C"
      ) + scale_x_continuous(limits = c(0, 1000)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + scale_fill_manual(values = "lightgrey")

acummulation_plot_zoom <-  data_plot %>% ggplot(aes(N, Chao)) +
  geom_line(aes(color = "black"), size = 2) + theme_minimal_grid() + scale_colour_manual(values = "darkgrey") + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  ) + scale_x_continuous(
    limits = c(0, 100)) + geom_ribbon(
        aes(
          ymin = Chao-Std.Dev,
          ymax = Chao+Std.Dev
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + labs(
        x = "",
        y = ""
        #subtitle = "B",
      ) + scale_fill_manual(values = "lightgrey") + guides(colour = guide_legend(nrow = 1)) + theme(panel.background = element_rect(fill = "white"), plot.margin = margin(0, 0, 0, 0, "cm"), panel.border = element_rect(colour = "darkgrey"))

acummulation_plot
acummulation_plot_zoom

b <- acummulation_plot + annotation_custom(ggplotGrob(acummulation_plot_zoom), xmin = 500, xmax = 1000, ymin = -600, ymax = 250)

#ggsave("plots/zoom_accumulation_chao1.png")

library(gtable)
g1 <- ggplotGrob(mapplot)
g2 <- ggplotGrob(small_leg_plot)
g3 <- ggplotGrob(b)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
grid.newpage()
grid.draw(g)

ggsave('plots/g_test.tiff', plot = g, width=300,height=600, units = "mm", dpi=300)

#figure1 <- ggarrange(mapplot, small_leg_plot, b, ncol = 1, align = "v")
#save_plot("plots/fig1_relabun_estimaterich.png", figure1, base_height = 10, base_width = 8)

### z-test  ##################
# how many unique observations in taxonomy
length(unique(taxonomy$order))
length(unique(taxonomy$family))
length(unique(taxonomy$genus))
length(unique(taxonomy$species))

# in the Dnaish species list
allearter %>% dplyr::filter(!Orden %in% remove) %>% group_by(Orden) %>% dplyr::summarise(value = n()) # count how many species are known for each order
dkarter_z <- allearter %>% dplyr::filter(!Orden %in% remove) # remove non-flying insect orders
length(unique(dkarter_z$Orden))
length(unique(dkarter_z$Familie))
length(unique(dkarter_z$Slægt))
length(unique(dkarter_z$`Videnskabeligt navn`))

# order
res <- prop.test(x = c(15, 19), n = c(19, 19))
# Printing the results
res

# family
res <- prop.test(x = c(215, 485), n = c(485, 485))
# Printing the results
res

# genus
res <- prop.test(x = c(998, 5467), n = c(5467, 5467))
# Printing the results
res

# species
res <- prop.test(x = c(1607, 18791), n = c(18791, 18791))
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
res <- prop.test(x = c(643, 3866), n = c(4546, 18882))
# Printing the results
res 

# Dermaptera
res <- prop.test(x = c(1, 6), n = c(4546, 18882))
# Printing the results
res 

# Diptera
res <- prop.test(x = c(2405, 5086), n = c(4546, 18882))
# Printing the results
res 

# Ephemeroptera
res <- prop.test(x = c(11, 43), n = c(4546, 18882))
# Printing the results
res 

# Hemiptera  
res <- prop.test(x = c(410, 1495), n = c(4546, 18882))
# Printing the results
res 

# Hymenoptera     
res <- prop.test(x = c(855, 5150), n = c(4546, 18882))
# Printing the results
res 

# Lepidoptera      
res <- prop.test(x = c(98, 2590), n = c(4546, 18882))
# Printing the results
res 
# Mecoptera         
res <- prop.test(x = c(2, 4), n = c(4546, 18882))
# Printing the results
res 
# Neuroptera        
res <- prop.test(x = c(2, 62), n = c(4546, 18882))
# Printing the results
res 
# Odonata           
res <- prop.test(x = c(5, 60), n = c(4546, 18882))
# Printing the results
res 
# Orthoptera        
res <- prop.test(x = c(4, 38), n = c(4546, 18882))
# Printing the results
res 
# Plecoptera        
res <- prop.test(x = c(2, 25), n = c(4546, 18882))
# Printing the results
res 
# Psocoptera (Psocodea)
res <- prop.test(x = c(40, 62), n = c(4546, 18882))
# Printing the results
res 
# Thysanoptera     
res <- prop.test(x = c(60, 113), n = c(4546, 18882))
# Printing the results
res 
# Trichoptera   
res <- prop.test(x = c(8, 172), n = c(4546, 18882))
# Printing the results
res 
