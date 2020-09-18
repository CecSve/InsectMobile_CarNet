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

### colour scheme ##################
landuseCols <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73") # colour friendly, ordered by land cover 
landuseOrder <- c("Urban","Farmland","Bog","Forest")

### load data #####################
asvs <- read.delim("cleaned-data/asvtable.txt",sep=" ")
data <- read.delim("cleaned-data/data.txt",sep="\t")
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

### renaming observations for standardized data input #########

# temperature
table(data$Temperature..celsius.) # none above 25 degrees celsius

data <- data %>% rename("Temperature" = Temperature..celsius.) # first rename column header
data$Temperature[data$Temperature %in% c("18", "19")] <- "15-20" # change observations to be a range
data$Temperature[data$Temperature %in% c("20", "23")] <- "20-25" # change observations to be a range - 20 will be placed in 20-25
table(data$Temperature) # two levels

# wind
table(data$Wind..m.s.)

data <- data %>% rename("Wind" = Wind..m.s.) # first rename column header
data$Wind[data$Wind %in% c("3", "Svag")] <- "Low" # change observations to be in English and three levels
data$Wind[data$Wind %in% c("Middel")] <- "Weak" 
data$Wind[data$Wind %in% c("7", "8", "6")] <- "Strong" 
table(data$Wind) # three levels, well balanced

#land use type
table(data$LandUseType)

data$LandUseType[data$LandUseType %in% c("AltBebyggelse")] <- "Urban" # change observations to be in English 
data$LandUseType[data$LandUseType %in% c("Mark")] <- "Farmland" 
data$LandUseType[data$LandUseType %in% c("Mose")] <- "Bog"
data$LandUseType[data$LandUseType %in% c("Skov")] <- "Forest"

# format date
data$yDay <- yday(data$Date)
#centering
data$cyDay <- data$yDay - median(data$yDay)

# format time
#sort time data to standard each around the time band
data$numberTime <- as.numeric(hms(data$TimeStart))#Denmark
#transform to minutes
data$numberTime <- data$numberTime/60 

table(data$TimeStart) # sampling occurred over 10 hrs
#order time band
table(data$Time_band) # had to manually add time_band since I couldn't get lubridate to work. Based on TimeStart it assigned in 4 categories: 12-15 = midday, 15-18 = afternoon, 18-20 = evening, 20-22+ = late evening, the shorter time periods in the evening is due to the clear shift in community composition we saw from the pilot study (lots and lots of mosquitoes come out after 20)

### preliminary data exploration #####
table(data$Date)
table(data$Route)
table(data$LandUseType)

### sampling map #######

#order habitat
data$Land_use <- factor(data$LandUseType,levels=c("Urban","Farmland",
                                               "Bog","Forest"))

# extract Latitude and Longitude variables and put them into simple data frame called lat.long.df.
lat.long.df <- data.frame(data$Lat, data$Long) 
str(lat.long.df)

coordinates(lat.long.df) <-  ~data.Lat + data.Long #Function coordinates() creates a spatial object

#str(lat.long.df) # at this point, this dataset doesn’t have CRS. Spatial data can’t exist without CRS(Coordinates Reference System). Spatial data requires, at least, coordinates and CRS.

proj4string(lat.long.df) <- CRS("+init=epsg:25832") #this is +proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
head(lat.long.df)

landuse.map <- 
  data.frame(Landuse = data$Land_use,
             lat = lat.long.df$data.Lat,
             long = lat.long.df$data.Long)

# DK map
worldmap <- ne_countries(scale = 'large', type = 'map_units',
                         returnclass = 'sf')

denmark <- worldmap[worldmap$name == 'Denmark',]
ggplot() + geom_sf(data = denmark) + theme_bw()

denmark_cropped <- st_crop(denmark, xmin = 8, xmax = 13,
                           ymin = 54.5, ymax = 58)
ggplot() + geom_sf(data = denmark_cropped) + theme_bw()

landuse.map %>% head()
main.landuse.map <- landuse.map[landuse.map$long < 13, ]
main.data <- data[data$Lat < 58, ]

# plot with whole of Denmark
denmark %>%
  ggplot() + 
  geom_sf(data = denmark, 
          fill="white", colour = "black") + 
  coord_sf() + 
  geom_point(data = landuse.map, 
             aes(x=long, y = lat, colour = data$Land_use), alpha = 0.5, shape = 19, size=5, show.legend = T) + theme_void() + scale_colour_manual("Predominant land cover", values = landuseCols) + scalebar(denmark, dist = 25, dist_unit = "km", transform = T, model = "WGS84", st.size = 3) + north(denmark, symbol = 4, scale = 0.07) + panel_border() + labs(subtitle = "Sampling areas in Denmark") + theme(plot.subtitle = element_text(face = "bold", size = 15))  

# plot without Bornholm
denmark_cropped %>%
  ggplot() + 
  geom_sf(data = denmark_cropped, 
          fill="white", colour = "black") + 
  coord_sf() + 
  geom_point(data = main.landuse.map, 
             aes(x=long, y = lat, colour = data$Land_use), shape = 19, size=6, show.legend = T) + theme_void() + scale_colour_manual("Predominant \nland cover", values = landuseCols) + scalebar(denmark_cropped, dist = 25, dist_unit = "km", transform = T, model = "WGS84", st.size = 3, location = "bottomleft") + north(denmark_cropped, symbol = 4, scale = 0.07) + panel_border() + labs(subtitle = "Sampling areas in Denmark") + theme(plot.subtitle = element_text(face = "bold", size = 15))  

ggsave("plots/Sampling_map_DK.png", height = 10) # remember to increase DPI for publication

### 

### analysis  #####

### insect order ############
taxon <- dplyr::select(taxonomy, "order") # or choose other taxonomic levels (but remember morphology only has order level for all individuals)
otus <- decostand(asvs, "pa") # make asvtable into presence absence

# To be able to merge the OTU table with the taxonomy table, they need to have a common column to call
otus <- rownames_to_column(otus, var = "otuid")
taxon <- rownames_to_column(taxon, var = "otuid") # remeber to remake the column into rownames for both datasets if you need to

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
1632/2350 # diptera
252/2350 # hymenoptera
224/2350 #coleoptera
140/2350 # hemiptera
38/2350 # Thysanoptera

plot_data <- merge(data, longdata, by.x = "PCRID", by.y = "variable")
plot_data <- as_tibble(plot_data)
plot_data <- plot_data %>% dplyr::rename(richness = value)

plot_data %>% dplyr::group_by(Time_band, LandUseType, order) %>% dplyr::summarise(richness = sum(richness)) # oturichness by time band, habitat and insect order

plot_data %>% dplyr::group_by(Time_band, order) %>% dplyr::summarise(richness = sum(richness))

# Psocodea needs to be renamed to Psocoptera
plot_data <- plot_data%>% 
  mutate(order = as.character(order)) %>% 
  mutate(order = replace(order, order == 'Psocodea', 'Psocoptera'))

custom_11 <- c("#771155", "#AA4488", "#EA6CC0", "#CC99BB", "#114477", "#4477AA","#1E78D2", "#77AADD", "#117777", "#44AAAA", "#3FE4E4")

custom_other = c("#77CCCC", "#117744","#44AA77", "#1ED278", "#88CCAA", "#771122", "#AA4455", "#D21E2C","#DD7788","#777711", "#AAAA44", "#D2D21E", "#DDDD77","#774411", "#AA7744", "#D2781E", "#DDAA77")

hist(plot_data$richness) # non-normal and log and sqrt does not resolve this

str(plot_data)

# Stacked bar chart oof OTU richness
g <- ggplot(plot_data, aes(LandUseType, richness))
l <- ggplot(plot_data, aes(Time_band, richness))
s <- ggplot(plot_data, aes(Route, richness))

p1 <- g + geom_bar(stat = "identity", aes(fill = order), show.legend = F, position = "fill") + labs(x = "Land cover", y = "", fill = "Insect order", title = "A") + theme_classic2() + scale_fill_manual(values = c("red", viridis::viridis(12))) + theme(title = element_text(face = "bold")) # position = "fill"in geom_bar gives relative

p2 <- l + geom_bar(stat = "identity", aes(fill = order), show.legend = T, position = "fill") + labs(x = "Time band", y = "", fill = "Insect order", title = "B") + theme_classic2() + scale_fill_manual(values = c("red", viridis::viridis(12))) + theme(title = element_text(face = "bold")) + guides(fill = guide_legend(ncol = 6))

p3 <- s + geom_bar(stat = "identity", aes(fill = order), show.legend = F, position = "fill") + labs(x = "Routes", y = "", fill = "Insect order", title = "C") + theme_classic2() + scale_fill_manual(values = c("red", viridis::viridis(12))) + theme(axis.text.x = element_blank(), title = element_text(face = "bold"), axis.ticks.x = element_blank())

figure <- ggarrange(p1,p2,p3, common.legend = T, legend = "bottom", nrow = 3, align = "hv")

# Annotate the figure by adding a common labels
fig1 <- annotate_figure(figure,
                left = text_grob("Relative richness", rot = 90, face = "bold")
)

save_plot("plots/fig1_richness_landcover_timeband_routes.png", fig1, base_height = 12, base_width = 10)

### insect family ######
taxon <- dplyr::select(taxonomy, "family") # or choose other taxonomic levels (but remember morphology only has order level for all individuals)
otus <- decostand(asvs, "pa") # make asvtable into presence absence

# To be able to merge the OTU table with the taxonomy table, they need to have a common column to call
otus <- rownames_to_column(otus, var = "otuid")
taxon <- rownames_to_column(taxon, var = "otuid") # remeber to remake the column into rownames for both datasets if you need to

taxon_data <- left_join(otus, taxon, by = "otuid")
taxon_data <- column_to_rownames(taxon_data, var = "otuid")

# working with reshape2 to wrangle data from wide to long format
longdata <- melt(taxon_data) 
#longdata <- longdata[longdata$value!=0,]
longdata <- aggregate(. ~ family + variable, data = longdata, FUN = sum) # if values in the insect order column and the sample (variable) column are identical, then sum up how many unique otus there were in the sample from that insect order (richness)

longdata2 <- longdata %>% group_by(family) %>% 
  mutate(value_mean = mean(value)) %>% 
  arrange(value_mean) %>% ungroup() 

longdata2$family <- factor(longdata2$family, levels=unique(longdata2$family))
str(longdata2)

p1 <- ggplot(data = longdata2, aes(x = family , y = variable , fill = value)) 
family_raster_plot <- p1 + geom_raster() + scale_fill_gradient(low="midnightblue", high="darkgoldenrod2") + 
  theme_bw() +
  labs(title = "B", fill = "Unique \nASVs") + xlab("Insect family") +
  theme_minimal_grid() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3, hjust = 1),
                               axis.text.y=element_blank(),
                               axis.ticks = element_blank(),
                               axis.title.y = element_blank(),
                               axis.title.x = element_text(face = "bold"),
                               plot.title=element_text(face = "bold")) + scale_x_discrete(limits = rev(levels(longdata2$family)))

# how frequent is the different orders detected
sum(longdata2$value) # 2116
aggregate(longdata2$value, by=list(family=longdata2$family), FUN=sum) # 126 families
40/126 # so many has only one unique ASV in the family
293/2116 # Chironomidae 
183/2116 # Cecidomyiidae 
155/2116 #Ceratopogonidae 
146/2116 # Staphylinidae 
143/2116 # Phoridae 
100/2116 # most have under hundred species, what does thta equal
119/126

p1 <- plot_grid(order_raster_plot, family_raster_plot, ncol = 1, align = "hv", rel_heights = 1, rel_widths = 2)

save_plot("plots/heatmap_order_family_unique_ASVs.png", p1, base_height = 16, base_width = 16)
