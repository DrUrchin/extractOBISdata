# R tool to extract data and plot data from the OBIS database (OBIS.org)
# Author: Narimane Dorey - narimane.dorey_at_gmail.com
# The initial idea is from Corentin Clerc - corentin.clerc_at_lmd.ens.fr

# 1. load the packages you need ####

list.of.packages<-c("ggplot2", #"ggExtra",
                    "ggspatial", "tidyverse",
                    "sf",               # Support for simple features, a standardized way to encode spatial vector data. 
                    #"rgdal",           #removed from CRAN
                    "rworldmap", "maps", "mapproj", # The legacy packages maptools, rgdal, and rgeos, underpinning the sp package, which was just loaded, will retire in October 2023.
                    "rnaturalearth","rnaturalearthdata", # Support for Spatial objects (`sp`) will be deprecated in {rnaturalearth}
                    "raster",           # Reading, writing, manipulating, analyzing and modeling of spatial data.
                    "robis",           # Client for the Ocean Biodiversity Information System
                    "viridis")         # Color palettes


# Do we need to install new packages?
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) {
  print("new packages to install")}else{print("no new packages needed")}

# Install new packages when needed
if(length(new.packages)) install.packages(new.packages)

# Library() for all the packages
lapply(list.of.packages, library, character.only=TRUE)


# 2. Download data from OBIS database ####

# Enter a taxon (species, family, order...)
taxon1="Lophelia pertusa"
taxon2="Madrepora oculata"
taxon3="Solenosmilia variabilis"
taxon4="Bathelia candida"
taxon5="Enallopsammia profunda"

# Load all occurences of your taxon from OBIS
# ?robis::occurence #finds occurence by species
# This is quite long (especially for Lophelia)
dataset1=occurrence(taxon1)
dataset2=occurrence(taxon2)
dataset3=occurrence(taxon3)
dataset4=occurrence(taxon4)
dataset5=occurrence(taxon5)


# 3. Reduce the dataset to bare minimum ####
# Make a function to reduce the data to only what I am interested in (X, Y, Z dimensions)
reduceDF<-function(dataset){
  data.frame(dataset$scientificName, 
             dataset$genus,dataset$decimalLatitude,
             dataset$decimalLongitude, dataset$depth, 
             dataset$maximumDepthInMeters, dataset$minimumDepthInMeters)
  }

# Use reduceDF to create smaller dataframes
Dfreduit1=reduceDF(dataset1)
Dfreduit2=reduceDF(dataset2)
Dfreduit3=reduceDF(dataset3)
Dfreduit4=reduceDF(dataset4)
Dfreduit5=reduceDF(dataset5)

# Make a single dataframe with all the 5 species
dataset=rbind(Dfreduit1, Dfreduit2, 
              Dfreduit3, Dfreduit4, Dfreduit5)

# Change the names of columns
names(dataset)=c("Sp_Name", "Genus", 
                 "DecLatitude", "DecLongitude",
                 "Depth_m", "maxDepth_m", "minDepth_m")

summary(dataset)

# ############################## #
# __ Load data for CWC         ####
# ############################## #

###### ____Offline option #####
# Use the file already written to do more plots:
# FiveCWC=read.csv(file="Five_CWC.csv")

# Remove the rows that do no have any depth
# FiveCWC2=FiveCWC[!is.na(FiveCWC$Depth_m),]

###### ____Online option #####
# Remove the rows that do not have any depth
FiveCWC2=dataset[!is.na(dataset$Depth_m),]

# How many data were removed ?
nrow(FiveCWC2)-nrow(dataset)
summary(FiveCWC2)

# Define number of observations for later use in graphs
B=nrow(FiveCWC2)
# How many observations are lower than 1000 m ?
A=nrow(FiveCWC2[FiveCWC2$Depth_m>=1000,])/B*100
A 
# 26 % of observations that are lower than 1000 m


# 4. Make means of observations by 1° X 1 ° grid ####

# Gather observations by 1°, round the coordinates
FiveCWC2$RoundLon=round(FiveCWC2$DecLongitude,0)
FiveCWC2$RoundLat=round(FiveCWC2$DecLatitude,0)

# Create coordinate couples for each 1° x 1° grid
FiveCWC2$Coord=as.factor(paste(FiveCWC2$RoundLon,FiveCWC2$RoundLat))

# Make one mean per species and 1° x 1° grid
MeanDepth2=data.frame(
  aggregate(Depth_m~RoundLat*RoundLon*Genus,
            data=FiveCWC2, mean),                #mean depth of occurences
  "sd"=aggregate(Depth_m~RoundLat*RoundLon*Genus,
                 data=FiveCWC2, sd)$Depth_m,     #standard deviation of the depth
  "n"=aggregate(Depth_m~RoundLat*RoundLon*Genus,
                data=FiveCWC2, length)$Depth_m)  # number of occurences in that 1° x 1° grid

# rename lat and long for consistency with shp file
names(MeanDepth2)[which(names(MeanDepth2) == 'RoundLon')] <- 'long'
names(MeanDepth2)[which(names(MeanDepth2) == 'RoundLat')] <- 'lat'

# import world map from ggplot2::map_data
world_map <- map_data("world")


# 5. Make the map for each coral ####

# Create a function for the map
map_species_obs_depth= function(G) {
  p4 = ggplot(world_map,
              aes(x = long, y = lat, group = group)) +
    geom_polygon(fill="lightgray", colour = "lightgray") +
    geom_point(data = MeanDepth2[MeanDepth2$Genus==G,],
               aes(x = long, y = lat,
                   size=n, fill=Depth_m, group = NULL),
               shape = 21, col="#ffffff40") +
    theme_minimal() +
    labs(title=bquote("CWC: "*.(G)*" observation map (OBIS dataset)"),
    subtitle="Location, number (n) and mean depth
    (m) of observations per grid of 1°x 1°",
    x ="longitude (°)", y = "latitude (°)") +
    scale_fill_viridis("Mean depth (m)
", option="H", # option = chosen palette
                       limits=c(0,3000),
                       direction=-1,
                       alpha=0.7, discrete = F) +
    scale_size_continuous("Observations (n)
", limits=c(0,10000),
                          breaks=c(10,100,1000,10000))

  p4 + theme(legend.key =
               element_rect(fill = "gray15",
                            color = NA),
             #legend.direction = "horizontal",
             #legend.position = c(0.2, 0.1)
             legend.direction = "vertical",
             legend.position = "right")
  }

# Plot the maps
map_species_obs_depth(G="Desmophyllum") + 
  coord_map(xlim=c(-180,180))  # allows the map to stay at the right dimension whatever the data
map_species_obs_depth(G="Madrepora") + 
  coord_map(xlim=c(-180,180))
map_species_obs_depth(G="Solenosmilia") + 
  coord_map(xlim=c(-180,180)) 
map_species_obs_depth(G="Enallopsammia") + 
  coord_map(xlim=c(-180,180)) 
map_species_obs_depth(G="Bathelia") + 
  coord_map(xlim=c(-180,180))



#Save the maps in your current repertory
G="Desmophyllum"          # Change accordingly
ggsave(filename=paste0(G,"points_map.png"),
       plot=map_species_obs_depth(G) + coord_map(xlim=c(-180,180)),
       device="png",
       width = 30, height = 15, units = "cm")





# II. Another alternative with leaflet ####
# Adapted from : https://www.r-graph-gallery.com/19-map-leafletr.html
# Library
library(leaflet)

G="Desmophyllum"          # Change genus accordingly
G="Solenosmilia"
G="Enallopsammia"
G="Madrepora"

# Create a color palette with handmade bins
max(MeanDepth2[MeanDepth2$Genus==G,]$Depth_m)

mybins=seq(0, 1500, by=500) # Change depending on species
mypalette=colorBin(palette="YlOrBr", 
                   domain=MeanDepth2[MeanDepth2$Genus==G,]$Depth_m,
                   na.color="transparent", 
                   bins=mybins)

# Prepare the text for the tooltip:
mytext=paste(
   "Depth = ", round(MeanDepth2[MeanDepth2$Genus==G,]$Depth_m), 
   " m, ",  
   "n = ", MeanDepth2[MeanDepth2$Genus==G,]$n,sep="",
   "observations")

# Reduce the extreme n for the plot
MeanDepth2$n2 = MeanDepth2$n; max(MeanDepth2$n2)

# This works for Desmophyllum at least
MeanDepth2$n2[MeanDepth2$n2 >= 6 & 
                MeanDepth2$n2 <= 29] <- 7.5

MeanDepth2$n2[MeanDepth2$n2 >= 30 & 
                MeanDepth2$n2 <= 99] <- 10

MeanDepth2$n2[MeanDepth2$n2 >= 100 & 
                MeanDepth2$n2 <= 199] <- 15

MeanDepth2$n2[MeanDepth2$n2 >= 200 & 
                MeanDepth2$n2 <= 699] <- 25

MeanDepth2$n2[MeanDepth2$n2 >= 700 & 
                MeanDepth2$n2 <= 3999] <- 50

MeanDepth2$n2[MeanDepth2$n2 >= 4000] <- 100 # Florida


# MeanDepth2[which.max(MeanDepth2$n),]

# Final Map
m <- leaflet(MeanDepth2[MeanDepth2$Genus==G,]) %>% 
  addTiles()  %>% 
  setView(lat=40, lng=0 , zoom=3) %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(~long, ~lat, 
    fillColor = ~mypalette(Depth_m), 
    fillOpacity = 0.9, color="white", 
    radius=~n2/2, stroke=FALSE,
    label = mytext,
    labelOptions = labelOptions( style = 
                                list("font-weight" ="normal", 
                                     padding = "3px 8px"), 
                                textsize = "13px", 
                                direction = "auto")
  ) %>%
  addLegend( pal=mypalette, 
             values=~Depth_m, 
             opacity=0.9, 
             title = "Depth (m)", 
             position = "bottomright" )

m 



# save the widget in a html file if needed. to
# library(htmlwidgets)
# saveWidget(m, file=paste0( getwd(), "/HtmlWidget/bubblemapQuakes.html"))
