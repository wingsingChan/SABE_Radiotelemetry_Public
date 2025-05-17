## Load Packages ----
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(labdsv)
library(sp)
library(sf)
library(sfnetworks)
library(rgeos)
library(geosphere)
library(stars)
library(adehabitatHS)
library(ggplot2)
library(ggmap)
library(tmap)
library(move)
library(ggspatial)
library(randomcoloR)

# Movement Metrics ----
## Load Movement Data ----
movement <- read.csv("data/SABE_movement_records.csv", header=TRUE, na.strings = c("", "NA"))
View(movement)
str(movement)

## Remove Specific Types of Records
## Inc. Release records, Acclimitization records, Missing records, Records of lost tag
movement$behaviour <- as.character(movement$behaviour)
movement <- movement[!(movement$behaviour == "Lost Tag" | movement$behaviour == "Lost Tag - Underwater" | 
                       movement$behaviour == "Missing" | movement$behaviour == "Release" |
                       movement$behaviour == "Acclimitization - Overland" | movement$behaviour == "Acclimitization - Underwater"),]
movement$behaviour <- as.factor(movement$behaviour)
levels(movement$behaviour)
#### 1679 Obs.

## Radio-Frequency
movement$frequency <- as.factor(movement$frequency)

## Turtle Id
movement$turtle.id <- as.factor(movement$turtle.id)
movement %>% 
  dplyr::count(turtle.id)

## Date, Time, Year, Month, Season
movement$date <- parse_date(movement$date, "%d-%b-%y")
movement$time <- paste(as.character(movement$date), as.character(movement$time))
movement$time <- as.POSIXct(movement$time, format = "%Y-%m-%d %H:%M", tz="GMT")
movement$month <- months(movement$time)
movement$month <- factor(movement$month, c("January", "February", "March", "April", "May", "June", 
                                           "July", "August", "September", "October", "November", "December"))
movement$year <- year(movement$time)
movement$year <- as.factor(movement$year)

## Define Season
## Overwintering Nov-Feb (Dec-Feb)
## Nesting Mar-Apr (Mar-May)
## Wet Season May-Aug (Jun-Aug)
## Mating Sept-Oct (Sept-Nov)
movement$season <- ifelse(movement$month %in% c("May", "June", "July", "August"), "Wet Season", 
                          ifelse(movement$month %in% c("September", "October"), "Mating Season", 
                                 ifelse(movement$month %in% c("November", "December"), "Dry Season", 
                                        ifelse(movement$month %in% c("January", "February"), "Dry Season (L)", 
                                               ifelse(movement$month %in% c("March", "April"), "Nesting Season",
                                                      "Unclassified")))))
movement$season <- paste(as.character(movement$year), as.character(movement$season))
movement$season[movement$season == "2021 Dry Season (L)"] <- "2020 Dry Season"
movement$season <- factor(movement$season, c("2020 Wet Season", "2020 Mating Season", "2020 Dry Season", "2021 Nesting Season", 
                                             "2021 Wet Season", "2021 Mating Season", "2021 Dry Season"))
levels(movement$season)

## Elevation
movement$stream.elev <- as.character(movement$stream.elev)
movement$stream.elev[movement$stream.elev == "*"] <- NA
movement$stream.elev <- as.integer(movement$stream.elev)

## Habitat Types
movement$location.habitat <- as.character(movement$location.habitat)
movement$location.habitat[movement$location.habitat=="Pool (Reservoir)"] <- "Pool"
movement$location.habitat <- as.factor(movement$location.habitat)
levels(movement$location.habitat)
movement$random.habitat <- as.character(movement$random.habitat)
movement$random.habitat[movement$random.habitat=="Pool (Reservoir)"] <- "Pool"
movement$random.habitat <- as.factor(movement$random.habitat)
levels(movement$random.habitat)
movement$stream.habitat <- as.character(movement$stream.habitat)
movement$stream.habitat[movement$stream.habitat=="Pool (Reservoir)"] <- "Pool"
movement$stream.habitat[movement$stream.habitat=="Run (Dried)"] <- "Run"
movement$stream.habitat <- as.factor(movement$stream.habitat)
levels(movement$stream.habitat)

## Water Depth
movement$location.water.depth <- as.character(movement$location.water.depth)
movement$location.water.depth[movement$location.water.depth == ">100"] <- 100
movement$location.water.depth <- as.integer(movement$location.water.depth)
movement$random.water.depth <- as.character(movement$random.water.depth)
movement$random.water.depth[movement$random.water.depth == ">100"] <- 100
movement$random.water.depth <- as.integer(movement$random.water.depth)
movement$stream.water.depth <- as.character(movement$stream.water.depth)
movement$stream.water.depth[movement$stream.water.depth == ">100"] <- 100
movement$stream.water.depth[movement$stream.water.depth == "<1"] <- 0
movement$stream.water.depth <- as.integer(movement$stream.water.depth) 

## Soil Types - Removed due to Great Homogeneity
movement <- movement[,!(names(movement) %in% c("location.soil", "random.soil", "stream.soil"))]

## Canopy Cover - Prepared for merging datasets
names(movement)[names(movement) == "location.canopy"] <- "location.canopyImg"
names(movement)[names(movement) == "random.canopy"] <- "random.canopyImg"
names(movement)[names(movement) == "stream.canopy"] <- "stream.canopyImg"

## Remarks
names(movement)[names(movement) == "Remark3"] <- "Remarks3"
movement$Remarks1 <- as.character(movement$Remarks1)
movement$Remarks2 <- as.character(movement$Remarks2)
movement$Remarks3 <- as.character(movement$Remarks3)

## Load Canopy Cover
relocationCanopy1 <- read.csv("data/SABE_Relocation_canopy01.csv", header = TRUE)
relocationCanopy2 <- read.csv("data/SABE_Relocation_canopy02.csv", header = TRUE)
relocationCanopy3 <- read.csv("data/SABE_Relocation_canopy03.csv", header = TRUE)
relocationCanopy <- rbind(relocationCanopy1, relocationCanopy2, relocationCanopy3)
View(relocationCanopy)
str(relocationCanopy)

randBearCanopy <- read.csv("data/SABE_RandomBearing_canopy.csv", header = TRUE)
View(randBearCanopy)
str(randBearCanopy)

randDirectCanopy <- read.csv("data/SABE_RandomDirection_canopy.csv", header = TRUE)
View(randDirectCanopy)
str(randDirectCanopy)

## Index
relocationCanopy <- relocationCanopy[,!names(relocationCanopy) %in% c("Index")]
randBearCanopy <- randBearCanopy[,!names(randBearCanopy) %in% c("Index")]
randDirectCanopy <- randDirectCanopy[,!names(randDirectCanopy) %in% c("Index")]

## Canopy Percentage
relocationCanopy$location.canopy <- 100 - relocationCanopy$gap_fraction
randBearCanopy$random.canopy <- 100 - randBearCanopy$gap_fraction
randDirectCanopy$stream.canopy <- 100 - randDirectCanopy$gap_fraction

## Photo File Name
## Remove file extension ".JPG"
relocationCanopy$photo <- as.character(relocationCanopy$photo)
relocationCanopy$photo <- gsub("\\..*", "", relocationCanopy$photo)
relocationCanopy$photo <- as.factor(relocationCanopy$photo)
randBearCanopy$photo <- as.character(randBearCanopy$photo)
randBearCanopy$photo <- gsub("\\..*", "", randBearCanopy$photo)
randBearCanopy$photo <- as.factor(randBearCanopy$photo)
randDirectCanopy$photo <- as.character(randDirectCanopy$photo)
randDirectCanopy$photo <- gsub("\\..*", "", randDirectCanopy$photo)
randDirectCanopy$photo <- as.factor(randDirectCanopy$photo)

## Combined Dataset - `movement` (relocation) & `canopy`
movement <- merge(movement, relocationCanopy[, c("photo", "location.canopy")], 
                  by.x = "location.canopyImg", by.y = "photo", all = TRUE)
movement <- merge(movement, randBearCanopy[, c("photo", "random.canopy")], 
                  by.x = "random.canopyImg", by.y = "photo", all = TRUE)
movement <- merge(movement, randDirectCanopy[, c("photo", "stream.canopy")], 
                  by.x = "stream.canopyImg", by.y = "photo", all = TRUE)

movement <- movement[order(movement$event.id),]
movement <- movement %>% drop_na(event.id)
movement <- movement[,!names(movement) %in% c("location.canopyImg", "random.canopyImg", "stream.canopyImg")]
View(movement)
str(movement)

## Load Turtle Biometric Data ----
turtlesBiometric <- read.csv("data/FWT_CMR_records.csv", header=TRUE, na.strings = c("", "NA"))
turtlesBiometric <- turtlesBiometric %>% drop_na(Species)
View(turtlesBiometric)
str(turtlesBiometric)

## Keep SABE Records Only
turtlesBiometric <- turtlesBiometric[turtlesBiometric$Species=="SABE",]
turtlesBiometric$Species <- as.factor(as.character(turtlesBiometric$Species))

## Turtle ID
turtlesBiometric$Turtle.ID <- as.character(turtlesBiometric$Turtle.ID)
turtlesBiometric$Turtle.ID[turtlesBiometric$Turtle.ID == "12/16"] <- "16"
turtlesBiometric$Turtle.ID <- as.factor(turtlesBiometric$Turtle.ID)
levels(turtlesBiometric$Turtle.ID)

## Radio-Frequency
turtlesBiometric$Transmitter.Freq. <- as.factor(turtlesBiometric$Transmitter.Freq.)
levels(turtlesBiometric$Transmitter.Freq.)

## Pit-Tag
turtlesBiometric$Pit.tag <- as.factor(as.character(turtlesBiometric$Pit.tag))

## Date
turtlesBiometric$Date <- parse_date(turtlesBiometric$Date, "%d-%b-%y")

## Location
turtlesBiometric$Location <- as.factor(as.character(turtlesBiometric$Location))

## Remarks
names(turtlesBiometric)[names(turtlesBiometric) == "Remarks"] <- "Remarks1"
names(turtlesBiometric)[names(turtlesBiometric) == "X"] <- "Remarks2"
names(turtlesBiometric)[names(turtlesBiometric) == "X.1"] <- "Remarks3"
turtlesBiometric$Remarks1 <- as.character(turtlesBiometric$Remarks1)
turtlesBiometric$Remarks2 <- as.character(turtlesBiometric$Remarks2)
turtlesBiometric$Remarks3 <- as.character(turtlesBiometric$Remarks3)

## Typo
turtlesBiometric$New.Recapture <- as.character(turtlesBiometric$New.Recapture)
turtlesBiometric$New.Recapture[turtlesBiometric$New.Recapture == "Reecapture"] <- "Recapture"
turtlesBiometric$New.Recapture <- as.factor(turtlesBiometric$New.Recapture)
levels(turtlesBiometric$New.Recapture)

## Removed Metrics Related to Other Studies
## Inc. TransDist
## Inc. HW, TL
## Inc. Accelerometer.Spec.
turtlesBiometric <- turtlesBiometric[,!names(turtlesBiometric) %in% c("TransDist", "HW", "TL", 
                                                                      "Accelerometer.Spec.")]

## Combine Dataset - `movement` (relocation) & `turtlesBiometric`
turtles <- merge(movement, turtlesBiometric, 
                 by.x = c("turtle.id", "frequency"), by.y = c("Turtle.ID", "Transmitter.Freq."))
turtles <- turtles[order(turtles$event.id),]
str(turtles, list.len=ncol(turtles))

## Check If SpatialPointsDataFrame Objects Include Missing Values
turtles[is.na(turtles$location.lat)]
turtles[is.na(turtles$location.lon)]

## Spatial Analysis -- Creating SpatialPointsDataFrame ----
## Creating Spatial Points
SABE.Sp <- turtles[, c("turtle.id", "location.lon", "location.lat")]
View(SABE.Sp)
str(SABE.Sp)

coordinates(SABE.Sp) <- c("location.lon", "location.lat")
str(SABE.Sp)

## Projection / Setting CRS
proj4string(SABE.Sp) <- CRS("+init=epsg:4326")

## Re-Projection
SABE.SpProj <- spTransform(SABE.Sp, 
                           CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs"))
str(SABE.SpProj)
plot(SABE.SpProj@coords)

## Creating Spatial Points by Seasons
SABE.season = list()
SABE.seasonSp = list()
SABE.seasonSpProj = list()
for(i in 1:length(levels(turtles$season))){
  SABE.season[[i]] <- turtles %>% filter(season == levels(turtles$season)[i])
  #### remove individuals with fewer than 5 records for MCP estimation
  tab <- table(SABE.season[[i]]$turtle.id)
  SABE.season[[i]]$turtle.id <- as.character(SABE.season[[i]]$turtle.id)
  SABE.season[[i]] <- SABE.season[[i]][SABE.season[[i]]$turtle.id %in% names(tab)[tab>4],]
  SABE.season[[i]]$turtle.id <- as.factor(SABE.season[[i]]$turtle.id)
  
  if(nrow(SABE.season[[i]]) != 0){
    ## Creat Spatial Points Data Frame
    SABE.seasonSp[[i]] <- SABE.season[[i]][, c("turtle.id", "location.lon", "location.lat")]
    coordinates(SABE.seasonSp[[i]]) <- c("location.lon", "location.lat")
    ## Projection
    proj4string(SABE.seasonSp[[i]]) <- CRS("+init=epsg:4326")
    ## Re-Projection
    SABE.seasonSpProj[[i]] <- spTransform(SABE.seasonSp[[i]], 
                                          CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs"))
  }
  i = i + 1
}

## Spatial Analysis -- Creating MCP ---- 
## Calculate MCPs for All Records
## Remove records with less than 5 relocations
tab <- table(SABE.SpProj$turtle.id)
SABE.SpProj$turtle.id <- as.character(SABE.SpProj$turtle.id)
SABE.SpProj <- SABE.SpProj[SABE.SpProj$turtle.id %in% names(tab)[tab>4],]
SABE.SpProj$turtle.id <- as.factor(SABE.SpProj$turtle.id)

SABE.MCP100 <- mcp(SABE.SpProj, percent = 100, unin = "m", unout="m2") 

SABE.PopuMCP100 <- mcp(SABE.SpProj[,!(names(SABE.SpProj) %in% c("turtle.id"))], percent = 100, unin = "m", unout="m2")

## Calculate MCPs by Seasons
SABE.seasonMCP100 = list()
SABE.seasonPopuMCP100 = list()
for(i in 1:length(SABE.seasonSpProj)){
  SABE.seasonMCP100[[i]] <- mcp(SABE.seasonSpProj[[i]], percent = 100, unin = "m", unout = "m2")
  SABE.seasonMCP100[[i]]$season <- levels(turtles$season)[i]
  SABE.seasonPopuMCP100[[i]] <- mcp(SABE.seasonSpProj[[i]][,!(names(SABE.seasonSpProj) %in% c("turtle.id"))], 
                                    percent = 100, unin = "m", unout = "m2")
  SABE.seasonPopuMCP100[[i]]$season <- levels(turtles$season)[i]
  i = i + 1
}

## Combine Summarized Datasets - Biometric & MCP Results
#### Number of Relocation Points
base <- ddply(turtles, .(turtle.id), nrow)
names(base)[names(base) == "V1"] <- "Pts"
#### Combine Biometric Datasets & Num of Relocation Points
base <- merge(base, turtlesBiometric, by.x = c("turtle.id"), by.y = c("Turtle.ID"))
base <- base %>% drop_na(Transmitter.Freq.)
#### Retain the First Measurement at Start Only
base <- base %>% 
  group_by(turtle.id) %>% 
  filter(Date == min(Date))
base <- as.data.frame(base)
base <- base[,!(names(base)) %in% c("Species", "New.Recapture", "Date", "Location", 
                                    "Lat.D", "Lat.M", "Lon.D", "Lon.M", "Lat.Dd", "Lon.Dd", 
                                    "Leech", "Pit.tag", "Transmitter.Freq.",
                                    "Remarks1", "Remarks2", "Remarks3")]
base <- with(base, data.frame(turtle.id, Pts, 'sex' = Sex, 'cl' = CA, 'pl' = PA, 'wt' = Wt))
str(base)

# write.csv(base, "data/SABE_biometry.csv", row.names = FALSE)

#### Combine Base & MCP100
SABE.MCP100results <- data.frame(ID = SABE.MCP100$id, 
                                 HR_Type = "MCP",
                                 MCP_Lv = 100,
                                 Area_m2 = SABE.MCP100$area)
SABE.MCP100results <- merge(SABE.MCP100results, base, by.x = "ID", by.y = "turtle.id")

SABE.MCP100results <- with(SABE.MCP100results, data.frame('turtle.id' = ID, sex, cl, pl, wt, 
                                                          HR_Type, MCP_Lv, Area_m2, Pts))
write.csv(SABE.MCP100results, "data/SABE_movement_mcp.csv", row.names = FALSE)

## Combine Summarised Datasets Based on Season - Biometric & MCP Results (Season)
#### Number of Relocation Points
baseSeason <- turtles %>% group_by(turtle.id, season) %>% dplyr::summarise(Pts = n()) %>% filter(Pts > 4)
baseSeason <- as.data.frame(baseSeason)
#### Update Base to Num of Relocation Points by Season
baseSeason <- merge(baseSeason, base[,!names(base) %in% c("Pts")], by = "turtle.id")
str(baseSeason)

#### Combine Base & MCP100 by Season
SABE.seasonMCP100results <- do.call("rbind", SABE.seasonMCP100)
SABE.seasonMCP100results <- SABE.seasonMCP100results@data
SABE.seasonMCP100results <- data.frame(ID = SABE.seasonMCP100results$id, 
                                       Season = SABE.seasonMCP100results$season,
                                       HR_Type = "MCP",
                                       MCP_Lv = 100,
                                       Area_m2 = SABE.seasonMCP100results$area)
SABE.seasonMCP100results <- merge(SABE.seasonMCP100results, baseSeason, 
                                  by.x = c("ID", "Season"), by.y = c("turtle.id", "season"))

SABE.seasonMCP100results <- with(SABE.seasonMCP100results, data.frame('turtle.id' = ID, 'season' = Season, 
                                                                      sex, cl, pl, wt, 
                                                                      HR_Type, MCP_Lv, Area_m2, Pts))
write.csv(SABE.seasonMCP100results, "data/SABE_movement_mcp_season.csv", row.names = FALSE)

## Spatial Analysis -- Dynamic Brownian Bridge Movement Models (dBBMMs) ----
## Set location error for dBBMM analyses
set_loc.error <- 10
set_grid.ext <- 50
set_dimsize <- 1000

## Set margin and window 
set_mrg <- 3
set_ws <- 7

## Inputting dataset
turtles <- turtles %>% arrange(turtle.id)
str(turtles)

dbbmm.list <- list()
for(i in 1:length(levels(turtles$turtle.id))){
  
  data <- turtles %>% filter(turtle.id == levels(turtles$turtle.id)[i])
  
  ## Calculating trajectory
  move <- move(x = data$location.lon, y = data$location.lat, 
               time = data$time, 
               proj = CRS("+init=epsg:4326"), 
               data = data)
  move <- spTransform(move, CRSobj = "+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs")
  
  ## Calculate the dynamic brownian motion variance
  dbbv <- brownian.motion.variance.dyn(object = move, 
                                       location.error = set_loc.error, 
                                       window.size = set_ws,
                                       margin = set_mrg)
  dbbv@interest[timeLag(move,"mins")>10*24*60] <- FALSE
  
  ## Calculating dBBMM
  dbbmm <- brownian.bridge.dyn(dbbv, 
                               dimSize = set_dimsize,
                               ext = set_grid.ext,
                               location.error = set_loc.error,
                               time.step = 768)
  
  ## Calculating and extracting dBBMM variance
  tempname <- paste0("outputs/SABE_ID", levels(turtles$turtle.id)[i], "_dBBMM.data.csv")
  data$var <- getMotionVariance(dbbmm)
  write.csv(data, tempname)
  
  ## Extracting Raster
  dbbmm.sp <- as(dbbmm, "SpatialPixelsDataFrame")
  dbbmm.sp.ud <- new("estUD", dbbmm.sp)
  dbbmm.sp.ud@vol = FALSE
  dbbmm.sp.ud@h$meth = "dBBMM"
  dbbmm.ud <- getvolumeUD(dbbmm.sp.ud, standardize = TRUE)
  
  r <- as(dbbmm.ud, "SpatialPixelsDataFrame") %>% raster()
  
  ## Extracting contours
  contour.099 <- raster2contour(dbbmm, levels = .99, maxpixels = 1000000)
  
  tempname.099 <- paste0("outputs/SABE_ID", levels(turtles$turtle.id)[i], "_dBBMM.099")
  write_sf(as(contour.099, "sf"), ".", tempname.099, driver = "ESRI shapefile")
  
  dbbmm.list[[i]] <- dbbmm
  
}

## Spatial Analysis -- Summarising and Plotting dBBMMs ----
## UD for Each Individuals on Plain Bkg
SABE_ud_shp_file <- list.files(path = ".", 
                               pattern = "outputs_SABE_ID.[0-9+]*._dBBMM.099.*shp",
                               full.names = TRUE, 
                               recursive = TRUE)
SABE_ud_shp <- lapply(SABE_ud_shp_file, function(x){
  dat <- st_read(x)
  dat <- st_as_sf(dat)
  dat <- st_polygonize(dat)
  dat <- st_transform(dat, crs = "+proj=longlat")
  dat$id <- as.character(x)
  dat$id <- gsub("^\\D+(\\d+).*", "\\1", dat$id)
  return(dat)
})
SABE_ud_shp <- bind_rows(SABE_ud_shp)

SABE.UDresults <- SABE_ud_shp
SABE.UDresults$area <- SABE.UDresults %>% 
  st_transform("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs") %>% 
  st_area()

SABE.UDresults <- data.frame(ID = SABE.UDresults$id, 
                             HR_Type = "dBBMM",
                             Ct_Vol = as.numeric(SABE.UDresults$level)*100,
                             Area_m2 = as.numeric(SABE.UDresults$area))
SABE.UDresults <- merge(SABE.UDresults, base, by.x = "ID", by.y = "turtle.id")

## Spatial Analysis -- Stream Distance ---- 
SABE.Sf <- st_as_sf(SABE.SpProj)
SABE.Sf <- merge(SABE.Sf, base, by.x = "turtle.id", by.y = "turtle.id")
names(SABE.Sf)[names(SABE.Sf)=="turtle.id"] <- "id"

SABE.Sf_group <- SABE.Sf %>% 
  group_by(id) %>% 
  group_split()

## MASK THE FILE NAME ----
hydroline_utm <- st_read("data/HydrographyLine.shp") %>%
  st_union() %>%  
  st_transform(crs = "+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs")

hydroline_sfn <- hydroline_utm %>%
  st_cast("LINESTRING") %>% 
  as_sfnetwork(directed = FALSE)  

distance_per_group <- map(SABE.Sf_group, function(x){
  
  id <- unique(x$id)
  
  distance_matrix <- st_distance(x)
  max_distance <- as.numeric(which(distance_matrix == max(distance_matrix), arr.ind = TRUE)[1,])
  x <- x[max_distance,]
  
  start_on_path <- st_nearest_points(x[1,], hydroline_utm) %>% 
    st_cast("POINT") %>% 
    st_difference(., st_geometry(x[1,]))
  end_on_path <- st_nearest_points(x[2,], hydroline_utm) %>%
    st_cast("POINT") %>%
    st_difference(., st_geometry(x[2,]))
  
  points_on_path <- st_nearest_points(x, hydroline_utm)
  
  dist1 <- st_distance(x[1,], start_on_path)
  dist2 <- st_distance(x[2,], end_on_path)
  
  joint_sfn <- st_network_blend(hydroline_sfn, st_geometry(x))
  
  stream_distances <- st_network_cost(x = joint_sfn,
                                      from = start_on_path,
                                      to = end_on_path,
                                      direction = "all") 
  
  dist <- as.numeric(dist1) + as.numeric(dist2) + as.numeric(stream_distances)
  
  data.frame(ID = id, 
             HR_Type = "Stream distance",
             Ct_Vol = 100,
             Area_m2 = dist)
  
})
distance_per_group <- bind_rows(distance_per_group)
SABE.streamDistResulsts <- merge(distance_per_group, base, by.x = "ID", by.y = "turtle.id")

## Spatial Analysis -- Displacement Distance ---- 
## Creating Spatial Point Data Frame
SABE.spData <- turtles[, c("event.id", "date", "time", "session", "month", "year", "season", 
                           "turtle.id", "frequency", "location.lon", "location.lat")]
View(SABE.spData)
str(SABE.spData)

SABE.spData$radio.id <- paste(SABE.spData$turtle.id, SABE.spData$frequency, sep = "_")
SABE.spData$radio.id <- as.factor(SABE.spData$radio.id)

coordinates(SABE.spData) <- c("location.lon", "location.lat")
str(SABE.spData)

## Projection
proj4string(SABE.spData) <- CRS(SRS_string = "EPSG:4326")

## Re-projection
SABE.spData <- spTransform(SABE.spData, 
                           CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs"))
proj4string(SABE.spData)

## Create a Trajectories for Individual Turtle
SABE.traj <- as.ltraj(coordinates(SABE.spData), date = SABE.spData$time, 
                      id = SABE.spData$radio.id, typeII = TRUE, 
                      proj4string = CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs"))
class(SABE.traj)
adehabitatLT::is.regular(SABE.traj)

SABE.traj
head(SABE.traj[[1]])

## Combine Lists of Individual Trajectories
SABE.trajDf <- lapply(1:length(SABE.traj), function(i){
  data.frame(id = attributes(SABE.traj[[i]])$id, SABE.traj[[i]])
})
SABE.trajDf <- do.call(rbind, SABE.trajDf)

SABE.trajDf$id <- as.character(SABE.trajDf$id)
SABE.trajDf$turtle.id <- sub("_.*", "", SABE.trajDf$id)
SABE.trajDf$turtle.id <- as.factor(SABE.trajDf$turtle.id)
SABE.trajDf$frequency <- sub(".*_", "", SABE.trajDf$id)
SABE.trajDf$frequency <- as.factor(SABE.trajDf$frequency)
SABE.trajDf <- dplyr::rename(SABE.trajDf, radio.id = id)
SABE.trajDf$radio.id <- as.factor(SABE.trajDf$radio.id)

SABE.trajDf$year <- year(SABE.trajDf$date)
SABE.trajDf$year <- as.factor(SABE.trajDf$year)
SABE.trajDf$month <- months(SABE.trajDf$date)
SABE.trajDf$month <- factor(SABE.trajDf$month, c("January", "February", "March", "April", "May", "June", 
                                                 "July", "August", "September", "October", "November", "December"))
SABE.trajDf$season <- ifelse(SABE.trajDf$month %in% c("May", "June", "July", "August"), "Wet Season", 
                             ifelse(SABE.trajDf$month %in% c("September", "October"), "Mating Season", 
                                    ifelse(SABE.trajDf$month %in% c("November", "December"), "Dry Season", 
                                           ifelse(SABE.trajDf$month %in% c("January", "February"), "Dry Season (L)", 
                                                  ifelse(SABE.trajDf$month %in% c("March", "April"), "Nesting Season",
                                                         "Unclassified")))))
SABE.trajDf$season <- paste(as.character(SABE.trajDf$year), as.character(SABE.trajDf$season))
SABE.trajDf$season[SABE.trajDf$season == "2021 Dry Season (L)"] <- "2020 Dry Season"
SABE.trajDf$season <- factor(SABE.trajDf$season, c("2020 Wet Season", "2020 Mating Season", "2020 Dry Season", "2021 Nesting Season", 
                                                   "2021 Wet Season", "2021 Mating Season", "2021 Dry Season"))
SABE.trajDf$year.2 <- sub("\\s.*", "", SABE.trajDf$season)
SABE.trajDf$year.2 <- as.factor(SABE.trajDf$year.2)
SABE.trajDf$season.4 <- sub("\\d+\\s", "", SABE.trajDf$season)
SABE.trajDf$season.4 <- factor(SABE.trajDf$season.4, c("Nesting Season", "Wet Season", "Mating Season", "Dry Season"))

SABE.trajDf <- merge(SABE.trajDf, base, by = "turtle.id")

View(SABE.trajDf)
str(SABE.trajDf)

SABE.trajDf <- SABE.trajDf %>% drop_na("dist") # Remove records with NA
SABE.trajDf <- SABE.trajDf[SABE.trajDf$dt <= 14*24*60*60, ] # Remove records that are more than 14 days apart

SABE.trajDf <- SABE.trajDf[order(SABE.trajDf$date),]

## Calculate Daily Displacement Distance for each Records
SABE.trajDf$dailyDist <- SABE.trajDf$dist/(SABE.trajDf$dt/(60*60*24)) 

## Save Dataset
SABE.trajDf <- with(SABE.trajDf, data.frame(turtle.id, radio.id, 
                                            sex, cl, pl, wt, 
                                            year, month, season, dx, dy, dist, dt, dailyDist))
write.csv(SABE.trajDf, "data/SABE_movement_dist.csv", row.names = FALSE)



## Habitat Selection ----
### Compositional Analysis -- Data Preparation ----
#### Study Stream: Habitat Data ----
#### Load Random Points Data
habitat <- read.csv("data/SABE_RandomPts.csv", header = TRUE, na.strings = c("NA", ""))
View(habitat)
str(habitat)

#### Date, Month and Year
habitat$Date <- parse_date(habitat$Date, "%d-%b-%y")
habitat$Month[habitat$Month == "Apirl"] <- "April"
habitat$Month <- factor(habitat$Month, c("January", "February", "March", "April", "May", "June", 
                                         "July", "August", "September", "October", "November", "December"))
habitat$Year <- year(habitat$Date)
habitat$Year <- as.factor(habitat$Year)
#### Season
habitat$season <- ifelse(habitat$Date >= "2020-05-01" & habitat$Date <= "2020-08-31", "2020 Wet Season", 
                         ifelse(habitat$Date >= "2020-09-01" & habitat$Date <= "2020-10-31", "2020 Mating Season", 
                                ifelse(habitat$Date >= "2020-11-01" & habitat$Date < "2021-03-01", "2020 Dry Season", 
                                       ifelse(habitat$Date >= "2021-03-01" & habitat$Date <= "2021-04-30", "2021 Nesting Season", 
                                              ifelse(habitat$Date >= "2021-05-01" & habitat$Date <= "2021-08-31", "2021 Wet Season", 
                                                     ifelse(habitat$Date >= "2021-09-01" & habitat$Date <= "2021-10-31", "2021 Mating Season", 
                                                            ifelse(habitat$Date >= "2021-11-01" & habitat$Date <= "2022-03-01", "2021 Dry Season", 
                                                                   "Unclassified")))))))
habitat$season <- as.factor(habitat$season)
#### Weather Info
habitat$Weather <- sub("(Cloudy).*", "\\1", habitat$Weather)
habitat$Weather[habitat$Weather == "Partly cloudy"] <- "Partly Cloudy"
habitat$Weather <- as.factor(habitat$Weather)
#### Habitat Types
habitat$Habitat[habitat$Habitat == "Pool (Temp.)"] <- "Pool"
habitat$Habitat <- as.factor(habitat$Habitat)
#### Water Depth
habitat$Water.Depth[habitat$Water.Depth == ">100"] <- 100
habitat$Water.Depth <- as.integer(habitat$Water.Depth)
#### Vegetation Types
habitat$Submerged <- factor(habitat$Submerged, levels = c("No", "Yes"))
habitat$Floating <- factor(habitat$Floating, levels = c("No", "Yes"))
habitat$Emergent <- sub("(Yes).*", "\\1", habitat$Emergent)
habitat$Emergent <- factor(habitat$Emergent, levels = c("No", "Yes"))
habitat$Tree.Root <- factor(habitat$Tree.Root, levels = c("No", "Yes"))
habitat$Filamentous.Algae <- factor(habitat$Filamentous.Algae, levels = c("No", "Yes"))
habitat$Iron.Bacteria <- factor(habitat$Iron.Bacteria, levels = c("No", "Yes"))
#### Soil (Discarded)
habitat <- habitat[!names(habitat) %in% c("Soil")]

## Load Canopy Cover
randomCanopy1 <- read.csv("data/SABE_Random_canopy01.csv", header = TRUE)
randomCanopy2 <- read.csv("data/SABE_Random_canopy02.csv", header = TRUE)
randomCanopy3 <- read.csv("data/SABE_Random_canopy03.csv", header = TRUE)
randomCanopy4 <- read.csv("data/SABE_Random_canopy04.csv", header = TRUE)
randomCanopy <- rbind(randomCanopy1, randomCanopy2, randomCanopy3, randomCanopy4)
View(randomCanopy)
str(randomCanopy)

#### Index
randomCanopy <- randomCanopy[,!names(randomCanopy) %in% c("Index")]
#### Photo File Name
#### Remove file extension ".JPG"
randomCanopy$photo <- as.character(randomCanopy$photo)
randomCanopy$photo <- gsub("\\..*", "", randomCanopy$photo)
#### Canopy Cover
randomCanopy$canopy <- 100 - randomCanopy$gap_fraction

#### Combined Dataset - `habitat` (random) & `canopy`
habitat <- merge(habitat, randomCanopy[, c("photo", "canopy")], 
                 by.x = "Canopy", by.y = "photo", all = TRUE)
habitat <- habitat[order(habitat$Date),]
habitat <- habitat %>% drop_na(canopy)
habitat <- habitat[,!names(habitat) %in% c("Canopy")]
View(habitat)
str(habitat)

#### Remove Records with Habitat Type as Pavement
habitat <- habitat[habitat$Habitat != "Pavement",]

#### Write Dataset for Compositional Analysis 
# write.csv(habitat, "data/SABE_habitat_studyStream.csv", row.names = FALSE)

#### Study Stream: Habitat Proportion ----
studyarea <- habitat %>% 
  drop_na("Habitat") %>%
  group_by(Habitat) %>%
  dplyr::summarise(n=length(Habitat)) %>%
  dplyr::mutate(prop = n/sum(n)*100)
studyarea <- as.data.frame(studyarea)
studyarea <- studyarea %>% dplyr::select(-n)

studyarea_coln <- studyarea$Habitat
studyarea <- as.data.frame(t(studyarea[,-1]))
colnames(studyarea) <- studyarea_coln
studyarea <- as.data.frame(lapply(studyarea, rep, length(SABE.MCP100$id)))
rownames(studyarea) <- SABE.MCP100$id
studyarea <- cbind(id = rownames(studyarea), studyarea)

#### Write Dataset for Compositional Analysis 
studyarea <- with(studyarea, data.frame('turtle.id' = id, Pool, Riffle, Run))
write.csv(studyarea, "data/SABE_habitatProp_studyarea.csv", row.names = FALSE)

#### Study Stream: Habitat Proportion Per Seasons ----
studyareaSeasonal = list()
for(i in 1:length(levels(habitat$season))){
  studyareaSeasonal[[i]] <- habitat %>%
    filter(season == levels(turtles$season)[i]) %>%
    drop_na("Habitat") %>%
    group_by(Habitat) %>%
    dplyr::summarise(n=length(Habitat)) %>%
    dplyr::mutate(prop = n/sum(n)*100)
  studyareaSeasonal[[i]] <- as.data.frame(studyareaSeasonal[[i]])
  studyareaSeasonal[[i]] <- studyareaSeasonal[[i]] %>% dplyr::select(-n)
  
  studyareaSeasonal_coln <- studyareaSeasonal[[i]]$Habitat
  studyareaSeasonal[[i]] <- as.data.frame(t(studyareaSeasonal[[i]][,-1]))
  colnames(studyareaSeasonal[[i]]) <- studyareaSeasonal_coln
  studyareaSeasonal[[i]] <- as.data.frame(lapply(studyareaSeasonal[[i]], rep, length(SABE.seasonMCP100[[i]]$id)))
  rownames(studyareaSeasonal[[i]]) <- SABE.seasonMCP100[[i]]$id
  studyareaSeasonal[[i]]$season <- unique(SABE.seasonMCP100[[i]]$season)
  studyareaSeasonal[[i]] <- cbind(id = rownames(studyareaSeasonal[[i]]), studyareaSeasonal[[i]])
  i = i + 1
}

#### Write Dataset for Compositional Analysis 
studyareaSeasonalDf <- do.call(rbind, studyareaSeasonal)
studyareaSeasonalDf <- with(studyareaSeasonalDf, data.frame('turtle.id' = id, season, Pool, Riffle, Run))
write.csv(studyareaSeasonalDf, "data/SABE_habitatProp_studyarea_season.csv", row.names = FALSE)

#### MCP: Habitat Data ----
#### Calculating Habitat Proportion within 100% MCP with 15m Buffer
#### Buffer was set as the the 75% quantile of the distance moved between consecutive relocations 

#### Create a SpatialPointDataFrame for the Habitat Dataset
habitat.Sp <- habitat %>% drop_na("Lon", "Lat")
coordinates(habitat.Sp) <- c("Lon", "Lat")
proj4string(habitat.Sp) <- CRS("+init=epsg:4326")
habitat.Sp <- spTransform(habitat.Sp, CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs"))

#### Draw a Buffer Zone for each Individual's MCP
mcp.Ply <- spTransform(SABE.MCP100, CRS("+proj=longlat"))
mcp.Ply <- spTransform(mcp.Ply, CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs"))
mcp.PlyBuff <- gBuffer(as(mcp.Ply, "SpatialPolygons"), byid = TRUE, width = ceiling(quantile(SABE.trajDf$dist, 0.75))) # 15m Buffer

#### Extract All Points in the Habitat Dataset that Intersects with the MCP Polygon for each Individual
PtsPly <- rgeos::gIntersects(habitat.Sp, as(mcp.PlyBuff, "SpatialPolygons"), byid = TRUE)

#### Collate Habitat Data for All Points Intersects with the MCP 
PtsMcp <- apply(PtsPly, 1, function(x){
  c <- which(x == TRUE)
  habitat.Sp[as.numeric(c), ] 
})
PtsMcp <- do.call(rbind, lapply(names(PtsMcp), function(x) data.frame(c(ID=x, PtsMcp[[x]]))))

#### Write Dataset for Compositional Analysis 
# write.csv(PtsMcp, "data/SABE_habitat_mcp.csv", row.names = FALSE)

#### MCP: Habitat Proportion ----
mcp <- PtsMcp %>%
  group_by(ID, Habitat) %>% 
  dplyr::summarise(n=length(Habitat)) %>%
  dplyr::mutate(prop = n/sum(n)*100)
mcp <- as.data.frame(mcp)
mcp <- mcp %>% dplyr::select(-n)
mcp <- matrify(mcp)
mcp <- cbind(id = rownames(mcp), mcp)

#### Write Dataset for Compositional Analysis 
mcp <- with(mcp, data.frame('turtle.id' = id, Pool, Riffle, Run))
write.csv(mcp, "data/SABE_habitatProp_mcp.csv", row.names = FALSE)

#### MCP: Habitat Data Per Season ----
habitat.SeasonalSp = list()
mcp.SeasonalPly = list()
mcp.SeasonalPlyBuff = list()
for(i in 1:length(levels(habitat$season))){
  habitat.SeasonalSp[[i]] <- habitat %>% filter(season == levels(turtles$season)[i])
  habitat.SeasonalSp[[i]] <- habitat.SeasonalSp[[i]] %>% drop_na("Lon", "Lat")
  coordinates(habitat.SeasonalSp[[i]]) <- c("Lon", "Lat")
  proj4string(habitat.SeasonalSp[[i]]) <- CRS("+init=epsg:4326")
  habitat.SeasonalSp[[i]] <- spTransform(habitat.SeasonalSp[[i]], CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs"))
  
  mcp.SeasonalPly[[i]] <- spTransform(SABE.seasonMCP100[[i]], CRS("+proj=longlat"))
  mcp.SeasonalPly[[i]] <- spTransform(mcp.SeasonalPly[[i]], CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs"))
  mcp.SeasonalPlyBuff[[i]] <- gBuffer(as(mcp.SeasonalPly[[i]], "SpatialPolygons"), byid = TRUE, width = ceiling(quantile(SABE.trajDf$dist, 0.75)))
  i = i + 1
}

PtsPlySeasonal = list()
PtsMcpSeasonal = list()
for(i in 1:length(mcp.SeasonalPlyBuff)){
  PtsPlySeasonal[[i]] <- rgeos::gIntersects(habitat.SeasonalSp[[i]], as(mcp.SeasonalPlyBuff[[i]], "SpatialPolygons"), byid = TRUE)
  
  PtsMcpSeasonal[[i]] <- apply(PtsPlySeasonal[[i]], 1, function(x){
    c <- which(x == TRUE)
    habitat.SeasonalSp[[i]][as.numeric(c), ]
  })
  PtsMcpSeasonal[[i]] <- do.call(rbind, lapply(names(PtsMcpSeasonal[[i]]), function(x) data.frame(c(ID=x, PtsMcpSeasonal[[i]][[x]]))))
  i = i + 1
}

#### Write Dataset for Compositional Analysis 
PtsMcpSeasonalDf <- do.call(rbind, PtsMcpSeasonal)
write.csv(PtsMcpSeasonalDf, "data/SABE_habitat_mcp_season.csv", row.names = FALSE)

#### MCP: Habitat Proportion Per Season ----
mcpSeasonal = list()
for(i in 1:length(PtsMcpSeasonal)){
  mcpSeasonal[[i]] <- PtsMcpSeasonal[[i]] %>% 
    group_by(ID, Habitat) %>% 
    dplyr::summarise(n=length(Habitat)) %>%
    dplyr::mutate(prop = n/sum(n)*100)
  mcpSeasonal[[i]] <- as.data.frame(mcpSeasonal[[i]])
  mcpSeasonal[[i]] <- mcpSeasonal[[i]] %>% dplyr::select(-n)
  mcpSeasonal[[i]] <- matrify(mcpSeasonal[[i]])
  mcpSeasonal[[i]]$season <- unique(PtsMcpSeasonal[[i]][["season"]])
  mcpSeasonal[[i]] <- cbind(id = rownames(mcpSeasonal[[i]]), mcpSeasonal[[i]])
  i = i + 1
}

#### Write Dataset for Compositional Analysis 
mcpSeasonalDf <- do.call(rbind, mcpSeasonal)
mcpSeasonalDf <- with(mcpSeasonalDf, data.frame('turtle.id' = id, season, Pool, Riffle, Run))
write.csv(mcpSeasonalDf, "data/SABE_habitatProp_mcp_season.csv", row.names = FALSE)



### Resource Selection Function -- Data Preparation ----
turtlesHS <- turtles
turtlesHS <- turtlesHS[!turtlesHS$behaviour %in% c("Surface (Trapped)", "Underwater (Trapped)", "Swimming"), ]
turtlesHS$behaviour <- as.character(turtlesHS$behaviour)
turtlesHS$behaviour <- as.factor(turtlesHS$behaviour)
levels(turtlesHS$behaviour)

#### Assigning Random Points ----
terrestrialHS <- turtlesHS[turtlesHS$location.habitat == "Forested Upland",] 
nrow(terrestrialHS) ## 19 obs.

aquaticHS <- turtlesHS[turtlesHS$location.habitat != "Forested Upland",]
used <- aquaticHS %>% dplyr::select(turtle.id : session, month : season, 
                                    location.lat.D : location.substrate, location.canopy)
avail <- aquaticHS %>% dplyr::select(turtle.id : session, month : season, 
                                     stream.section : stream.substrate, stream.canopy)

#### Extract Relocation Data Already Paired
#### during the initial period (May 2020 and June 2020) and 
usedInit <- used %>% filter(date <= as.Date("2020-06-30")) 
availInit <- avail %>% filter(date <= as.Date("2020-06-30"))
#### during the ending period (Sept 2021 and Oct 2021)
usedEnd <- used %>% filter(date >= as.Date("2021-09-01")) 
availEnd <- avail %>% filter(date >= as.Date("2021-09-01"))
#### combine the two periods
usedPaired <- rbind(usedInit, usedEnd)
usedPaired <- usedPaired %>% 
  dplyr::rename(
    lat.D = location.lat.D,
    lat.M = location.lat.M,
    lon.D = location.lon.D,
    lon.M = location.lon.M,
    lat = location.lat,
    lon = location.lon,
    elev = location.elev,
    habitat = location.habitat,
    airTemp = location.air.temp,
    waterTemp = location.water.temp,
    waterDepth = location.water.depth,
    streamWidth = location.water.width,
    gravel = location.substrate.gravel,
    pebble = location.substrate.pebble,
    cobble = location.substrate.cobble,
    boulder = location.substrate.boulder,
    submerged = location.submerged,
    floating = location.floating,
    emergent = location.emergent,
    root = location.tree.root,
    algae = location.filamentous.algae,
    ironBacteria = location.iron.bacteria,
    leafLitter = location.leaf.litter,
    canopy = location.canopy,
    substrate = location.substrate
  )

availPaired <- rbind(availInit, availEnd)
availPaired <- availPaired %>% 
  dplyr::rename(
    direction = stream.section,
    dist = stream.dist,
    perc = stream.width,
    lat.D = stream.lat.D,
    lat.M = stream.lat.M,
    lon.D = stream.lon.D,
    lon.M = stream.lon.M,
    lat = stream.lat,
    lon = stream.lon,
    elev = stream.elev,
    habitat = stream.habitat,
    airTemp = stream.air.temp,
    waterTemp = stream.water.temp,
    waterDepth = stream.water.depth,
    streamWidth = stream.water.width,
    gravel = stream.substrate.gravel,
    pebble = stream.substrate.pebble,
    cobble = stream.substrate.cobble,
    boulder = stream.substrate.boulder,
    submerged = stream.submerged,
    floating = stream.floating,
    emergent = stream.emergent,
    root = stream.tree.root,
    algae = stream.filamentous.algae,
    ironBacteria = stream.iron.bacteria,
    leafLitter = stream.leaf.litter,
    canopy = stream.canopy,
    substrate = stream.substrate
  )

#### Extract Relocation Data that are Un-paired
usedUnpaired <- used %>% filter(date > as.Date("2020-06-30") & date < as.Date("2021-09-01")) 
availUnpaired <- avail %>% filter(date > as.Date("2020-06-30") & date < as.Date("2021-09-01"))

#### Match Pairs Automatically
set.seed(1234)

usedMonth = list()
availMonth = list()
habitatMonth = list()

for(i in 1:length(levels(aquaticHS$year))){
  usedMonth[[i]] <- list(i)
  availMonth[[i]] <- list(i)
  habitatMonth[[i]] <- list(i)
  set <- usedUnpaired %>% filter(year == levels(year)[i]) %>% dplyr::select(month)
  set$month <- as.factor(as.character(set$month))
  for(j in 1:length(levels(set$month))){
    usedMonth[[i]][[j]] <- list(j)
    availMonth[[i]][[j]] <- list(j)
    habitatMonth[[i]][[j]] <- list(j)
    
    usedMonth[[i]][[j]] <- usedUnpaired %>% filter(year == levels(usedUnpaired$year)[i] & month == levels(set$month)[j])
    habitatMonth[[i]][[j]] <- habitat %>% filter(Year == levels(usedUnpaired$year)[i] & Month == levels(set$month)[j])
    
    gridMat <- expand.grid(1:nrow(usedMonth[[i]][[j]]), 1:nrow(habitatMonth[[i]][[j]]))
    distMat <- matrix(distGeo(usedMonth[[i]][[j]][gridMat[,1], c("location.lon", "location.lat")], 
                              habitatMonth[[i]][[j]][gridMat[,2], c("Lon","Lat")]), 
                      nrow=nrow(usedMonth[[i]][[j]]))
    
    availMeta <- usedMonth[[i]][[j]] %>% dplyr::select(turtle.id:season)
    availMonth[[i]][[j]] <- habitatMonth[[i]][[j]][apply(distMat, 1, function(x) {
      c <- which(x >= round(as.numeric(quantile(SABE.trajDf$dist, 0.25))) 
                 & x <= round(as.numeric(quantile(SABE.trajDf$dist, 0.75)))+10) # GPS Fix Error of 10M
      if(length(c)==0) NA
      else c[sample(length(c), 1)]
      # Not Run: if(length(c)<=2) c <- c(c, rep("NA", 2-length(c)))
      # Not Run: c[sample(length(c), 2)]
    }),]
    availMonth[[i]][[j]] <- availMonth[[i]][[j]] %>% dplyr::select(Trans:Leaf.Litter, canopy, Substrate)
    availMeta$dist <- distGeo(usedMonth[[i]][[j]][,c("location.lon", "location.lat")], 
                              availMonth[[i]][[j]][,c("Lon","Lat")])
    availMonth[[i]][[j]] <- cbind(availMeta, availMonth[[i]][[j]])
    
    j = j + 1
  }
  
  i = i + 1
}
#### Random location was based on the inter-quartile range of distances moved between locations
#### within an activity area as defined by Compton et al. (2002)

#### Combine Lists of Available Data
availData = list()
for(i in 1:length(availMonth)){
  availData[[i]] <- do.call(rbind, availMonth[[i]])
  i = i + 1
}
availData <- do.call(rbind, availData)
availData <- availData[order(availData$event.id),]

#### Return Available Data Collected on the Days of Tracking 
#### beyond the initial period and ending period 
#### as a result of sampling away from the main stream
availSupp <- availUnpaired %>% drop_na(c(stream.lon, stream.lat))
for(i in 1:nrow(availSupp)){
  dataName = c("dist", "Width", "Lat.Deg", "Lat.Min", "Lon.Deg", "Lon.Min", "Lat", "Lon", 
               "Elev", "Habitat", "Water.Depth", "Water.Width", 
               "Substrate.Gravel", "Substrate.Pebble", "Substrate.Cobble", "Substrate.Boulder", 
               "Submerged", "Floating", "Emergent", "Tree.Root", "Filamentous.Algae", "Iron.Bacteria", 
               "Leaf.Litter", "canopy", "Substrate")
  suppName = c("stream.dist", "stream.width", "stream.lat.D", "stream.lat.M", "stream.lon.D", "stream.lon.M", "stream.lat", "stream.lon",
               "stream.elev", "stream.habitat", "stream.water.depth", "stream.water.width", 
               "stream.substrate.gravel", "stream.substrate.pebble", "stream.substrate.cobble", "stream.substrate.boulder", 
               "stream.submerged", "stream.floating", "stream.emergent", "stream.tree.root", "stream.filamentous.algae", "stream.iron.bacteria",
               "stream.leaf.litter", "stream.canopy", "stream.substrate")
  availData[availData$event.id == availSupp$event.id[i], dataName] <- availSupp[availSupp$event.id == availSupp$event.id[i], suppName]
}

#### Re-structured the Matched Data Points
usedData <- usedUnpaired %>% 
  dplyr::rename(
    lat.D = location.lat.D,
    lat.M = location.lat.M,
    lon.D = location.lon.D,
    lon.M = location.lon.M,
    lat = location.lat,
    lon = location.lon,
    elev = location.elev,
    habitat = location.habitat,
    airTemp = location.air.temp,
    waterTemp = location.water.temp,
    waterDepth = location.water.depth,
    streamWidth = location.water.width,
    gravel = location.substrate.gravel,
    pebble = location.substrate.pebble,
    cobble = location.substrate.cobble,
    boulder = location.substrate.boulder,
    submerged = location.submerged,
    floating = location.floating,
    emergent = location.emergent,
    root = location.tree.root,
    algae = location.filamentous.algae,
    ironBacteria = location.iron.bacteria,
    leafLitter = location.leaf.litter,
    canopy = location.canopy,
    substrate = location.substrate
  )

availData <- availData %>% 
  dplyr::rename(
    trans = Trans,
    dist = dist,
    perc = Width,
    lat.D = Lat.Deg,
    lat.M = Lat.Min,
    lon.D = Lon.Deg,
    lon.M = Lon.Min,
    lat = Lat,
    lon = Lon,
    elev = Elev,
    habitat = Habitat,
    waterDepth = Water.Depth,
    streamWidth = Water.Width,
    gravel = Substrate.Gravel,
    pebble = Substrate.Pebble,
    cobble = Substrate.Cobble,
    boulder = Substrate.Boulder,
    submerged = Submerged,
    floating = Floating,
    emergent = Emergent,
    root = Tree.Root,
    algae = Filamentous.Algae,
    ironBacteria = Iron.Bacteria,
    leafLitter = Leaf.Litter,
    canopy = canopy,
    substrate = Substrate
  )

#### Combine ALL Data 
usedAll <- rbind(usedPaired, usedData)
usedAll <- usedAll %>% dplyr::select(-c(airTemp, waterTemp, substrate))
usedAll$status <- 1

availAll <- rbind(availPaired[,!names(availPaired) %in% c("direction","airTemp", "waterTemp", "substrate")], 
                  availData[,!names(availData) %in% c("trans", "substrate")])
availAll$status <- 0
randomDist <- availAll$dist
plot(stats::density(randomDist))
randomPerc <- availAll$perc
plot(stats::density(randomPerc))

turtlesRSF <- rbind(usedAll, availAll[,!names(availAll) %in% c("dist", "perc")])
str(turtlesRSF)



#### Relocation Points: Habitat Proportion ----
locs <- turtlesRSF %>% 
  filter(status == 1) %>%
  group_by(turtle.id, habitat) %>% 
  dplyr::summarise(n=length(habitat)) %>%
  dplyr::mutate(prop = n/sum(n)*100)
locs <- as.data.frame(locs)
locs <- locs %>% dplyr::select(-n)
locs <- matrify(locs)
locs <- cbind(id = rownames(locs), locs)

#### Write Dataset for Compositional Analysis 
locs <- with(locs, data.frame('turtle.id' = id, Pool, Riffle, Run))
write.csv(locs, "data/SABE_habitatProp_loc.csv", row.names = FALSE)

#### Relocation Points: Habitat Proportion per Seasons ----
locsProp <- turtlesRSF %>% 
  filter(status == 1) %>% 
  group_by(turtle.id, habitat, season) %>%
  dplyr::summarise(n=length(habitat)) %>%
  dplyr::mutate(prop = n/sum(n)*100)
locsProp <- as.data.frame(locsProp)

locsSeasonal = list()
for(i in 1:length(mcpSeasonal)){
  locsSeasonal[[i]] <- locsProp %>% filter(season == levels(locsProp$season)[i])
  locsSeasonal[[i]] <- locsSeasonal[[i]] %>% dplyr::select(-n, -season)
  locsSeasonal[[i]] <- matrify(locsSeasonal[[i]])
  locsSeasonal[[i]] <- locsSeasonal[[i]][which(rownames(locsSeasonal[[i]]) %in% rownames(mcpSeasonal[[i]])),]
  locsSeasonal[[i]]$season <- levels(locsProp$season)[i]
  locsSeasonal[[i]] <- cbind(id = rownames(locsSeasonal[[i]]), locsSeasonal[[i]])
  i = i + 1
}

#### Write Dataset for Compositional Analysis 
locsSeasonalDf <- do.call(rbind, locsSeasonal)
locsSeasonalDf <- with(locsSeasonalDf, data.frame('turtle.id' = id, season, Pool, Riffle, Run))
write.csv(locsSeasonalDf, "data/SABE_habitatProp_loc_season.csv", row.names = FALSE)

#### Write Dataset for RSF 
turtlesRSF <- with(turtlesRSF, data.frame(event.id, turtle.id, frequency, date, time, session, month, year, season, 
                                          visible, behaviour, elev, habitat, waterDepth, streamWidth, 
                                          gravel, pebble, cobble, boulder, 'litter' = leafLitter, canopy, status))
write.csv(turtlesRSF, "data/SABE_habitatselection.csv", row.names = FALSE)