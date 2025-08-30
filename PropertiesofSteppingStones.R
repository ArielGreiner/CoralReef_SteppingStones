library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(maps)
library(ggplot2)
library(tidyverse)
library(sf)

#function below written by Cole Baird Brookson
plot.map<- function(database,center,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1]
  Obj[[2]] <- polygons[,2]
  
  map(Obj,...)
}


#set working directory to your personal working directory
setwd('/Users/akg6325/Dropbox/Github')

##Load background stuff
#load thresholded/censored connectivity matrix
connmat_reduced <- readRDS('CoralReef_SteppingStones/datafiles/censoredconnectivitymatrix.rds') #connmat_reduced, #from Greiner et al. (2022a) in Global Ecology and Biogeography

#load reefdata dataframe
reefdata <- readRDS('CoralReef_SteppingStones/datafiles/reefdata_colsaftersteppingstones_round1_4.2.2024.rds')
fullreefdata <- readRDS("CoralReef_SteppingStones/datafiles/fullreefdata_4.2.2024.rds")


world <- readOGR(dsn = "CoralReef_SteppingStones/datafiles/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")


##Load in the stepping stones
#only looking at the round 1 stepping stones because the round 2 stepping stones don't add much
steppingstones <- readRDS('CoralReef_SteppingStones/datafiles/steppingstones_round1_4.2.2024.rds')

steppingstones_reefdata <- fullreefdata[fullreefdata$reefnum %in% steppingstones,]
#they're all in Marine Realm 13
#2624, 2632 are partially protected? (see Prtcdr and pct_Prt columns)
#EEZ: Madagascar + France, Indonesia, Malaysia, Mozambique

##Adding in some pressures data as well

#going back to cells12292
cells12292 <- st_read("CoralReef_SteppingStones/datafiles/cells12292.shp")
#cells12292

#extract the rows from cells12292 that correspond to the steppingstones, using PolyNo
cells12292_ssonly <- cells12292[cells12292$PolyNo %in% steppingstones_reefdata$PolyNo,]

load("CoralReef_SteppingStones/datafiles/allreefs.RData")
#allreefs

cells12292_ssonly %>% st_transform(st_crs(allreefs)) -> cells12292_ssonly_Mercator_2SP
plot(st_geometry(cells12292_ssonly_Mercator_2SP)) # looks OK
allreefs %>% select(OBJECTID, grav_NC, pop_count, num_ports, reef_value, sediment, nutrient, cumul_score, top_threat, Region, grav_NC_raw, pop_count_raw, num_ports_raw, reef_value_raw, sediment_raw, nutrient_raw) -> allreefs_slim

#st_join with join = st_intersects leads to duplicate rows if a big grid cell contains more than one value of a small grid cell, but that's fine
#calculates the small grid cell values that are contained within the 
st_join(cells12292_ssonly_Mercator_2SP, allreefs_slim, join = st_intersects) -> cells12292_ssonly_Mercator_2SP_all 
st_write(cells12292_ssonly_Mercator_2SP_all,dsn="CoralReef_SteppingStones/datafiles/cells12292_ssonly_Mercator_2SP_all.shp")
#open in QGIS and extract .csv file and then port back into R
reefdata_ssonly_full <- read.csv(file = "CoralReef_SteppingStones/datafiles/cells12292_ssonly_Mercator_2SP_all.csv")

#extract the OBJECTIDs (some duplicates bc some small grid cells in multiple large grid cells)
steppingstones_smallgridcellIDs <- unique(reefdata_ssonly_full$OBJECTI)

allreefs_ssonly <- allreefs[allreefs$OBJECTID %in% steppingstones_smallgridcellIDs,]
#top threats: fishing (13), nitrogen (3), coastal pop (2), industrial dev (1)
#regions: western indian ocean (5), southeast asia (14)
#cumulative threat scores vary from 0.052-0.808 

#ADD PROTECTED AREA STATUS
#Downloaded Protected Area Data from protectedplanet.net on April 3, 2024
#WDPA_WDOECM_Apr2024_Public_all_csv/WDPA_WDOECM_Apr2024_Public_all_csv.csv
protectedarea_data <- read.csv(file = "CoralReef_SteppingStones/datafiles/WDPA_WDOECM_Apr2024_Public_all_csv.csv")

#just used the shape files downloaded from protectedplanet.net, and then went into QGIS -> 'ProtectedArea_Stressor_Figuring' file and 'Steppingstones_protectedareastatus_thoughts.pptx'
#the Madagascar ones don't overlap anymore? which is odd
#PolyNo 3975 overlaps with Pantai Panjang dan P.Baai nature recreation park
protectedarea_data[protectedarea_data$WDPA_PID == 555571198,]
#partial no take, managed by the government, IUCN level V "Protected Landscape/Seascape", status = designated, owned by the state, not a UNESCO world heritage site

