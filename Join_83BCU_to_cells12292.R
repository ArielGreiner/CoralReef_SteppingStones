library(tidyverse)
library(sf)

#set working directory to your personal working directory
setwd('/Users/akg6325/Dropbox/Github')

cells12292 <- st_read("CoralReef_SteppingStones/datafiles/cells12292.shp") #Wood et al. 2014 grid cells
cells12292

load("CoralReef_SteppingStones/datafiles/allreefs.RData") #Andrello et al. 2022 grid cells
allreefs

# allreefs %>% st_transform(st_crs(cells12292)) -> allreefs_Mollweide # creates some weird polygons stretching all over the world
cells12292 %>% st_transform(st_crs(allreefs)) -> cells12292_Mercator_2SP
plot(st_geometry(cells12292_Mercator_2SP)) # looks OK
allreefs %>% select(is.bcu, BCU_name) -> allreefs_slim #NOTE: BCU = Bioclimatic Unit is another term for what we called 'Climate Refugia' in the manuscript
                    
#Determining which Wood et al. (2014) grid cells can be classified as 'refugia' from the Andrello et al. (2022) grid cells 
#if more than 1/2 of the smaller grid cells within the larger grid cell are BCU then it's BCU, otherwise it's not
#just classify the BCU grid cells properly, then the rest of the grid cells will have no ID...can just manually add that in -> extracted csv from QGIS from file above, worked with it in MakingNewReefData_Dataframe.R
allreefs_rslim <- allreefs_slim[allreefs_slim$is.bcu == "refugia",]
st_join(cells12292_Mercator_2SP, allreefs_rslim, largest = T) -> cells12292_Mercator_2SP_83BCU_BCUonly
st_write(cells12292_Mercator_2SP_83BCU_BCUonly,dsn="CoralReef_SteppingStones/datafiles/cells12292_Mercator_2SP_83BCU_BCUonly.shp")
