library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(maps)
library(dplyr)

#set working directory to your personal working directory
setwd('/Users/akg6325/Dropbox/Github')

#load reefdata dataframe

#new reef data
reefdata <- read.csv(file = "CoralReef_SteppingStones/datafiles/cells12292_Mercator_2SP_83BCU_BCUonly.csv")
reefdata$is_bcu[reefdata$is_bcu != "refugia"] <- "non-refugia"

full_reefdata <- reefdata
reefdata <- full_reefdata[,c(1,2,3,4,5,27,28)] #removing the rows that contain information that i'm not using

#correct the longitude (for plotting in R)
reefdata$Longitude <- reefdata$Longitd
reefdata$Longitude[reefdata$Longitude > 180] <- reefdata$Longitude[reefdata$Longitude > 180] - 360
reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre
reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre
#remove the Longitd column just because it's confusing
reefdata$Longitd <- NULL

#make a new BCU_ID column
BCUnames <- unique(reefdata$BCU_name)
#unique(reefdata$is_bcu[reefdata$BCU_name == BCUnames[1]]) #all are non-refugia or blank so that's good
reefdata$BCU_ID <- NA
for(i in 1:length(BCUnames)){
  reefdata$BCU_ID[reefdata$BCU_name == BCUnames[i]] <- (i-1)
}

#checking things, it looks good (only one name per number)
reefdata %>%
  select(BCU_name, BCU_ID) %>%
  unique()


saveRDS(reefdata, file = "CoralReef_SteppingStones/datafiles/reefdata_11.10.2023.rds")  #old: 11.6.2023

##4.2.2024: Recovering some of the old columns with information that I need
fullreefdata <- read.csv(file = "CoralReef_SteppingStones/datafiles/cells12292_Mercator_2SP_83BCU_BCUonly.csv")
fullreefdata$is_bcu[fullreefdata$is_bcu != "refugia"] <- "non-refugia"

#correct the longitude (for plotting in R)
fullreefdata$Longitude <- fullreefdata$Longitd
fullreefdata$Longitude[fullreefdata$Longitude > 180] <- fullreefdata$Longitude[fullreefdata$Longitude > 180] - 360
fullreefdata$Longitude_corrected <- fullreefdata$Longitude
newcentre <- 180
range(fullreefdata$Longitude_corrected)
fullreefdata$Longitude_corrected[fullreefdata$Longitude > 0 & fullreefdata$Longitude  <  180] <- fullreefdata$Longitude_corrected[fullreefdata$Longitude > 0 & fullreefdata$Longitude  <  180] - newcentre
fullreefdata$Longitude_corrected[fullreefdata$Longitude < 0 & fullreefdata$Longitude  >  -180] <- fullreefdata$Longitude_corrected[fullreefdata$Longitude < 0 & fullreefdata$Longitude  >  -180] + newcentre
#remove the Longitd column just because it's confusing
fullreefdata$Longitd <- NULL

#make a new BCU_ID column
BCUnames <- unique(fullreefdata$BCU_name)
#unique(reefdata$is_bcu[reefdata$BCU_name == BCUnames[1]]) #all are non-refugia or blank so that's good
fullreefdata$BCU_ID <- NA
for(i in 1:length(BCUnames)){
  fullreefdata$BCU_ID[fullreefdata$BCU_name == BCUnames[i]] <- (i-1)
}

#checking things, it looks good (only one name per number)
fullreefdata %>%
  select(BCU_name, BCU_ID) %>%
  unique()

#add reefnum column
fullreefdata$reefnum <- seq(1,dim(fullreefdata)[1],1)


saveRDS(fullreefdata, file = "CoralReef_SteppingStones/datafiles/fullreefdata_4.2.2024.rds")

