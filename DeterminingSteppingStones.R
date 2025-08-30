library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(maps)
library(ggplot2)
library(dplyr)

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

#load thresholded/censored connectivity matrix
connmat_reduced <- readRDS('CoralReef_SteppingStones/datafiles/censoredconnectivitymatrix.rds') #connmat_reduced, #from Greiner et al. (2022a) in Global Ecology and Biogeography

#load reefdata dataframe
reefdata <- readRDS('CoralReef_SteppingStones/datafiles/reefdata_11.10.2023.rds')


world <- readOGR(dsn = "CoralReef_SteppingStones/datafiles/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")

#which reefs are in 50 reef BCUs?

fiftyreefs <- which(reefdata$BCU_ID > 0)
nonfiftyreefs <- which(reefdata$BCU_ID == 0)
#make a column in the reefdata dataset that has NAs for the reefs that aren't in BCUs at all, but otherwise looks like BCU_ID
reefdata$fiftyreef <- reefdata$BCU_ID
reefdata$fiftyreef[reefdata$BCU_ID == 0] <- NA #this seems to have worked but is throwing a bunch of error messages for some reason, anyways

##which reefs are in the same networks as the 50 reef BCUs? basically for BCU x, which networks contain at least one reef from BCU x
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
og_networks <- clusters(g_orig,mode = "weak")
reefdata$networks <- og_networks$membership #making a column that says which network each reef is in

#step 1: go through the BCUs and see which networks their member reef cells are in
numBCUs <- range(reefdata$BCU_ID)[2]
BCU_networks <- list()
allBCUntwks <- unique(reefdata$networks[reefdata$BCU_ID == 1])
for(i in 1:numBCUs){
  BCU_networks[[i]] <- unique(reefdata$networks[reefdata$BCU_ID == i])
  if(i > 1){
    allBCUntwks <- append(allBCUntwks, BCU_networks[[i]]) 
  }
}
allBCUntwks <- unique(allBCUntwks)
#^ can see that there is some overlap in BCU_networks[[i]]

#plot a map showing all of the BCUs and all of the reefs in networks touched by those BCUs
#colouring some of the networks the same specific colours as in Greiner et al. GEB 2022
set.seed(3) 
cols <- sample(rainbow(clusters(g_orig, mode = "weak")$no))
#colouring some of the networks specific colours
cols[166] <- "purple"
cols[30] <- "pink"
cols[1] <- "blue"
cols[53] <- "green"
cols[54] <- "red"
cols[19] <- "black"
cols[540] <- "brown"
cols[428] <- "Aquamarine"
cols[40] <- "SpringGreen"
cols[454] <- "orange"
cols[453] <- "violetred2"
cols[326] <- "salmon"
cols[60] <- "darkgoldenrod1"
cols[156] <- "chocolate4"
cols[170] <- "darkgreen"
cols[138] <- "mediumpurple2"
cols[162] <- "indianred2"

#FIGURE 2A
png("CoralReef_SteppingStones/NetworksofFiftyReefs_50Randothers_fixedcolours.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected, reefdata$Latitud, pch = 20, cex = 0.2, col = alpha("lightblue",0.5))
for(i in 1:length(allBCUntwks)){
  points(reefdata$Longitude_corrected[reefdata$is_bcu == "refugia" & reefdata$networks == allBCUntwks[i]], reefdata$Latitud[reefdata$is_bcu == "refugia" & reefdata$networks  == allBCUntwks[i]], pch = 20, cex = 0.2, col = cols[allBCUntwks[i]])
}
dev.off()

#how many of the networks that have BCU reefs in them have BCU reefs from multiple networks
#if a network only has BCU reefs from one BCU, then the other reefs in it aren't going to be useful as stepping stones because they can't have connections to multiple BCUs
BCUntwkoverlap <- data.frame(BCUntwk = allBCUntwks, num_BCUs = 0)
for(i in 1:length(allBCUntwks)){
  BCUntwkoverlap$num_BCUs[BCUntwkoverlap$BCUntwk == allBCUntwks[i]] <- length(unique(reefdata$fiftyreef[reefdata$networks == allBCUntwks[i] & reefdata$BCU_ID > 0]))
}
morethanoneBCU <- BCUntwkoverlap$BCUntwk[BCUntwkoverlap$num_BCUs > 1]
#BCUntwkoverlap[BCUntwkoverlap$num_BCUs > 1,]

#which BCUs are in multi-BCU networks initially? most of them, apparently
BCUsinmultiBCUntwks <- unique(reefdata$BCU_ID[reefdata$BCU_ID > 0 & reefdata$networks %in% morethanoneBCU])

#SUPPLEMENTARY MATERIALS 2 - FIGURE 1
cols_adjusted <- cols
cols_adjusted[96] <- "deepskyblue" #8.29.2025: this colour was too similar to another one
#plot the reefs in the networks that cross multiple BCUs
png("CoralReef_SteppingStones/NetworksofFiftyReefs_multiBCUnetworks.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0)) 
points(reefdata$Longitude_corrected, reefdata$Latitud, pch = 20, cex = 0.2, col = alpha("lightblue",0.5))
points(reefdata$Longitude_corrected[fiftyreefs],reefdata$Latitud[fiftyreefs],col="black",pch=20,cex=0.2)
points(reefdata$Longitude_corrected[reefdata$BCU_ID %in% BCUsinmultiBCUntwks],reefdata$Latitud[reefdata$BCU_ID %in% BCUsinmultiBCUntwks],col="darkgrey",pch=20,cex=0.2) #8.29.2025: for the reef cells that are in refugia in multi-reef networks that are not themselves in those multi-refugia networks
for(i in 1:length(morethanoneBCU)){
  points(reefdata$Longitude_corrected[reefdata$networks == morethanoneBCU[i]], reefdata$Latitud[reefdata$networks  == morethanoneBCU[i]], pch = 20, cex = 0.2, col = cols_adjusted[morethanoneBCU[i]])
}
dev.off()

#how many reefs in those candidate networks are not in BCUs
length(reefdata$PolyNo[(reefdata$networks %in% morethanoneBCU) & reefdata$BCU_ID == 0]) #4165

#####BCU-only networks
#if only the 50 reefs are left, which networks are those 50 reefs in and how many BCUs are still in the same networks as each other
reefdata$reefnum <- seq(1,dim(reefdata)[1],1) #because TARGET_FID is off by 1 and PolyNo seems to skip a number also 

connmat_fiftyr <- connmat_reduced[fiftyreefs,fiftyreefs]
g_fiftyreef <- graph.adjacency(as.matrix(connmat_fiftyr), weighted = TRUE)
fiftyr_networks <- clusters(g_fiftyreef,mode = "weak")
reefdata$fiftyr_ntwks <- NA #to fill in the non-BCU rows
reefdata$fiftyr_ntwks[reefdata$is_bcu == "refugia"] <- fiftyr_networks$membership

#step 1: go through the BCUs and see which networks their member reef cells are in
numBCUs <- range(reefdata$BCU_ID)[2]
BCU_fiftyr_networks <- list()
allBCU_fiftyrntwks <- unique(reefdata$fiftyr_ntwks[reefdata$BCU_ID == 1])
for(i in 1:numBCUs){
  BCU_fiftyr_networks[[i]] <- unique(reefdata$fiftyr_ntwks[reefdata$BCU_ID == i])
  if(i > 1){
    allBCU_fiftyrntwks <- append(allBCU_fiftyrntwks, BCU_fiftyr_networks[[i]]) 
  }
}
allBCU_fiftyrntwks <- unique(allBCU_fiftyrntwks)

#how many of the networks that have BCU reefs in them have BCU reefs from multiple networks
#if a network only has BCU reefs from one BCU, then the other reefs in it aren't going to be useful as stepping stones because they can't have connections to multiple BCUs
BCU_fiftyr_ntwkoverlap <- data.frame(BCUntwk = allBCU_fiftyrntwks, num_BCUs = 0)
for(i in 1:length(allBCU_fiftyrntwks)){
  BCU_fiftyr_ntwkoverlap$num_BCUs[BCU_fiftyr_ntwkoverlap$BCUntwk == allBCU_fiftyrntwks[i]] <- length(unique(reefdata$fiftyreef[reefdata$fiftyr_ntwks == allBCU_fiftyrntwks[i] & reefdata$BCU_ID > 0]))
}
morethanonefiftyrBCU <- BCU_fiftyr_ntwkoverlap$BCUntwk[BCU_fiftyr_ntwkoverlap$num_BCUs > 1]
#BCU_fiftyr_ntwkoverlap[BCU_fiftyr_ntwkoverlap$num_BCUs > 1,] #many more networks with multiple BCUs in them, i think because that big network got split up

cols_two <- sample(rainbow(clusters(g_fiftyreef, mode = "weak")$no)) #make a new one of these bc the colours were blending too much

#FIGURE 2B
png("CoralReef_SteppingStones/Only50Reefs_50RNetworks_fixedcolours.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
#points(reefdata$Longitude_corrected, reefdata$Latitud, pch = 20, cex = 0.2, col = alpha("lightblue",0.5))
for(i in 1:length(allBCU_fiftyrntwks)){
  points(reefdata$Longitude_corrected[reefdata$fiftyr_ntwks == allBCU_fiftyrntwks[i]], reefdata$Latitud[reefdata$fiftyr_ntwks == allBCU_fiftyrntwks[i]], pch = 20, cex = 0.2, col = cols_two[allBCU_fiftyrntwks[i]]) #cols[allBCU_fiftyrntwks[i]]
}
#points(reefdata$Longitude_corrected[fiftyreefs],reefdata$Latitud[fiftyreefs],col="black",pch=20,cex=0.2)
dev.off()

##Want to figure out which of the 50R BCUs that used to be in the same networks, are still in the same networks when only the 50R remain
#make a matrix of 1s and 0s, 1s saying that they're in the same network as each other (in the full conn mat)
fiftyr_fullconnmat_networkoverlap <- matrix(data = 0, nrow = numBCUs, ncol = numBCUs)

for(i in 1:numBCUs){
  ntwks <- BCU_networks[[i]]
  bcus <- unique(reefdata$BCU_ID[reefdata$is_bcu == "refugia" & reefdata$networks %in% ntwks]) #which BCUs also have reefs in those ntwks, will include numBCUs[i] but that's fine 
  fiftyr_fullconnmat_networkoverlap[i,bcus] <- fiftyr_fullconnmat_networkoverlap[bcus,i] <- 1 #may as well include both the upper and lower diagonal
  ntwks <- bcus <- NA
}

#make a matrix of 1s and 0s, 1s saying that they're in the same network as each other (in the 50r conn mat)
fiftyr_50rconnmat_networkoverlap <- matrix(data = 0, nrow = numBCUs, ncol = numBCUs)

for(i in 1:numBCUs){
  ntwks <- BCU_fiftyr_networks[[i]]
  bcus <- unique(reefdata$BCU_ID[reefdata$is_bcu == "refugia" & reefdata$fiftyr_ntwks %in% ntwks]) #which BCUs also have reefs in those ntwks, will include numBCUs[i] but that's fine 
  fiftyr_50rconnmat_networkoverlap[i,bcus] <- fiftyr_50rconnmat_networkoverlap[bcus,i] <- 1 #may as well include both the upper and lower diagonal
  ntwks <- bcus <- NA
}

#then see how the two matrices compare
#then turn those 1s into -1s if those two BCUs are no longer in the same network and into 2s if they're still in the same network
#if any two 50Rs are in the same network now but weren't in the initial matrix, make those 1s (0s -> 1s)
fiftyr_ba_networkoverlap <- matrix(data = 0, nrow = numBCUs, ncol = numBCUs)
BCU_overlap_networks <- list()
BCU_nonoverlap_networks <- list()

for(i in 1:numBCUs){
  for(j in 1:numBCUs){
    #still in the same network as each other 
    if(fiftyr_50rconnmat_networkoverlap[i,j] == 1 & fiftyr_fullconnmat_networkoverlap[i,j] == 1){
      fiftyr_ba_networkoverlap[i,j] <- 2
    }
    #no longer in the same network as each other
    if(fiftyr_50rconnmat_networkoverlap[i,j] == 0 & fiftyr_fullconnmat_networkoverlap[i,j] == 1){
      fiftyr_ba_networkoverlap[i,j] <- -1
    }
    #in the same network, but only when only the 50 reefs are the only reefs left (i don't think that this should be possible)
    if(fiftyr_50rconnmat_networkoverlap[i,j] == 1 & fiftyr_fullconnmat_networkoverlap[i,j] == 0){
      fiftyr_ba_networkoverlap[i,j] <- 1
    }
  }
  BCU_overlap_networks[[i]] <- which(fiftyr_ba_networkoverlap[i,] == 2) #only need to do [i,] because a symmetric matrix
  BCU_nonoverlap_networks[[i]] <- which(fiftyr_ba_networkoverlap[i,] == -1)
}
#which(fiftyr_ba_networkoverlap == 1) #none, good

#removing the self-loops
BCU_overlap_networks_nosl <- list()
for(i in 1:numBCUs){
  BCU_overlap_networks_nosl[[i]] <- BCU_overlap_networks[[i]][!(BCU_overlap_networks[[i]] == i)]
}

#plot a map showing the reefs that are in BCUs and still in the same networks as other BCU reefs when only the BCU reefs remain
png("CoralReefSteppingStones/SteppingStones_v2/Only50Reefs_overlapnetworks.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
 for(i in 1:numBCUs){
  if(length(BCU_overlap_networks[[i]]) > 1){
  points(reefdata$Longitude_corrected[reefdata$BCU_ID %in% BCU_overlap_networks[[i]]], reefdata$Latitud[reefdata$BCU_ID %in% BCU_overlap_networks[[i]]], pch = 20, cex = 0.2, col = cols[i])
  }
  if(length(BCU_overlap_networks[[i]]) == 1){
    points(reefdata$Longitude_corrected[reefdata$is_bcu == "refugia" & reefdata$BCU_ID == i], reefdata$Latitud[reefdata$is_bcu == "refugia" & reefdata$BCU_ID == i], pch = 20, cex = 0.2, col = alpha("lightblue",0.5))
  }
}
dev.off()

#which BCUs are no longer in ANY multiBCUnetworks? i.e. in nonoverlapBCUs but not in overlapBCUs
#some can be in both lists because they lost network connections with some BCUs but not all of them
overlapBCUs <- sort(unique(unlist(BCU_overlap_networks_nosl))) #still in the same networks as some other BCUs
nonoverlapBCUs <- sort(unique(unlist(BCU_nonoverlap_networks))) #no longer in the same networks as all of the other BCUs that they used to be in the same networks as
lostBCUs <- nonoverlapBCUs[which(!(nonoverlapBCUs %in% overlapBCUs))] #no longer in the same networks as any other BCUs (i.e. used to be in the same networks)
notlostBCUs <- overlapBCUs[which(!(overlapBCUs %in% nonoverlapBCUs))] #BCUs that didn't lose any of their other in-network-BCUs NOT DONE
partiallostBCUs <- overlapBCUs[which(!(overlapBCUs %in% notlostBCUs))] #BCUs that lost some in-network-BCUs but not all of them
allBCUs <- seq(1,83,1)
reallylostBCUs <- allBCUs[which(!(allBCUs %in% nonoverlapBCUs))]
reallylostBCUs <- reallylostBCUs[which(!(reallylostBCUs %in% overlapBCUs))] #these are the ones that were never in a multiBCU network, [1] 29 32 33 34 35 36 39 41 78 79

#give every BCU a colour that is easy enough to distinguish from each other and then change the plots below to use this instead #NOT DONE
BCUcols <- sample(rainbow((numBCUs*2)))

#SUPPLEMENTARY MATERIALS 3 - FIGURE 1
png("CoralReef_SteppingStones/Only50Reefs_notlostBCUs.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
for(i in 1:length(notlostBCUs)){
  points(reefdata$Longitude_corrected[reefdata$BCU_ID == notlostBCUs[i]], reefdata$Latitud[reefdata$BCU_ID == notlostBCUs[i]], pch = 20, cex = 0.2, col = cols[i])
}
dev.off()  

png("CoralReef_SteppingStones/Only50Reefs_partiallostBCUs.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
for(i in 1:length(partiallostBCUs)){
  points(reefdata$Longitude_corrected[reefdata$BCU_ID == partiallostBCUs[i]], reefdata$Latitud[reefdata$BCU_ID == partiallostBCUs[i]], pch = 20, cex = 0.2, col = cols[i]) #BCUcols[partiallostBCUs[i]]
}
dev.off()  


png("CoralReef_SteppingStones/Only50Reefs_lostBCUs.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
for(i in 1:length(lostBCUs)){
  points(reefdata$Longitude_corrected[reefdata$BCU_ID == lostBCUs[i]], reefdata$Latitud[reefdata$BCU_ID == lostBCUs[i]], pch = 20, cex = 0.2, col = cols[i])
}
dev.off()


#okay so which multiBCU networks are the 'lostBCUs' in?
lostBCUnetworks <- list()
lostBCUconnections <- list()
for(i in 1:length(lostBCUs)){
  lostBCUnetworks[[i]] <- BCU_networks[[lostBCUs[i]]][which(BCU_networks[[lostBCUs[i]]] %in% morethanoneBCU)] #because only care which multiBCU networks it was initially in
  lostBCUconnections[[i]] <- which(fiftyr_ba_networkoverlap[lostBCUs[i],] == -1) #the BCUs they lost connection with
}
lostBCUntwks <- sort(unique(unlist(lostBCUnetworks))) # 54 166

#how many reefs are in those networks and are not in BCUs? 3421
length(reefdata$reefnum[reefdata$networks %in% lostBCUntwks & reefdata$is_bcu == "non-refugia"])


#for the partiallostBCUs, which BCUs did they lose connections with? and which networks were those connections through?
#which BCUs did they lose connections with
lostconnectionBCUs <- list()
lostconnectionnetworks <- list()
for(i in 1:length(partiallostBCUs)){
lostconnectionBCUs[[i]] <- which(fiftyr_ba_networkoverlap[partiallostBCUs[i],] == -1)
for(j in 1:length(lostconnectionBCUs[[i]])){ #which networks were the partiallostBCUs and the lostconnectionBCUs both in
  lostconnectionnetworks[[i]] <- BCU_networks[[partiallostBCUs[i]]][which(BCU_networks[[partiallostBCUs[i]]] %in% BCU_networks[[lostconnectionBCUs[[i]][j]]])]
  lostconnectionnetworks[[i]] <- lostconnectionnetworks[[i]][which(lostconnectionnetworks[[i]] %in% morethanoneBCU)] #bc also has to be one of the original multiBCUnetworks
}
}
lostconnectionntwks <- sort(unique(unlist(lostconnectionnetworks))) #53  54 166

#so the lostBCUnetworks are the priority for finding joiner reefs but the lostconnectionntwks are also relevant to look in
#lostBCUntwks # 54 166 #lostconnectionntwks #53 54 166 #overlap: 54,166 #non-overlap: 53
#need to figure out which of the reefs in
#the lostBCUnetworks have connections between the lostBCUs + lostBCUconnections and 
#the lostconnectionnetworks have connections between the partiallostBCUs + lostconnectionBCUs
#need to look at the neighbourhoods of each of the non-BCU-reefs in these networks and see which ones have reefs from the relevant BCUs in them

finallostBCU_candidates <- list() 
for(i in 1:length(lostBCUs)){
  nonBCUreefs_inntwk <- reefdata$reefnum[reefdata$is_bcu == "non-refugia" & reefdata$networks %in% lostBCUnetworks[[i]]]
  #which of the reefs in lostBCUs[i] are in the reef ntwks that used to also include other BCU reefs (i.e. the lostBCUnetworks)
  lostBCUreefs_inntwk <- reefdata$reefnum[reefdata$is_bcu == "refugia" & reefdata$networks %in% lostBCUnetworks[[i]] & reefdata$BCU_ID == lostBCUs[i]]
  #which of the reefs in lostBCUconnections[[i]] are in the reef ntwks that used to also include the lostBCUreefs_inntwk
  lostBCUreefcnctns_inntwk <- reefdata$reefnum[reefdata$is_bcu == "refugia" & reefdata$networks %in% lostBCUnetworks[[i]] & reefdata$BCU_ID %in% lostBCUconnections[[i]]]
  
  #which of the nonBCUreefs_inntwk are connected with lostBCUreefs_inntwk 
  lostBCU_candidates <- -1 #need to remove this later
  for(j in 1:length(nonBCUreefs_inntwk)){
    if(sum(ego(g_orig,order = 1, nodes = nonBCUreefs_inntwk[j], mode = "all")[[1]] %in% lostBCUreefs_inntwk) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
      lostBCU_candidates <- append(lostBCU_candidates, nonBCUreefs_inntwk[j])
      if(lostBCU_candidates[1] == -1){
        lostBCU_candidates <- lostBCU_candidates[lostBCU_candidates >0] #to remove the -1
      }
    }
  }
  finallostBCU_candidates[[i]] <- -1
  #and also with at least one reef from lostBCUreefcnctns_inntwk
  for(j in 1:length(lostBCU_candidates)){
    if(sum(ego(g_orig,order = 1, nodes = lostBCU_candidates[j], mode = "all")[[1]] %in% lostBCUreefcnctns_inntwk) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
      finallostBCU_candidates[[i]] <- append(finallostBCU_candidates[[i]], lostBCU_candidates[j])
      if(finallostBCU_candidates[[i]][1] == -1){
        finallostBCU_candidates[[i]] <- finallostBCU_candidates[[i]][finallostBCU_candidates[[i]] >0] #to remove the -1
      }
    }
  }
  #don't think i care which networks these finallostBCU_candidates are connecting at this juncture, just that their neighbourhoods contain reefs from multiple BCUs and thus if they were brought back in...they would reconnect some BCUs together
}
lostBCU_candidates <- NA

#now for the partiallostBCUs
finalpartiallostBCU_candidates <- list() 
for(i in 1:length(partiallostBCUs)){
  nonBCUreefs_inntwk <- reefdata$reefnum[reefdata$is_bcu == "non-refugia" & reefdata$networks %in% lostconnectionnetworks[[i]]]
  #which of the reefs in lostBCUs[i] are in the reef ntwks that used to also include other BCU reefs (i.e. the lostBCUnetworks)
  lostBCUreefs_inntwk <- reefdata$reefnum[reefdata$is_bcu == "refugia" & reefdata$networks %in% lostconnectionnetworks[[i]] & reefdata$BCU_ID == partiallostBCUs[i]]
  #which of the reefs in lostconnectionBCUs[[i]] are in the reef ntwks that used to also include the lostBCUreefs_inntwk
  lostBCUreefcnctns_inntwk <- reefdata$reefnum[reefdata$is_bcu == "refugia" & reefdata$networks %in% lostconnectionnetworks[[i]] & reefdata$BCU_ID %in% lostconnectionBCUs[[i]]]
  
  #which of the nonBCUreefs_inntwk are connected with lostBCUreefs_inntwk 
  lostBCU_candidates <- -1 #need to remove this later
  for(j in 1:length(nonBCUreefs_inntwk)){
    if(sum(ego(g_orig,order = 1, nodes = nonBCUreefs_inntwk[j], mode = "all")[[1]] %in% lostBCUreefs_inntwk) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
      lostBCU_candidates <- append(lostBCU_candidates, nonBCUreefs_inntwk[j])
      if(lostBCU_candidates[1] == -1){
        lostBCU_candidates <- lostBCU_candidates[lostBCU_candidates >0] #to remove the -1
      }
    }
  }
  finalpartiallostBCU_candidates[[i]] <- -1
  #and also with at least one reef from lostBCUreefcnctns_inntwk
  for(j in 1:length(lostBCU_candidates)){
    if(sum(ego(g_orig,order = 1, nodes = lostBCU_candidates[j], mode = "all")[[1]] %in% lostBCUreefcnctns_inntwk) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
      finalpartiallostBCU_candidates[[i]] <- append(finalpartiallostBCU_candidates[[i]], lostBCU_candidates[j])
      if(finalpartiallostBCU_candidates[[i]][1] == -1){
        finalpartiallostBCU_candidates[[i]] <- finalpartiallostBCU_candidates[[i]][finalpartiallostBCU_candidates[[i]] >0] #to remove the -1
      }
    }
  }
  #don't think i care which networks these finalpartiallostBCU_candidates are connecting at this juncture, just that their neighbourhoods contain reefs from multiple BCUs and thus if they were brought back in...they would reconnect some BCUs together
}
lostBCU_candidates <- NA
#okay so some of them aren't connectable by just one reef, lots of -1s in finalpartiallostBCU_candidates and finallostBCU_candidates
steppingstonecandidates_1stpriority <- sort(unique(unlist(finallostBCU_candidates)))
steppingstonecandidates_1stpriority <- steppingstonecandidates_1stpriority[steppingstonecandidates_1stpriority > 0] #removing the -1

steppingstonecandidates_2ndpriority <- sort(unique(unlist(finalpartiallostBCU_candidates))) #interesting how many repeat with the 1st priority ones
steppingstonecandidates_2ndpriority <- steppingstonecandidates_2ndpriority[steppingstonecandidates_2ndpriority > 0] #removing the -1

steppingstonecandidates_roundone <- steppingstonecandidates_1stpriority
steppingstonecandidates_roundone <- sort(unique(append(steppingstonecandidates_roundone, steppingstonecandidates_2ndpriority)))

saveRDS(steppingstonecandidates_roundone, file = 'CoralReef_SteppingStones/datafiles/steppingstones_round1_4.2.2024.rds') 
saveRDS(reefdata, file = 'CoralReef_SteppingStones/datafiles/reefdata_colsaftersteppingstones_round1_4.2.2024.rds')

##Add those reefs back in and make a map showing where they are and then plot the new networks

#FIGURE 4A
#fixed colours stepping stone candidate plot
cols_fix <- cols_two
cols_fix[8] <- "#cc47e9"
png("CoralReef_SteppingStones/SteppingStoneCandidates_fixedcolours.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002, 
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0)) #"gainsboro"
for(i in 1:length(steppingstonecandidates_roundone)){
  points(reefdata$Longitude_corrected[reefdata$reefnum == steppingstonecandidates_roundone[i]], reefdata$Latitud[reefdata$reefnum == steppingstonecandidates_roundone[i]], pch = 8, col = cols_fix[i]) 
}  
dev.off()

#8.19.2025: map of the stepping stone candidates, fixed colours
worldmap <- map_data ("world", wrap = c(0, 360))

#add only the stepping stone candidates back in
fiftyreefs_andss <- sort(append(fiftyreefs,steppingstonecandidates_roundone))
connmat_fiftyr_ss <- connmat_reduced[fiftyreefs_andss,fiftyreefs_andss]
g_fiftyreef_ss <- graph.adjacency(as.matrix(connmat_fiftyr_ss), weighted = TRUE)
fiftyr_ss_networks <- clusters(g_fiftyreef_ss,mode = "weak")
reefdata$fiftyr_ss_ntwks <- NA #to fill in the non-BCU rows
reefdata$fiftyr_ss_ntwks[reefdata$reefnum %in% fiftyreefs_andss] <- fiftyr_ss_networks$membership

#FIGURE 2C
png("CoralReef_SteppingStones/FiftyR_SS1_Networks_fixedcolours.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000 #remade on 4.9.2024 bc the one before was mislabeled as this i think
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
for(i in 1:fiftyr_ss_networks$no){
  points(reefdata$Longitude_corrected[reefdata$fiftyr_ss_ntwks == i], reefdata$Latitud[reefdata$fiftyr_ss_ntwks == i], pch = 20, cex = 0.2, col = cols[i])
}
points(reefdata$Longitude_corrected[reefdata$reefnum %in% steppingstonecandidates_roundone], reefdata$Latitud[reefdata$reefnum %in% steppingstonecandidates_roundone], pch = 8, col = "black")
dev.off()

#comparing the networks before and after, 123 networks vs 113 networks
png("CoralReefSteppingStones/SteppingStones_v2/FiftyR_Networks_Hist.png")
hist(fiftyr_networks$csize)
dev.off()

png("CoralReefSteppingStones/SteppingStones_v2/FiftyR_ss1_Networks_Hist.png")
hist(fiftyr_ss_networks$csize)
dev.off()

quantile(fiftyr_networks$csize)
#0%  25%  50%  75% 100% 
#1    1    2   34  455 
quantile(fiftyr_ss_networks$csize)
#  0%  25%  50%  75% 100% 
#1    1    2   36  455 

#how many multiBCU networks are there now?

#step 1: go through the BCUs and see which networks their member reef cells are in
numBCUs <- range(reefdata$BCU_ID)[2]
BCU_fiftyr_ss_networks <- list()
allBCU_fiftyr_ss_ntwks <- unique(reefdata$fiftyr_ss_ntwks[reefdata$BCU_ID == 1])
for(i in 1:numBCUs){
  BCU_fiftyr_ss_networks[[i]] <- unique(reefdata$fiftyr_ss_ntwks[reefdata$BCU_ID == i])
  if(i > 1){
    allBCU_fiftyr_ss_ntwks <- append(allBCU_fiftyr_ss_ntwks, BCU_fiftyr_ss_networks[[i]]) 
  }
}
allBCU_fiftyr_ss_ntwks <- unique(allBCU_fiftyr_ss_ntwks)

#how many of the networks that have BCU reefs in them have BCU reefs from multiple networks
#if a network only has BCU reefs from one BCU, then the other reefs in it aren't going to be useful as stepping stones because they can't have connections to multiple BCUs
BCU_fiftyr_ss_ntwkoverlap <- data.frame(BCUntwk = allBCU_fiftyr_ss_ntwks, num_BCUs = 0)
for(i in 1:length(allBCU_fiftyr_ss_ntwks)){
  BCU_fiftyr_ss_ntwkoverlap$num_BCUs[BCU_fiftyr_ss_ntwkoverlap$BCUntwk == allBCU_fiftyr_ss_ntwks[i]] <- length(unique(reefdata$fiftyreef[reefdata$fiftyr_ss_ntwks == allBCU_fiftyr_ss_ntwks[i] & reefdata$BCU_ID > 0]))
}
morethanonefiftyrssBCU <- BCU_fiftyr_ss_ntwkoverlap$BCUntwk[BCU_fiftyr_ss_ntwkoverlap$num_BCUs > 1]
#BCU_fiftyr_ss_ntwkoverlap[BCU_fiftyr_ss_ntwkoverlap$num_BCUs > 1,] 


#make a matrix of 1s and 0s, 1s saying that they're in the same network as each other (in the 50r conn mat)
fiftyr_50rssconnmat_networkoverlap <- matrix(data = 0, nrow = numBCUs, ncol = numBCUs)

for(i in 1:numBCUs){
  ntwks <- BCU_fiftyr_ss_networks[[i]]
  bcus <- unique(reefdata$BCU_ID[reefdata$is_bcu == "refugia" & reefdata$fiftyr_ss_ntwks %in% ntwks]) #which BCUs also have reefs in those ntwks, will include numBCUs[i] but that's fine 
  fiftyr_50rssconnmat_networkoverlap[i,bcus] <- fiftyr_50rssconnmat_networkoverlap[bcus,i] <- 1 #may as well include both the upper and lower diagonal
  ntwks <- bcus <- NA
}

#then see how the two matrices compare
#then turn those 1s into -1s if those two BCUs are no longer in the same network and into 2s if they're still in the same network
#if any two 50Rs are in the same network now but weren't in the initial matrix, make those 1s (0s -> 1s)
fiftyr_ba_ss_networkoverlap <- matrix(data = 0, nrow = numBCUs, ncol = numBCUs)
BCU_ss_overlap_networks <- list()
BCU_ss_nonoverlap_networks <- list()

for(i in 1:numBCUs){
  for(j in 1:numBCUs){
    #still in the same network as each other 
    if(fiftyr_50rssconnmat_networkoverlap[i,j] == 1 & fiftyr_fullconnmat_networkoverlap[i,j] == 1){
      fiftyr_ba_ss_networkoverlap[i,j] <- 2
    }
    #no longer in the same network as each other
    if(fiftyr_50rssconnmat_networkoverlap[i,j] == 0 & fiftyr_fullconnmat_networkoverlap[i,j] == 1){
      fiftyr_ba_ss_networkoverlap[i,j] <- -1
    }
    #in the same network, but only when only the 50 reefs are the only reefs left (i don't think that this should be possible)
    if(fiftyr_50rssconnmat_networkoverlap[i,j] == 1 & fiftyr_fullconnmat_networkoverlap[i,j] == 0){
      fiftyr_ba_ss_networkoverlap[i,j] <- 1
    }
  }
  BCU_ss_overlap_networks[[i]] <- which(fiftyr_ba_ss_networkoverlap[i,] == 2) #only need to do [i,] because a symmetric matrix
  BCU_ss_nonoverlap_networks[[i]] <- which(fiftyr_ba_ss_networkoverlap[i,] == -1)
}
#which(fiftyr_ba_ss_networkoverlap == 1) #none, good

#removing the self-loops
BCU_ss_overlap_networks_nosl <- list()
for(i in 1:numBCUs){
  BCU_ss_overlap_networks_nosl[[i]] <- BCU_ss_overlap_networks[[i]][!(BCU_ss_overlap_networks[[i]] == i)]
}

#which BCUs are no longer in ANY multiBCUnetworks? i.e. in nonoverlap_ss_BCUs but not in overlap_ss_BCUs
#some can be in both lists because they lost network connections with some BCUs but not all of them
overlap_ss_BCUs <- sort(unique(unlist(BCU_ss_overlap_networks_nosl))) #still in the same networks as some other BCUs
nonoverlap_ss_BCUs <- sort(unique(unlist(BCU_ss_nonoverlap_networks))) #no longer in the same networks as all of the other BCUs that they used to be in the same networks as
lost_ss_BCUs <- nonoverlap_ss_BCUs[which(!(nonoverlap_ss_BCUs %in% overlap_ss_BCUs))] #no longer in the same networks as any other BCUs (i.e. used to be in the same networks)
notlost_ss_BCUs <- overlap_ss_BCUs[which(!(overlap_ss_BCUs %in% nonoverlap_ss_BCUs))] #BCUs that didn't lose any of their other in-network-BCUs NOT DONE
partiallost_ss_BCUs <- overlap_ss_BCUs[which(!(overlap_ss_BCUs %in% notlost_ss_BCUs))] #BCUs that lost some in-network-BCUs but not all of them
#^ better than before (more overlapping, less lost)

#SUPPLEMENTARY MATERIALS 3 - FIGURE 2
png("CoralReef_SteppingStones/50Reefsandss_notlost_ss_BCUs.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
for(i in 1:length(notlost_ss_BCUs)){
  points(reefdata$Longitude_corrected[reefdata$BCU_ID == notlost_ss_BCUs[i]], reefdata$Latitud[reefdata$BCU_ID == notlost_ss_BCUs[i]], pch = 20, cex = 0.2, col = cols[i])
}
dev.off()  

png("CoralReef_SteppingStones/50Reefsandss_partiallost_ss_BCUs.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
for(i in 1:length(partiallost_ss_BCUs)){
  points(reefdata$Longitude_corrected[reefdata$BCU_ID == partiallost_ss_BCUs[i]], reefdata$Latitud[reefdata$BCU_ID == partiallost_ss_BCUs[i]], pch = 20, cex = 0.2, col = cols[i])
}
dev.off()  


png("CoralReef_SteppingStones/50Reefsandss_lost_ss_BCUs.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
for(i in 1:length(lost_ss_BCUs)){
  points(reefdata$Longitude_corrected[reefdata$BCU_ID == lost_ss_BCUs[i]], reefdata$Latitud[reefdata$BCU_ID == lost_ss_BCUs[i]], pch = 20, cex = 0.2, col = cols[i])
}
dev.off()

#okay so which multiBCU networks are the 'lost_ss_BCUs' in?
lostssBCUnetworks <- list()
lostssBCUconnections <- list()
for(i in 1:length(lost_ss_BCUs)){
  lostssBCUnetworks[[i]] <- BCU_networks[[lost_ss_BCUs[i]]][which(BCU_networks[[lost_ss_BCUs[i]]] %in% morethanoneBCU)]
  lostssBCUconnections[[i]] <- which(fiftyr_ba_ss_networkoverlap[lost_ss_BCUs[i],] == -1) #the BCUs they lost connection with
}
lostssBCUntwks <- sort(unique(unlist(lostssBCUnetworks))) #166
#basically they're all completely in 166 except for BCU 65 which has two reef cells in network 316

#how many reefs are in those networks and are not in BCUs? 3162
length(reefdata$reefnum[reefdata$networks %in% lostssBCUntwks & reefdata$is_bcu == "non-refugia"])


#for the partiallost_ss_BCUs, which BCUs did they lose connections with? and which networks were those connections through?
lostconnection_ss_BCUs <- list()
lostconnection_ss_networks <- list()
for(i in 1:length(partiallost_ss_BCUs)){
  lostconnection_ss_BCUs[[i]] <- which(fiftyr_ba_ss_networkoverlap[partiallost_ss_BCUs[i],] == -1)
  for(j in 1:length(lostconnection_ss_BCUs[[i]])){ #which networks were the partiallostBCUs and the lostconnectionBCUs both in
    lostconnection_ss_networks[[i]] <- BCU_networks[[partiallost_ss_BCUs[i]]][which(BCU_networks[[partiallost_ss_BCUs[i]]] %in% BCU_networks[[lostconnection_ss_BCUs[[i]][j]]])]
    lostconnection_ss_networks[[i]] <- lostconnection_ss_networks[[i]][which(lostconnection_ss_networks[[i]] %in% morethanoneBCU)] #bc also has to be one of the original multiBCUnetworks
  }
}
lostconnection_ss_ntwks <- sort(unique(unlist(lostconnection_ss_networks))) #53 166

#9.26.2024: calculate some grid cell level stats about all of this
#how many grid cells are in each BCU?
BCU_size <- reefdata %>%
  count(BCU_ID)

#which BCUs are initially in multi-BCU networks, and how many grid cells is that?
#which BCUs are in multi-BCU networks initially? most of them, apparently
#BCUsinmultiBCUntwks <- unique(reefdata$BCU_ID[reefdata$BCU_ID > 0 & reefdata$networks %in% morethanoneBCU])
BCU_size$originmultiBCUntwk <- 0
BCU_size$originmultiBCUntwk[BCU_size$BCU_ID %in% BCUsinmultiBCUntwks] <- 1
sum(BCU_size$n*BCU_size$originmultiBCUntwk) #3760 grid cells

#which BCUs are still in some multi-BCU networks if only BCU cells left, and how many grid cells is that?
BCUsinmultiBCUntwks_only50Rleft <- unique(reefdata$BCU_ID[reefdata$BCU_ID > 0 & reefdata$fiftyr_ntwks %in% morethanonefiftyrBCU])
#overlapBCUs <- sort(unique(unlist(BCU_overlap_networks_nosl))) #still in the same networks as some other BCUs; same as the above
BCU_size$fiftyRonly_inmultiBCUntwk <- 0
BCU_size$fiftyRonly_inmultiBCUntwk[BCU_size$BCU_ID %in% BCUsinmultiBCUntwks_only50Rleft] <- 1
sum(BCU_size$n*BCU_size$fiftyRonly_inmultiBCUntwk) #3388 grid cells

#which BCUs are in some multi-BCU networks if BCU + stepping stone cells left, and how many grid cells is that?
#overlap_ss_BCUs <- sort(unique(unlist(BCU_ss_overlap_networks_nosl))) #still in the same networks as some other BCUs
BCU_size$fiftyRss_inmultiBCUntwk <- 0
BCU_size$fiftyRss_inmultiBCUntwk[BCU_size$BCU_ID %in% overlap_ss_BCUs] <- 1
sum(BCU_size$n*BCU_size$fiftyRss_inmultiBCUntwk) #3501 grid cells

#What about the original red (54) and purple (166) networks? how many BCU grid cells were initially connected in those two networks?
#unique(reefdata$BCU_name[reefdata$networks == 54])
dim(reefdata[reefdata$networks == 166 & reefdata$is_bcu == "refugia",]) #2332 BCU grid cells in 166 #purple
dim(reefdata[reefdata$networks == 54 & reefdata$is_bcu == "refugia",]) #249 BCU grid cells in 54 #red
#which BCUs were in those?
redBCUs <- unique(reefdata$BCU_ID[reefdata$networks == 54 & reefdata$is_bcu == "refugia"]) #8 BCUs: 17 18 19 21 23 26 27 28 
purpleBCUs <- unique(reefdata$BCU_ID[reefdata$networks == 166 & reefdata$is_bcu == "refugia"]) #39 BCUs: 37 38 40 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77
#10.24.2024: how many grid cells are in those BCUs in total? 
dim(reefdata[reefdata$BCU_ID %in% redBCUs & reefdata$is_bcu == "refugia",]) #261 grid cells
dim(reefdata[reefdata$BCU_ID %in% purpleBCUs & reefdata$is_bcu == "refugia",]) #2379 grid cells

#which of the red and purple BCUs are still in multi-BCU ntwks if only BCU cells left? 
redBCUs[which(redBCUs %in% BCUsinmultiBCUntwks_only50Rleft)] #6 BCUs: 17 18 19 21 23 27
purpleBCUs[which(purpleBCUs %in% BCUsinmultiBCUntwks_only50Rleft)] #34 BCUs: 37 38 40 42 43 44 46 47 48 49 50 51 52 53 54 56 57 58 59 60 61 62 63 64 66 67 68 69 70 71 73 74 75 76
#how many grid cells is that
dim(reefdata[reefdata$BCU_ID %in% redBCUs[which(redBCUs %in% BCUsinmultiBCUntwks_only50Rleft)] & reefdata$is_bcu == "refugia",]) #181 grid cells
dim(reefdata[reefdata$BCU_ID %in% purpleBCUs[which(purpleBCUs %in% BCUsinmultiBCUntwks_only50Rleft)] & reefdata$is_bcu == "refugia",]) #2087 grid cells

#which of the red and purple BCUs are still in multi-BCU ntwks if BCU + ss cells left? 
#how many grid cells is that 
dim(reefdata[reefdata$BCU_ID %in% redBCUs[which(redBCUs %in% overlap_ss_BCUs)] & reefdata$is_bcu == "refugia",]) #261 grid cells
dim(reefdata[reefdata$BCU_ID %in% purpleBCUs[which(purpleBCUs %in% overlap_ss_BCUs)] & reefdata$is_bcu == "refugia",]) #2120 grid cells

#which of the red and purple BCUs are still in multi-BCU ntwks if only BCU cells + stepping stones left? 
redBCUs[which(redBCUs %in% overlap_ss_BCUs)] #8 BCUs: 17 18 19 21 23 26 27 28
purpleBCUs[which(purpleBCUs %in% overlap_ss_BCUs)] #35 BCUs: 37 38 40 42 43 44 45 46 47 48 49 50 51 52 53 54 56 57 58 59 60 61 62 63 64 66 67 68 69 70 71 73 74 75 76

#8.20.2025: How did average network size change as non-refugia reef cells were lost and then stepping stones added back in?
#each grid cell = 324km^2
#original, but only considering the networks that the 50 reefs are in; median = 2
quantile(as.data.frame(table(reefdata$networks[!is.na(reefdata$fiftyreef)]))$Freq) #1,1,2,53.25,2332
mean(as.data.frame(table(reefdata$networks[!is.na(reefdata$fiftyreef)]))$Freq) #78.28571
#when only 50 reef networks left: median = 2
quantile(as.data.frame(table(reefdata$fiftyr_ntwks))$Freq) #1,1,2,34,455
mean(as.data.frame(table(reefdata$fiftyr_ntwks))$Freq) #35.64228
#50 reef networks + ss: median = 2
quantile(as.data.frame(table(reefdata$fiftyr_ss_ntwks))$Freq) #1,1,2,36,455
mean(as.data.frame(table(reefdata$fiftyr_ss_ntwks))$Freq) #38.88496

#8.20.2025: 9 refugia become more connected to other Refugia, 3 Refugia become connected to at least one other Refugia again
#which are these 9 and these 3?
#3: 2 lost refugia to not lost, 1 lost to partially lost
#2 lost to not lost: 26,28 
#1 lost to partially lost: 45
#26, 28, 45: 36+44+33 grid cells => 113 grid cells

#9: 6 partially lost to not lost (17,18,19,21,23,27), 2 lost to not lost (26,28), 1 lost to partially lost (45)
#removed from lost: 26,28,45
#added to notlost: 17,18,19,21,23,26,27,28
#removed from notlost: NA
#removed from partially lost: 17,18,19,21,23,27
#added to partially lost: 45
#17,18,19,21,23,26,27,28,45: 47+44+27+28+21+36+14+44+33 => 294 grid cells




#SUPPLEMENTARY MATERIALS 5 - PAIRED STEPPING STONE ANALYSIS
##then go in and try to figure out which sets of 2 reef cells could connect the ones still left unconnected? for lost_ss and partiallost_ss
#need to find non-BCU reef cells in network 166 that link the lost_ss_BCUs back to the lostssBCUconnections BCUs (they all used to be connected to each other, and more)
#need to find non-BCU reef cells in the lostconnection_ss_networks that link the partiallost_ss_BCUs back to the lostconnection_ss_BCUs
#figure out which non-BCU reef cells have a connection with a lostBCU or a former-in-network-BCU and then see if any of those are themselves linked together?

#lost_ss_BCUs first
finallostBCUss_candidates <- list() 
for(i in 1:length(lost_ss_BCUs)){
  nonBCUreefs_inntwk <- reefdata$reefnum[reefdata$is_bcu == "non-refugia" & reefdata$networks %in% lostssBCUnetworks[[i]]]
  #which of the reefs in lostBCUs[i] are in the reef ntwks that used to also include other BCU reefs (i.e. the lostBCUnetworks)
  lostBCUreefs_inntwk <- reefdata$reefnum[reefdata$is_bcu == "refugia" & reefdata$networks %in% lostssBCUnetworks[[i]] & reefdata$BCU_ID == lost_ss_BCUs[i]]
  #which of the reefs in lostBCUconnections[[i]] are in the reef ntwks that used to also include the lostBCUreefs_inntwk
  lostBCUreefcnctns_inntwk <- reefdata$reefnum[reefdata$is_bcu == "refugia" & reefdata$networks %in% lostssBCUnetworks[[i]] & reefdata$BCU_ID %in% lostssBCUconnections[[i]]]
  
  #which of the nonBCUreefs_inntwk are connected with lostBCUreefs_inntwk 
  lostBCU_candidates <- -1 #need to remove this later
  for(j in 1:length(nonBCUreefs_inntwk)){
    if(sum(ego(g_orig,order = 1, nodes = nonBCUreefs_inntwk[j], mode = "all")[[1]] %in% lostBCUreefs_inntwk) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
      lostBCU_candidates <- append(lostBCU_candidates, nonBCUreefs_inntwk[j])
      if(lostBCU_candidates[1] == -1){
        lostBCU_candidates <- lostBCU_candidates[lostBCU_candidates >0] #to remove the -1
      }
    }
  }
  
  #which of the nonBCUreefs_inntwk are connected with lostBCUreefcnctns_inntwk 
  alsolostBCU_candidates <- -1 #need to remove this later
  for(j in 1:length(nonBCUreefs_inntwk)){
    if(sum(ego(g_orig,order = 1, nodes = nonBCUreefs_inntwk[j], mode = "all")[[1]] %in% lostBCUreefcnctns_inntwk) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
      alsolostBCU_candidates <- append(alsolostBCU_candidates, nonBCUreefs_inntwk[j])
      if(alsolostBCU_candidates[1] == -1){
        alsolostBCU_candidates <- alsolostBCU_candidates[alsolostBCU_candidates >0] #to remove the -1
      }
    }
  }
  
  #which of the lostBCU_candidates are linked together with at least one alsolostBCU_candidates, those two reef cells could act (in a pair) as a stepping stone
  #i.e. which of the lostBCU_candidates have at least one alsolostBCU_candidates in their neighbourhood (don't have to search the alsolostBCU_candidates neighbourhoods bc doing a non-directional search)
  finallostBCUss_candidates[[i]] <- -1
  #and also with at least one reef from lostBCUreefcnctns_inntwk
  for(j in 1:length(lostBCU_candidates)){
    for(k in 1:length(alsolostBCU_candidates)){
      if(sum(ego(g_orig,order = 1, nodes = lostBCU_candidates[j], mode = "all")[[1]] %in% alsolostBCU_candidates[k]) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
        finallostBCUss_candidates[[i]] <- append(finallostBCUss_candidates[[i]], c(lostBCU_candidates[j],alsolostBCU_candidates[k]))
        if(finallostBCUss_candidates[[i]][1] == -1){
          finallostBCUss_candidates[[i]] <- finallostBCUss_candidates[[i]][finallostBCUss_candidates[[i]] >0] #to remove the -1
        }
      }
    }
  }
  #don't think i care which networks these finallostBCU_candidates are connecting at this juncture, just that their neighbourhoods contain reefs from multiple BCUs and thus if they were brought back in...they would reconnect some BCUs together
}
#finallostBCUss_candidates are all -1 :(

#now for the partiallostBCUs
#need to find non-BCU reef cells in the lostconnection_ss_networks that link the partiallost_ss_BCUs back to the lostconnection_ss_BCUs
finalpartiallostBCUss_candidates_a <- list() 
finalpartiallostBCUss_candidates_b <- list() 
for(i in 1:length(partiallost_ss_BCUs)){
  nonBCUreefs_inntwk <- reefdata$reefnum[reefdata$is_bcu == "non-refugia" & reefdata$networks %in% lostconnection_ss_networks[[i]]]
  #which of the reefs in lostBCUs[i] are in the reef ntwks that used to also include other BCU reefs (i.e. the lostBCUnetworks)
  lostBCUreefs_inntwk <- reefdata$reefnum[reefdata$is_bcu == "refugia" & reefdata$networks %in% lostconnection_ss_networks[[i]] & reefdata$BCU_ID == partiallost_ss_BCUs[i]]
  #which of the reefs in lostconnectionBCUs[[i]] are in the reef ntwks that used to also include the lostBCUreefs_inntwk
  lostBCUreefcnctns_inntwk <- reefdata$reefnum[reefdata$is_bcu == "refugia" & reefdata$networks %in% lostconnection_ss_networks[[i]] & reefdata$BCU_ID %in% lostconnection_ss_BCUs[[i]]]
  
  #which of the nonBCUreefs_inntwk are connected with lostBCUreefs_inntwk 
  lostBCU_candidates <- -1 #need to remove this later
  for(j in 1:length(nonBCUreefs_inntwk)){
    if(sum(ego(g_orig,order = 1, nodes = nonBCUreefs_inntwk[j], mode = "all")[[1]] %in% lostBCUreefs_inntwk) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
      lostBCU_candidates <- append(lostBCU_candidates, nonBCUreefs_inntwk[j])
      if(lostBCU_candidates[1] == -1){
        lostBCU_candidates <- lostBCU_candidates[lostBCU_candidates >0] #to remove the -1
      }
    }
  }
  
  #which of the nonBCUreefs_inntwk are connected with lostBCUreefs_inntwk 
  alsolostBCU_candidates <- -1 #need to remove this later
  for(j in 1:length(nonBCUreefs_inntwk)){
    if(sum(ego(g_orig,order = 1, nodes = nonBCUreefs_inntwk[j], mode = "all")[[1]] %in% lostBCUreefcnctns_inntwk) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
      alsolostBCU_candidates <- append(alsolostBCU_candidates, nonBCUreefs_inntwk[j])
      if(alsolostBCU_candidates[1] == -1){
        alsolostBCU_candidates <- alsolostBCU_candidates[alsolostBCU_candidates >0] #to remove the -1
      }
    }
  }
  
  finalpartiallostBCUss_candidates_a[[i]] <- -1
  finalpartiallostBCUss_candidates_b[[i]] <- -1
  #which of the lostBCU_candidates are linked together with at least one alsolostBCU_candidates, those two reef cells could act (in a pair) as a stepping stone
  #i.e. which of the lostBCU_candidates have at least one alsolostBCU_candidates in their neighbourhood (don't have to search the alsolostBCU_candidates neighbourhoods bc doing a non-directional search)
  for(j in 1:length(lostBCU_candidates)){
    for(k in 1:length(alsolostBCU_candidates)){
      if(sum(ego(g_orig,order = 1, nodes = lostBCU_candidates[j], mode = "all")[[1]] %in% alsolostBCU_candidates[k]) > 0){ #don't need to do the as.character() in this case bc they went into the graph structure in order
        finalpartiallostBCUss_candidates_a[[i]] <- append(finalpartiallostBCUss_candidates_a[[i]], lostBCU_candidates[j])
        finalpartiallostBCUss_candidates_b[[i]] <- append(finalpartiallostBCUss_candidates_b[[i]],alsolostBCU_candidates[k])
        if(finalpartiallostBCUss_candidates_a[[i]][1] == -1){
          finalpartiallostBCUss_candidates_a[[i]] <- finalpartiallostBCUss_candidates_a[[i]][finalpartiallostBCUss_candidates_a[[i]] >0] #to remove the -1
        }
        if(finalpartiallostBCUss_candidates_b[[i]][1] == -1){
          finalpartiallostBCUss_candidates_b[[i]] <- finalpartiallostBCUss_candidates_b[[i]][finalpartiallostBCUss_candidates_b[[i]] >0] #to remove the -1
        }
      }
    }
  }
  
  #don't think i care which networks these finalpartiallostBCU_candidates are connecting at this juncture, just that their neighbourhoods contain reefs from multiple BCUs and thus if they were brought back in...they would reconnect some BCUs together
}
#put the a's and b's together, there shouldn't be any duplicates (bc then they should've been found in round 1) #there are not, good
pairedss_partiallostBCUs <- list() 
for(i in 1:length(partiallost_ss_BCUs)){
  if(finalpartiallostBCUss_candidates_a[[i]][1] == -1){next}
  pairedss_partiallostBCUs[[i]] <- data.frame(pairnum = seq(1,length(finalpartiallostBCUss_candidates_a[[i]])), a = NA, b = NA, same = NA)
  pairedss_partiallostBCUs[[i]]$a <- finalpartiallostBCUss_candidates_a[[i]]
  pairedss_partiallostBCUs[[i]]$b <- finalpartiallostBCUss_candidates_b[[i]]
  pairedss_partiallostBCUs[[i]]$same[which(pairedss_partiallostBCUs[[i]]$a == pairedss_partiallostBCUs[[i]]$b)] <- pairedss_partiallostBCUs[[i]]$a[which(pairedss_partiallostBCUs[[i]]$a == pairedss_partiallostBCUs[[i]]$b)]
}
#ultimately im just going to add all of these in at one time (even tho they're only effective as pairs), so just collate _a and _b into one set of values
pairedsteppingstones <- sort(unique(unlist(finalpartiallostBCUss_candidates_a)))
dummy <- sort(unique(unlist(finalpartiallostBCUss_candidates_b)))
#^ pairedsteppingstones and dummy are the same, which makes sense since all of the partiallost_ss_BCUs are contained in the lostconnection_ss_BCUs  
pairedsteppingstones <- unique(append(pairedsteppingstones,dummy))
pairedsteppingstones <- pairedsteppingstones[which(pairedsteppingstones > 0)]
#but, just to have it, these are the pairs...adding manually bc it's just easier
pairedsteppingstones_pairs <- matrix(data = c(4736,4772,11976,4756,11976,4772,11976,4773,5406,12006,5291,5827,5291,5898,5293,5827,12070,7099), nrow = 9, ncol = 2, byrow = T)
#but there are so many instances where one particular reef shows up in multiple pairs...so plotting this doesn't rly work

##Add those reefs back in and make a map showing where they are and then plot the new networks

#map of all of the stepping stone candidates after round 2, e.g. round 1 + round 2
#cols is very repetitive for this one, so going to use something else
rnbw <- rainbow(length(pairedsteppingstones))


#SUPPLEMENTARY MATERIALS 5 - FIGURE 1
png("CoralReefSteppingStones/SteppingStones_v2/FiftyR_SS2_Networks_fixedcolours.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=180, col="#e8e8e8",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
for(i in 1:fiftyr_sstwo_networks$no){
  points(reefdata$Longitude_corrected[reefdata$fiftyr_sstwo_ntwks == i], reefdata$Latitud[reefdata$fiftyr_sstwo_ntwks == i], pch = 20, cex = 0.2, col = cols[i])
}
points(reefdata$Longitude_corrected[reefdata$reefnum %in% steppingstonesafterroundtwo], reefdata$Latitud[reefdata$reefnum %in% steppingstonesafterroundtwo], pch = 8, col = "black")
dev.off()

#How this has changed the number of lost and partiallost BCUs, shouldn't have changed the former but should've changed the latter

#how many multiBCU networks are there now?

#step 1: go through the BCUs and see which networks their member reef cells are in
numBCUs <- range(reefdata$BCU_ID)[2]
BCU_fiftyr_sstwo_networks <- list()
allBCU_fiftyr_sstwo_ntwks <- unique(reefdata$fiftyr_sstwo_ntwks[reefdata$BCU_ID == 1])
for(i in 1:numBCUs){
  BCU_fiftyr_sstwo_networks[[i]] <- unique(reefdata$fiftyr_sstwo_ntwks[reefdata$BCU_ID == i])
  if(i > 1){
    allBCU_fiftyr_sstwo_ntwks <- append(allBCU_fiftyr_sstwo_ntwks, BCU_fiftyr_sstwo_networks[[i]]) 
  }
}
allBCU_fiftyr_sstwo_ntwks <- unique(allBCU_fiftyr_sstwo_ntwks)

#how many of the networks that have BCU reefs in them have BCU reefs from multiple networks
#if a network only has BCU reefs from one BCU, then the other reefs in it aren't going to be useful as stepping stones because they can't have connections to multiple BCUs
BCU_fiftyr_sstwo_ntwkoverlap <- data.frame(BCUntwk = allBCU_fiftyr_sstwo_ntwks, num_BCUs = 0)
for(i in 1:length(allBCU_fiftyr_sstwo_ntwks)){
  BCU_fiftyr_sstwo_ntwkoverlap$num_BCUs[BCU_fiftyr_sstwo_ntwkoverlap$BCUntwk == allBCU_fiftyr_sstwo_ntwks[i]] <- length(unique(reefdata$fiftyreef[reefdata$fiftyr_sstwo_ntwks == allBCU_fiftyr_sstwo_ntwks[i] & reefdata$BCU_ID > 0]))
}
morethanonefiftyrsstwoBCU <- BCU_fiftyr_sstwo_ntwkoverlap$BCUntwk[BCU_fiftyr_sstwo_ntwkoverlap$num_BCUs > 1]
#BCU_fiftyr_ss_ntwkoverlap[BCU_fiftyr_ss_ntwkoverlap$num_BCUs > 1,] 


#make a matrix of 1s and 0s, 1s saying that they're in the same network as each other (in the 50r conn mat)
fiftyr_50rsstwoconnmat_networkoverlap <- matrix(data = 0, nrow = numBCUs, ncol = numBCUs)

for(i in 1:numBCUs){
  ntwks <- BCU_fiftyr_sstwo_networks[[i]]
  bcus <- unique(reefdata$BCU_ID[reefdata$is_bcu == "refugia" & reefdata$fiftyr_sstwo_ntwks %in% ntwks]) #which BCUs also have reefs in those ntwks, will include numBCUs[i] but that's fine 
  fiftyr_50rsstwoconnmat_networkoverlap[i,bcus] <- fiftyr_50rsstwoconnmat_networkoverlap[bcus,i] <- 1 #may as well include both the upper and lower diagonal
  ntwks <- bcus <- NA
}

#then see how the two matrices compare
#then turn those 1s into -1s if those two BCUs are no longer in the same network and into 2s if they're still in the same network
#if any two 50Rs are in the same network now but weren't in the initial matrix, make those 1s (0s -> 1s)
fiftyr_ba_sstwo_networkoverlap <- matrix(data = 0, nrow = numBCUs, ncol = numBCUs)
BCU_sstwo_overlap_networks <- list()
BCU_sstwo_nonoverlap_networks <- list()

for(i in 1:numBCUs){
  for(j in 1:numBCUs){
    #still in the same network as each other 
    if(fiftyr_50rsstwoconnmat_networkoverlap[i,j] == 1 & fiftyr_fullconnmat_networkoverlap[i,j] == 1){
      fiftyr_ba_sstwo_networkoverlap[i,j] <- 2
    }
    #no longer in the same network as each other
    if(fiftyr_50rsstwoconnmat_networkoverlap[i,j] == 0 & fiftyr_fullconnmat_networkoverlap[i,j] == 1){
      fiftyr_ba_sstwo_networkoverlap[i,j] <- -1
    }
    #in the same network, but only when only the 50 reefs are the only reefs left (i don't think that this should be possible)
    if(fiftyr_50rsstwoconnmat_networkoverlap[i,j] == 1 & fiftyr_fullconnmat_networkoverlap[i,j] == 0){
      fiftyr_ba_sstwo_networkoverlap[i,j] <- 1
    }
  }
  BCU_sstwo_overlap_networks[[i]] <- which(fiftyr_ba_sstwo_networkoverlap[i,] == 2) #only need to do [i,] because a symmetric matrix
  BCU_sstwo_nonoverlap_networks[[i]] <- which(fiftyr_ba_sstwo_networkoverlap[i,] == -1)
}
#which(fiftyr_ba_sstwo_networkoverlap == 1) #none, good

#removing the self-loops
BCU_sstwo_overlap_networks_nosl <- list()
for(i in 1:numBCUs){
  BCU_sstwo_overlap_networks_nosl[[i]] <- BCU_sstwo_overlap_networks[[i]][!(BCU_sstwo_overlap_networks[[i]] == i)]
}

#which BCUs are no longer in ANY multiBCUnetworks? i.e. in nonoverlap_sstwo_BCUs but not in overlap_sstwo_BCUs
#some can be in both lists because they lost network connections with some BCUs but not all of them
overlap_sstwo_BCUs <- sort(unique(unlist(BCU_sstwo_overlap_networks_nosl))) #still in the same networks as some other BCUs
nonoverlap_sstwo_BCUs <- sort(unique(unlist(BCU_sstwo_nonoverlap_networks))) #no longer in the same networks as all of the other BCUs that they used to be in the same networks as
lost_sstwo_BCUs <- nonoverlap_sstwo_BCUs[which(!(nonoverlap_sstwo_BCUs %in% overlap_sstwo_BCUs))] #no longer in the same networks as any other BCUs (i.e. used to be in the same networks)
notlost_sstwo_BCUs <- overlap_sstwo_BCUs[which(!(overlap_sstwo_BCUs %in% nonoverlap_sstwo_BCUs))] #BCUs that didn't lose any of their other in-network-BCUs 
partiallost_sstwo_BCUs <- overlap_sstwo_BCUs[which(!(overlap_sstwo_BCUs %in% notlost_sstwo_BCUs))] #BCUs that lost some in-network-BCUs but not all of them
#lost_sstwo_BCUs = lost_ss_BCUs as expected, but also partiallost_sstwo_BCUs = partiallost_ss_BCUs so...that's not ideal

##some reefs still not connected, maybe need to go through and manually add 3 reefs from ntwk 166 and ntwk 53 and try to fill in the gaps
#only 192 non-BCU reef cells in ntwk 53, maybe even less if we remove the ones we've already designated as steppingstones

