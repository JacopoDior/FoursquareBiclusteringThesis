library(cluster)
library(ggplot2)
library(ggdendro)
library(maptools)
library(ggmap)

####################
# Prima Iterazione #
####################


################
# Preparazione #
################

load(file='../Dati/Data3D.Rdata')

# Trasformo gli NA in 0
data3d[which(is.na(data3d))] <- 0

# Carico i dati
data = read.csv("~/Desktop/Tesi/BiClustering/R-studio/Dati/4sq.csv")

setwd('~/Desktop/Tesi/BiClustering/R-studio/Script')

# carico shp
area <- readShapePoly("../Dati/NIL/NILZone_wgs84.shp")

# scarica la mappa e la salva in mapImage
mapImage <- get_map(location = c(lon = 9.18592, lat = 45.46542),
                    color = "color",
                    source = "osm",
                    # maptype = "terrain",
                    zoom = 11)
# prepara area in formato leggibile
area.points <- fortify(area)

# Riordinamento dei NIL causa file corrotto
map_index <- cbind(0:87, area[[1]])

new_id <- array(NA, length(area.points$id))
for (i in 0:dim(map_index)[1]-1){
  new_id[which(area.points$id == i)] <- map_index[i+1,2]
}
area.points$id <- as.character(new_id)


############################
# Hclust su Nil x Category #
############################

#clust_eucl<-clara(log(data3d[,,15]+0.001), 3, metric="euclidean")

# Distanza euclide e metodo di agglomerazione ward.D
for (w in (1:15))
{
diss <- dist(data3d[,,w], method='euclidean') # distanza euclidea
clust_hier <- hclust(diss, method = "ward.D") 
ggdendrogram(clust_hier,segments=TRUE,leaf_labels=TRUE) 

# Taglio in cluster
cluster.ec <- cutree(clust_hier, k=6)
#which(cluster.ec==2)

############
# Plotting #
############

colori <- array(NA, dim = length(area.points$group))
numclust <- length(unique(cluster.ec))

for (i in 1:numclust){
  idx_i <- unique(data$nil_id[which(data$nil_name %in% names(which(cluster.ec==i)))])
  colori[which(area.points$id %in% idx_i)] <- i
}

filename <- paste0('ClusterNIL_',w,'.png')
mypath <- file.path("~/Desktop/CANE",filename)

ggmap(mapImage) +
  geom_polygon(aes(x = long,
                   y = lat,
                   group = group),
               data = area.points,
               color = 'black',
               fill = colori,
               alpha = 0.5) +
  labs(x = "Longitudine",
       y = "Latitudine")

ggsave(file = mypath, width=15.875, height=15.875)
}

############################
# Kmeans su Nil x Category #
############################

K <- kmeans(data3d[,,6], 4)


######################################
# Hclust su Nil x [Category + Month] #
######################################

###########################
# Month VS Nil X Category # corretto
###########################

# Ottengo dimensione dati
dataDm <- dim(data3d)

# DataFrame
dataRGB <- data.frame(
  x = rep(colnames(data3d),each=dataDm[1]), # crea un array di 274*87 dove un contatore scatta ogni 87
  y = rep(sort(rownames(data3d)),dataDm[2]),  # array di 274*87  
  lug2014 = as.vector(data3d[,,1]), 
  ago2014 = as.vector(data3d[,,2]), 
  set2014 = as.vector(data3d[,,3]), 
  ott2014 = as.vector(data3d[,,4]), 
  nov2014 = as.vector(data3d[,,5]), 
  dic2014 = as.vector(data3d[,,6]), 
  gen2015 = as.vector(data3d[,,7]),
  feb2015 = as.vector(data3d[,,8]), 
  mar2015 = as.vector(data3d[,,9]), 
  apr2015 = as.vector(data3d[,,10]), 
  mag2015 = as.vector(data3d[,,11]), 
  giu2015 = as.vector(data3d[,,12]), 
  lug2015 = as.vector(data3d[,,13]),
  ago2015 = as.vector(data3d[,,14]), 
  set2015 = as.vector(data3d[,,15]) 
)

# Mergo le prime due colonne e aggiungo il risultato al dataset
dataRGB <-within(dataRGB, Nil_Category <- paste(y, x, sep='_'))
dataRGB <- dataRGB[,3:18]
dataRGB <- cbind(dataRGB$Nil_Category,dataRGB)[-(dim(dataRGB)[2]+1)]
colnames(dataRGB)[1]<- "Nil X Category"

rownames(dataRGB) <- dataRGB[,1]
dataRGB[,1] <- NULL

# Hclust

diss <- dist(t(dataRGB), method='euclidean') # distanza euclidea
clust_hier <- hclust(diss, method = "ward.D") 
ggdendrogram(clust_hier,segments=TRUE,leaf_labels=TRUE) 

###########################
# Nil VS Category X Month # corretto
###########################

# Ottengo dimensione dati
dataDm <- dim(data3d) # 87 274 15

dataRGB <- data.frame(
  x = rep(colnames(data3d),each=dataDm[3]), # crea un array di categorie di 274*15 dove un contatore scatta ogni 15
  y = rep(sort(dimnames(data3d)[[3]]),dataDm[2]) #
  )
  
for (i in 1:dataDm[1]){
  dataRGB <- cbind(dataRGB,as.vector(t(data3d[i,,])))
}
 
dataRGB <- cbind(dataRGB, Category_name <- paste(dataRGB$x, dataRGB$y, sep='_'))
dataRGB <- dataRGB[,3:dim(dataRGB)[2]]
dataRGB <- cbind(dataRGB[dim(dataRGB)[2]],dataRGB)[-(dim(dataRGB)[2]+1)]
colnames(dataRGB)<- c("Nil X Category",rownames(data3d))

rownames(dataRGB) <- dataRGB[,1]
dataRGB[,1] <- NULL

# hclust Accessories Vs Nil X Month
diss <- dist(t(dataRGB), method='euclidean') # distanza euclidea
clust_hier <- hclust(diss, method = "ward.D") 
ggdendrogram(clust_hier,segments=TRUE,leaf_labels=TRUE) 


###########################
# Category VS Nil X Month #  corretto
###########################

# Ottengo dimensione dati
dataDm <- dim(data3d) # 87 274 15

dataRGB <- data.frame(
  x = rep(rownames(data3d),each=dataDm[3]), # crea un array di categorie di 87*15 dove un contatore scatta ogni 15
  y = rep(sort(dimnames(data3d)[[3]]),dataDm[1]) #
)

for (i in 1:dataDm[2]){
  dataRGB <- cbind(dataRGB,as.vector(t(data3d[,i,]))) # devo traslare i dati per non commettere errori
}

dataRGB <- cbind(dataRGB,  paste(dataRGB$x, dataRGB$y, sep='_'))
dataRGB <- dataRGB[,3:dim(dataRGB)[2]]
dataRGB <- cbind(dataRGB[dim(dataRGB)[2]],dataRGB)[-(dim(dataRGB)[2]+1)]
colnames(dataRGB)<- c("Nil X Month",colnames(data3d))


rownames(dataRGB) <- dataRGB[,1]
dataRGB[,1] <- NULL

# hclust Accessories Vs Nil X Month
diss <- dist(t(dataRGB), method='euclidean') # distanza euclidea
clust_hier <- hclust(diss, method = "ward.D") 
ggdendrogram(clust_hier,segments=TRUE,leaf_labels=TRUE) 

######################
##### SECONDA ITERAZIONE #####
######################

load(file='../Dati/NewData3D.Rdata')

# Trasformo gli NA in 0
new_data3d[which(is.na(new_data3d))] <- 0

for (w in 1:15){
# Distanza euclide e metodo di agglomerazione ward.D
  
filename <- paste0('45NILdendro',w,'.png')
mypathdendro <- file.path("~/Desktop/CANE",filename)
  
diss <- dist(new_data3d[,,w], method='euclidean') # distanza euclidea
clust_hier <- hclust(diss, method = "ward.D") 
ggdendrogram(clust_hier,segments=TRUE,leaf_labels=TRUE) 
ggsave(mypathdendro, width=15, height=8)

# Taglio in cluster
cluster.ec <- cutree(clust_hier, k=8)
#which(cluster.ec==2)

############
# Plotting #
############

colori <- array(NA, dim = length(area.points$group))
numclust <- length(unique(cluster.ec))

for (i in 1:numclust){
  idx_i <- unique(data$nil_id[which(data$nil_name %in% names(which(cluster.ec==i)))])
  colori[which(area.points$id %in% idx_i)] <- i
}


filename <- paste0('NILsummed','.png')
mypathnew <- file.path("~/Desktop/CANE",filename)


p <- ggmap(mapImage) +
  geom_polygon(aes(x = long,
                   y = lat,
                   group = group),
               data = area.points,
               color = 'black',
               fill = colori,
               alpha = 0.5) +
  labs(x = "Longitudine",
       y = "Latitudine")

ggsave(file = mypathnew, plot=p, dpi=120)

}



########################
# APPROCCIO VETTORIALE #
########################

load(file='../Dati/NewData3D.Rdata')

# Trasformo gli NA in 0
new_data3d[which(is.na(new_data3d))] <- 0

# Ottengo dimensione dati
dataDm <- dim(new_data3d) # 87 274 15

dataRGB <- data.frame(
  x = rep(colnames(new_data3d),each=dataDm[3]), # crea un array di categorie di 274*15 dove un contatore scatta ogni 15
  y = rep(sort(dimnames(new_data3d)[[3]]),dataDm[2]) #
)

for (i in 1:dataDm[1]){
  dataRGB <- cbind(dataRGB,as.vector(t(new_data3d[i,,])))
}

dataRGB <- cbind(dataRGB, Category_name <- paste(dataRGB$x, dataRGB$y, sep='_'))
dataRGB <- dataRGB[,3:dim(dataRGB)[2]]
dataRGB <- cbind(dataRGB[dim(dataRGB)[2]],dataRGB)[-(dim(dataRGB)[2]+1)]
colnames(dataRGB)<- c("Nil X Category",rownames(new_data3d))

rownames(dataRGB) <- dataRGB[,1]
dataRGB[,1] <- NULL

filename <- paste0('NILvectorialdendro','.png')
mypath <- file.path("~/Desktop/CANE",filename)

# hclust Accessories Vs Nil X Month
diss <- dist(t(dataRGB), method='euclidean') # distanza euclidea
clust_hier <- hclust(diss, method = "ward.D") 
ggdendrogram(clust_hier,segments=TRUE,leaf_labels=TRUE) 

ggsave(mypath, width=15, height=8)

#### SUM ####

sumdata <- new_data3d[,,1]

for (i in 2:15)
  sumdata <- sumdata + new_data3d[,,i]

filename <- paste0('NILsumDendro','.png')
mypath <- file.path("~/Desktop/CANE",filename)

diss <- dist((sumdata), method='euclidean') # distanza euclidea
clust_hier <- hclust(diss, method = "ward.D") 
ggdendrogram(clust_hier,segments=TRUE,leaf_labels=TRUE) 

ggsave(mypath, width=15, height=8)
