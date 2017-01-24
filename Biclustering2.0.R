library(biclust)
library(xtable)
library(ggplot2)
library(ggdendro)
library(maptools)
library(ggmap)
library(RColorBrewer)
#### NEWDATA3d ####
setwd('~/Desktop/Tesi/BiClustering/R-studio/Dati/')

# carico e pulisco i dati originali [247 categorie]
load(file='../Dati/Data3D.Rdata')
data3d[which(is.na(data3d))] <- 0

# carico e pulisco i dati new [45 categorie invece che 247]
load(file='../Dati/NewData3D.Rdata')
new_data3d[which(is.na(new_data3d))] <- 0

# carico i dati originali per plottare la mappa
data = read.csv("~/Desktop/Tesi/BiClustering/R-studio/Dati/4sq.csv")

#### DATASUM ####
# Matrice somma dei 15 mesi
data_sum <- data3d[,,1]
for (i in 2:15)
  data_sum <- data_sum + data3d[,,i]


#### NEWDATASUM ####
# Matrice somma dei 15 mesi
new_data_sum <- new_data3d[,,1]
for (i in 2:15)
  new_data_sum <- new_data_sum + new_data3d[,,i]

#### PREPARAZIONE MAPPA ####
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

#### DATA CHOICE ####

# prepara i dati sommati con tutte le categorie
bi_data <- data_sum

# prepara i dati sommati con 45 categorie
bi_data <- new_data_sum

#### BICLUSTERING TEXT ####

# NOTA BENE: bisogna cambiare il Metodo manualmente
# tieni presente che sui SUM BCPLAID non da risultati

res <- biclust(bi_data, method=BCCC())
#summary(res)
#names(attributes(res))
nomefile <- paste0("BCCC_Sum",".html")
sink(nomefile)
for (k in 1:res@Number){
  print(paste("cluster numero",k))
  print(xtable(bi_data[res@RowxNumber[,k], res@NumberxCol[k,]]),type='html')
  
}
sink()


#### BICLUSTERING MAP ####

res <- biclust(bi_data, method=BCCC())

colori <- array(NA, dim = res@Number)
numclust <- res@Number

for (i in 1:numclust){
  idx_i <- unique(data$nil_id[which(data$nil_name %in% rownames(bi_data[res@RowxNumber[,i],]))])
  colori[which(area.points$id %in% idx_i)] <- i
}

filename <- paste0('BCC_SuMMap_','.png')
mypath <- file.path("~/Desktop/SUM_BICLUST",filename)

map <- ggmap(mapImage) +
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


#### BICLUSTERING BARPLOT ####

res <- biclust(bi_data, method=BCCC())

colori <- array(NA, dim = res@Number)
numclust <- res@Number

# plottiamo ora il barplot
filename <- paste0('BCC_SUMBar_', '.png')
mypath <- file.path("~/Desktop/SUM_BICLUST",filename)
png(file= mypath, width=800, height=1000)
biclustbarchart(bi_data,res)
dev.off()
