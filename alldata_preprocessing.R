library(plyr)
data = read.csv("~/Desktop/Tesi/BiClustering/R-studio/Dati/export_milan_city_venues_stats.csv", sep=';')

# elimino colonne inutili
drops <- c("id", "name", "NIL_id", "category_id", "checkins")
data <- data[, !(names(data) %in% drops)]

# raggruppo per NIL_name e category_name
prova <- ddply(data, c("NIL_name","category_name"), colwise(sum))

# elimino le righe senza category_name
drops <- which(prova$category_name=='')
prova <- prova[-drops,]

# modifico i levels
prova$category_name <- factor(prova$category_name)

#creo la struttura 3d
nil_num <- length(unique(prova$NIL_name)) #88
cat_num <- length(unique(prova$category_name)) #590
time_num <- 20


alldata3d <- array(rep(NA,nil_num*cat_num*time_num), c(nil_num, cat_num, time_num),
                dimnames=list(sort(levels(unique(prova$NIL_name))),
                              sort(levels(unique(prova$category_name))),
                              names(prova)[3:length(names(prova))]))


bar<-0
pb <- txtProgressBar(style=3)

for (i in sort(levels(unique(prova$NIL_name)))){
  bar<-bar+1
  setTxtProgressBar(pb, value = bar/nil_num)
  for (j in sort(levels(unique(data$category_name)))){
    if (length(which(prova$NIL_name==i & prova$category_name==j))>0){
      for(k in 1:20)
      alldata3d[i,j,k] <- prova[which(prova$NIL_name==i & prova$category_name==j),k+2]
    }
    
  }
}


alldata3d

save(alldata3d, file = 'AllData3D.Rdata')
