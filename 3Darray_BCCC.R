# World of biclustering --------
library(biclust)

# matrice test --------------------------------------------------------------------------------------
test1 <- matrix(rbinom(400, 50, 0.4), 120, 120)
test2 <- matrix(rbinom(400, 50, 0.4), 120, 120)

# 3d array test
lista3d <- list(test1, test2)
library(abind)
mat <- do.call(abind, c(lista3d, along = 3))

# FourSquare: Urbanscope ----------------------------------------------------------------------------
load("/Users/Jacopo/Desktop/Tesi/BiClustering/R-studio/Dati/AllNewData.Rdata")
dim(allnewdata)
mat <- allnewdata

# Cheng & Church for 3d Array

# Algoritmo per trovare bicluster -------------------------------------------------------------------

# Find biggest 3D Bicluster for 3D Array

bigcc<-function(mat,delta,alpha=1.5)
{
  # crea due vettori della dimensione delle righe (r) e delle colonne (c) con TRUE
  logr<-rep(TRUE,nrow(mat))
  logc<-rep(TRUE,ncol(mat))
  
  # Multiple Node deletion
  step1<-cc2(mat,logr,logc,delta,alpha)
  # Single Node deletion
  step2<-cc1(mat,step1[[1]],step1[[2]],delta)
  # controllo se ne trova qualcuna
  if(sum(step2[[1]])==0)
  {ret<-list(0,warning(paste('No matrix with score smaller than', delta,'found')))
  }
  else{
    # se ne ha trovate: Node Addition
    ret<-cc3(mat,step2[[1]],step2[[2]])
  }
  ret
}

# Trova numero di bicluster ------------------------------------------------------------------------

ccbiclust<-function(mat,delta,alpha=1.5,number=100, rand=FALSE)
{
  MYCALL <- match.call()
  ma<-max(mat)
  mi<-min(mat)
  x<-matrix(FALSE,nrow=nrow(mat),ncol=number)
  y<-matrix(FALSE,nrow=number,ncol=ncol(mat))
  logr<-rep(TRUE,nrow(mat))
  Stop <- FALSE
  logr<-rep(TRUE,nrow(mat))
  for(i in 1:number)
  {
    if(sum(logr)<2)
    {
      Stop <- TRUE
      break
    }
    erg<-bigcc(mat[logr,,],delta,alpha)
    if(sum(erg[[1]])==0)
    {
      Stop <- TRUE
      break
    }
    else
    {
      x[logr,i]<-erg[[1]]
      y[i,]<-erg[[2]]
      if(rand)
      {
        mat[erg[[1]],erg[[2]]]<-runif(sum(erg[[1]])*sum(erg[[2]]),mi,ma)
      }
      else
      {
        logr[logr][erg[[1]]] <- FALSE
      }
    }
  }
  
  if(Stop)
  {
    return(BiclustResult(as.list(MYCALL),as.matrix(x[,1:(i-1)]),as.matrix(y[1:(i-1),]),(i-1),list(0)))
  }
  else
  {
    return(BiclustResult(as.list(MYCALL),as.matrix(x),as.matrix(y),i,list(0)))
  }
  
  
}

BiclustResult <- function(mypara, a, b, c, d) {
  return(new('Biclust', Parameters=mypara, RowxNumber=a, NumberxCol=b, Number=c, info=d))
}


# Multiple Node deletion ---------------------------------------------------------------------------

# CRITICITA' scelta del delta
delta <-1.0

cc2<-function(mat,logr,logc,delta,alpha=1.5)
{
  # imposto a 1 il numero di multiple row e col da eliminare per entrare nel ciclo
  mdi<-1
  mdj<-1
  
  # Finché H, residual score di ccscore per la sottomatrice selezionata [logr,logc] è maggiore del 
  # delta ricercato e finché faccio delle eliminazioni multiple allora non ho il bicluster che cerco
  while((h<-ccscore(mat[logr,logc,]))>delta & (sum(mdi)+sum(mdj))>0)
  {
    
    # se le righe sono più di 100 procedo con la multiple node deletion
    if(sum(logr)>100)
    {
      # rowscore calcola lo score riga per tutte le righe
      di<-rowscore(mat[logr,logc,])
      
      # trovo tutte le righe con di maggiore di alpha volte lo score attuale della matrice
      # restituisce un vettore con tanti elementi TRUE/FALSE quanti le righe
      mdi<-di>(alpha*h)
      
      # se non sto cancellando tutte le righe tranne una in un colpo -> Multiple Node Deletion
      if(sum(mdi) < (sum(logr)-1))
      {
        # tolgo la riga e aggiorno logr
        logr[logr][mdi]<-FALSE
        #calcolo il nuovo H
        h<-ccscore(mat[logr,logc,])
      }
      else
      {
        # altrimenti... sto cancellendo tutto perché ho preso un alpha troppo piccolo
        print(warning(paste('Alpha', alpha,'to small!')))
        # non faccio la multiple node deletion 
        mdi <- 0
      }
    }
    # se le righe invece sono poche non posso fare una multiple deletion
    else{mdi<-0}
    
    # passo alle colonne con lo score H eventualmente aggiornato dalla multiple row deletion
    if(sum(logc)>100)
    {
      dj<-colscore(mat[logr,logc,])
      mdj<-dj>(alpha*h)
      if(sum(mdj) < (sum(logc)-1))
      {
        logc[logc][mdj]<-FALSE
      }
      else
      {
        
        print(warning(paste('Alpha', alpha,'to small!')))
        mdi <- 0
      }
    }
    else{mdj<-0}
  }
  
  # esco dal WHILE e restituisco lo storico delle eliminazioni per righe e colonne come lista 
  ret<-list(logr,logc)
  ret
}

# Single Node Deletion for 3D array ---------------------------------------------------------------

cc1<-function(mat,logr,logc,delta=1.5)
{
  # finché H, residual score di ccscore, è maggiore del delta richiesto allora non il bicluster
  while(ccscore(mat[logr,logc,])>delta)
  {
    # calcolo lo score riga per ogni riga
    di<-rowscore(mat[logr,logc,])
    # calcolo lo score colonna per ogni colonna
    dj<-colscore(mat[logr,logc,])
    # identifico la riga e la colonna con gli scores maggiori
    mdi<-which.max(di)
    mdj<-which.max(dj)
    
    # se score riga più grande supera quello colonna allora elimino la riga, altrimenti viceversa
    ifelse(di[mdi]>dj[mdj] ,logr[logr][mdi]<-FALSE ,logc[logc][mdj]<-FALSE)
    
    # se non c'è più una colonna e una riga disponibile esco dal WHILE
    if (!(sum(logr)>1 & sum(logc)>1))
      break
  }
  
  # uscito da WHILE restituisco se possibile la lista delle deletion altrimenti il warning
  ifelse(sum(logr)>1 & sum(logc)>1,ret<-list(logr,logc),ret<-list(0,warning(paste('No matirx with score smaller', delta,'found'))))
  ret
}

# Node Addition for 3D Array-----------------------------------------------------------------------------------

cc3<-function(mat,logr,logc)
{
  # impongo br<-1 per fare la prima iterazione
  br<-1
  # vettore con tanti elementi quante righe tutti FALSE
  ilogr<-rep(FALSE,length(logr))
  
  # finché si può aggiungere (br>0)
  while(br>0)
  {
    # impongo br1 e br2 come il numero di righe non eliminate
    br1<-sum(logc)
    br2<-sum(logr)
    
    # calcolo H
    h<-ccscore(mat[logr,logc,])
    
    # calcolo lo score per ogni colonna non presente
    dj<-addcolscore(mat,logr,logc)
    
    # trovo tutte le colonne che si possono aggiungere
    mdj<-dj<=h
    #le aggiungo
    logc[mdj]<-TRUE
    
    #calcolo H nuovo
    h<-ccscore(mat[logr,logc,])
    # calcolo lo scoro per ogni riga non presente
    di<-addrowscore(mat,logr,logc)
    # calcolo lo score INVERSO per ogni riga non presente
    idi<-iaddrowscore(mat,logr,logc)
    
    # multiple row addition
    mdi<-di<=h
    logr[mdi]<-TRUE
    
    # multiple inverse row addition
    imdi<-idi<=h
    # CONTROLLA MEGLIO!!!
    mat[!(logr==imdi)&imdi]<- -mat[!(logr==imdi)&imdi]
    logr[imdi]<-TRUE
    
    br<-sum(logc)+sum(logr)-br1-br2
  }
  ret<-list(logr,logc)
  ret
  
}

# Alcuni Scores 3D --------------------------------------------------------------------------------

# Mean squared residual score per array 3d di qualsiasi dimensione

ccscore<-function(mat)
{
  scorevector <- c()
  for( i in 1:dim(mat)[3]){
    partial_score <- sum((mat[,,i]-rowMeans(mat[,,i])-matrix(colMeans(mat[,,i]),nrow=nrow(mat[,,i]),ncol=ncol(mat[,,i]),byrow=TRUE)+mean(mat[,,i]))^2)/(nrow(mat[,,i])*ncol(mat[,,i]))
    scorevector <- c(scorevector, partial_score)
  }
  score <- sqrt(sum(scorevector^2))
  return(score)
}

# d(i) usato nella node deletion
rowscore<-function(mat)
{
  scorevector <- c()
  for( i in 1:dim(mat)[3]){
    partial_score<-rowSums((mat[,,i]-rowMeans(mat[,,i])-matrix(colMeans(mat[,,i]),nrow=nrow(mat[,,i]),ncol=ncol(mat[,,i]),byrow=TRUE)+mean(mat[,,i]))^2)/ncol(mat[,,i])
    scorevector <- rbind(scorevector, partial_score)
  }
  score <- sqrt(colSums(scorevector^2))
  return(score)
}

# d(j) usato nella node deletion
colscore<-function(mat,logr,logc)
{
  scorevector <- c()
  for( i in 1:dim(mat)[3]){
    partial_score <- colSums((mat[,,i]-rowMeans(mat[,,i])-matrix(colMeans(mat[,,i]),nrow=nrow(mat[,,i]),ncol=ncol(mat[,,i]),byrow=TRUE)+mean(mat[,,i]))^2)/nrow(mat[,,i])
    scorevector <- rbind(scorevector, partial_score)
  }
  score <- sqrt(colSums(scorevector^2))
  return(score)
}

# Scores per Node Addition: identici agli altri tranne che considera le righe o le colonne eliminate

addrowscore<-function(mat,logr,logc)
{
  scorevector <- c()
  for( i in 1:dim(mat)[3]){
    partial_score<-rowSums((mat[,,i]-rowMeans(mat[,logc,i])-matrix(colMeans(mat[logr,,i]),nrow=nrow(mat[,,i]),ncol=ncol(mat[,,i]),byrow=TRUE)+mean(mat[logr,logc,i]))^2)/ncol(mat[logr,logc,i])
    scorevector <- rbind(scorevector, partial_score)
  }
  score <- sqrt(colSums(scorevector^2))
  return(score)
}

iaddrowscore <- function(mat,logr,logc)
{
  scorevector <- c()
  for (i in 1:dim(mat)[3]){
    partial_score <- rowSums((-mat[,,i]+rowMeans(mat[,logc,i])-matrix(colMeans(mat[logr,,i]),nrow=nrow(mat[,,i]),ncol=ncol(mat[,,i]),byrow=TRUE)+mean(mat[logr,logc,i]))^2)/ncol(mat[logr,logc,i])
    scorevector <- rbind(scorevector, partial_score)
  }
  score <- sqrt(colSums(scorevector^2))
  return(score)
}

addcolscore <- function(mat,logr,logc)
{
  scorevector <- c()
  for (i in 1:dim(mat)[3]){
    partial_score<-colSums((mat[,,i]-rowMeans(mat[,logc,i])-matrix(colMeans(mat[logr,,i]),nrow=nrow(mat[,,i]),ncol=ncol(mat[,,i]),byrow=TRUE)+mean(mat[logr,logc,i]))^2)/nrow(mat[logr,logc,i])
    scorevector <- rbind(scorevector, partial_score)
  }
  score <- sqrt(colSums(scorevector^2))
  return(score)
}


# IDEE PER IL FUNZIONALE

# Creare biclustering CC vettoriale che funzioni con i 3d array  (FATTO!!!)
#   - scelta da norma euclidea per ricondurre tutti i layers ad un unico score
#   - capire l'utilità
#   - capire se con lat e lon ci sono problemi perché sono diveri (UTC?)


# Lanciare e visualizzare -------------------------------------------------------------------------

res <- ccbiclust(mat,delta = 1,alpha = 1.5,rand = FALSE)
mat[res@RowxNumber[,8], res@NumberxCol[8,],]

# 3d Array Telecom
load("~/Desktop/Telecom/Traiettorie/user_tracking_export/newsolvedsecondsTraj_unocento.Rda")

library(rgdal)
cord.dec <- SpatialPoints(cbind(as.numeric(secondsTraj$lon_most), as.numeric(secondsTraj$lon_most)), proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32632"))
secondsTraj <- cbind(secondsTraj, as.data.frame(cord.UTM))

# rifare la procedura per la matrice in cord.UTM dato che hanno la stessa scala fra di loro

