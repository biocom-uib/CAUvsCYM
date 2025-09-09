#' Funció EasyCoda::WARD 
#' 1. per aplicar directament al resultat de EasyCoda::CLR
#' 2. eliminats tots els warnings,
#' fa més via a les simulacions

WARD.noprints <- function(LRdata, weight = TRUE) {
  LRfoo <- LRdata$LR
  weights <- LRdata$LR.wt
  row.wt <- rep(1/nrow(LRfoo), nrow(LRfoo))        
  ### make equivalent to hierclust
  
  a <- LRfoo
  wt <- row.wt
  
  n <- nrow(a)    
  m <- ncol(a)                           
  diss <- matrix(0, n, n)
  rownames(diss) <- colnames(diss) <- rownames(a)
  
  for (i1 in 2:n) {
    for (i2 in 1:(i1-1)) {
      temp <- 0.0
      for (j in 1:m) {
        # We use the squared Euclidean distance, weighted row & column wise
        temp <- temp + (wt[i1]*wt[i2])/(wt[i1]+wt[i2]) *
          weights[j] * (a[i1,j]-a[i2,j])^2 
      }
      diss[i2,i1] <- diss[i1,i2] <- temp
    }
  }
  
  flag <- rep(1, n)                          # active/dead indicator
  a <- rep(0, n-1)                           # left subnode on clustering
  b <- rep(0, n-1)                           # right subnode on clustering
  ia <- rep(0, n-1)                          # R-compatible version of a
  ib <- rep(0, n-1)                          # R-compatible version of b
  lev <- rep(0, n-1)                         # level or criterion values
  card <- rep(1, n)                          # cardinalities
  mass <- wt
  order <- rep(0, n)                         # R-compatible order for plotting
  

  nn <- rep(0, nrow(diss))
  nndiss <- rep(0.0, nrow(diss))
  MAXVAL <- 1.0e12
  for (i1 in 1:nrow(diss)) {
    minobs <- -1
    mindis <- MAXVAL
    for (i2 in 1:ncol(diss)) {
      if ( (diss[i1,i2] < mindis) && (i1 != i2) ) {
        mindis <- diss[i1,i2]
        minobs <- i2
      }
    }
    nn[i1] <- minobs
    nndiss[i1] <- mindis
  }
  
### compatibility with hierclust after "function" call
  
  nnsnnsdiss <- list(nn=nn, nndiss=nndiss)
  
  
  
  clusmat <- matrix(0, n, n)                 # cluster memberships
  for (i in 1:n) clusmat[i,n] <- i           # init. trivial partition
  
  MAXVAL <- 1.0e12
  
  for (ncl in (n-1):1) {                      # main loop 
    # check for agglomerable pair
    minobs <- -1;  
    mindis <- MAXVAL;
    for (i in 1:n) {
      if (flag[i] == 1) {
        if (nnsnnsdiss$nndiss[i] < mindis) {
          mindis <- nnsnnsdiss$nndiss[i]
          minobs <- i
        }
      }
    }
    # find agglomerands clus1 and clus2, with former < latter
    if (minobs < nnsnnsdiss$nn[minobs]) {
      clus1 <- minobs
      clus2 <- nnsnnsdiss$nn[minobs]
    }
    if (minobs > nnsnnsdiss$nn[minobs]) {
      clus2 <- minobs
      clus1 <- nnsnnsdiss$nn[minobs]
    }
    # So, agglomeration of pair clus1 < clus2 defines cluster ncl
    
    ### Block for subnode labels 
    a[ncl] <- clus1                       # aine, or left child node
    b[ncl] <- clus2                       # benjamin, or right child node
    # Now build up ia, ib as version of a, b which is R-compliant
    if (card[clus1] == 1) ia[ncl] <- (-clus1)     # singleton
    if (card[clus2] == 1) ib[ncl] <- (-clus2)     # singleton
    if (card[clus1] > 1) {                # left child is non-singleton
      lastind <- 0
      for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
        if (a[i2] == clus1) lastind <- i2    # Only concerns a[i2]
      }
      ia[ncl] <- n - lastind             # label of non-singleton
    }
    if (card[clus2] > 1) {                # right child is non-singleton
      lastind <- 0
      for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
        if (a[i2] == clus2) lastind <- i2    # Can only concern a[i2]
      }
      ib[ncl] <- n - lastind             # label of non-singleton
    }
    if (ia[ncl] > 0 || ib[ncl] > 0) {     # Check that left < right
      left <- min(ia[ncl],ib[ncl])
      right <- max(ia[ncl],ib[ncl])
      ia[ncl] <- left                    # Just get left < right
      ib[ncl] <- right
    }
    #### End of Block for subnode labels
    
    lev[ncl] <- mindis
    for (i in 1:n) {
      clusmat[i,ncl] <- clusmat[i,ncl+1]
      if (clusmat[i,ncl] == clus2) clusmat[i,ncl] <- clus1
    }
    # Next we need to update diss array
    for (i in 1:n) {
      if ( (i != clus1) && (i != clus2) && (flag[i] == 1) ) {
        diss[clus1,i] <- 
          ((mass[clus1]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,i] +
          ((mass[clus2]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus2,i] -
          (mass[i]/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,clus2] 
        diss[i,clus1] <- diss[clus1,i]
      }
    }
    mass[clus1] <- mass[clus1] + mass[clus2]    # Update mass of new cluster
    card[clus1] <- card[clus1] + card[clus2]    # Update card of new cluster
    # Cluster label clus2 is knocked out; following not nec. but no harm
    flag[clus2] <- 0
    nnsnnsdiss$nndiss[clus2] <- MAXVAL
    mass[clus2] <- 0.0
    for (i in 1:n) {
      diss[clus2,i] <- MAXVAL
      diss[i,clus2] <- diss[clus2,i]
    }
    # Finally update nnsnnsdiss$nn and nnsnnsdiss$nndiss
    # i.e. nearest neighbors and the nearest neigh. dissimilarity
    #       nnsnnsdiss <- getnns(diss, flag)
    
    ### replaced with code to find smallest distance, no error checking needed
    ### flag all 1s, no need to check
    
    nn <- rep(0, nrow(diss))
    nndiss <- rep(0.0, nrow(diss))
    MAXVAL <- 1.0e12
    for (i1 in 1:nrow(diss)) {
      minobs <- -1
      mindis <- MAXVAL
      for (i2 in 1:ncol(diss)) {
        if ( (diss[i1,i2] < mindis) && (i1 != i2) ) {
          mindis <- diss[i1,i2]
          minobs <- i2
        }
      }
      nn[i1] <- minobs
      nndiss[i1] <- mindis
    }
    
    ### compatibility with hierclust after "function" call
    
    nnsnnsdiss <- list(nn=nn, nndiss=nndiss)
    
  }
  
  temp <- cbind(a,b)
  merge2 <- temp[nrow(temp):1, ]
  temp <- cbind(ia,ib)
  merge <- temp[nrow(temp):1,]
  dimnames(merge) <- NULL
  # merge is R-compliant; later suppress merge2
  
  #### Build R-compatible order from ia, ib
  orderlist <- c(merge[n-1,1], merge[n-1,2])
  norderlist <- 2
  for (i in 1:(n-2)) {           # For precisely n-2 further node expansions
    for (i2 in 1:norderlist) {       # Scan orderlist
      if (orderlist[i2] > 0) {     # Non-singleton to be expanded
        tobeexp <- orderlist[i2]
        if (i2 == 1) {
          orderlist <- c(merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[2:norderlist])
        }
        if (i2 == norderlist) {
          orderlist <- c(orderlist[1:(norderlist-1)],
                         merge[tobeexp,1],merge[tobeexp,2])
        }
        if (i2 > 1 && i2 < norderlist) {
          orderlist <- c(orderlist[1:(i2-1)], 
                         merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[(i2+1):norderlist])
        }
        norderlist <- length(orderlist)
      }
    }
  }
  orderlist <- (-orderlist)
  class(orderlist) <- "integer"
  
  xcall <- "hierclust(a,wt)"
  class(xcall) <- "call"
  
  #  output square roots of heights, so that sum of squares gives TotVar
  retlist <- list(merge=merge,height=sqrt(lev[(n-1):1]),order=orderlist)
  class(retlist) <- "hclust"
  retlist
}

##---------------------------------------------##
##---------------------------------------------##

#' Encerts: Compara la classificació obtinguda amb la real

Encerts=function(A,tt){
  if(diff(dim(A))!=0){stop("A no és quadrada")}
  n=dim(A)[1]
  if (n==1){optim=0} 
  else {
    perms=gtools::permutations(n,n)
    valors=c(sum(diag(A)),sum((A-diag(tt))^2))
    for (j in 2:dim(perms)[1]){
      x=perms[j,]  
      valors=rbind(valors,c(sum(diag(A[,x])),sum((A[,x]-diag(tt))^2)))
      optim=min(valors[valors[1,]==max(valors[1,]),2])
    }
  }
  return(optim)
}

##---------------------------------------------##
##---------------------------------------------##

#' Potpourri de funcions que calcula clustering jeràrquic a partir de
#' distància d'Aitchinson, el dibuixa (o no), i 
#' dibuixa (o no) barplot de composició de taxa
#' Retorna el clustering, el dendrograma, i la taula de la classificació 
#' de les classes de Grups en les classes del  clustering 
#' 
ClusterHC=function(X,dendrograma=TRUE,barplot=FALSE,Grups,colores=colors){
  M=length(levels(Grups))
  hc=WARD.noprints(easyCODA::CLR(X),weight=TRUE)
  dend=as.dendrogram(hc)
  clustHC=data.frame(Orig=Grups[hc$order],
                     clust=cutree(dend, k = M)[order.dendrogram(dend)])
  labels_colors(dend)=colores[hc$order]
  labels(dend)=row.names(X)
  #'
  if(dendrograma==TRUE){
      XOr=X[hc$order,]
    XOr.CLR=acomp(XOr)
    if (barplot==TRUE){
      layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(6,2), height=c(4,4))
      par(mar=c(2,1,1,1)+0.1,cex=0.75)
      plot(dend, main = "")
      par(cex=1)}
    if (barplot==FALSE){plot(dend, main = "")}
  }  
  #'
  if(barplot==TRUE){ 
    nb.cols=dim(X)[2]
    d.names=colnames(X)[order(apply(X, 2, sum), decreasing=T) ]
    colors.OTU=colorRampPalette(brewer.pal(length(d.names),"Spectral"))(nb.cols)
    barplot(XOr.CLR, legend.text=F, col=colors.OTU, axisnames=F, border=NA, xpd=T,)
    par(mar=c(0,1,1,1)+0.1,cex=1)
    plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
    legend(x="center", legend=d.names, col=colors.OTU, lwd=5, cex=.6, border=NULL)}
  #'
  list(clustering=hc,dendr=dend,tabla=table(clustHC))
}


##---------------------------------------------##
##---------------------------------------------##
#' Calcula el Compositional Biplot guapo
#' No es pot aplicar si hi ha més OTU que mostres
#' 

Biplot=function(DF){
  pcx=princomp(compositions::acomp(DF),scale.=TRUE)  
  coloredBiplot(pcx, cex=0.5, col="red",  arrow.len=0, scale=1,var.axes=TRUE,
                xlab=paste("PC1", round(pcx$sdev[1]^2 / sum(pcx$sdev^2),3), sep=": "),
                ylab=paste("PC2", round(pcx$sdev[2]^2 / sum(pcx$sdev^2),3), sep=": "),
                xlabs.col=colors, main="Form biplot")
}

##---------------------------------------------##
##---------------------------------------------##
#' Llegeix l'output d'un Aldex KW binari guardat en un Fitxer
#' j=1 mira el p-valor ajustat del test KW i j=2 el p-valor ajustatdel test glm del KW
#' ordena els p-valors i per a cada un pren tots els OTU que donen p-valor menor
#' i guarda com classifiquen
#' Finalment es queda amb el p-valor més gran que dóna millor valor d'Encerts
#' i ho guarda tot en un fitxer

Optim.HC=function(Fitxer,DF,j,Grups=tipo,print.i=TRUE){
  Dades.eBH=readRDS(Fitxer)[,c(2,4)] 
  DFTemp=DF
  GrupsTemp=Grups
  
  ps=sort(unique(Dades.eBH[,j]))
  
  enc.p=c()
  for (i in 3:length(ps)){
    if (print.i==TRUE){print(i) }
    signif.p1=which(Dades.eBH[,j]<=ps[i])  
    DF.s=DFTemp[,signif.p1] 
    CC=ClusterHC(DF.s,dendrograma=FALSE,barplot=FALSE,Grups=GrupsTemp,colores=colors)
    enc.p=rbind(enc.p,c(ps[i],Encerts(CC$tabla,as.vector(table(GrupsTemp)))))
  }
  signifHC.Ald=which(Dades.eBH[,j] <= enc.p[max(which(enc.p[,2]==min(enc.p[,2]))),1])
  taxa.signifHC.Ald=taxa[signifHC.Ald]
  DFTemp.Aldex=DFTemp[,signifHC.Ald]
  RR=list(encerts=enc.p,taxa=taxa.signifHC.Ald,Frecs=DFTemp.Aldex,Grups=GrupsTemp)
  #saveRDS(RR, file=paste("RRenc_p",j,Fitxer,sep="_"))
  if (j==1){
    saveRDS(RR, file=paste("RRenc_p",Fitxer,sep="_"))}
  else{
    saveRDS(RR, file=paste("RRenc_p",j, Fitxer,sep="_"))
  }
  return(RR)
}

##---------------------------------------------##
##---------------------------------------------##
#' Com l'anterior, però per a Aldex GLM binari

Optim.HC.glm=function(Fitxer,DF,Grups=tipo,print.i=TRUE){
  Dades.eBH=readRDS(Fitxer)[,10] 
  DFTemp=DF
  GrupsTemp=tipo
  
  ps=sort(unique(Dades.eBH))
  
  enc.p=c()
  for (i in 3:length(ps)){
    if (print.i==TRUE){print(i) }
    signif.p1=which(Dades.eBH<=ps[i])  
    DF.s=DFTemp[,signif.p1] 
    CC=ClusterHC(DF.s,dendrograma=FALSE,barplot=FALSE,Grups=GrupsTemp)
    enc.p=rbind(enc.p,c(ps[i],Encerts(CC$tabla,as.vector(table(GrupsTemp)))))
  }
  signifHC.Ald=which(Dades.eBH <= enc.p[max(which(enc.p[,2]==min(enc.p[,2]))),1])
  taxa.signifHC.Ald=taxa[signifHC.Ald]
  DFTemp.Aldex=DFTemp[,signifHC.Ald]
  RR=list(encerts=enc.p,taxa=taxa.signifHC.Ald,Frecs=DFTemp.Aldex,Grups=GrupsTemp)
  saveRDS(RR, file=paste("RRenc_p",Fitxer,sep="_"))
  return(RR)
}


##---------------------------------------------##
##---------------------------------------------##
#' Llegeix el resultat d'una funció Optim.HC(.glm), dibuixa els encerts,
#' sóna els taxa òptims, dibuixa el biplot si ha lugar i fa el clustering jeràrquic

I.Despres=function(Fitxer){
  RR=readRDS(Fitxer)
  plot(RR$encerts,type="b",pch=20)
  RR$taxa
  if(length(RR$taxa)<length(RR$Grups)){Biplot(RR$Frecs)}
  HC=ClusterHC(RR$Frecs,dendrograma=TRUE,barplot=TRUE, Grups=RR$Grups)
}


##---------------------------------------------##
##---------------------------------------------##
#' A partir del resultat d'un Aldex, pren els bitxos amb p-valor (triat amb j)
#' per davall de q i els dóna, fa biplot, clustering etc.
#' QQ.HC pren p-valors ajustats, QQ.HC.no adjust... eso
#' Bichos.Signif.Aldex dóna els bitxos, es fa apart perquè no sempre ho voldrem fer al mateix lloc

QQ.HC=function(Fitxer,DF,Grups,taxa,q=0.05,j,BP=FALSE){
  Dades=readRDS(Fitxer)[,c(2*j)] 
  signif.p1=which(Dades<q) 
  DF.s=DF[,signif.p1] 
  taxa.signif=taxa[signif.p1]
  if(length(taxa.signif)<length(Grups)){Biplot(DF.s)}
  HC=ClusterHC(DF.s,dendrograma=TRUE,barplot=BP, Grups=Grups)
  
  return(c(length(taxa.signif),Sep(HC)))
}

QQ.HC.noadjust=function(p.valores,DF,Grups,q=0.05,j,BP=FALSE){
  signif.p1=which(p.valores[,j]<q) 
  HC=ClusterHC(DF[,signif.p1] ,dendrograma=TRUE,barplot=BP, Grups=Grups)
  return(c(length(signif.p1),Sep(HC)))
}


Bichos.Signif.Aldex=function(Fitxer,q=0.05,j){
  signif.p1=which(readRDS(Fitxer)[,c(2*j)]<=q) 
  taxa.signif=taxa[signif.p1]
  data.frame("OTUS"=signif.p1,"Taxa"=taxa.signif)
}


##---------------------------------------------##
##---------------------------------------------##
#' Calcula la separació entre els dos nivells superiros d'un clustering jeràrquic

Sep=function(HC){
  HC$clustering$height[length(HC$clustering$height)]/sqrt(HC$clustering$height[length(HC$clustering$height)-1]^2+HC$clustering$height[length(HC$clustering$height)-2]^2)
}

##---------------------------------------------##
##---------------------------------------------##
#'  Calcula la matriu de logratios a partir de la matriu CLR
#'  
logratios_matrix_clr=function (x) 
{
  if (is.null(colnames(x))) 
    colnames(x) <- (1:ncol(x))
  k <- ncol(x)
  m <- nrow(x)
  lrcolnames <- NULL
  lrX <- matrix(0, m, k * (k - 1)/2)
  idlrX <- matrix(0, k * (k - 1)/2, 2)
  nameslrX <- matrix(0, k * (k - 1)/2, 2)
  colnamesx <- colnames(x)
  lloc <- 0
  for (i in (1:(k - 1))) {
    print(i)
    for (j in ((i + 1):k)) {
      lloc = lloc + 1
      idlrX[lloc, ] <- c(i, j)
      nameslrX[lloc, ] <- c(colnamesx[i], colnamesx[j])
      lrX[, lloc] <- x[, i] - x[, j]
      lrcolnames <- c(lrcolnames, paste(paste("lr", i, 
                                              sep = ""), j, sep = "."))
    }
  }
  colnames(lrX) <- lrcolnames
  results <- list(`matrix of log-ratios` = lrX, `pairs of variables in the logratio` = idlrX, 
                  `names of the variables in the logratio` = nameslrX)
  return(results)
}


##---------------------------------------------##
##---------------------------------------------##
#' Fa coda_glmnet binari a partir de la matriu de LR; 
# tunejat per emprar mínima memòria

coda_glmnet_lr_bin=function (lrx,y, taxa,
                             alpha = 0.9, nfolds = 10, 
                             showPlots = FALSE, 
                             coef_threshold = 0) 
{
  k <- length(taxa)
  lassocv <- glmnet::cv.glmnet(lrx[[1]], y, family = "binomial", 
                               alpha = alpha, type.measure = "auc", nfolds = nfolds, 
                               keep = TRUE)
  
  if (showPlots == TRUE) {
    plot(lassocv)
  }
  
  lambdavalue <- lassocv$lambda.1se
  #idrow <- max(which(lassocv$glmnet.fit$lambda >= lambdavalue))
  coeflr <- as.vector(coef(lassocv, s = lambdavalue))[-1]
  #    lrselect <- which(coeflr != 0)
  coeflogcontrast <- rep(0, k)
  for (i in (1:length(coeflr))) {
    coeflogcontrast[lrx[[2]][i, 1]] <- coeflogcontrast[lrx[[2]][i, 
                                                                1]] + coeflr[i]
    coeflogcontrast[lrx[[2]][i, 2]] <- coeflogcontrast[lrx[[2]][i, 
                                                                2]] - coeflr[i]
  }
  varlogcontrast <- which(abs(coeflogcontrast) > coef_threshold)
  coeflogcontrast <- coeflogcontrast[varlogcontrast]
  names.select <- colnames(DF)[varlogcontrast]
  positive <- ifelse(coeflogcontrast > 0, 1, 0)
  positive <- factor(positive, levels = c(0, 1), labels = c("negative", 
                                                            "positive"))
  
  plot2 <- plot_signature(names.select, coeflogcontrast, showPlots = TRUE)
  
  results <- list(taxa.num = varlogcontrast, 
                  taxa.name = names.select, 
                  `log-contrast coefficients` = coeflogcontrast, 
                  `signature plot` = plot2)   
}    


##---------------------------------------------##
##---------------------------------------------##
#' Fan clustering kmeans i fuzzy kmeans els dibuixen i donen els encerts
#' 
ClusterKM=function(X){
  DF.FF=data.frame(compositions::acomp(X))
  km.res=kmeans(DF.FF, centers=N, nstart=N)
  return(list(clustering=km.res,tabla=t(table(km.res$cluster,Grups)),dades=DF.FF))
}

ClusterFKM=function(X){
  DF.FF=data.frame(compositions::acomp(X))
  res.fcm=fcm(DF.FF, centers=N) 
  km.res=ppclust2(res.fcm, "kmeans")
  tabla=t(table(km.res$cluster,Grups))
  while(dim(tabla)[1]<dim(tabla)[2]){tabla=rbind(tabla,rep(0,dim(tabla)[2]))}
  while(dim(tabla)[1]>dim(tabla)[2]){tabla=cbind(tabla,rep(0,dim(tabla)[1]))}
  return(list(clustering=km.res,tabla=tabla,dades=DF.FF))
}
#


##---------------------------------------------##
##---------------------------------------------##
#' Adaptada de EasyCODA, que al començament em donava error
#' ara ja no ho sé
#' 
codaSeq.outlier=function(y, plot.me=TRUE, col=rgb(1,0,0,0.3),p){
  pcx=prcomp(y)
  mv=sum(pcx$sdev^2)
  sample.var= apply(pcx$x,1,function(y){sum(y^2/mv)})
  cut=quantile(apply(pcx$x,1,function(x){sum(x^2/mv)}),0.75) + p * IQR(apply(pcx$x,1,function(x){sum(x^2/mv)}))
  bad=which(apply(pcx$x,1,function(x){sum(x^2/mv)}) > cut) ##
  good=which(apply(pcx$x,1,function(x){sum(x^2/mv)}) <= cut) ##
  if(plot.me == TRUE){
    hist(sample.var, breaks=100)
    boxplot(sample.var, horizontal=TRUE, col=col, add=TRUE)
    abline(v=cut, lty=2)
  }
  return(list(sample.var=sample.var, bad=bad, good=good) )
}


##---------------------------------------------##
##---------------------------------------------##
#' 
#' Adaptada de CODA4microbiome
#' Fa servir la AUC de predicció multinomial de 
#' https://link.springer.com/article/10.1023/A:1010920819831   
#'           
ELRmulti.HT=function (x, y){
  k <- ncol(x)
  lrmatrix <- coda4microbiome::logratios_matrix(x)
  lrX <- lrmatrix[[1]]
  idlrX <- lrmatrix[[2]]
  logratio_cor <- matrix(rep(0, k * k), nrow = k)
  s = 0
  for (i in (1:(k - 1))) {
    for (j in ((i + 1):k)) {
      s <- s + 1
      lr_ij <- lrX[, s]
      model1=nnet::multinom(y ~ lr_ij, trace = F)
      result1= predict(model1, lr_ij, type='probs')
      logratio_cor[i, j]=
        HandTill2001::auc(HandTill2001::multcap(response = y,
                                                predicted = as.matrix(result1)))
      logratio_cor[j, i] <- logratio_cor[i, j]
    }
  }
  o <- order(colSums(abs(logratio_cor)), decreasing = T)
  M <- logratio_cor[o, o]
  colnames(M) <- o
  rownames(M) <- colnames(M)
  #
  results <- list(
    `max log-ratio` = colnames(M)[which(M == max(abs(M)), arr.ind = TRUE)[(2:1)]],
    `names max log-ratio` = colnames(x)[as.numeric(colnames(M)[which(M ==max(abs(M)), arr.ind = TRUE)[(2:1)]])], 
    `order of importance` = o, 
    `most important variables` = colnames(x)[o], 
    `association matrix` = M
  )
  return(results)
}


##---------------------------------------------##
##---------------------------------------------##
#' 
#' Separat el gràfic del càlcul de l'anterior
#' X el resultat de ELRmulti.HT()
#' 
Gràfic.CP=function(X,shownames = FALSE, maxrow = 15, maxcol = 15){
  title <- "AUC multinomial regression y~log(xi/xj)"
  GnBu8 <- c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", 
             "#4eb3d3", "#2b8cbe", "#08589e")
  Blues4 <- c("#eff3ff", "#bdd7e7", "#6baed6", "#2171b5")
  col2 <- grDevices::colorRampPalette(GnBu8, space = "Lab")
  k <- ncol(X$`association matrix`)
  if (maxrow > k) maxrow <- k
  if (maxcol > k) maxcol <- k
  o=X$`order of importance`
  M=X$`association matrix`
  if (shownames == TRUE) {rownames(M) <- colnames(x)[o]}
  corrplot::corrplot(M[(1:maxrow), (1:maxcol)],
                     tl.pos = "lt", tl.col="black", title = title, mar =c(0, 0, 1, 0), 
                     col.lim = c(min(M), max(M)), col = col2(200), 
                     is.corr = FALSE)
}
#
#'


##---------------------------------------------##
##---------------------------------------------##
#' 
#
#' Adaptacions d'ALDEx
#' 

#'rdirichlet usual, però no pot donar 0
#'

rdirichlet.nn <- function (n, alpha)
{
  n <- as.integer(n)
  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  
  x[which(x==0)]=10^(-100)
  return(x / rowSums(x))
}


##---------------------------------------------##
##---------------------------------------------##
#' 

#' aldexCesc.clr.function d'ALDEx, però (a) empra la rdirichlet anterior;
#' (b) s'aplica a matriu amb zeros ja imputats
#'

aldexCesc.clr <- function( reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE) {
  
  # INPUT
  # The 'reads' data.frame MUST have row
  # and column names that are unique, and
  # looks like the following:
  #
  #              T1a T1b  T2  T3  N1  N2
  #   Gene_00001   0   0   2   0   0   1
  #   Gene_00002  20   8  12   5  19  26
  #   Gene_00003   3   0   2   0   0   0
  #       ... many more rows ...
  #
  # ---------------------------------------------------------------------
  
  # OUTPUT
  # The output returned is a list (x) that contains Monte-Carlo instances of
  # the centre log-ratio transformed values for each sample
  # Access to values
  # sample IDs: names(x)
  # number of features (genes, OTUs): length(x[[1]][,1])
  # number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
  # feature names: rownames(x[[1]])
  
  
  # Fully validate and coerce the data into required formats
  # make sure that the multicore package is in scope and return if available
  has.BiocParallel <- FALSE
  if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
    print("multicore environment is is OK -- using the BiocParallel package")
    #require(BiocParallel)
    has.BiocParallel <- TRUE
  }
  else {
    print("operating in serial mode")
  }
  
  # make sure that mc.samples is an integer, despite it being a numeric type value
  as.numeric(as.integer(mc.samples))
  
  #  remove all rows with reads less than the minimum set by minsum
  minsum <- 0
  
  # remove any row in which the sum of the row is 0
  z <- as.numeric(apply(reads, 1, sum))
  reads <- as.data.frame( reads[(which(z > minsum)),]  )
  
  if (verbose) print("removed rows with sums equal to zero")
  
  
  #  SANITY CHECKS ON THE DATA INPUT
  # if ( any( round(reads) != reads ) ) stop("not all reads are integers")
  if ( any( reads < 0 ) )             stop("one or more reads are negative")
  
  for ( col in names(reads) ) {
    if ( any( ! is.finite( reads[[col]] ) ) )  stop("one or more reads are not finite")
  }
  
  if ( length(rownames(reads)) == 0 ) stop("rownames(reads) cannot be empty")
  if ( length(colnames(reads)) == 0 ) stop("colnames(reads) cannot be empty")
  
  if ( length(rownames(reads)) != length(unique(rownames(reads))) ) stop ("row names are not unique")
  if ( length(colnames(reads)) != length(unique(colnames(reads))) ) stop ("col names are not unique")
  if ( mc.samples < 128 ) warning("values are unreliable when estimated with so few MC smps")
  
  # add a prior expection to all remaining reads that are 0
  # this should be by a Count Zero Multiplicative approach, but in practice
  # this is not necessary because of the large number of features
  prior <- 0
  
  # This extracts the set of features to be used in the geometric mean computation
  feature.subset <- aldex.set.mode(reads, conds, denom)
  
  
  reads <- reads + prior
  
  if (verbose == TRUE) print("data format is OK")
  
  # ---------------------------------------------------------------------
  # Generate a Monte Carlo instance of the frequencies of each sample via the Dirichlet distribution,
  # returns frequencies for each feature in each sample that are consistent with the
  # feature count observed as a proportion of the total counts per sample given
  # technical variation (i.e. proportions consistent with error observed when resequencing the same library)
  
  nr <- nrow( reads )
  rn <- rownames( reads )
  
  #this returns a list of proportions that are consistent with the number of reads per feature and the
  #total number of reads per sample
  
  # environment test, runs in multicore if possible
  if (has.BiocParallel){
    p <- bplapply( reads ,
                   function(col) {
                     q <- t( rdirichlet.nn( mc.samples, col ) ) ;
                     rownames(q) <- rn ;
                     q })
    names(p) <- names(reads)
  }
  else{
    p <- lapply( reads ,
                 function(col) {
                   q <- t( rdirichlet.nn( mc.samples, col ) ) ;
                   q[q==0]=10^(-100);
                   rownames(q) <- rn ; q } )
  }
  
  # sanity check on the data, should never fail
  for ( i in 1:length(p) ) {
    if ( any( ! is.finite( p[[i]] ) ) ) stop("non-finite frequencies estimated")
  }
  
  if (verbose == TRUE) print("dirichlet samples complete")
  
  # ---------------------------------------------------------------------
  # Take the log2 of the frequency and subtract the geometric mean log2 frequency per sample
  # i.e., do a centered logratio transformation as per Aitchison
  
  # apply the function over elements in a list, that contains an array
  
  # DEFAULT
  if(length(feature.subset) == nr)
  {
    # Default ALDEx2
    if (has.BiocParallel){
      l2p <- bplapply( p, function(m) {
        apply( log2(m), 2, function(col) { col - mean(col) } )
      })
      names(l2p) <- names(p)
    }
    else{
      l2p <- lapply( p, function(m) {
        apply( log2(m), 2, function(col) { col - mean(col) } )
      })
    }
  } else {
    ## IQLR or ZERO
    feat.result <- vector("list", length(unique(conds))) # Feature Gmeans
    condition.list <- vector("list", length(unique(conds)))    # list to store conditions
    
    for (i in 1:length(unique(conds)))
    {
      condition.list[[i]] <- which(conds == unique(conds)[i]) # Condition list
      feat.result[[i]] <- lapply( p[condition.list[[i]]], function(m) {
        apply(log2(m), 2, function(x){mean(x[feature.subset[[i]]])})
      })
    }
    set.rev <- unlist(feat.result, recursive=FALSE) # Unlist once to aggregate samples
    p.copy <- p
    for (i in 1:length(set.rev))
    {
      p.copy[[i]] <- as.data.frame(p.copy[[i]])
      p[[i]] <- apply(log2(p.copy[[i]]),1, function(x){ x - (set.rev[[i]])})
      p[[i]] <- t(p[[i]])
    }
    l2p <- p    # Save the set in order to generate the aldexCesc.clr variable
  }
  
  
  # sanity check on data
  for ( i in 1:length(l2p) ) {
    if ( any( ! is.finite( l2p[[i]] ) ) ) stop("non-finite log-frequencies were unexpectedly computed")
  }
  if (verbose == TRUE) print("clr transformation complete")
  
  return(new("aldex.clr",reads=reads,conds=conds,mc.samples=mc.samples,verbose=verbose,useMC=useMC,analysisData=l2p))
}

Indicador=function(x,k=n.OTUs.original){
  xx=rep(0,k)
  y=as.numeric(gsub("\\D+", "", OTUs[x]))
  xx[y]=1
  
  return(xx)
}