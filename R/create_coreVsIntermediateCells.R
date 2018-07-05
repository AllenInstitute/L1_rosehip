##################################################################

print("Read in the data and load relevant libraries")

outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/human/"

library(dplyr);
library(feather)
library(ggplot2)
library(parallel)
library(Biobase)
library(gplots)

Samp.dat <- read_feather(paste(inputFolder,"anno.feather",sep="")) 
Expr.dat <- feather(paste(inputFolder,"data.feather",sep=""))   # FPKM
Expr.dat <- Expr.dat[match(Samp.dat$sample_id,Expr.dat$sample_id),] # Make sure the expression matches the sample information
datIn    <- as.matrix(Expr.dat[,names(Expr.dat)!="sample_id"])
rownames(datIn) <- Expr.dat$sample_id
datIn    <- t(datIn)

##################################################################

print("Load required functions.")

map.by.cor <- function(train.dat, train.cl, test.dat){
  cl.med= do.call("cbind",tapply(names(train.cl), train.cl, function(x){
    rowMedians(train.dat[,x,drop=F])
  }))
  row.names(cl.med)=row.names(train.dat)
  test.cl.cor = cor(as.matrix(test.dat), cl.med)
  pred.cl= setNames(colnames(test.cl.cor)[apply(test.cl.cor, 1, which.max)], row.names(test.cl.cor))
  if(is.factor(train.cl)){
    pred.cl = setNames(factor(pred.cl, levels=levels(train.cl)), names(pred.cl))
  }
  return(list(pred.cl=pred.cl,pred.score=test.cl.cor))
}

test.cv.cor <- function(norm.dat, cl, markers, n.bin=5,g.perc=1){
  bins=unlist(tapply(names(cl), cl, function(x){
    if(length(x) > n.bin){
      tmp=rep_len(1:n.bin, length(x))
    }else{
      tmp = sample(1:n.bin, length(x))
    }
    setNames(tmp[sample(length(tmp))], x)
  }))
  names(bins) = gsub(".*\\.", "", names(bins))
  bins= bins[names(cl)]
  pred.cl = setNames(rep(NA, length(cl)), names(cl))
  for(i in 1:n.bin){
    train.cells = names(cl)[bins!=i]
    test.cells =names(cl)[bins==i]
    select.markers=sample(markers, round(length(markers)*g.perc))
    pred.cl[test.cells] <- as.character(map.by.cor(norm.dat[select.markers,], cl[train.cells], norm.dat[select.markers, test.cells])$pred.cl)
  }
  return(pred.cl)
}

getBetaScore <- function(propExpr,returnScore=TRUE,spec.exp = 2){
  # Marker score is combination of specificity and sparsity
  # propExpr = proportions of cells in a given cluster with CPM/FPKM > 1 (or 0, HCT uses 1)
  calc_beta <- function(y, spec.exp = 2) {
    d1 <- as.matrix(dist(y))
    eps1 <- 1e-10
    score1 <- sum(d1^spec.exp) / (sum(d1) + eps1)
    return(score1)
  }
  betaScore <- apply(propExpr, 1, calc_beta)
  betaScore[is.na(betaScore)] <- 0
  if(returnScore) return(betaScore)
  scoreRank = rank(-betaScore)
  return(scoreRank)
}

##################################################################

print("Subset the cells to interneurons and find the top 1200 marker genes")

kpSamp     <- Samp.dat$cluster_type_label=="inh"
cl3        <- factor(Samp.dat$cluster_label,levels = sort(unique(Samp.dat$cluster_label)))
names(cl3) <- colnames(datIn)
clustersF  <- droplevels(cl3[kpSamp])

norm.dat   <- log2(datIn[,kpSamp]+1)
propExpr   <- do.call("cbind", tapply(names(clustersF), clustersF, function(x) rowMeans(norm.dat[,x]>1))) 
propExpr   <- propExpr[,levels(clustersF)]
rownames(propExpr) <- rownames(norm.dat) 

topBeta    <- getBetaScore(propExpr,FALSE)
kpGene     <- topBeta<=1200
markers    <- sort(names(kpGene)[kpGene])

colnames(norm.dat) <- names(clustersF) <- gsub("\\.","_",names(clustersF))

##################################################################

print("Find the core vs. intermediate cells.")

cl.cv <-mclapply(1:100, function(i){
  set.seed(i)
  tmp=test.cv.cor(norm.dat, clustersF, markers, n.bin=5)
}, mc.cores=1)

corMap = do.call("rbind",sapply(cl.cv, function(x){
  df = data.frame(cell=names(x),cl=x)
},simplify=F))
corMap = table(corMap[,1],corMap[,2])
corMap = corMap / rowSums(corMap)

isCoreCell = apply(corMap,1,function(x) max(x)>=0.9)
clusterMatrix = matrix(NA,nrow=dim(corMap)[1],ncol=4)
clusterMatrix[,1] = colnames(corMap)[apply(corMap,1,which.max)]
rownames(clusterMatrix) = rownames(corMap)
clusterMatrix[,4] = as.character(clustersF[rownames(clusterMatrix)])
colnames(clusterMatrix) = c("PrimaryCluster","SecondaryCluster","TertiaryCluster","AssignedCluster")
for (i in which(!isCoreCell)) {
  out = -sort(-corMap[i,])
  out = out[out>0]
  l   = min(3,length(out))
  clusterMatrix[i,1:l] = names(out)[1:l]
}
isCoreCell = isCoreCell&(clusterMatrix[,1]==clusterMatrix[,4])

print(table(apply(clusterMatrix,1,function(x) length(unique(x[(!is.na(x))])))))
#   1   2   4 
# 443  26   1
sum(!isCoreCell)
# 27 intermediate cells and only 1 maps to >2 types


##################################################################

print("Find cluster confidences.")

clusterOverlap = matrix(0,nrow=length(levels(clustersF)),ncol = length(levels(clustersF)))
rownames(clusterOverlap) <- colnames(clusterOverlap) <- levels(clustersF)
clMat = clusterMatrix[!isCoreCell,]

for (i in 1:dim(clMat)[1]) for (j in 2:4) if(!is.na(clMat[i,j])) {
  clusterOverlap[clMat[i,1],clMat[i,j]] = clusterOverlap[clMat[i,1],clMat[i,j]] + 1
  clusterOverlap[clMat[i,j],clMat[i,1]] = clusterOverlap[clMat[i,j],clMat[i,1]] + 1
}

clusterMatrix = cbind(clusterMatrix,isCoreCell)
save(corMap,clusterMatrix,clusterOverlap,file=paste0(outputFolder,"coreCells_clusterConfidences.RData"))

