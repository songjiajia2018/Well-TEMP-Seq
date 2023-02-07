library(MASS)
library(nloptr)
args <- commandArgs(T)
sample_name <- args[1]
pqdata <- read.table(paste0(sample_name,'_sampling.txt'),row.names = 1)

LL <- function(params,data){
  t1 <- dbinom(data[,2],data[,1],params[1])
  t2 <- dbinom(data[,2],data[,1],params[2])
  l <- params[3]*t1+(1-params[3])*t2
  l <- sum(log(l))
  return(-l)
}

r_mtx <- data.frame()
i <- 1
while(i < 101){
  init <- runif(3)
  result <- neldermead(init,LL,data=pqdata)
  par <- result$par
  value <- result$value
  if(par[3] >=0 && par[3] <= 1){
    r_mtx[i,1] <- par[1]
    r_mtx[i,2] <- par[2]
    r_mtx[i,3] <- value
    i <- i + 1
  }
}

result <- r_mtx[which.min(r_mtx$V3),]
p <- as.numeric(result[1])
q <- as.numeric(result[2])

write.table(data.frame(p=max(p,q),q=min(p,q)), file = paste0(sample_name,'_pq.txt'),row.names=F)

genedata <- read.table(paste0(sample_name,'_count.txt'), row.names = 1)
genedata <- split(genedata,genedata$V4)

LL_2 <- function(params,data){
  t1 <- dbinom(data[,2],data[,1],p)
  t2 <- dbinom(data[,2],data[,1],q)
  l <- params*t1+(1-params)*t2
  l <- sum(log(l))
  return(-l)
}


theta_for_gene <- lapply(genedata, function(x){
  umi_num <- dim(x)[1]
  if(umi_num >= 100){
    result=optim(par = runif(1),fn = LL_2,data = x, method = 'Brent',lower = 1e-8, upper = 1)
    theta <- result$par
  }
  else{theta <- NULL}
  return(theta)
})



theta_for_gene <- do.call(rbind, theta_for_gene)
if(p < q){
    theta_for_gene <- 1-theta_for_gene
}
theta_for_gene <- as.data.frame(theta_for_gene)
write.table(theta_for_gene, file = paste0(sample_name,'_theta_for_gene.txt'),col.names=F)
theta_for_gene$gene <- rownames(theta_for_gene)

mtx <- readRDS(paste0(sample_name,'_TC_matrix.rds'))
mtx <- as.matrix(mtx)

mtx.n <- mtx[grep('*--C',rownames(mtx)),]
mtx.o <- mtx[grep('*--T',rownames(mtx)),]
rownames(mtx.n) <- strsplit(rownames(mtx.n),'--C')
rownames(mtx.o) <- strsplit(rownames(mtx.o),'--T')
allgene <- unique(c(rownames(mtx.n),rownames(mtx.o)))
complete_mtx <- matrix(data = 0, nrow = length(allgene), ncol = dim(mtx)[2])
colnames(complete_mtx) <- colnames(mtx)
rownames(complete_mtx) <- allgene
complete_mtx_new <- complete_mtx
complete_mtx_old <- complete_mtx
tm.m <- complete_mtx_new[-which(rownames(complete_mtx_new) %in% rownames(mtx.n)),]
complete_mtx_new <- rbind(tm.m, mtx.n)
tm.m <- complete_mtx_old[-which(rownames(complete_mtx_old) %in% rownames(mtx.o)),]
complete_mtx_old <- rbind(tm.m, mtx.o)
complete_mtx_new <- complete_mtx_new[allgene,]
complete_mtx_old <- complete_mtx_old[allgene,]
complete_mtx <- complete_mtx_old + complete_mtx_new

theta_for_gene <- theta_for_gene[which(rownames(theta_for_gene) %in% rownames(complete_mtx)),]
mtx_filterby_theta <- complete_mtx[rownames(theta_for_gene),]
mtx_new_filterby_theta <- complete_mtx_new[rownames(theta_for_gene),]

add_theta <- cbind(mtx_filterby_theta,theta_for_gene$V1)
estN_limited_gene <- apply(add_theta,1,function(x){
  x <- x*x[length(x)]
})
estN_limited_gene <- estN_limited_gene[-dim(estN_limited_gene)[1],]


estcellN <- colSums(t(estN_limited_gene))
labelcellN <- colSums(mtx_new_filterby_theta)
labelcellN <- labelcellN[names(estcellN)]
alphacell <- labelcellN/estcellN

discard_a <- alphacell[which(alphacell > 1)]
alpha <- alphacell[which(alphacell <= 1)]
write.table(as.data.frame(alpha),file = paste0(sample_name,'_alpha.txt'),col.names=F)

complete_mtx <- complete_mtx[allgene,names(alpha)]
complete_mtx_new <- complete_mtx_new[allgene,names(alpha)]
estAllN <- rbind(complete_mtx_new,alpha)
estAllN <- apply(estAllN, 2, function(x){
  a <- x[length(x)]
  x <- x[-length(x)]
  x <- round(x/a)
})

estNew <- pmin(estAllN,complete_mtx)
estOld <- complete_mtx - estNew

write.csv(estNew,file = paste0(sample_name,'_pre_new.csv'),row.names = T)
write.csv(estOld,file = paste0(sample_name,'_pre_old.csv'),row.names = T)



