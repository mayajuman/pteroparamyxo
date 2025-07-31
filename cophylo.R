###cophylogenetic analysis

#####cophylogeny stuff
library(paco)
library(ape)
library(vegan)
library(ggtree)
library(tangler)
library(phytools)
library(knitr)
library(stringr)
library(plyr)
library(Hmisc)
library(dplyr)
library(patchwork)

###read in data
setwd("~/Documents/PhD/1. pteropodidae/pteroparamyxo submission")
rev <- read.csv("Supplementary Data 1.csv")

#trees
setwd("~/Documents/PhD/1. pteropodidae/phylo")

#host
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
## trim to family
taxa=taxa[taxa$fam=="PTEROPODIDAE",]
taxa$tip=taxa$Species_Name
## trim tree
tree=keep.tip(tree,taxa$tiplabel)
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))

#virus
vtree=read.tree('ONSR.Paramyxo.Fig3.v6.tree.nwk')
vtree$tip.label=sapply(strsplit(vtree$tip.label,'_'),function(x) paste(x[2],x[3],x[4],sep='_'))
vtree$tip.label <- str_replace(vtree$tip.label,"_NA", "")

####create distance matrix
matrix <- matrix(nrow = length(tree[["tip.label"]]), ncol = length(vtree[["tip.label"]]))
colnames(matrix) <- vtree[["tip.label"]]
rownames(matrix) <- tree[["tip.label"]]
matrix[,] <- 0

rev$Host.species <- plyr::revalue(rev$Host.species,
                            c("magna"="moluccensis",
                              "medius"="giganteus",
                              "melanotus natalis"="melanotus",
                              "natalis"="melanotus",
                              "seychellensis comorensis"="seychellensis",
                              "labiatus minor"="labiatus"))

#calculate missing prevalence
rev$Prevalence <- ifelse(is.na(rev$Prevalence), rev$nPositive/rev$nSamples, rev$Prevalence)

knowns <- rev[which(rev$Prevalence > 0 | rev$Isolation.successful > 0),]
knowns <- knowns[which(!is.na(knowns$Virus.strain.or.sequence)),]
knowns <- knowns[which(!is.na(knowns$Host.species)),]
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "BatPV/Eid.hel/GH10/2008, BatPV/Eid.hel/GH45/2008, BatPV/Eid.hel/GH48/2008")] <- "Ghana_henipavirus"
knowns <- knowns[!grepl("-|\\.|/", knowns$Virus.strain.or.sequence),]
knowns <- knowns[!grepl(",|\\.|/", knowns$Host.genus),]
knowns <- knowns[!grepl(",|\\.|/", knowns$Host.species),]
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "achimota virus 1")] <- "Achimota_pararubulavirus_1"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "achimota virus 2")] <- "Achimota_pararubulavirus_2"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "alston virus")] <- "Alston_orthorubulavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "cedar virus")] <- "Cedar_henipavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "hendra virus")] <- "Hendra_henipavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "hervey virus")] <- "Hervey_pararubulavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "human mumps virus")] <- "Mumps_orthorubulavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "menangle virus")] <- "Menangle_pararubulavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "nipah virus")] <- "Nipah_henipavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "sosuga virus")] <- "Sosuga_pararubulavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "teviot virus")] <- "Teviot_pararubulavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "tioman virus")] <- "Tioman_pararubulavirus"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "tuhoko virus 1")] <- "Tuhoko_pararubulavirus_1"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "tuhoko virus 2")] <- "Tuhoko_pararubulavirus_2"
knowns$Virus.strain.or.sequence[which(knowns$Virus.strain.or.sequence == "tuhoko virus 3")] <- "Tuhoko_pararubulavirus_3"

###cut serological detections -- not specific enough
knowns <- knowns[-which(knowns$Detection.method == "serology"),]

knowns$host <- with(knowns, paste0(Host.genus," ",Host.species)) #create host name variable
knowns$host <- gsub(" ","_",capitalize(knowns$host))

knowns <- knowns %>% dplyr::select(host, Virus.strain.or.sequence, Prevalence)

for (i in 1:nrow(knowns)) {
  matrix[rownames(matrix) == knowns$host[i],colnames(matrix) == knowns$Virus.strain.or.sequence[i]] <- 1
}

###paco

NLinks = sum(matrix)
host.D <- cophenetic(tree)
vir.D <- cophenetic(vtree)

host.D <- host.D[rownames(matrix),rownames(matrix)]
vir.D <- vir.D[colnames(matrix),colnames(matrix)] 

PACo <- function (H.dist, P.dist, HP.bin) {
  HP.bin <- which(HP.bin > 0, arr.in=TRUE)
  H.PCo <- pcoa(H.dist, correction="cailliez")$vectors
  P.PCo <- pcoa(P.dist, correction="cailliez")$vectors
  H.PCo <- H.PCo[HP.bin[,1],]
  P.PCo <- P.PCo[HP.bin[,2],]
  list (H.PCo = H.PCo, P.PCo = P.PCo)
}

PACo.fit <- PACo(host.D, vir.D, matrix)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo)

##goodness of fit test
m2.obs <- HP.proc$ss
N.perm = 5000
P.value = 0
set.seed(.Random.seed[trunc(runif(1,1,626))]) 

for (n in c(1:N.perm)) {
  if (NLinks <= nrow(matrix) | NLinks <= ncol(matrix)) {
    flag2 <- TRUE
    while (flag2 == TRUE) {
      HP.perm <- t(apply(matrix,1,sample))
      if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE
      else flag2 <- FALSE
    } 
  } 
  else { 
    HP.perm <- t(apply(matrix,1,sample))}
  PACo.perm <- PACo(host.D, vir.D, HP.perm)
  m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss
  #write(m2.perm, file="PACo/example/m2_perm.txt", sep="\t", append=TRUE)
  if (m2.perm <= m2.obs){P.value = P.value + 1}
}

P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", P.value, " based on ", N.perm,"
permutations.") 

###extract residuals

HP.ones <- which(matrix > 0, arr.in=TRUE) 
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)

HostX <- HP.proc$X
VirY <- HP.proc$Yrot 

colnames(SQres.jackn) <- paste(rownames(HostX),rownames(VirY), sep="-")
t.critical = qt(0.975,NLinks-1)

for(i in c(1:NLinks)) {
  HP.ind <- matrix
  HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
  PACo.ind <- PACo(host.D, vir.D, HP.ind)
  Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo)
  res.Proc.ind <- c(residuals(Proc.ind))
  res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
  SQres.jackn [i, ] <- res.Proc.ind
}

SQres.jackn <- SQres.jackn**2 
SQres <- residuals(HP.proc)**2
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres))
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE)
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE)
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks)

pat.bar <- barplot(phi.mean, names.arg = " ", space = 0.25, col="white", 
                   xlab="Host-parasite link", ylab= "Squared residuals", 
                   ylim=c(0, max(phi.UCI)), cex.lab=1.2)
text(pat.bar, par("usr")[3] - 0.001, srt = 330, adj = 0, labels = colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.6)
arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
abline(a=median(phi.mean), b=0, lty=2)

###visualizing

knowns <- knowns[,1:2]
knowns <- knowns[-which(knowns$Virus.strain.or.sequence == "unclassified EPMV"),]
knowns <- knowns[-which(knowns$Virus.strain.or.sequence == "dawn bat paramyxovirus"),]
knowns <- knowns[-which(knowns$Virus.strain.or.sequence == "bat parainfluenza"),]
knowns <- knowns[-which(knowns$Virus.strain.or.sequence == "geelong paramyxovirus"),]
knowns <- knowns[-which(knowns$Virus.strain.or.sequence == "yarra bend paramyxovirus"),]
knowns <- knowns[-which(knowns$Virus.strain.or.sequence == "angavokely virus"),]
knowns <- knowns[-which(knowns$Virus.strain.or.sequence == "achimota virus 3"),]
knowns <- knowns %>% distinct()

knowns$residuals <- NA
phi.mean <- as.data.frame(phi.mean)
phi.UCI <- as.data.frame(phi.UCI)
phi <- cbind(phi.mean,phi.UCI)

for (i in 1:nrow(knowns)) {
  knowns$residuals[i] <- phi.mean$phi.mean[which(rownames(phi.mean) == paste(c(knowns$host[i],"-",knowns$Virus.strain.or.sequence[i]), collapse=''))]
}

#smaller residuals = thicker lines = more support for co-evolution
knowns$residualsrev <- 3.5-log10(knowns$residuals)

vtree2=keep.tip(vtree,knowns$Virus.strain.or.sequence)
tree2=keep.tip(tree,knowns$host)

###label offset
knowns$host <- paste(knowns$host," ")
knowns$Virus.strain.or.sequence <- paste(" ",knowns$Virus.strain.or.sequence)
tree2$tip.label <- paste(tree2$tip.label," ")
vtree2$tip.label <- paste(" ",vtree2$tip.label)

#tanglegram
tanglegram2 <- cophylo(tree2,vtree2,assoc=knowns)

setwd("~/Documents/PhD/pteropodidae")
png("Figure 6.png",width=8,height=5,units="in",res=600)
plot(tanglegram2, fsize=c(1,1), link.lty="solid", 
     link.lwd=knowns$residualsrev, ftype=c("i","i"))
dev.off()
