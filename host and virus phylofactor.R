####host phylofactor

library(phylofactor)
library(dplyr)
library(tidyr)
library(stringr)
library(caper)
library(treeio)
library(ggtree)
library(scales)
library(knitr)
library(patchwork)

## load in data
setwd("~/Documents/PhD/pteropodidae")
rev <- read.csv("pteroparamyxo.csv")

####host phylofactor

##aligning with upham et al. taxonomy
rev$Host.species <- revalue(rev$Host.species,
                            c("magna"="moluccensis",
                              "medius"="giganteus",
                              "melanotus natalis"="melanotus",
                              "natalis"="melanotus",
                              "seychellensis comorensis"="seychellensis",
                              "labiatus minor"="labiatus"))

rev$Host.genus <- revalue(rev$Host.genus,
                          c("lissonycteris"="myonycteris"))

#creating df with # of samples and studies per species

sp <- rev[!grepl(",|\\.", rev$Host.species),]
sp <- sp[!is.na(sp$Host.species),]
sp$host <- with(sp, paste0(Host.genus," ",Host.species)) #create host name variable
sp$host <- gsub(" ","_",capitalize(sp$host))
sp$Host.species[which(sp$Host.genus == "macroglossus" & sp$Host.species == "minimus")] <- "Mminimus"
sp$Host.species[which(sp$Host.genus == "rousettus" & sp$Host.species == "celebensis")] <- "Rcelebensis"

sp <- sp %>% group_by(host, Host.species) %>% dplyr::summarise(nSamp = sum(nSamples, na.rm=TRUE))

sp_studies <- rev %>% dplyr::select(Title, Host.genus, Host.species) %>% drop_na(Host.species)
sp_studies$Host.species[which(sp_studies$Host.genus == "macroglossus" & sp_studies$Host.species == "minimus")] <- "Mminimus"
sp_studies$Host.species[which(sp_studies$Host.genus == "rousettus" & sp_studies$Host.species == "celebensis")] <- "Rcelebensis"
sp_studies <- sp_studies[!grepl("\\.", sp_studies$Host.species),]
sp_studies <- sp_studies %>%
  mutate(Host.species = strsplit(as.character(Host.species), ",")) %>%
  unnest(Host.species) %>% 
  mutate(Host.species = str_trim(as.character(Host.species))) %>%
  mutate(Host.species=dplyr::recode(Host.species, "magna"="moluccensis")) %>%
  group_by(Host.species) %>% 
  dplyr::summarise(Studies = n_distinct(Title))

sp <- join_all(list(sp, sp_studies), by='Host.species', type='left')
sp <- sp %>% dplyr::select(-Host.species)
rm(sp_studies)

## load in bat phylogeny and taxonomy
setwd("~/Documents/PhD/pteropodidae/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)

## trim to family
taxa=taxa[taxa$fam=="PTEROPODIDAE",]
taxa$tip=taxa$Species_Name

## trim tree
tree=keep.tip(tree,taxa$tiplabel)

tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$host=taxa$Species_Name

## all species
adata=data.frame(host=tree$tip.label)

## merge
sdata=merge(adata,sp,by="host",all=T)
rm(adata)

## no studies
sdata$Studies=ifelse(is.na(sdata$Studies),0,sdata$Studies)

## no bats
sdata$nSamp=ifelse(is.na(sdata$nSamp),0,sdata$nSamp)

## binary
sdata$binstudy=ifelse(sdata$Studies==0,0,1)

## merge with taxa
taxa=taxa[c("host","tiplabel","gen","fam","clade","tip")]
sdata=merge(sdata,taxa,by="host")

sdata$tip=sdata$host
cdata=comparative.data(phy=tree,data=sdata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)
cdata$data$label=cdata$data$host
cdata$data$Species=cdata$data$host
cdata$data$tip=cdata$data$host

## taxonomy
cdata$data$taxonomy=with(cdata$data,paste(fam,gen,tip,sep='; '))

## just sampled bats
sdata=cdata[cdata$data$binstudy==1,]

## D statistic on sampled/not sampled
set.seed(1)
dstat=phylo.d(data=cdata,binvar=binstudy,permut=1000)
dstat

## number of studies and number of bats
hist(log10(sdata$data$Studies))
hist(log10(sdata$data$nSamp))

## transform
sdata$data$lstudies=log10(sdata$data$Studies)
sdata$data$ltested=log10(sdata$data$nSamp)

## range
range(sdata$data$Studies)
range(sdata$data$nSamp)

## pagel's lambda
pmod1=pgls(lstudies~1,data=sdata,lambda="ML") ## lambda = 0.00
pmod2=pgls(ltested~1,data=sdata,lambda="ML") ## lambda = 0.588

## summarize
summary(pmod1)
summary(pmod2)

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
pfsum=function(pf,letter){
  
  ## get formula
  chars=as.character(pf$frmla.phylo)[2]
  #chars=strsplit(chars," ")[[1]]
  
  ## response
  resp=chars[1]
  
  ## fix
  #resp=ifelse(resp=='cbind(pos, neg)','prevalence',resp)
  resp=ifelse(str_detect(resp,"cbind"),"prevalence",resp)
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$tip%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=paste(letter,i,sep = "")
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}


## GPF for study binary
set.seed(1)
study_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=binstudy~phylo,
             family=binomial,
             algorithm='phylo',nfactors=2,
             min.group.size = 3)

## summarize
HolmProcedure(study_pf)
study_res=pfsum(study_pf,letter="A")$results

## lower/greater
study_res$check=ifelse(study_res$clade>study_res$other,"more","less")
table(study_res$check)

## number of studies
set.seed(1)
nstudies_pf=gpf(Data=sdata$data,tree=sdata$phy,
                frmla.phylo=Studies~phylo,
                family=poisson,
                algorithm='phylo',nfactors=5,
                min.group.size = 3)

## summarize
HolmProcedure(nstudies_pf)
nstudies_res=pfsum(nstudies_pf, letter="B")$results

## lower/greater
nstudies_res$check=ifelse(nstudies_res$clade>nstudies_res$other,"more","less")
table(nstudies_res$check)

## number of samples
set.seed(1)

sdata$data$lnSamp=log1p(sdata$data$nSamp)

nsamples_pf=gpf(Data=sdata$data,tree=sdata$phy,
                frmla.phylo=lnSamp~phylo,
                family=gaussian,
                algorithm='phylo',nfactors=25,
                min.group.size = 3)

## summarize
HolmProcedure(nsamples_pf)
nsamples_res=pfsum(nsamples_pf, letter="C")$results

## lower/greater
nsamples_res$check=ifelse(nsamples_res$clade>nsamples_res$other,"more","less")
table(nsamples_res$check)

nsamples_res <- nsamples_res[which(nsamples_res$tips > 2),]

## save trees
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")
stree=treeio::full_join(as.treedata(sdata$phy),sdata$data,by="label")


#######viral phylofactor
v <- rev
v$Virus.strain.or.sequence[which(v$Virus.strain.or.sequence == "BatPV/Eid.hel/GH10/2008, BatPV/Eid.hel/GH45/2008, BatPV/Eid.hel/GH48/2008")] <- "Ghana_henipavirus"

v <- v[!grepl("-|\\.|/", v$Virus.strain.or.sequence),]

v <- v %>%
  mutate(Virus.strain.or.sequence = strsplit(as.character(Virus.strain.or.sequence), ",")) %>%
  unnest(Virus.strain.or.sequence) %>% 
  mutate(Virus.strain.or.sequence = str_trim(as.character(Virus.strain.or.sequence))) %>%
  group_by(Virus.strain.or.sequence) %>% 
  dplyr::summarise(Studies = n_distinct(Title),
                   nSamp = sum(nSamples, na.rm=TRUE))

v$nSamp[which(v$Virus.strain.or.sequence == "hervey virus")] <- 9 #we don't know sample size, but there were 9 positives, so it's a minimum of 9
v <- v %>% dplyr::rename(virus = Virus.strain.or.sequence)

##cleaning and matching
v <- v[-which(v$virus == "YN2017A"),]
v$virus[which(v$virus == "achimota virus 1")] <- "Achimota_pararubulavirus_1"
v$virus[which(v$virus == "achimota virus 2")] <- "Achimota_pararubulavirus_2"
v <- v[-which(v$virus == "achimota virus 3"),]
v$virus[which(v$virus == "alston virus")] <- "Alston_orthorubulavirus"
v <- v[-which(v$virus == "angavokely virus"),]
v <- v[-which(v$virus == "bat parainfluenza"),]
v$virus[which(v$virus == "cedar virus")] <- "Cedar_henipavirus"
v <- v[-which(v$virus == "dawn bat paramyxovirus"),]
v <- v[-which(v$virus == "geelong paramyxovirus"),]
v <- v[-which(v$virus == "grove virus"),]
v$virus[which(v$virus == "hendra virus")] <- "Hendra_henipavirus"
v$virus[which(v$virus == "hervey virus")] <- "Hervey_pararubulavirus"
v$virus[which(v$virus == "human mumps virus")] <- "Mumps_orthorubulavirus"
v$virus[which(v$virus == "menangle virus")] <- "Menangle_pararubulavirus"
v$virus[which(v$virus == "nipah virus")] <- "Nipah_henipavirus"
v$virus[which(v$virus == "sosuga virus")] <- "Sosuga_pararubulavirus"
v$virus[which(v$virus == "teviot virus")] <- "Teviot_pararubulavirus"
v$virus[which(v$virus == "tioman virus")] <- "Tioman_pararubulavirus"
v$virus[which(v$virus == "tuhoko virus 1")] <- "Tuhoko_pararubulavirus_1"
v$virus[which(v$virus == "tuhoko virus 2")] <- "Tuhoko_pararubulavirus_2"
v$virus[which(v$virus == "tuhoko virus 3")] <- "Tuhoko_pararubulavirus_3"
v <- v[-which(v$virus == "unclassified EPMV"),]
v <- v[-which(v$virus == "yarra bend paramyxovirus"),]
v <- v[-which(v$virus == "yeppoon virus"),]
v <- v[!is.na(v$virus),]

## load in viral phylogeny and taxonomy
setwd("~/Documents/PhD/pteropodidae/phylo")
vtree=read.tree('ONSR.Paramyxo.Fig3.v6.tree.nwk')

vtree$tip.label=sapply(strsplit(vtree$tip.label,'_'),function(x) paste(x[2],x[3],x[4],sep='_'))
vtree$tip.label <- str_replace(vtree$tip.label,"_NA", "")

## all species
vdata=data.frame(virus=vtree$tip.label)

## merge
vdata=merge(vdata,v,by="virus",all=T)

## no studies
vdata$Studies=ifelse(is.na(vdata$Studies),0,vdata$Studies)

## no bats
vdata$nSamp=ifelse(is.na(vdata$nSamp),0,vdata$nSamp)

## binary
vdata$binstudy=ifelse(vdata$Studies==0,0,1)

vdata$tip=vdata$virus
vdata$taxonomy=paste0("PARAMYXOVIRIDAE; ",vdata$virus)

vtree$node.label<-NULL
vtree <- di2multi(vtree)
vtree <- root(vtree, "009094051.1_Sunshine_Coast", resolve.root = TRUE) ##root tree

vdata=comparative.data(phy=vtree,data=vdata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)

vdata$data$label=vdata$data$virus
vdata$data$Species=vdata$data$virus
vdata$data$tip=vdata$data$virus

## just sampled bats
svdata=vdata[vdata$data$binstudy==1,]

## D statistic on sampled/not sampled
set.seed(1)
vdata$phy <- di2multi(vdata$phy)
dstat=phylo.d(data=vdata,binvar=binstudy,permut=1000)
dstat ###Brownian phylogenetic structure

## number of studies and number of bats
hist(log10(svdata$data$Studies))
hist(log10(svdata$data$nSamp))

## transform
svdata$data$lstudies=log10(svdata$data$Studies)
svdata$data$ltested=log10(svdata$data$nSamp)

## range
range(svdata$data$Studies)
range(svdata$data$nSamp)

## pagel's lambda
pmod1=pgls(lstudies~1,data=svdata,lambda="ML") ## lambda = 0.903; departs from randomness but not BM, suggesting phylogenetic clustering
pmod2=pgls(ltested~1,data=svdata,lambda="ML") ## lambda = 1.000; departs from randomness but not BM, suggesting phylogenetic clustering

## summarize
summary(pmod1)
summary(pmod2)

## GPF for study binary
set.seed(1)
study_pf_v=gpf(Data=vdata$data,tree=vdata$phy,
               frmla.phylo=binstudy~phylo,
               family=binomial,
               algorithm='phylo',nfactors=2,
               min.group.size = 3)

## summarize
HolmProcedure(study_pf_v)

## number of studies
set.seed(1)
nstudies_pf_v=gpf(Data=svdata$data,tree=svdata$phy,
                  frmla.phylo=Studies~phylo,
                  family=poisson,
                  algorithm='phylo',nfactors=5,
                  min.group.size = 3)

## summarize
HolmProcedure(nstudies_pf_v)
nstudies_res_v=pfsum(nstudies_pf_v, letter="E")$results

## number of samples
set.seed(1)
svdata$data$snSamp=sqrt(svdata$data$nSamp)
svdata$data$lnSamp=log1p(svdata$data$nSamp)

nsamples_pf_v=gpf(Data=svdata$data,tree=svdata$phy,
                  frmla.phylo=lnSamp~phylo,
                  family=gaussian,
                  algorithm='phylo',nfactors=10,
                  min.group.size = 3)

## summarize
HolmProcedure(nsamples_pf_v)
nsamples_res_v=pfsum(nsamples_pf_v, letter="F")$results

## lower/greater
nstudies_res_v$check=ifelse(nstudies_res_v$clade>nstudies_res_v$other,"more","less")
nsamples_res_v$check=ifelse(nsamples_res_v$clade>nsamples_res_v$other,"more","less")

###figures
knit_hooks$set(crop = hook_pdfcrop)
opts_chunk$set(crop = TRUE)

## save trees
vtree=treeio::full_join(as.treedata(vdata$phy),vdata$data,by="label")
svtree=treeio::full_join(as.treedata(svdata$phy),svdata$data,by="label")

##########FIGURE

## set x max
plus=1
pplus=plus+0.75

## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## make low and high
pcols=afun(2)

## function to loop and add clades
cadd=function(gg,pf,pmax,label="yes",let="A"){
  
  ## ifelse
  if(HolmProcedure(pf)==0){
    gg=gg
  }else{
    
    ## make result
    result=pfsum(pf, letter=let)$results
    
    ## ifelse 
    if(nrow(result)>pmax){
      result=result[1:pmax,]
    }else{
      result=result
    }
    
    ## set tree
    for(i in 1:nrow(result)){
      
      ## highlight clade
      gg=gg+
        geom_hilight(node=result$node[i],
                     alpha=0.25,
                     fill=ifelse(result$clade>
                                   result$other,pcols[2],pcols[1])[i])+
        
        ## add label
        if(label=="yes"){
          geom_cladelabel(node = result$node[i], 
                          label = result$factor[i], 
                          offset = 10,
                          offset.text = 0.2,
                          fontsize=4)
        }    
      
    }
  }
  return(gg)
}

## state pmax
pmax=10

###host trees

## make base
base=ggtree(dtree,size=0.05)
base2=ggtree(stree,size=0.2,branch.length='none',layout="circular")
base3=ggtree(stree,size=0.1)

## binary tree (first fig)
gg=cadd(base,study_pf,pmax,label="no")
dtree@extraInfo$binstudy <- (as.factor(dtree@extraInfo$binstudy))
binstudy <- dtree@extraInfo$binstudy
plot1 = gg + geom_tippoint(aes(color=dtree@extraInfo$binstudy), size=0.3) +
  scale_color_manual(values=c("white","black")) +
  theme(legend.position = "none") + ggtitle("(A) binary studied bat species") +
  geom_cladelabel(node = 328, 
                  label = "A1", 
                  offset = 0.5, 
                  offset.text = 0.5,
                  fontsize=4)

#number of studies (second fig)
plot2=cadd(base3,nstudies_pf,pmax,label="no")

## get tree data
tdata=base3$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend=scales::rescale(tdata$Studies,c(max(tdata$x),xmax)),
                species=tdata$Species)

## fix gg
plot2=plot2+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.25,alpha=0.5) +
  ggtitle("(B) studies per bat species") +
    geom_cladelabel(node = nstudies_res$node[1], 
                    label = nstudies_res$factor[1], 
                    offset = 10,
                    offset.text = 0.2,
                    vjust = -1,
                    fontsize=4) +
    geom_cladelabel(node = nstudies_res$node[2], 
                  label = nstudies_res$factor[2], 
                  offset = 10,
                  offset.text = 0.2,
                  fontsize=4) +
  geom_cladelabel(node = nstudies_res$node[3], 
                  label = nstudies_res$factor[3], 
                  offset = 10,
                  offset.text = 0.2,
                  vjust = 1,
                  fontsize=4) +
  geom_cladelabel(node = nstudies_res$node[4], 
                  label = nstudies_res$factor[4], 
                  offset = 10,
                  offset.text = 0.2,
                  fontsize=4)
  

plot3=cadd(base3,nsamples_pf,pmax,let="C")

## get tree data
tdata=base3$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend=scales::rescale(log10(tdata$nSamp),c(max(tdata$x),xmax)),
                species=tdata$Species)

## fix gg
plot3=plot3+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.25,alpha=0.5)
plot3=plot3+ggtitle("(C) log(samples) per bat species")

####viral trees

## make base
vbase=ggtree(vtree,size=0.05, branch.length='none')
vbase2=ggtree(svtree,size=0.4,branch.length='none',layout="circular")
vbase3=ggtree(svtree,size=0.1, branch.length='none')

## binary tree (first fig)
gg=cadd(vbase,study_pf_v,pmax)
vtree@extraInfo$binstudy <- (as.factor(vtree@extraInfo$binstudy))
binstudy <- vtree@extraInfo$binstudy
vplot1 = gg + geom_tippoint(aes(color=vtree@extraInfo$binstudy), size=0.3) +
  scale_color_manual(values=c("white","black")) +
  theme(legend.position = "none") + ggtitle("(D) binary studied viral species")

#number of studies (second fig)
vplot2=cadd(vbase3,nstudies_pf_v,pmax,let="E")

## get tree data
tdata=vbase3$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend=scales::rescale(tdata$Studies,c(max(tdata$x),xmax)),
                species=tdata$Species)

## fix gg
vplot2=vplot2+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.4,alpha=0.5) +
  ggtitle("(E) studies per viral species")

## number of samples (third fig)
vplot3=cadd(vbase3,nsamples_pf_v,pmax,let="F")

## get tree data
tdata=vbase3$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend=scales::rescale(tdata$Studies,c(max(tdata$x),xmax)),
                species=tdata$Species)

## fix gg
vplot3=vplot3+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.4,alpha=0.5) +
  ggtitle("(F) log(samples) per viral species")

#export
setwd("~/Documents/PhD/pteropodidae")
png("Figure 3.png",width=10,height=12,units="in",res=600)
(plot1|(plot2/plot3))/(vplot1|(vplot2/vplot3))+plot_layout(widths=c(2,1))
dev.off()
