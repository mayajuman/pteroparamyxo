####rma.mv paramyxo models

library(metafor)
library(MetBrewer)
library(scales)
library(dplyr)
library(ape)
library(plyr)
library(Hmisc)
library(tidyr)
library(ggrepel)
library(pals)

##set wd
setwd("~/Documents/PhD/pteropodidae")
#load dataset
rev <- read.csv("pteroparamyxo.csv")

#calculate missing prevalence
rev$Prevalence <- ifelse(is.na(rev$Prevalence), rev$nPositive/rev$nSamples, rev$Prevalence)

###data cleaning
rev$tissue_simple <- rev$Tissue.tested.
rev$tissue_simple <- revalue(rev$tissue_simple,
                             c("alimentary, respiratory"="pooled swabs/samples",
                               "anus"="fecal, rectal, or anal",
                               "body"="skin swab",
                               "body, head, oral, rectal"="pooled swabs/samples",
                               "brain, heart, lung, liver, spleen, intestine, kidney, and reproductive organ homogenate"="pooled tissue",
                               "brain, salivary gland, spleen"="pooled tissue",
                               "fecal"="pooled swabs/samples",
                               "fecal, intestine, rectum, spleen"="pooled swabs/samples",
                               "fecal, intestine, spleen"="pooled swabs/samples",
                               "fecal, kidney"="pooled swabs/samples",
                               "fecal, spleen"="pooled swabs/samples",
                               "fecal, urine"="pooled swabs/samples",
                               "feces"="fecal, rectal, or anal",
                               "feces, rectal"="fecal, rectal, or anal",
                               "head"="skin swab",
                               "intestine, spleen"="pooled tissue",
                               "kidney, liver, lung"="pooled tissue",
                               "kidney, liver, spleen"="pooled tissue",
                               "kidney, lung, spleen"="pooled tissue",
                               "kidney, spleen"="pooled tissue",
                               "liver, kidney, lung, spleen, reproductive organs"="pooled tissue",
                               "liver, spleen"="pooled tissue",
                               "oral, rectal"="pooled swabs/samples",
                               "organs, feces, oropharyngeal"="pooled swabs/samples",
                               "oropharyngeal"="oral",
                               "pharynx/anus"="pooled swabs/samples",
                               "placenta"="reproductive",
                               "pooled fetal tisses (brain, lung, heart, kidney and/or spleen)"="pooled tissue",
                               "rectal"="fecal, rectal, or anal",
                               "rectum"="fecal, rectal, or anal",
                               "reproductive tract"="reproductive",
                               "saliva"="oral",
                               "salivary gland"="oral",
                               "small intestine"="intestine",
                               "testes, uterus"="reproductive",
                               "throat"="oral",
                               "throat, rectal"="pooled swabs/samples",
                               "tonsil"="oral",
                               "trachea"="oral",
                               "urine"="urinary",
                               "urine, throat"="pooled swabs/samples",
                               "urine, urogenital"="urinary",
                               "urine, urogenital/prepucial"="urinary",
                               "urogenital"="urinary",
                               "vaginal/preputial"="reproductive",
                               "various"="unspecified"))

rev[which(rev$tissue_simple == "urinary" & rev$Detection.type == "pooled"),]$tissue_simple <- "urinary (pooled)"
rev[which(rev$tissue_simple == "urinary" & rev$Detection.type == "single"),]$tissue_simple <- "urinary (single)"
rev[which(rev$tissue_simple == "urinary" & rev$Detection.type == "single, pooled"),]$tissue_simple <- "urinary (unspecified)"
rev[which(rev$Tissue.tested. == "feces" & rev$Detection.type == "pooled"),]$tissue_simple <- "feces (pooled)"
rev[which(rev$Tissue.tested. == "feces" & rev$Detection.type == "single"),]$tissue_simple <- "feces (single)"

prev <- rev[which(!is.na(rev$Prevalence)),]
prev$Host.genus <- as.factor(prev$Host.genus)
prev$Host.species <- as.factor(prev$Host.species)
prev$Host.subfamily <- as.factor(prev$Host.subfamily)
prev$Virus.subfamily <- as.factor(prev$Virus.subfamily)
prev$Virus.genus <- as.factor(prev$Virus.genus)
prev$tissue_simple <- as.factor(prev$tissue_simple)
prev$Detection.type <- as.factor(prev$Detection.type)
prev$Study.type <- as.factor(prev$Study.type)
prev$Detection.method <- as.factor(prev$Detection.method)
prev$Detection.method.specific <- as.factor(prev$Detection.method.specific)
prev$Gene.Antigen.Target <- as.factor(prev$Gene.Antigen.Target)
prev$Sample.type <- as.factor(prev$Sample.type)
prev$Country <- as.factor(prev$Country)
####averaging year range
prev$SamplingYear <- ifelse(is.na(prev$SamplingYear), (prev$StartYear+prev$EndYear)/2, prev$SamplingYear)
prev$SamplingYear <- as.vector(scale(as.numeric(prev$SamplingYear), center=TRUE, scale=TRUE)) #scaling and centering year

prev$region <- prev$Country

prev$region <- revalue(prev$region, c("australia" = "oceania",
                                      "bangladesh" = "asia",
                                      "cambodia" = "asia",
                                      "cameroon" = "africa",
                                      "china" = "asia",
                                      "democratic republic of the congo" = "africa",
                                      "democratic republic of the congo, central african republic, gabon, republic of the congo, ghana" = "africa",
                                      "east timor" = "asia",
                                      "equatorial guinea" = "africa",
                                      "fiji" = "oceania",
                                      "gabon" = "africa",
                                      "gabon, republic of the congo" = "africa",
                                      "gabon, republic of the congo, democratic republic of the congo" = "africa",
                                      "gabon, republic of the congo, ghana, central african republic" = "africa",
                                      "ghana" = "africa",
                                      "ghana, central african republic, gabon, democratic republic of the congo" = "africa",
                                      "ghana, tanzania, uganda, malawi, zambia, sao tome and principe, equatorial guinea" = "africa",
                                      "india" = "asia",
                                      "indonesia" = "asia",
                                      "kenya" = "africa",
                                      "madagascar" = "africa",
                                      "malawi" = "africa",
                                      "malaysia" = "asia",
                                      "mayotte" = "africa",
                                      "myanmar" = "asia",
                                      "nigeria" = "africa",
                                      "papua new guinea" = "oceania",
                                      "republic of the congo" = "africa",
                                      "rwanda" = "africa",
                                      "sao tome and principe" = "africa",
                                      "saudi arabia" = "asia",
                                      "singapore" = "asia",
                                      "south africa" = "africa",
                                      "tanzania" = "africa",
                                      "thailand" = "asia",
                                      "the comoros" = "africa",
                                      "uganda" = "africa",
                                      "vietnam" = "asia",
                                      "zambia" = "africa"))

prev$Sample.type <- revalue(prev$Sample.type, c("tissue, swab" = "various",
                                                "tissue, feces" = "various",
                                                "tissue, feces, swab" = "various",
                                                "urine, swab" = "various",
                                                "urine, feces" = "various",
                                                "feces, swab" = "various",
                                                "specimen" = "various",
                                                "packed haemocytes" = "serum/haemocytes",
                                                "serum" = "serum/haemocytes"))
prev$Detection.type <- factor(prev$Detection.type, levels=c("single","pooled","single, pooled"))
prev$tissue_simple <- revalue(prev$tissue_simple, c("pooled tissue" = "pooled samples",
                                                    "pooled swabs/samples" = "pooled samples"))
prev$Gene.Antigen.Target <- revalue(prev$Gene.Antigen.Target, c("L (AR)" = "L",
                                                                "L (PAR and RMH)" = "L",
                                                                "L (PAR)" = "L",
                                                                "L (RMH)" = "L",
                                                                "N, M" = "various",
                                                                "pol" = "L",
                                                                "CedV sF" = "F",
                                                                "CedV sG" = "G",
                                                                "HeV sF" = "F",
                                                                "HeV sG" = "G",
                                                                "MenV sHN" = "HN",
                                                                "NiV sG" = "G"))
prev$Detection.method.specific <- revalue(prev$Detection.method.specific, 
                                          c("serum neutralization test" = "neutralization test",
                                            "virus neutralization test" = "neutralization test"))

prev$Host.species <- revalue(prev$Host.species,
                             c("magna"="moluccensis",
                               "medius"="giganteus",
                               "melanotus natalis"="melanotus",
                               "natalis"="melanotus",
                               "seychellensis comorensis"="seychellensis",
                               "labiatus minor"="labiatus"))

prev$Host.genus <- revalue(prev$Host.genus,
                           c("lissonycteris"="myonycteris"))

prev <- prev[!grepl(",|\\.|/", prev$Host.genus),] #exclude rows with multiple genera
prev <- prev[!is.na(prev$Host.genus),] #exclude rows with missing genus
prev <- prev[!grepl(",|\\.|/", prev$Detection.method.specific),] #exclude rows with multiple detection methods

prev$host <- with(prev, paste0(Host.genus," ",Host.species)) #create host name variable
prev$host <- gsub(" ","_",capitalize(prev$host))
prev$host <- as.factor(prev$host)

## load in bat phylogeny and taxonomy
setwd("~/Documents/PhD/pteropodidae/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)

PCRprev <- prev %>% filter(Detection.method == "PCR")
seroprev <- prev %>% filter(Detection.method == "serology")

## trim to family
taxa=taxa[taxa$fam=="PTEROPODIDAE",]
taxa$tip=taxa$gen #genus level phylo and random effects

## trim tree to genus level phylogeny
tree=keep.tip(tree,taxa$tiplabel)
tree$tip.label=tolower(sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1])))
taxa$host=taxa$gen

###genus level
stree1a=keep.tip(tree,as.character(unique(PCRprev$Host.genus)))
stree4d=keep.tip(tree,as.character(unique(seroprev$Host.genus)))

## convert tree to correlation matrix
#all PMVs PCR
cmatrix1a=vcv.phylo(stree1a,cor=T)

#all PMVs sero
cmatrix4d=vcv.phylo(stree4d,cor=T)

## make observation and study-level random effect
#all PMV PCR
PCRprev$observation=factor(1:nrow(PCRprev))
PCRprev$study=factor(PCRprev$Title)

#all PMV sero
seroprev$observation=factor(1:nrow(seroprev))
seroprev$study=factor(seroprev$Title)

## pft in escalc for yi and vi 
#all PMV PCR
PCRprev=data.frame(PCRprev,escalc(xi=PCRprev$nPositive,ni=PCRprev$nSamples,measure="PFT"))
#all PMV sero
seroprev=data.frame(seroprev,escalc(xi=seroprev$nPositive,ni=seroprev$nSamples,measure="PFT"))

## back transform
#all PMV PCR
PCRprev$backtrans=transf.ipft(PCRprev$yi,PCRprev$nSamples)
#all PMV sero
seroprev$backtrans=transf.ipft(seroprev$yi,seroprev$nSamples)

## species and (genus-level) phylo effect
#all PMV PCR
PCRprev$phylo=PCRprev$Host.genus
PCRprev$host=PCRprev$phylo

#all PMV sero
seroprev$phylo=seroprev$Host.genus
seroprev$host=seroprev$phylo

####PCR model
mod1=rma.mv(yi=yi,V=vi,
            random=list(~1|study/observation,~1|host,~1|phylo,~1|Country),
            R=list(phylo=cmatrix1a),method="REML",
            mods=~region + Study.type + Sample.type + Host.subfamily + SamplingYear, data=PCRprev,
            control=list(optimizer="optim", optmethod="BFGS"))

###serological model
mod2=rma.mv(yi=yi,V=vi,
            random=list(~1|study/observation,~1|host,~1|phylo,~1|Country),
            R=list(phylo=cmatrix4d),method="REML",
            mods=~Host.subfamily + region + SamplingYear, data=seroprev,
            control=list(optimizer="optim", optmethod="BFGS"))

#figures

global_betas_PCR=data.frame(beta=mod1$beta, ci.lb=mod1$ci.lb, ci.ub=mod1$ci.ub, pval=mod1$pval)
global_betas_sero=data.frame(beta=mod2$beta, ci.lb=mod2$ci.lb, ci.ub=mod2$ci.ub, pval=mod2$pval)

global_betas_PCR$coefficient <- rownames(global_betas_PCR)
global_betas_sero$coefficient <- rownames(global_betas_sero)

global_betas_PCR$type <- "PCR"
global_betas_sero$type <- "Serology"
global_betas <- rbind(global_betas_PCR,global_betas_sero)

global_betas$type <- as.factor(global_betas$type)
global_betas$type <- factor(global_betas$type, levels=c("PCR","Serology"))

global_betas$coefficient[which(global_betas$coefficient=="intrcpt")] <- "intercept"
global_betas$coefficient[which(global_betas$coefficient=="Study.typecross-sectional, longitudinal")] <- "both cross-sectional and longitudinal"
global_betas$coefficient[which(global_betas$coefficient=="Study.typelongitudinal")] <- "longitudinal"
global_betas$coefficient[which(global_betas$coefficient=="Sample.typevarious")] <- "various samples"
global_betas$coefficient[which(global_betas$coefficient=="Sample.typeswab")] <- "swab"
global_betas$coefficient[which(global_betas$coefficient=="Sample.typetissue")] <- "tissue"
global_betas$coefficient[which(global_betas$coefficient=="Sample.typeurine")] <- "urine"
global_betas$coefficient[which(global_betas$coefficient=="Sample.typefeces")] <- "feces"
global_betas$coefficient[which(global_betas$coefficient=="Sample.typesaliva")] <- "saliva"
global_betas$coefficient[which(global_betas$coefficient=="Sample.typeserum/haemocytes")] <- "serum"
global_betas$coefficient[which(global_betas$coefficient=="Host.subfamilyeidolinae")] <- "Eidolinae"
global_betas$coefficient[which(global_betas$coefficient=="Host.subfamilyepomophorinae")] <- "Epomophorinae"
global_betas$coefficient[which(global_betas$coefficient=="Host.subfamilypteropodinae")] <- "Pteropodinae"
global_betas$coefficient[which(global_betas$coefficient=="Host.subfamilyrousettinae")] <- "Rousettinae"
global_betas$coefficient[which(global_betas$coefficient=="Host.subfamilymacroglossusinae")] <- "Macroglossusinae"
global_betas$coefficient[which(global_betas$coefficient=="SamplingYear")] <- "year"
global_betas$coefficient[which(global_betas$coefficient=="regionasia")] <- "Asia"
global_betas$coefficient[which(global_betas$coefficient=="regionafrica")] <- "Africa"

global_betas$coefficient <- as.factor(global_betas$coefficient)
global_betas$coefficient <- factor(global_betas$coefficient, 
                                   levels=rev(c("intercept",
                                                "both cross-sectional and longitudinal",
                                                "longitudinal",
                                                "feces",
                                                "various samples",
                                                "serum",
                                                "saliva",
                                                "swab",
                                                "tissue",
                                                "urine",
                                                "Asia",
                                                "Africa",
                                                "Eidolinae",
                                                "Epomophorinae",
                                                "Macroglossusinae",
                                                "Pteropodinae",
                                                "Rousettinae",
                                                "year")))

global_betas$Variable <- as.factor(c("intercept",
                                     "region",
                                     "region",
                                     "study type",
                                     "study type",
                                     "sample type",
                                     "sample type",
                                     "sample type",
                                     "sample type",
                                     "sample type",
                                     "sample type",
                                     "sample type",
                                     "host subfamily",
                                     "host subfamily",
                                     "host subfamily",
                                     "host subfamily",
                                     "host subfamily",
                                     "year",
                                     "intercept",
                                     "host subfamily",
                                     "host subfamily",
                                     "host subfamily",
                                     "host subfamily",
                                     "host subfamily",
                                     "region",
                                     "region",
                                     "year"))

global_betas$Variable <- factor(global_betas$Variable, levels=c("intercept",
                                                                "study type",
                                                                "sample type",
                                                                "region",
                                                                "host subfamily",
                                                                "year"))

list <- c()
for(i in 1:27){
  if(global_betas$ci.lb[i] < 0 & global_betas$ci.ub[i] < 0){
    new_element <- "no"
    list <- c(list, new_element)
  }
  else if(global_betas$ci.lb[i] > 0 & global_betas$ci.ub[i] > 0){
    new_element <- "no"
    list <- c(list, new_element)
  }
  else{
    new_element <- "yes"
    list <- c(list, new_element)
  }
}
global_betas$cicrosseszero <- list
global_betas$cicrosseszero <- as.factor(global_betas$cicrosseszero)

library(MetBrewer)
library("scales")
colors=palette.colors(palette="R4",n=6)

fig4 <- ggplot(global_betas) + geom_hline(yintercept=0,linetype="dashed",colour="gray") + 
  geom_errorbar(aes(x=coefficient,ymin=ci.lb,ymax=ci.ub,color=Variable,alpha=cicrosseszero),width=0,size=1.5) +
  scale_alpha_manual(values=c(1, 0.25),guide="none") +
  geom_point(aes(x=coefficient,y=beta,color=Variable),size=3.5) + coord_flip() + 
  facet_grid(.~type, scales = "free_x") + theme_bw(base_size = 14) + labs(x=NULL) + 
  ylab("Meta-analysis model coefficient and 95% confidence interval") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  scale_color_manual(values=colors)

setwd("~/Documents/PhD/pteropodidae")
ggsave("Fig4.jpg", width = 8, height = 8)

####heterogeneity in prevalence

## function for I2 for rma.mv
i2=function(model){
  
  ## metafor site code for I2
  W=diag(1/model$vi)
  X=model.matrix(model)
  P=W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2=100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  I2=round(I2,2)
  
  ## summarize by each variance component
  allI2=100 * model$sigma2 / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  allI2=round(allI2,3)
  return(list(I2=I2,allI2=allI2))
}

####intercept-only models

mod1i=rma.mv(yi=yi,V=vi,
            random=list(~1|study/observation,~1|host,~1|phylo,~1|Country),
            R=list(phylo=cmatrix1a),method="REML",
            #mods=~region + Study.type + Sample.type, 
            data=PCRprev,
            control=list(optimizer="optim", optmethod="BFGS"))

###sero all
mod2i=rma.mv(yi=yi,V=vi,
            random=list(~1|study/observation,~1|host,~1|phylo,~1|Country),
            R=list(phylo=cmatrix4d),method="REML",
            #mods=~Host.subfamily + region + SamplingYear, 
            data=seroprev,
            control=list(optimizer="optim", optmethod="BFGS"))

#~1|study/observation,~1|host,~1|phylo,~1|Country
i2(mod1i)
i2(mod2i)


######tissue tropism

sampletable <- rev %>% group_by(tissue_simple, Detection.method) %>% 
  dplyr::summarise(n=n(), nsamp=sum(nSamples, na.rm=TRUE))

sampletable_percent <- rev %>% group_by(tissue_simple, Detection.method) %>% tidyr::drop_na(nPositive) %>%
  dplyr::summarise(percentpos=round(sum(nPositive>0)/n()*100,2),avg=round(mean(Prevalence, na.rm=TRUE)*100,2))

PMVtable <- merge(sampletable, sampletable_percent)

ggplot(PMVtable %>% filter(Detection.method == "PCR", tissue_simple != "unspecified"), 
       aes(x=log10(nsamp), y=avg, label=tissue_simple)) + 
  geom_point(aes(color=tissue_simple)) + 
  xlab("Logged number of samples") +
  ylab("Average prevalence (%)") +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank()) + xlim(c(1.5,4.7)) +
  geom_text_repel(data = PMVtable %>% filter(Detection.method == "PCR", tissue_simple != "unspecified"), 
                  aes(label = tissue_simple, color=tissue_simple)) +
  scale_colour_manual(values=glasbey(17))

ggsave("Fig5.jpg", width=7, height=5)