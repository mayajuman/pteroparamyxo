#####spatial bias

library(tidyr)
library(maptools)
library(car)
library(performance)
library(visreg)
library(emmeans)
library(ggpubr)
library(letsR)
library(raster)

#load data
setwd("~/Documents/PhD/pteropodidae")
rev <- read.csv("pteroparamyxo.csv")

#calculate missing prevalence
rev$Prevalence <- ifelse(is.na(rev$Prevalence), rev$nPositive/rev$nSamples, rev$Prevalence)

#drop rows with >1 country and then calculate samples using summarize() and sum()
countries <- rev %>% dplyr::select(Title, nSamples, Country, Prevalence)
uniq_countries<- rev[!grepl(',', rev$Country),]
uniq_countries$Country <- str_to_title(uniq_countries$Country)
unique(uniq_countries$Country)

uniq_countries$country <- revalue(uniq_countries$Country, 
                                  c("Sao Tome And Principe"="Sao Tome and Principe",
                                    "Republic Of The Congo"="Republic of Congo",
                                    "The Comoros"="Comoros",
                                    "Democratic Republic Of The Congo"="Democratic Republic of the Congo",
                                    "East Timor"="Timor-Leste",
                                    "Republic of the Congo"="Congo"))

percountry <- uniq_countries %>% 
  group_by(country) %>% 
  dplyr::summarise(nSamp = sum(nSamples, na.rm=TRUE),
                   PCRPrev = mean(Prevalence[Detection.method=="PCR"], na.rm=TRUE)*100,
                   seroprev = mean(Prevalence[Detection.method=="serology"], na.rm=TRUE)*100)

#create second df with study name and country, 
#turn >1 country rows into one row per country, 
#then calculate # of studies per country using summarize and ndistinct()

studies <- countries %>% dplyr::select(Title, Country) %>%
  mutate(country = strsplit(as.character(Country), ",")) %>%
  unnest(country) %>% 
  mutate(country = str_trim(as.character(country))) %>%
  group_by(country) %>% 
  dplyr::summarise(Studies = n_distinct(Title))

studies$country <- str_to_title(studies$country)
studies$country <- revalue(studies$country, 
                           c("Sao Tome And Principe"="Sao Tome and Principe",
                             "Republic Of The Congo"="Republic of Congo",
                             "The Comoros"="Comoros",
                             "Democratic Republic Of The Congo"="Democratic Republic of the Congo",
                             "East Timor"="Timor-Leste"))

percountry <- merge(percountry, studies)
percountry <- percountry %>% drop_na(country)

#loading global map
wdata=map_data("world")
wdata=wdata[-which(wdata$region=='Antarctica'),]
wdata$country=wdata$region

setdiff(unique(percountry$country),unique(wdata$country))

#load georegion
geo=read.csv("georegion.csv")
geo$country=geo$name
geo=geo[c("country","region","sub.region")]
setdiff(percountry$country,geo$country)

geo$country=revalue(geo$country,
                    c("Congo"="Republic of Congo",
                      "Tanzania, United Republic of"="Tanzania",
                      "Congo, Democratic Republic of the"="Democratic Republic of the Congo",
                      "Viet Nam"="Vietnam"))

setdiff(percountry$country,geo$country)

#merge
gdata=merge(geo,percountry,by="country",all.x=T)

#binary sampling
gdata$binstudy=ifelse(is.na(gdata$Studies),0,1)

#clean
gdata=gdata[!gdata$region=="",]
gdata <- gdata %>% filter(region %in% c("Africa","Asia","Oceania"))

## GLM for binary sampling
mod1=glm(binstudy~region,data=gdata,family=binomial)

## within sampled GLMs
gdata2=gdata[gdata$binstudy==1,]
mod2=glm(Studies~region,data=gdata2,family=poisson)
mod3=glm(nSamp~region,data=gdata2,family=poisson)

#range
range(gdata2$Studies)
range(gdata2$nSamp)

## Anova
Anova(mod1) #insig for binary effort
Anova(mod2) #vsig for # of studies
Anova(mod3) #vsig for # of samples; prob australia

## R2
r2_mcfadden(mod1)
r2_mcfadden(mod2)
r2_mcfadden(mod3)

## visreg
visreg(mod1,"region",scale="response")
visreg(mod2,"region",scale="response")
visreg(mod3,"region",scale="response") #australia

## posthoc
s1=emmeans(mod1,list(pairwise~region),level=0.95,adjust="bonferroni",type="response") %>% confint()
s2=emmeans(mod2,list(pairwise~region),level=0.95,adjust="bonferroni",type="response") %>% confint()
s3=emmeans(mod3,list(pairwise~region),level=0.95,adjust="bonferroni",type="response") %>% confint()

#merge w wdata
cdata=dplyr::left_join(wdata,percountry,by="country",copy=T)

##load IUCN pteropodidae ranges
setwd("~/Documents/PhD/pteropodidae/databases/IUCN")
IUCN <- vect('data_0.shp')
maps <- letsR::lets.presab(IUCN)

#aggregate so you're only looking at presence/absence
maps[["Richness_Raster"]][which(values(maps[["Richness_Raster"]]) > 0)] <- 1
spatraster <- maps[["Richness_Raster"]]
raster <- as(spatraster, "Raster")

## Extract polygons
pp <- rasterToPolygons(raster, dissolve=TRUE)

## Convert SpatialPolygons to a format usable by ggplot2
outline <- fortify(pp)
outline <- outline[-which(outline$long < -15),]
outline <- outline[-which(outline$long > 180),]
outline <- outline[-which(outline$lat < -39),]
outline <- outline[-which(outline$lat > 37),]
outline <- outline[-which(outline$group == 1.1),]

## visualize studies
studiesplot <- ggplot(cdata,aes(long.x,lat.x))+
  scale_fill_viridis_c(na.value="grey90",option="rocket", direction=-1)+
  geom_polygon(aes(group=group,fill=Studies),size=0.1,colour='white')+
  theme_void()+
  coord_sf(xlim=c(-15,185), ylim=c(-45,50))+
  theme(legend.position="bottom")+
  guides(fill=guide_colorbar(title="Studies",
                             barwidth = 15)) +
  geom_path(aes(x = long, y = lat, group = group), data = outline, size=0.2, col="blue")

## visualize samples
samplesplot <- ggplot(cdata,aes(long.x,lat.x))+
  geom_polygon(aes(group=group,fill=log1p(nSamp)),size=0.1,colour='white')+
  theme_void()+
  scale_fill_viridis_c(na.value="grey90",option="rocket", direction=-1) +
  coord_sf(xlim=c(-15,185), ylim=c(-45,50))+
  theme(legend.position="bottom")+
  guides(fill=guide_colorbar(title="log(samples + 1)"),
         barwidth = 15) +
  geom_path(aes(x = long, y = lat, group = group), data = outline, size=0.2, col="blue")

fig2 <- ggarrange(studiesplot, samplesplot,
                  labels = c("A", "B"),
                  ncol = 1, nrow = 2)

ggsave("Fig2.jpg", width = 8, height = 8)


####timeline figure

years <- rev %>% dplyr::select(Title, Year) %>% distinct()
samples <- rev %>% dplyr::select(Year, nSamples) %>% 
  group_by(Year) %>% summarise(nSamples = sum(nSamples, na.rm=TRUE))

yearplot <- years %>% 
  ggplot(aes(x = Year)) +
  geom_bar() +
  labs(y = "Number of studies", x = "Year of publication") +
  theme_classic2() +
  geom_vline(xintercept=2000, colour="red", linetype="longdash", size=0.3) +
  geom_vline(xintercept=2002, colour="blue", linetype="longdash", size=0.3) +
  scale_y_continuous(breaks = round(seq(0, 12, by = 3),2))

sampleyearplot <- samples %>% 
  ggplot(aes(x = Year, y = nSamples)) +
  geom_bar(stat="identity") +
  labs(y = "Number of samples", x = "Year of publication") +
  theme_classic2() +
  geom_vline(xintercept=2000, colour="red", linetype="longdash", size=0.3) +
  geom_vline(xintercept=2002, colour="blue", linetype="longdash", size=0.3)

fig1 <- ggarrange(yearplot, sampleyearplot,
                  labels = c("A", "B"),
                  ncol = 2, nrow = 1)

ggsave("Fig1.jpg", width = 8, height = 4)
