
library(maps)
library(mapproj)
library(mapdata)
library(maptools)

setwd("C:/Users/Francy/Dropbox/MCC_CMMID_COVID19/Programs/Analysis/Paper_V2")
load("C:/Users/Francy/Dropbox/MCC_CMMID_COVID19/Programs/Analysis/Build_analysis_dataset/mcc_covid_data_analysis.RData")

countries<-as.character(unique(data$countryname))

table(data$countryname)

table(data$kgclzone)

table(data$emisphere)

t.test(data$r.mean~data$emisphere)

summary(data$r.mean)
sd(data$r.mean)
hist(data$r.mean)
table(data$r.mean)
data[data$r.mean<1,]

# MAPS

# MAP CASES

sum(data$casesw)

summary(data$casesw)

labels <- paste0(c("<100","100-200","200-400","400-800",">=800")," Cases")
data$casescat <- cut(data$casesw,c(-1,100,200,400,800,1000000),labels=labels)
funcol <- colorRampPalette(c("blue","skyblue","yellow","red"))
col <- funcol(length(labels))[data$casescat]
size <- rep(1.3,nrow(data))

pdf("FigureS1.pdf",width=4.75,height=2.75)
layout(1)
map("world",xlim=c(-170,180),ylim=c(-60,80),mar=c(2,1,1,1)+0.1,
    col=grey(0.6),myborder=0)
symbols(data$long,data$lat,circles=size,inches=0.015,fg=NULL,bg=col,add=T)

map.scale(2,-40,ratio=F,cex=0.6,relwidth=0.1)

legend(-125,-20,labels,xjust=0.5,yjust=0.5,pt.cex=0.4,
       pch=21,lwd=0,pt.bg=funcol(length(labels)),bty="n",cex=0.4,
       title=expression(paste("Covid cases",sep="")))

dev.off()

# MAP MEAN TEMPERATURE

summary(data$tmean)
labels <- paste0(c("<5","5-10","10-15","15-20",">20")," C")
avgtmeancat <- cut(data$tmean,c(-10,5,10,15,20,200),labels=labels)
funcol <- colorRampPalette(c("blue","skyblue","yellow","red"))
col <- funcol(length(labels))[avgtmeancat]
size <- rep(1.3,nrow(data))

pdf("FigureS2.pdf",width=4.75,height=2.75)
layout(1)

map("world",xlim=c(-170,180),ylim=c(-60,80),mar=c(2,1,1,1)+0.1,
    col=grey(0.6),myborder=0)
symbols(data$long,data$lat,circles=size,inches=0.015,fg=NULL,bg=col,add=T)

map.scale(2,-40,ratio=F,cex=0.6,relwidth=0.1)

legend(-125,-20,labels,xjust=0.5,yjust=0.5,pt.cex=0.4,
       pch=21,lwd=0,pt.bg=funcol(length(labels)),bty="n",cex=0.5,
       title=expression(paste("T mean (C)",sep="")))

dev.off()

# MAP R

summary(data$r.mean)
labels <- paste0(c("<1.0","1.0-1.2","1.2-1.4","1.4-1.6",">1.6"),"")
rcat <- cut(data$r.mean,c(-10,1.0,1.2,1.4,1.6,200),labels=labels)
funcol <- colorRampPalette(c("blue","skyblue","yellow","red"))
col <- funcol(length(labels))[rcat]
size <- rep(1.3,nrow(data))

pdf("Figure1.pdf",width=4.75,height=2.75)
layout(1)

map("worldHires",xlim=c(-170,180),ylim=c(-60,80),mar=c(2,1,1,1)+0.1,
    col=grey(0.6),myborder=0)
symbols(data$long,data$lat,circles=size,inches=0.015,fg=NULL,bg=col,add=T)

map.scale(2,-40,ratio=F,cex=0.6,relwidth=0.1)

legend(-125,-20,labels,xjust=0.5,yjust=0.5,pt.cex=0.4,
       pch=21,lwd=0,pt.bg=funcol(length(labels)),bty="n",cex=0.5,
       title=expression(paste("R",sep="")))

dev.off()

supplementary.table1<-data[,c("countryname","cityname","casesw","ndays","kgclzone","tmean",
  "pm25","population","density","oldpop","gdp", "oxgov10r","r.mean")]

write.csv(supplementary.table1,file="TableS2.csv")

# DESCRIPTIVE EXPOSURES

weather<-data[,c("tmean","rh","ah","uv","ws","prectot")]

names(weather)<-c("Ta","RH", "AH","UV","WS", "Prec")

res <- cor(weather)
round(res, 2)

#write.csv(res,file="TableS3.csv")

#### Figure S2

library(corrplot)
library(ragg)
agg_png("corrplot1.png", width = 4.75, height = 4.75, units = 'in', res = 300)
corrplot(res, type = "upper", 
  tl.col = "black", tl.srt = 45,cl.cex=0.6)
invisible(dev.off())

desweather<-matrix(nrow=ncol(weather),ncol=5)

for (i in 1:ncol(weather)) {
  
  desweather[i,1]<-colnames(weather)[i]
  desweather[i,2]<-mean(weather[,i])
  desweather[i,3]<-sd(weather[,i])
  desweather[i,4]<-min(weather[,i])
  desweather[i,5]<-max(weather[,i])
  
}

write.csv(desweather,file="Table2b.csv")

# DESCRIPTIVE COVARIATES

covariates<-data[,c("pm25", "population", "density","oldpop", "gdp" , 
                    "oxgov10r")]

names(covariates)<-c("PM25", "Population", "Density","%Pop>65years", "GDP" , "OxCGRT")

res <- cor(covariates, use="pairwise.complete.obs")
round(res, 2)

#write.csv(res,file="TableS2.csv")


#### Figure S3

library(corrplot)
agg_png("corrplot2.png", width = 6, height = 6, units = 'in', res = 300)
corrplot(res, type = "upper", 
  tl.col = "black", tl.srt = 45,cl.cex=0.6)
invisible(dev.off())

descovariates<-matrix(nrow=ncol(covariates),ncol=5)

for (i in 1:ncol(covariates)) {
  
  descovariates[i,1]<-colnames(covariates)[i]
  descovariates[i,2]<-mean(covariates[,i],na.rm=TRUE)
  descovariates[i,3]<-sd(covariates[,i],na.rm=TRUE)
  descovariates[i,4]<-min(covariates[,i],na.rm=TRUE)
  descovariates[i,5]<-max(covariates[,i],na.rm=TRUE)
  
}

write.csv(descovariates,file="Table2a.csv")

# TABLE BY COUNTRY

data$calyearo<-as.Date(data$calyear,origin = "1970-01-01")
data$const<-rep(1,nrow(data))

table(data$calyearo)

ncity<-aggregate(data$const,by=list(data$countryname),FUN=sum)
ncases<-aggregate(data$casesw,by=list(data$countryname),FUN=sum)
period<-aggregate(data$calyearo,by=list(data$countryname),FUN=mean)
reprnumb<-aggregate(data$r.mean,by=list(data$countryname),FUN=mean)
oxgovn<-aggregate(data$oxgov,by=list(data$countryname),FUN=mean)
oxgovb<-aggregate(data$oxgovb,by=list(data$countryname),FUN=mean)
oxgov10<-aggregate(data$oxgov10,by=list(data$countryname),FUN=mean)
oxgov10r<-aggregate(data$oxgov10r,by=list(data$countryname),FUN=mean)

table.country<-data.frame(country=ncity[,1], cites=ncity[,2],covid=ncases[,2], midperiod=period[,2],r=reprnumb[,2],
                          oxgov=oxgov10r[,2])

plot(table.country$oxgov,table.country$r)

write.csv(table.country,file="Table1.csv")

# observation period

periodmin<-aggregate(data$calyearo,by=list(data$countryname),FUN=min)
periodmean<-aggregate(data$calyearo,by=list(data$countryname),FUN=mean)
periodmax<-aggregate(data$calyearo,by=list(data$countryname),FUN=max)

period.country<-data.frame(country=periodmin[,1], start=periodmin[,2],mid=periodmean[,2],max=periodmax[,2])

relevel<-as.integer(period.country$mid)
relevelid<-data.frame(level=relevel,id=seq(1:length(relevel)))
relevelid<-relevelid[order(relevelid$level),]
relevelid$id2<-seq(1:length(relevel))
relevelid<-relevelid[order(relevelid$id2),]

period.country<-period.country[relevelid$id,]



#### Figure S1

library(tidyverse)

df_id <- filter(period.country) %>% mutate(idord = -relevelid$id2)

## join order and data
period.country <- filter(df_id) %>% mutate(country = fct_reorder(country, idord))

period.country$country<-droplevels(period.country$country)

Sys.setlocale('LC_TIME', "English")

ggplot(period.country, aes(mid, country)) + 
  geom_linerange(aes(xmin=start, xmax=max), colour = "#ef3b2c", size = .8) +
  geom_point(size=2, shape = 21, colour = "#ef3b2c", fill = "white") +
  theme_bw()+
  scale_y_discrete(limits = (levels(period.country$country))) + 
  scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
  labs(title = "", y = "", x = "")

ggsave("FigureS1.png", width = 5.47, height = 5.47) 


