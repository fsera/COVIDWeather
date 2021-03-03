
library(maps)
library(mapproj)
library(mapdata)
library(maptools)
library(tidyverse)
library(classInt)
library(patchwork)

setwd("C:/Users/Francy/Dropbox/MCC_CMMID_COVID19/Programs/Analysis/GitHub")
load("mcc_covid_data_analysis.RData")

# TABLE BY COUNTRY

data$calyearo<-as.Date(data$calyear,origin = "1970-01-01")
data$const<-rep(1,nrow(data))

table(data$calyearo)

ncity<-aggregate(data$const,by=list(data$countryname),FUN=sum)
ncases<-aggregate(data$casesw,by=list(data$countryname),FUN=sum)
period<-aggregate(data$calyearo,by=list(data$countryname),FUN=mean)
reprnumb<-aggregate(data$r.mean,by=list(data$countryname),FUN=mean)
oxgov10r<-aggregate(data$oxgov10r,by=list(data$countryname),FUN=mean)

table.country<-data.frame(country=ncity[,1], cites=ncity[,2],covid=ncases[,2], midperiod=period[,2],r=reprnumb[,2],
  oxgov=oxgov10r[,2])

write.csv(table.country,file="Table1.csv")

