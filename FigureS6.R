
setwd("C:/Users/Francy/Dropbox/MCC_CMMID_COVID19/Programs/Analysis/Paper_V2")
load("C:/Users/Francy/Dropbox/MCC_CMMID_COVID19/Programs/Analysis/Build_analysis_dataset/mcc_covid_data_analysis.RData")

data$abslat<-abs(data$lat)
abs
conf<-data[,c("tmean","rh","ah","uv","oxgov","oxgov10","oxgov10r","calyear","r.mean","abslat")]

res <- cor(conf)
round(res, 2)

library(corrplot)
corrplot(res, type = "upper", 
         tl.col = "black", tl.srt = 45,cl.cex=0.6)


exc<-!(data$r.mean<1)

datar<-data[exc,]

# I CALCULATE THE VARIANCE
datar$r.var<-datar$r.sd^2

datar$sh<-datar$sh*100000

datar$kgclzone <- factor(datar$kgclzone, levels = c("C","A","B","D"))

library(splines)

library(mixmeta)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}

# PLOTS

library(dlnm)

# MODEL BASE

fixed <- formula(r.mean~1)

model0 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="reml", na.action=na.exclude)
summary(model0)

# MODEL 1

fixed <- formula(r.mean~1+ns(oxgov10r,df=2))

model1 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="reml", na.action=na.exclude)
summary(model1)
fwald(model1,"oxgov10r")

# MODEL 1

fixed <- formula(r.mean~1+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))

model1 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="reml", na.action=na.exclude)
summary(model1)
fwald(model1,"oxgov10r")

# MEAN TEMPERATURE

fixed <- formula(r.mean~1+ns(tmean,df=4)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2)+ns(oxgov10r,df=2))

model3 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="reml", na.action=na.exclude)
summary(model3)
fwald(model3,"tmean")


coef<-model3$coefficients[2:5]
vcov<-model3$vcov[2:5,2:5]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$tmean,df=4)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

predper <- c(1:100)

tm_meancountry <- quantile(jitter(data$tmean), predper/100,na.rm=T)

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(tm_meancountry,fun="ns",
                 knots=knots)

cp1 <- crosspred(bvar,coef=coef,vcov=vcov,
                 model.link="lin", at=tm_meancountry,cen=20)

xused<-datar$tmean

xbas=onebasis(xused, "ns",df=4, intercept=F)

mn<-mean(xused)

quants<-quantile(xused,probs=seq(0.05,0.95,.05))

cp1<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=20)

plot(cp1)

c(100*sd(cp1$allfit)/0.15)
c(100*var(cp1$allfit)/0.2434^2)

#Relative humidity

fixed <- formula(r.mean~1+ns(rh,df=3)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))

model3 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="reml", na.action=na.exclude)
summary(model3)
fwald(model3,"rh")

coef<-model3$coefficients[2:4]
vcov<-model3$vcov[2:4,2:4]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$rh,df=3)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

predper <- c(1:100)

tm_meancountry <- quantile(jitter(data$rh), predper/100,na.rm=T)

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(tm_meancountry,fun="ns",
                 knots=knots)

cp3 <- crosspred(bvar,coef=coef,vcov=vcov,
                 model.link="lin", at=tm_meancountry,cen=70)

xused<-datar$rh

xbas=onebasis(xused, "ns",df=3, intercept=F)

mn<-mean(xused)

quants<-quantile(xused,probs=seq(0.05,0.95,.05))

cp3<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=75)

plot(cp3)

c(100*sd(cp3$allfit)/0.15)
c(100*var(cp3$allfit)/0.24^2)


#Absolute humidity


fixed <- formula(r.mean~1+ns(ah,df=4)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))

model3 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="reml", na.action=na.exclude)
summary(model3)
fwald(model3,"ah")

coef<-model3$coefficients[2:5]
vcov<-model3$vcov[2:5,2:5]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$ah,df=4)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

predper <- c(1:100)

tm_meancountry <- quantile(jitter(data$ah), predper/100,na.rm=T)

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(tm_meancountry,fun="ns",
                 knots=knots)

cp4 <- crosspred(bvar,coef=coef,vcov=vcov,
                 model.link="lin", at=tm_meancountry,cen=11)

xused<-datar$ah

xbas=onebasis(xused, "ns",df=4, intercept=F)

mn<-mean(xused)

quants<-quantile(xused,probs=seq(0.05,0.95,.05))

cp4<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=11)

plot(cp4)

c(100*sd(cp4$allfit)/0.15)
c(100*var(cp4$allfit)/0.24^2)


#UV


fixed <- formula(r.mean~1+ns(uv,df=2)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))

model3 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="reml", na.action=na.exclude)
summary(model3)
fwald(model3,"uv")

coef<-model3$coefficients[2:3]
vcov<-model3$vcov[2:3,2:3]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$uv,df=2)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

predper <- c(1:100)

tm_meancountry <- quantile(jitter(data$uv), predper/100,na.rm=T)

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(tm_meancountry,fun="ns",
                 knots=knots)

cp6 <- crosspred(bvar,coef=coef,vcov=vcov,
                 model.link="lin", at=tm_meancountry,cen=200)

xused<-datar$uv

xbas=onebasis(xused, "ns",df=2, intercept=F)

mn<-mean(xused)

quants<-quantile(xused,probs=seq(0.05,0.95,.05))

cp6<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=180)

plot(cp6)

c(100*sd(cp6$allfit)/0.15)
c(100*var(cp6$allfit)/0.24^2)

#WS


fixed <- formula(r.mean~1+ns(ws,df=2)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))

model3 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="reml", na.action=na.exclude)
summary(model3)
fwald(model3,"ws")

coef<-model3$coefficients[2:3]
vcov<-model3$vcov[2:3,2:3]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$ws,df=2)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

predper <- c(1:100)

tm_meancountry <- quantile(jitter(data$ws), predper/100,na.rm=T)

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(tm_meancountry,fun="ns",
                 knots=knots)

cp7 <- crosspred(bvar,coef=coef,vcov=vcov,
                 model.link="lin", at=tm_meancountry,cen=3)

xused<-datar$ws

xbas=onebasis(xused, "ns",df=2, intercept=F)

mn<-mean(xused)

quants<-quantile(xused,probs=seq(0.05,0.95,.05))

cp7<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=3)

plot(cp7)

c(100*sd(cp7$allfit)/0.15)
c(100*var(cp7$allfit)/0.24^2)

# PRECIPITATION

fixed <- formula(r.mean~1+prectot+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))

model0 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T), datar ,method="reml", na.action=na.exclude)
summary(model0)

coef<-model0$coefficients[2:2]
vcov<-model0$vcov[2:2,2:2]

library(dlnm)

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns((datar$prectot),df=1)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

predper <- c(1:100)

meancountry <- quantile(jitter((datar$prectot)), predper/100,na.rm=T)

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(meancountry,fun="lin")

vcov<-as.matrix(vcov)

cp8 <- crosspred(bvar, model=NULL, coef=coef,vcov=vcov,
                 model.link="lin", meancountry,cen=5)

xused<-datar$prectot

xbas=onebasis(xused, fun="lin")

mn<-mean(xused)

quants<-quantile(xused,probs=seq(0.05,0.95,.05))

cp8<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=5)

plot(cp8)

c(100*sd(cp8$allfit)/0.15)
c(100*var(cp8$allfit)/0.24^2)


# oxford government index 

fixed <- formula(r.mean~1+ns(oxgov10r,df=2)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25))

model0 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T), datar ,method="reml", na.action=na.exclude)
summary(model0)
fwald(model0,"oxgov10r")

coef<-model0$coefficients[2:3]
vcov<-model0$vcov[2:3,2:3]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$oxgov10r,df=2)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

predper <- c(1:100)

tm_meancountry <- quantile(jitter(data$oxgov10r), predper/100,na.rm=T)

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(tm_meancountry,fun="ns",
                 knots=knots)

cp9 <- crosspred(bvar,coef=coef,vcov=vcov,
                 model.link="lin", at=tm_meancountry,cen=50)

xused<-datar$oxgov10r

xbas=onebasis(xused, "ns",df=2, intercept=F)

mn<-mean(xused)

quants<-quantile(xused,probs=seq(0.05,0.95,.05))

cp9<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=50)

plot(cp9)

c(100*sd(cp9$allfit)/0.15)
c(100*var(cp9$allfit)/0.24^2)

library(dlnm)


# MULTIVARIABLE MODEL

fixed <- formula(r.mean~1+ns(tmean,df=4)+ns(rh,df=3)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))

model3 <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
                  control=list(showiter=T),datar,method="ml", na.action=na.exclude)
summary(model3)
fwald(model3,"tmean")
fwald(model3,"rh")

drop1(model3,test="Chisq")

##################################################
###### Figure 
##################################################

# packages 
library(tidyverse)
library(patchwork)

# functions

rrchange_extract <- function(x){
  
  tibble(var = as.numeric(names(x$allfit)), fit = x$allfit, low = x$alllow, high = x$allhigh)
  
} 

# plot preparation

change <- list(A = cp1, B = cp3, C = cp4, D = cp6)

df <- map_df(change, rrchange_extract, .id = "panel")
df <- mutate(df, panel = factor(panel, LETTERS[1:4]))


range_val <- group_by(df, panel) %>% summarise(mx = max(var),  min = min(var))

vardistr <- dplyr::select(datar, tmean, rh, ah, uv) %>% 
  pivot_longer(1:4, names_to = "panel", values_to = "var") %>% 
  mutate(panel = factor(panel, c("tmean", "rh", "ah", "uv"),
                        LETTERS[1:4]))

vardistr <- left_join(vardistr, range_val)
vardistr <- filter(vardistr, var <= mx, var >= min)

# plots

p1 <- filter(df, panel == "A") %>% ggplot() + 
  geom_hline(yintercept = 0, size = .5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(var, ymin = low, ymax = high), fill = "grey90", alpha = .6) +
  geom_line(aes(var, fit, colour = panel), show.legend = FALSE, size = .8) +
  geom_rug(data = filter(vardistr, panel == "A"), aes(var, colour = panel), inherit.aes = TRUE, alpha = .5, show.legend = FALSE) +
  scale_color_manual(values = c("#cb181d", "#08519c", "#08519c", "#ae017e")) +
  scale_y_continuous(limits = c(-0.2, .2)) +
  labs(y = expression(R[e]~"difference"), x = "Mean temperature (ºC)", title = "p = 0.001") +
  
  theme_bw() +
  theme(strip.text = element_text(hjust = 0.08, size = 10),
        plot.title = element_text(face = "plain", size = 10))


p2 <- filter(df, panel == "B") %>% ggplot() + 
  geom_hline(yintercept = 0, size = .5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(var, ymin = low, ymax = high), fill = "grey90", alpha = .6) +
  geom_line(aes(var, fit), show.legend = FALSE, size = .8, colour = "#08519c") +
  geom_rug(data = filter(vardistr, panel == "B"), aes(var), colour = "#08519c", inherit.aes = TRUE, alpha = .5, show.legend = FALSE) +
  scale_y_continuous(limits = c(-0.21, .2)) +
  labs(y = expression(R[e]~"difference"), x = "Relative humidity (%)", title = "p = 0.009") +
  theme_bw() +
  theme(strip.text = element_text(hjust = 0.08, size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "plain", size = 10))


p3 <- filter(df, panel == "C") %>% ggplot() + 
  geom_hline(yintercept = 0, size = .5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(var, ymin = low, ymax = high), fill = "grey90", alpha = .6) +
  geom_line(aes(var, fit), show.legend = FALSE, size = .8, colour = "#08519c") +
  geom_rug(data = filter(vardistr, panel == "C"), aes(var), colour = "#08519c", inherit.aes = TRUE, alpha = .5, show.legend = FALSE) +
  #scale_color_manual(values = c("#cb181d", "#08519c", "#08519c", "#ae017e")) +
  scale_y_continuous(limits = c(-0.2, .2)) +
  labs(y = expression(R[e]~"difference"), x = expression("Absolute humidity (g/"*m^3*")"), title = "p = 0.003") +
  theme_bw() +
  theme(strip.text = element_text(hjust = 0.08, size = 10),
        plot.title = element_text(face = "plain", size = 10))



p4 <- filter(df, panel == "D") %>% ggplot() + 
  geom_hline(yintercept = 0, size = .5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(var, ymin = low, ymax = high), fill = "grey90", alpha = .6) +
  geom_line(aes(var, fit), show.legend = FALSE, size = .8, colour = "#fd8d3c") +
  geom_rug(data = filter(vardistr, panel == "D"), aes(var), colour = "#fd8d3c", inherit.aes = TRUE, alpha = .5, show.legend = FALSE) +
  scale_y_continuous(limits = c(-0.2, .2)) +
  labs(y = expression(R[e]~"difference"), x = expression("Solar surface radiation (J/"*m^2*")"), title = "p = 0.047") +
  theme_bw() +
  theme(strip.text = element_text(hjust = 0.08, size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "plain", size = 10))



((p1 | p2) / (p3 | p4) ) +   
  plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = 6))

ggsave("FigureS6.png", width = 5, height = 5)



