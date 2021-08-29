# PACKAGES
library(tidyverse)
library(classInt)
library(patchwork)
library(sf)

library(splines)
library(mixmeta)
library(dlnm)

# LOAD THE DATA
load("mcc_covid_data_analysis.RData")

# CALCULATE THE VARIANCE
data$r.var<-data$r.sd^2

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}

# MODEL A REML
fixed <- formula(r.mean~1)
modelAreml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelAreml)

# MODEL A ML
fixed <- formula(r.mean~1)
modelAml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelAml)

# MODEL B REML
fixed <- formula(r.mean~1+ns(oxgov10r,df=2))
modelBreml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelBreml)
fwald(modelBreml,"oxgov10r")

# MODEL B ML
fixed <- formula(r.mean~1+ns(oxgov10r,df=2))
modelBml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelBml)
fwald(modelBml,"oxgov10r")

# MODEL C REML
fixed <- formula(r.mean~1+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelCREML <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelCREML)

# MODEL C ML
fixed <- formula(r.mean~1+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelCML <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelCML)

# MODEL D1 REML
fixed <- formula(r.mean~1+ns(tmean,df=4)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD1reml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelD1reml)
fwald(modelD1reml,"tmean")

coef<-modelD1reml$coefficients[2:5]
vcov<-modelD1reml$vcov[2:5,2:5]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$tmean,df=4)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATION
xused<-data$tmean
xbas=onebasis(xused, "ns",df=4, intercept=F)
quants<-quantile(xused,probs=seq(0.05,0.95,.05))
cp1<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=20)
plot(cp1)

# MODEL D1 ML
fixed <- formula(r.mean~1+ns(tmean,df=4)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD1ml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
   control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelD1ml)

# MODEL D2 REML
fixed <- formula(r.mean~1+ns(rh,df=3)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD2reml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelD2reml)
fwald(modelD2reml,"rh")

coef<-modelD2reml$coefficients[2:4]
vcov<-modelD2reml$vcov[2:4,2:4]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$rh,df=3)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS
xused<-data$rh
xbas=onebasis(xused, "ns",df=3, intercept=F)
quants<-quantile(xused,probs=seq(0.05,0.95,.05))
cp3<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=75)
plot(cp3)

# MODEL D2 ML
fixed <- formula(r.mean~1+ns(rh,df=3)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD2ml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelD2ml)

# MODEL D3 REML
fixed <- formula(r.mean~1+ns(ah,df=4)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD3reml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelD3reml)
fwald(modelD3reml,"ah")

coef<-modelD3reml$coefficients[2:5]
vcov<-modelD3reml$vcov[2:5,2:5]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$ah,df=4)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS
xused<-data$ah
xbas=onebasis(xused, "ns",df=4, intercept=F)
quants<-quantile(xused,probs=seq(0.05,0.95,.05))
cp4<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=11)
plot(cp4)

# MODEL D3 ML
fixed <- formula(r.mean~1+ns(ah,df=4)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD3ml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelD3ml)
fwald(modelD3ml,"ah")

# MODEL D4 REML
fixed <- formula(r.mean~1+ns(uv,df=2)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD4reml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelD4reml)
fwald(modelD4reml,"uv")

coef<-modelD4reml$coefficients[2:3]
vcov<-modelD4reml$vcov[2:3,2:3]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$uv,df=2)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS
xused<-data$uv
xbas=onebasis(xused, "ns",df=2, intercept=F)
quants<-quantile(xused,probs=seq(0.05,0.95,.05))
cp6<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=180)
plot(cp6)

# MODEL D4 ML
fixed <- formula(r.mean~1+ns(uv,df=2)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD4ml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelD4ml)

# MODEL D5 REML
fixed <- formula(r.mean~1+ns(ws,df=2)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD5reml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelD5reml)
fwald(modelD5reml,"ws")

coef<-modelD5reml$coefficients[2:3]
vcov<-modelD5reml$vcov[2:3,2:3]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$ws,df=2)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS
xused<-data$ws
xbas=onebasis(xused, "ns",df=2, intercept=F)
quants<-quantile(xused,probs=seq(0.05,0.95,.05))
cp7<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=3)
plot(cp7)

# MODEL D5 ML
fixed <- formula(r.mean~1+ns(ws,df=2)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD5ml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelD5ml)

# MODEL D6 REML
fixed <- formula(r.mean~1+prectot+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD6reml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T), data ,method="reml", na.action=na.exclude)
summary(modelD6reml)

coef<-modelD6reml$coefficients[2:2]
vcov<-as.matrix(modelD6reml$vcov[2:2,2:2])

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns((data$prectot),df=1)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS
xused<-data$prectot
xbas=onebasis(xused, fun="lin")
quants<-quantile(xused,probs=seq(0.05,0.95,.05))
cp8<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=5)
plot(cp8)

# MODEL D6 ML
fixed <- formula(r.mean~1+prectot+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25)+ns(oxgov10r,df=2))
modelD6ml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T), data ,method="ml", na.action=na.exclude)
summary(modelD6ml)

# MODEL C oxford government index 
fixed <- formula(r.mean~1+ns(oxgov10r,df=2)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25))
modelCremla <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T), data ,method="reml", na.action=na.exclude)
summary(modelCremla)
fwald(modelCremla,"oxgov10r")

coef<-modelCremla$coefficients[2:3]
vcov<-modelCremla$vcov[2:3,2:3]

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
nsobj<-ns(data$oxgov10r,df=2)
knots<-attr(nsobj,"knots")

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS
xused<-data$oxgov10r
xbas=onebasis(xused, "ns",df=2, intercept=F)
quants<-quantile(xused,probs=seq(0.05,0.95,.05))
cp9<- crosspred(xbas,coef=coef,vcov=vcov,at=quants , cen=50)
plot(cp9)

# MODEL D7 REML
fixed <- formula(r.mean~1+ns(tmean,df=4)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25))
modelD7reml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelD7reml)

# MODEL D7 ML
fixed <- formula(r.mean~1+ns(tmean,df=4)+log(population)+log(density)+(oldpop)+(log(gdp))+log(pm25))
modelD7ml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelD7ml)

# MODEL D8 REML
fixed <- formula(r.mean~1+ns(tmean,df=4)+ns(oxgov10r,df=2))
modelD8reml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="reml", na.action=na.exclude)
summary(modelD8reml)

# MODEL D8 ML
fixed <- formula(r.mean~1+ns(tmean,df=4)+ns(oxgov10r,df=2))
modelD8ml <- mixmeta(fixed,r.var, random=~1|countrycode/mcccitiescode,
  control=list(showiter=T),data,method="ml", na.action=na.exclude)
summary(modelD8ml)

##################################################
###### Figure 3
##################################################

# RR extract function

rrchange_extract <- function(x){
  
  tibble(var = as.numeric(names(x$allfit)), fit = x$allfit, low = x$alllow, high = x$allhigh)
  
} 


# plot preparation

change <- list(A = cp1, B = cp3, C = cp4, D = cp9)

df <- map_df(change, rrchange_extract, .id = "panel")
df <- mutate(df, panel = factor(panel, LETTERS[1:4]))

range_val <- group_by(df, panel) %>% summarise(mx = max(var),  min = min(var))

vardistr <- dplyr::select(data, tmean, rh, ah, oxgov10r) %>% 
  pivot_longer(1:4, names_to = "panel", values_to = "var") %>% 
  mutate(panel = factor(panel, c("tmean", "rh", "ah", "oxgov10r"),
  LETTERS[1:4]))

vardistr <- left_join(vardistr, range_val)
vardistr <- filter(vardistr, var <= mx, var >= min)

# plots

p1 <- filter(df, panel == "A") %>% 
       ggplot() + 
        geom_hline(yintercept = 0, 
                   size = .5, 
                   linetype = "dashed", 
                   colour = "grey50") +
        geom_ribbon(aes(var, 
                        ymin = low, 
                        ymax = high), 
                    fill = "grey90", 
                    alpha = .6) +
        geom_line(aes(var, 
                      fit, 
                      colour = panel),
                  size = .8,
                  show.legend = FALSE) +
        geom_rug(data = filter(vardistr, panel == "A"), 
                 aes(var, colour = panel), 
                 alpha = .5, 
                 inherit.aes = TRUE, 
                 show.legend = FALSE) +
       scale_color_manual(values = c("#cb181d", "#08519c", "#08519c", "#ae017e")) +
       scale_y_continuous(limits = c(-0.21, .2)) +
       labs(y = expression(R[e]~"difference"), 
            x = "Mean temperature (ºC)", 
            title = "p = 0.014") +
      theme_bw() +
      theme(strip.text = element_text(hjust = 0.08, size = 10),
            plot.title = element_text(face = "plain", size = 10))


p2 <- filter(df, panel == "B") %>% 
        ggplot() + 
        geom_hline(yintercept = 0, 
                   size = .5, 
                   linetype = "dashed", 
                   colour = "grey50") +
        geom_ribbon(aes(var, 
                        ymin = low, 
                        ymax = high), 
                    fill = "grey90", 
                    alpha = .6) +
        geom_line(aes(var, 
                      fit), 
                  size = .8,
                  colour = "#08519c",
                  show.legend = FALSE) +
        geom_rug(data = filter(vardistr, panel == "B"), 
                 aes(var), 
                 colour = "#08519c", 
                 alpha = .5, 
                 inherit.aes = TRUE,
                 show.legend = FALSE) +
        scale_y_continuous(limits = c(-0.21, .2)) +
        labs(y = expression(R[e]~"difference"), 
             x = "Relative humidity (%)", 
             title = "p = 0.058") +
        theme_bw() +
        theme(strip.text = element_text(hjust = 0.08, size = 10),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(face = "plain", size = 10))


p3 <- filter(df, panel == "C") %>% 
        ggplot() + 
        geom_hline(yintercept = 0, 
                   size = .5, 
                   linetype = "dashed", 
                   colour = "grey50") +
        geom_ribbon(aes(var, 
                        ymin = low, 
                        ymax = high), 
                    fill = "grey90", 
                    alpha = .6) +
        geom_line(aes(var, 
                     fit),
                 size = .8, 
                 colour = "#08519c",
                 show.legend = FALSE) +
        geom_rug(data = filter(vardistr, panel == "C"), aes(var), colour = "#08519c", inherit.aes = TRUE, alpha = .5, show.legend = FALSE) +
        scale_y_continuous(limits = c(-0.21, .2)) +
        labs(y = expression(R[e]~"difference"), 
             x = expression("Absolute humidity (g/"*m^3*")"), 
             title = "p = 0.036") +
        theme_bw() +
        theme(strip.text = element_text(hjust = 0.08, size = 10),
              plot.title = element_text(face = "plain", size = 10))



p4 <- filter(df, panel == "D") %>%  
        ggplot() + 
        geom_hline(yintercept = 0, 
                   size = .5, 
                   linetype = "dashed", 
                   colour = "grey50") +
        geom_ribbon(aes(var, 
                        ymin = low, 
                        ymax = high), 
                    fill = "grey90", 
                    alpha = .6) +
        geom_line(aes(var, 
                      fit), 
                  size = .8, 
                  colour = "#ae017e",
                  show.legend = FALSE) +
        geom_rug(data = filter(vardistr, panel == "D"), 
                 aes(var), 
                 colour = "#ae017e",
                 alpha = .5, 
                 inherit.aes = TRUE,
                 show.legend = FALSE) +
        scale_y_continuous(limits = c(-0.21, .2)) +
        labs(y = expression(R[e]~"difference"), 
             x = "Oxford Government Index", 
             title = "p < 0.001") +
        theme_bw() +
        theme(strip.text = element_text(hjust = 0.08, size = 10),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(face = "plain", size = 10))



((p1 | p2) / (p3 | p4) ) +   
  plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = 6))

ggsave("Figure3.png", width = 5, height = 5)


