library(tidyverse)
library(classInt)
library(patchwork)

#setwd("C:/Users/Francy/Dropbox/MCC_CMMID_COVID19/Programs/Analysis/Paper_V2")
load("C:/Users/Francy/Dropbox/MCC_CMMID_COVID19/Programs/Analysis/Build_analysis_dataset/mcc_covid_data_analysis.RData")

covid <- select(data, long, lat, casesw, tmean, r.mean, ah, rh, uv, kgclzone, oxgov10r,
                r.sd) %>% mutate(prec = 1/(r.sd^2)) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)


### scatterplots

# precision classes
br <- classIntervals(covid$prec, 5)

## Temperature


p1 <- mutate(covid, kgclzone = factor(kgclzone, LETTERS[1:4], 
                                      str_c("(", LETTERS[1:4],")") %>%
                                        str_c(c("Tropical", "Arid", "Temperate",
                                                "Continental"), sep = " "))) %>% 
  arrange(-prec) %>%
  ggplot(aes(tmean, r.mean, 
             size = prec, fill = kgclzone)) + 
  geom_hline(yintercept = 1, size = .7, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = .6, shape = 21, 
             colour = alpha("black", .5), 
             stroke = .2) +
  scale_y_continuous(breaks = seq(0, 3, .2), limits = c(0.6, 2.2)) +
  scale_x_continuous(breaks = seq(-10, 30, 5)) +
  scale_fill_manual(values = c("#08519c", "#cb181d", "#41ab5d", "#ae017e")) +
  scale_size(breaks = round(br$brks), range = c(2, 7)) +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4))) +
  labs(size = "Precision",
       fill = "Köppen climate zones",
       x = "Mean temperature (ºC)",
       y = expression(R[e])) +
  theme_bw()

#ggsave("scatterplot_TvsR.png", type = "cairo")


## AH
p2 <- mutate(covid, kgclzone = factor(kgclzone, LETTERS[1:4], 
                                      str_c("(", LETTERS[1:4],")") %>%
                                        str_c(c("Tropical", "Arid", "Temperate",
                                                "Continental"), sep = " "))) %>% 
  arrange(-prec) %>%
  ggplot(aes(ah, r.mean, 
             size = prec, fill = kgclzone)) + 
  geom_hline(yintercept = 1, size = .7, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = .6, shape = 21, 
             colour = alpha("black", .5), 
             stroke = .2) +
  scale_y_continuous(breaks = seq(0, 3, .2), limits = c(0.6, 2.2)) +
  scale_x_continuous(breaks = seq(0, 22, 2)) +
  scale_fill_manual(values = c("#08519c", "#cb181d", "#41ab5d", "#ae017e")) +
  scale_size(breaks = round(br$brks), range = c(2, 7)) +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4))) +
  labs(size = "Precision",
       fill = "Köppen climate zones",
       x = expression("Absolute"~"humidity"~"("~"g/"~m^3~")"),
       y = expression(R[e])) +
  theme_bw()

#ggsave("scatterplot_AHvsR.png", type = "cairo")


## RH
p3 <- mutate(covid, kgclzone = factor(kgclzone, LETTERS[1:4], 
                                      str_c("(", LETTERS[1:4],")") %>%
                                        str_c(c("Tropical", "Arid", "Temperate",
                                                "Continental"), sep = " "))) %>% 
  arrange(-prec) %>%
  ggplot(aes(rh, r.mean, 
             size = prec, fill = kgclzone)) + 
  geom_hline(yintercept = 1, size = .7, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = .6, shape = 21, 
             colour = alpha("black", .5), 
             stroke = .2) +
  scale_y_continuous(breaks = seq(0, 3, .2), limits = c(0.6, 2.2)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("#08519c", "#cb181d", "#41ab5d", "#ae017e")) +
  scale_size(breaks = round(br$brks), range = c(2, 7)) +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4))) +
  labs(size = "Precision",
       fill = "Köppen climate zones",
       x = "Relative humidity (%)",
       y = expression(R[e])) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



## Radiation
p4 <- mutate(covid, kgclzone = factor(kgclzone, LETTERS[1:4], 
                                      str_c("(", LETTERS[1:4],")") %>%
                                        str_c(c("Tropical", "Arid", "Temperate",
                                                "Continental"), sep = " "))) %>% 
  arrange(-prec) %>%
  ggplot(aes(uv, r.mean, 
             size = prec, fill = kgclzone)) + 
  geom_hline(yintercept = 1, size = .7, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = .6, shape = 21, 
             colour = alpha("black", .5), 
             stroke = .2) +
  scale_y_continuous(breaks = seq(0, 3, .2), limits = c(0.6, 2.2)) +
  scale_x_continuous(breaks = seq(80, 320, 20)) +
  scale_fill_manual(values = c("#08519c", "#cb181d", "#41ab5d", "#ae017e")) +
  scale_size(breaks = round(br$brks), range = c(2, 7)) +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4))) +
  labs(size = "Precision",
       fill = "Köppen climate zones",
       x = expression("Solar surface radiation"~"("~"J/"~m^2~")"),
       y = expression(R[e])) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


((p1 | p3) / (p2 | p4) ) +   
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = 6))


ggsave("Fig2.png", width = 7, height = 5)