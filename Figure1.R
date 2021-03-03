library(tidyverse)
library(sf)
library(rnaturalearth)
library(classInt)
library(biscale)
library(patchwork)
library(sf)

setwd("C:/Users/Francy/Dropbox/MCC_CMMID_COVID19/Programs/Analysis/GitHub")
load("mcc_covid_data_analysis.RData")

# MAPS
wld <- ne_countries(scale = 50, returnclass = "sf")
covid <- select(data, long, lat, casesw, tmean, r.mean, ah, rh, uv, kgclzone, oxgov10r,
  r.sd) %>% mutate(prec = 1/(r.sd^2)) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)


# graticulate grid
grid <- st_graticule(lon = seq(-180, 180, length.out = 7), 
  lat = seq(-90, 90, length.out = 5),
  ndiscr=1000, crs = 4326, margin = 0.00001) 

### BIVARIATE MAP

## Temperature vs Re GLOBAL
databi2 <- bi_class(covid, tmean, r.mean, style = "quantile", dim = 3)

legend2 <- bi_legend(pal = "DkViolet",
                     dim = 3,
                     xlab = "Higher Ta",
                     ylab = "Higher R",
                     size = 8)

ggplot(wld) + 
  geom_sf(fill = "grey93", colour = "white", size = .2) +
  geom_sf(data = databi2, aes(fill = bi_class), size = 1.8, stroke = .3,
          shape = 21, alpha = .7, colour = alpha("white", .5),
          show.legend = FALSE) +
  geom_sf(data = eu_ext_sf, fill = NA, colour = "red") +
  annotation_custom(ggplotGrob(legend2), xmin = -15616528, xmax = -9616528,
                    ymin = -3478079, ymax = 0) +
  geom_sf(data = grid, size = .2, linetype = "dashed", colour = "grey50", alpha = 0.5) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  bi_theme() +
  coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 

ggsave("world_map.pdf", width = 10.07, height = 5.46)


# inset map EUROPE 
eu_ext <- st_bbox(c(xmin = -20, xmax = 27, ymax = 61, ymin = 26), crs = st_crs(4326))

eu_ext_sf <- st_bbox(c(xmin = -1834289, xmax = 2076290, ymax = 6433621, ymin = 2980741), crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" )
eu_ext_sf <- st_as_sfc(eu_ext_sf)


eu_circ <- st_centroid(st_as_sfc(eu_ext)) %>% 
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
  st_buffer(2500*1000) 

eu_wld <- st_crop(wld, st_transform(eu_circ, 4326))
eu_databi2 <- st_crop(databi2, st_transform(eu_circ, 4326))

ggplot(eu_wld) + 
  geom_sf(fill = "grey93", colour = "white", size = .2) +
  geom_sf(data = eu_databi2, aes(fill = bi_class), size = 2, stroke = .3,
          shape = 21, alpha = .7, colour = alpha("white", .5),
          show.legend = FALSE) +
  #geom_sf(data = grid, size = .2, linetype = "dashed", colour = "grey50", alpha = 0.5) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  bi_theme() +
  theme(panel.border = element_rect(colour = "grey80", fill=NA, size=.8)) +
  coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
           xlim = c(-1834289, 2076290), ylim = c(2980741, 6433621)) 

ggsave("eu_map.pdf", width = 5, height = 2.5)


# inset map USA
usa_ext <- st_bbox(c(xmin = -130, xmax = -66, ymax = 56, ymin = 16), crs = st_crs(4326))

usa_circ <- st_centroid(st_as_sfc(usa_ext)) %>% 
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
  st_buffer(2600*1000) 

usa_wld <- st_crop(wld, st_transform(usa_circ, 4326))
usa_databi2 <- st_crop(databi2, st_transform(usa_circ, 4326))

usa_ext %>% st_as_sfc() %>% st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

ggplot(wld) + 
  geom_sf(fill = "grey93", colour = "white", size = .2) +
  geom_sf(data = usa_databi2, aes(fill = bi_class), size = 2, stroke = .3,
          shape = 21, alpha = .7, colour = alpha("white", .5),
          show.legend = FALSE) +
  geom_sf(data = grid, size = .2, linetype = "dashed", colour = "grey50", alpha = 0.5) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  bi_theme() +
  theme(panel.border = element_rect(colour = "grey80", fill=NA, size=.8)) +
  coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
           xlim = c(-11041860, -5162682), ylim = c(1711230, 5939113)) 

ggsave("usa_map.pdf", width = 5, height = 2.5)


# inset map Asia

jap_ext <- st_bbox(c(xmin = 100, xmax = 148, ymax = 46, ymin = 19), crs = st_crs(4326))

jap_circ <- st_centroid(st_as_sfc(jap_ext)) %>% 
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
  st_buffer(2000*1000) 

jap_wld <- st_crop(wld, st_transform(jap_circ, 4326))
jap_databi2 <- st_crop(databi2, st_transform(jap_circ, 4326))

jap_ext %>% st_as_sfc() %>% st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

ggplot(wld) + 
  geom_sf(fill = "grey93", colour = "white", size = .2) +
  geom_sf(data = jap_databi2, aes(fill = bi_class), size = 2, stroke = .3,
          shape = 21, alpha = .7, colour = alpha("white", .5),
          show.legend = FALSE) +
  #geom_sf(data = grid, size = .2, linetype = "dashed", colour = "grey50", alpha = 0.5) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  bi_theme() +
  theme(panel.border = element_rect(colour = "grey80", fill=NA, size=.8)) +
  coord_sf(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
           xlim = c(9015998, 12564670), ylim = c(2352084, 4910090)) 

ggsave("jap_map.pdf", width = 5, height = 2.5)


# distribution plot 

dkviol <- c("#3F2949", "#435786", "#4885C1", "#77324C", "#806A8A", "#89A1C8", "#AE3A4E", "#BC7C8F", "#CABED0")
names(dkviol) <- c("3-3", "2-3", "1-3", "3-2", "2-2", "1-2", "3-1", "2-1", "1-1")

dkviol <- dkviol[order(names(dkviol))]

tb <- table(databi2$bi_class) %>% as.data.frame()
names(tb) <- c("Class", "Count")

mutate(databi2, dummy = str_extract(bi_class, "[1-3]$") %>% as.numeric() %>% magrittr::multiply_by(-1) ) %>% 
  ggplot(aes(bi_class)) + geom_bar(aes(fill = bi_class), show.legend = FALSE) +
  facet_wrap(dummy~., nrow = 3, scales = "free_x") +
  scale_fill_manual(values = dkviol) +
  ggthemes::theme_hc() +
  theme(strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", size = .2),
        panel.ontop = TRUE)

ggsave("distr_class_bivariate.pdf")

