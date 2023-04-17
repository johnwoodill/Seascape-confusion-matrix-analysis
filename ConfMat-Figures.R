library(tidyverse)
library(plyr)


setwd("~/Projects/Seascape-confusion-matrix-analysis/")

dat = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/processed/conf_mats/ca_coast_2017.csv")

df_cm <- tidyr::gather(dat, key = "Predicted", value = "Count", -seascape, -acc, -region, -year)
df_cm = filter(df_cm, seascape != 0 | Predicted != 0)
seascapes_ordered = arrange(dat, -acc)$seascape

df_cm$seascape = factor(df_cm$seascape, levels=seascapes_ordered)
df_cm$Predicted = factor(df_cm$Predicted, levels=seascapes_ordered)

df_cm$Predicted
df_cm$seascape

labels = dplyr::select(dat, seascape, acc)
labels$seascape2 = labels$seascape
labels$acc = round(labels$acc, 2)
labels$seascape = factor(labels$seascape, levels=seascapes_ordered)
labels$seascape2 = factor(labels$seascape2, levels=seascapes_ordered)



# create a ggplot2 heatmap of the confusion matrix
ggplot(df_cm, aes(x = factor(seascape), y = Predicted, fill=Count)) +
  geom_tile() +
  geom_text(aes(label=Count), color='white') +
  ggtitle("2017 California Seascapes Confusion Matrix") + 
  labs(x = "True Class", y = "Predicted Class", fill = "Count") +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ggsave("figures/ca_coast_2017_conf_mat.png", width=10, height=10)
#



dat1 = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/processed/conf_mats/ca_coast_2017.csv")
dat2 = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/processed/conf_mats/ca_coast_2018.csv")
dat3 = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/processed/conf_mats/ca_coast_2019.csv")
dat4 = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/processed/conf_mats/ca_coast_2020.csv")
dat5 = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/processed/conf_mats/ca_coast_2021.csv")
dat6 = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/processed/conf_mats/ca_coast_2022.csv")

mdat = rbind.fill(dat1, dat2, dat3, dat4, dat5, dat6)
mdat <- tidyr::gather(mdat, key = "Predicted", value = "Count", -seascape, -acc, -region, -year)
mdat = filter(mdat, seascape != 0)

ggplot(mdat, aes(year, seascape, fill=acc)) + 
  geom_tile() +
  geom_text(aes(label=round(acc, 2)), color='white') +
  scale_y_continuous(breaks=seq(0, 33)) +
  scale_x_continuous(breaks=seq(2017, 2022, 1)) +
  theme_minimal(14) +
  ggtitle("California Seascapes") + 
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1), 
      legend.position = "none") +
  labs(x=NULL, y="Seascape Classification \n (% Accurate)") +
  NULL

ggsave("figures/ca_coast_2017_seascape_accuracy.png", width=8, height=8)

#





# Create map
ssdat = read_sf("~/Projects/test/data/cb_2018_us_states_500k/cb_2018_us_state_500k.shp")
st_crs(ssdat) <- "EPSG:4326" 

ssdat = filter(ssdat, NAME %in% c("California", "Nevada"))

cdat = read_csv("data/processed/oceans/ca_coast_2017.csv")
cdat$month = month(cdat$time)
cdat$month = month.abb[cdat$month]
cdat$month = factor(cdat$month, levels=month.abb)

cdat = filter(cdat, MCLASS != 0)
cdat$true = ifelse(cdat$MCLASS == cdat$VCLASS, "True", "False")

ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y, fill=factor(true))) +
  ggtitle("2017 California Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  theme(legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal(14) +
  facet_wrap(~month)

ggsave("figures/ca_seascape_map_2017.png", width=8, height=8)


