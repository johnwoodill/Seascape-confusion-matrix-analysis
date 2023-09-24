library(tidyverse)
library(cowplot)
library(scales)
library(sf)
library(lubridate)

setwd("~/Projects/Seascape-confusion-matrix-analysis/")


loc = "pnw_coast"
loc = "ca_coast"
loc = "gulf_of_mexico"  
loc = "north_atlantic"
loc = "chukchi_sea"

plot_title = "test"

time_step = "monthly"

# Create tile map for most dominant seascapes in all years
create_acc_tile_map = function(loc, plot_title, time_step){
  # Get data for all years given location
  # Get a list of all CSV files that contain "my_string" in their name
  if (time_step == "monthly"){
  
    files <- list.files("data/processed/conf_mats/", pattern = paste0(loc, ".*\\.csv$"), full.names = TRUE)
    files <- files[!grepl("8day", files)]
  }
  
  if (time_step == "8day"){
  
    files <- list.files("data/processed/conf_mats/", pattern = paste0(loc, ".*8day.*\\.csv$"), full.names = TRUE)
  
  }
  
  # Loop through each file, read the CSV, and bind it to the data frame
  indat <- data.frame()
  for (file in files) {
    print(file)
    temp_data <- read_csv(file)
    indat <- bind_rows(indat, temp_data)
  }
  
  # Get most dominant seascapes in all years (90% percentile)
  indat$totals <- rowSums(indat[, c(5:length(names(indat)))], na.rm=TRUE)
  
  # Filter seascapes less than 23 
  indat = filter(indat, seascape < 23 & seascape != 0)
  indat = filter(indat, acc > 0)
  
  seascape_totals <- indat %>% 
    group_by(seascape) %>% 
    summarize(seascape_totals = sum(acc),
              n = n()) %>% 
    filter(seascape_totals >= .10 & n == 6)
  
  indat = filter(indat, seascape %in% seascape_totals$seascape)
  
  # Clean up data
  # mdat <- tidyr::gather(indat, key = "Seascape_Predicted", value = "Total", -seascape, -acc, -region, -year, -month)
  
  indat$seascape <- factor(as.numeric(as.character(indat$seascape)), levels = sort(unique(as.numeric(as.character(indat$seascape)))))
  
  # PLot
  plot_ = ggplot(indat, aes(year, seascape, fill=acc)) + 
    geom_tile() +
    geom_text(aes(label=paste0(round(acc * 100, 1), "%")), color='white', size=3.5) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_continuous(breaks=seq(2017, 2022, 1), expand=c(0,0)) +
    theme_minimal(14) +
    labs(title=plot_title) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=1), 
           legend.position = "none") +
    labs(x=NULL, y="Seascape (% Accurate)") +
    NULL
  
  # Get filename and save
  ggsave(plot = plot_, filename=paste0("figures/May2023_poster_figures/", loc, "_conf_matrix_year.png"), width=6, height=5)
  plot_
  # (NOT REALLY INFORMATIVE)
  # indat2 = indat %>% group_by(seascape, month) %>% summarise(acc = mean(acc))
  # ggplot(indat2, aes(month, seascape, fill=acc)) + 
  #   geom_tile()  +
  #   geom_text(aes(label=paste0(round(acc * 100, 1), "%")), color='white', size=3.5) +
  #   scale_y_discrete(expand=c(0,0)) +
  #   scale_x_continuous(breaks=seq(1, 12, 1), expand=c(0,0)) +
  #   theme_minimal(14) +
  #   labs(title=plot_title) +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1), 
  #         legend.position = "none") +
  #   labs(x=NULL, y="Seascape (% Accurate)") +
  #   # facet_wrap(~year)
  #   NULL
  #   
  # plot2_

  indat3 = indat %>% group_by(year, month) %>% summarise(acc = mean(acc, na.rm=TRUE))
  indat3$day = 15
  indat3$date = paste0(indat3$year, "-", indat3$month, "-", indat3$day)
  indat3$date <- as.Date(indat3$date, format = "%Y-%m-%d")

  indat3$acc = ifelse(is.na(indat3$acc), 0, indat3$acc)
  
  plot2_ = ggplot(indat3, aes(date, acc)) + 
    geom_bar(stat='identity') +
    theme_minimal(15) +
    labs(x=NULL, y="Monthly Accuracy (%)", title=plot_title) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ylim(0, 1) +
    theme(plot.title = element_text(hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1), 
            legend.position = "none")
    NULL

    ggsave(plot = plot2_, filename=paste0("figures/May2023_poster_figures/", loc, "_conf_matrix_year_month.png"), width=8, height=5)

  
}




# Create tile map for most dominant seascapes in all years
create_f1_tile_map = function(loc, plot_title){
  # Get data for all years given location
  # Get a list of all CSV files that contain "my_string" in their name
  files <- list.files("data/processed/conf_mats/", pattern = paste0(loc, ".*\\.csv$"), full.names = TRUE)
  
  # Loop through each file, read the CSV, and bind it to the data frame
  indat <- data.frame()
  for (file in files) {
    print(file)
    temp_data <- read_csv(file)
    indat <- bind_rows(indat, temp_data)
  }
  
  # Get most dominant seascapes in all years (90% percentile)
  # indat$totals <- rowSums(indat[, c(5:length(names(indat)))], na.rm=TRUE)
  
  # Filter seascapes less than 23 
  indat = filter(indat, seascape < 23 & seascape != 0)
  indat = filter(indat, f1 > 0)
  
  seascape_totals <- indat %>%
    group_by(seascape) %>%
    summarize(seascape_totals = sum(f1),
              n = n()) %>%
    filter(seascape_totals >= .10 & n == 6)

  indat = filter(indat, seascape %in% seascape_totals$seascape)
  
  # Clean up data
  # mdat <- tidyr::gather(indat, key = "Seascape_Predicted", value = "Total", -seascape, -acc, -region, -year, -month)
  
  indat$seascape <- factor(as.numeric(as.character(indat$seascape)), levels = sort(unique(as.numeric(as.character(indat$seascape)))))
  
  # PLot
  ggplot(indat, aes(year, seascape, fill=f1)) + 
    geom_tile() +
    geom_text(aes(label=paste0(round(f1 * 100, 1), "%")), color='white', size=3.5) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_continuous(breaks=seq(2017, 2022, 1), expand=c(0,0)) +
    theme_minimal(14) +
    labs(title=plot_title) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=1), 
           legend.position = "none") +
    labs(x=NULL, y="Seascape (F1-Score)") +
    NULL
  
  # Get filename and save
  ggsave(plot = plot_, filename=paste0("figures/May2023_poster_figures/", loc, "_conf_matrix_f1score_year.png"), width=6, height=5)
  
  # (NOT REALLY INFORMATIVE)
  # indat2 = indat %>% group_by(seascape, month) %>% summarise(acc = mean(acc))
  # ggplot(indat2, aes(month, seascape, fill=acc)) + 
  #   geom_tile()  +
  #   geom_text(aes(label=paste0(round(acc * 100, 1), "%")), color='white', size=3.5) +
  #   scale_y_discrete(expand=c(0,0)) +
  #   scale_x_continuous(breaks=seq(1, 12, 1), expand=c(0,0)) +
  #   theme_minimal(14) +
  #   labs(title=plot_title) +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1), 
  #         legend.position = "none") +
  #   labs(x=NULL, y="Seascape (% Accurate)") +
  #   # facet_wrap(~year)
  #   NULL
  #   
  # plot2_

  indat3 = indat %>% group_by(year, month) %>% summarise(f1 = mean(f1, na.rm=TRUE))
  indat3$day = 15
  indat3$date = paste0(indat3$year, "-", indat3$month, "-", indat3$day)
  indat3$date <- as.Date(indat3$date, format = "%Y-%m-%d")

  indat3$acc = ifelse(is.na(indat3$f1), 0, indat3$f1)
  
  plot2_ = ggplot(indat3, aes(date, f1)) + 
    geom_bar(stat='identity') +
    theme_minimal(15) +
    labs(x=NULL, y="Monthly Accuracy (%)", title=plot_title) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ylim(0, 1) +
    theme(plot.title = element_text(hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1), 
            legend.position = "none")
    NULL

    ggsave(plot = plot2_, filename=paste0("figures/May2023_poster_figures/", loc, "_conf_matrix_year_month.png"), width=8, height=5)

  
}





create_tile_map("pnw_coast", "Pacific Northwest")
create_tile_map("ca_coast", "California Coast")
create_tile_map("gulf_of_mexico", "Gulf of Mexico")
create_tile_map("north_atlantic", "North Atlantic")
create_tile_map("chukchi_sea", "Chukchi Sea")




loc = "pnw_coast"
loc = "ca_coast"
loc = "gulf_of_mexico"  
loc = "north_atlantic"
loc = "chukchi_sea"


# Create tile map for most dominant (90-95%) seascapes in all months in a year
create_tile_map = function(loc, plot_title){
  # Get data for all years given location
  # Get a list of all CSV files that contain "my_string" in their name
  files <- list.files("data/processed/conf_mats/", pattern = paste0(loc, ".*\\.csv$"), full.names = TRUE)

  # Loop through each file, read the CSV, and bind it to the data frame
  indat <- data.frame()
  for (file in files) {
    print(file)
    temp_data <- read_csv(file)
    indat <- bind_rows(indat, temp_data)
  }
  
  # Get most dominant seascapes in all years (90% percentile)
  indat$totals <- rowSums(indat[, c(5:length(names(indat)))], na.rm=TRUE)
  
  # Filter seascapes less than 23 
  indat = filter(indat, seascape < 23 & seascape != 0)
  indat = filter(indat, acc > 0)
  
  seascape_totals <- indat %>% 
    group_by(seascape) %>% 
    summarize(seascape_totals = sum(acc),
              n = n()) %>% 
    filter(seascape_totals >= .10 & n == 6)
  
  indat = filter(indat, seascape %in% seascape_totals$seascape)
  
  # Clean up data
  mdat <- tidyr::gather(indat, key = "Predicted", value = "Count", -seascape, -acc, -region, -year)
  
  mdat$seascape <- factor(as.numeric(as.character(mdat$seascape)), levels = sort(unique(as.numeric(as.character(mdat$seascape)))))
  
  # PLot
  plot_ = ggplot(mdat, aes(year, seascape, fill=acc)) + 
    geom_tile() +
    geom_text(aes(label=paste0(round(acc * 100, 1), "%")), color='white', size=3.5) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_continuous(breaks=seq(2017, 2022, 1), expand=c(0,0)) +
    labs(title=plot_title)
    theme_minimal(14) +
    theme(plot.title = element_text(hjust = 0.5) +  # Center the title, 
          panel.border = element_rect(colour = "black", fill=NA, size=1), 
          legend.position = "none") +
    labs(x=NULL, y="Seascape (% Accurate)") +
    NULL
  
  # Get filename and save
  ggsave(plot = plot_, filename=paste0("figures/May2023_poster_figures/", loc, "_conf_matrix_year.png"), width=6, height=5)
  return(plot_)
}

# -------------------------------------------
# California Coast

# Create map
ssdat = read_sf("~/Projects/test/data/cb_2018_us_states_500k/cb_2018_us_state_500k.shp")
st_crs(ssdat) <- "EPSG:4326" 

ssdat = filter(ssdat, NAME %in% c("California", "Nevada"))

cdat = read_csv("data/processed/oceans/ca_coast_2017.csv")
cdat$month = month(cdat$time)
cdat$month = month.abb[cdat$month]
cdat$month = factor(cdat$month, levels=month.abb)

cdat = filter(cdat, MCLASS != 0 & month == "Jun")
cdat$true = ifelse(cdat$MCLASS == cdat$VCLASS, "True", "False")

cdat$MCLASS = as.factor(cdat$MCLASS)

ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y, fill=factor(MCLASS))) +
  ggtitle("June 2017 California Coast Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  scale_fill_viridis_d(option = "G") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow=1))


ggsave("figures/May2023_poster_figures/ca_seascape_map_2017.png", width=8, height=8)


ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y)) +
  ggtitle("June 2017 California Coast Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  # scale_fill_viridis_d(option = "G") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal(14) +
  # facet_wrap(~month) +
  NULL

ggsave("figures/May2023_poster_figures/ca_seascape_map_2017_blank.png", width=8, height=8)



# -------------------------------------------
# Chukchi Sea


# Create map
ssdat = read_sf("~/Projects/test/data/cb_2018_us_states_500k/cb_2018_us_state_500k.shp")
st_crs(ssdat) <- "EPSG:4326" 

ssdat = filter(ssdat, NAME %in% c("Alaska"))

cdat = read_csv("data/processed/oceans/chukchi_beaufort_sea_2017.csv")
cdat$month = month(cdat$time)
cdat$month = month.abb[cdat$month]
cdat$month = factor(cdat$month, levels=month.abb)

cdat = filter(cdat, MCLASS != 0 & month == "Jun")
cdat$true = ifelse(cdat$MCLASS == cdat$VCLASS, "True", "False")

cdat$MCLASS = as.factor(cdat$MCLASS)

ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y, fill=factor(MCLASS))) +
  ggtitle("June 2017 Chukchi Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  scale_fill_viridis_d(option = "G") +
  theme_minimal() +
  xlim(-180, -120) +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow=1))


ggsave("figures/May2023_poster_figures/chukchi_seascape_map_2017.png", width=8, height=8)


ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y)) +
  ggtitle("June 2017 Chukchi Sea Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  # scale_fill_viridis_d(option = "G") +
  xlim(-180, -120) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal(14) +
  # facet_wrap(~month) +
  NULL

ggsave("figures/May2023_poster_figures/chukchi_seascape_map_2017_blank.png", width=8, height=8)



# -------------------------------------------
# Gulf of Mexico

# Create map
ssdat = read_sf("~/Projects/test/data/cb_2018_us_states_500k/cb_2018_us_state_500k.shp")
st_crs(ssdat) <- "EPSG:4326" 

ssdat = filter(ssdat, NAME %in% c("Alabama", "Mississippi", "Florida", "Texas", "Louisiana", "Georgia", "Arkansas", "Oklahoma"))

cdat = read_csv("data/processed/oceans/gulf_of_mexico_2017.csv")
cdat$month = month(cdat$time)
cdat$month = month.abb[cdat$month]
cdat$month = factor(cdat$month, levels=month.abb)

cdat = filter(cdat, MCLASS != 0 & month == "Jun")
cdat$true = ifelse(cdat$MCLASS == cdat$VCLASS, "True", "False")

cdat$MCLASS = as.factor(cdat$MCLASS)

ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y, fill=factor(MCLASS))) +
  ggtitle("June 2017 Gulf of Mexico Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  scale_fill_viridis_d(option = "G") +
  theme_minimal() +
  ylim(22, 35) +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow=1))


ggsave("figures/May2023_poster_figures/gulf_of_mexico_seascape_map_2017.png", width=8, height=8)


ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y)) +
  ggtitle("June 2017 Gulf of Mexico Sea Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  # scale_fill_viridis_d(option = "G") +
  theme_minimal() +
  ylim(22, 35) +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal(14) +
  # facet_wrap(~month) +
  NULL

ggsave("figures/May2023_poster_figures/gulf_of_mexico_seascape_map_2017_blank.png", width=8, height=8)




# -------------------------------------------
# North Atlantic

ssdat = read_sf("~/Projects/test/data/cb_2018_us_states_500k/cb_2018_us_state_500k.shp")
st_crs(ssdat) <- "EPSG:4326" 

ssdat = filter(ssdat, NAME %in% c("Maine", "New Hampshire", "Massachusettes", "Deleware", "Maryland", "Virginia",
                                  "West Virginia", "North Carolina",
                                  "New Jersey", "Vermont", "New York", "Pennsylvania", "Connecticut"))

cdat = read_csv("data/processed/oceans/north_atlantic_2017.csv")
cdat$month = month(cdat$time)
cdat$month = month.abb[cdat$month]
cdat$month = factor(cdat$month, levels=month.abb)

cdat = filter(cdat, MCLASS != 0 & month == "Jun")
cdat$true = ifelse(cdat$MCLASS == cdat$VCLASS, "True", "False")

cdat$MCLASS = as.factor(cdat$MCLASS)

ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y, fill=factor(MCLASS))) +
  ggtitle("June 2017 North Atlantic Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  scale_fill_viridis_d(option = "G") +
  theme_minimal() +
  xlim(-80, -60) +
  ylim(35, 45) +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow=2))


ggsave("figures/May2023_poster_figures/north_atlantic_seascape_map_2017.png", width=8, height=8)


ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=ccdat, aes(x, y)) +
  ggtitle("June 2017 North Atlantic Sea Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  # scale_fill_viridis_d(option = "G") +
  xlim(-80, -60) +
  ylim(35, 45) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal(14) +
  # facet_wrap(~month) +
  NULL

ggsave("figures/May2023_poster_figures/north_atlantic_seascape_map_2017_blank.png", width=8, height=8)



# -------------------------------------------
# Pacific Northwest

ssdat = read_sf("~/Projects/test/data/cb_2018_us_states_500k/cb_2018_us_state_500k.shp")
st_crs(ssdat) <- "EPSG:4326" 

ssdat = filter(ssdat, NAME %in% c("Washington", "Oregon"))

cdat = read_csv("data/processed/oceans/pnw_coast_2017.csv")
cdat$month = month(cdat$time)
cdat$month = month.abb[cdat$month]
cdat$month = factor(cdat$month, levels=month.abb)

cdat = filter(cdat, MCLASS != 0 & month == "Jun")
cdat$true = ifelse(cdat$MCLASS == cdat$VCLASS, "True", "False")

cdat$MCLASS = as.factor(cdat$MCLASS)

ccdat = filter(cdat, x <= -60)
ccdat = filter(ccdat, y >= 35)

ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y, fill=factor(MCLASS))) +
  ggtitle("June 2017 Pacific Northwest Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  scale_fill_viridis_d(option = "G") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.direction = "horizontal",
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow=1))


ggsave("figures/May2023_poster_figures/pnw_seascape_map_2017.png", width=8, height=8)


ggplot(NULL) + 
  geom_sf(data=ssdat, fill=NA) + 
  geom_tile(data=cdat, aes(x, y)) +
  ggtitle("June 2017 Pacific Northwest Sea Seascapes") + 
  labs(x=NULL, y=NULL, fill=NULL) + 
  # scale_fill_viridis_d(option = "G") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal(14) +
  # facet_wrap(~month) +
  NULL

ggsave("figures/May2023_poster_figures/pnw_seascape_map_2017_blank.png", width=8, height=8)






