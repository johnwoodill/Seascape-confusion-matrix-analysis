dat = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/test.csv")
dat$time = as.Date(dat$time)
pdat = dplyr::filter(dat, time == as.Date("2017-06-15"))

pdat$diff = ifelse(pdat$diff == 0, 1, 0)

ggplot(pdat, aes(lon, lat, color=diff)) + geom_tile()


[-100, 20, -80, 30]

test = filter(pdat, lon >= -100 & lon <= -81.5)
test = filter(test, lat >= 18 & lat <= 30)

ggplot(test, aes(lon, lat, fill=factor(diff))) + geom_sf(data=states, inherit.aes = FALSE) + geom_tile()


coasts = read_sf("~/Projects/Seascape-confusion-matrix-analysis/data/shapefiles/gulf_of_mexico/World_Seas_IHO_v3.shp")
ggplot(coasts) + geom_sf()



test = read_csv("~/Projects/Seascape-confusion-matrix-analysis/data/test_gom.csv")
test$time = as.Date(test$time)

test = filter(test, VCLASS != 0)

test = filter(test, time == "2017-06-15")
test$diff = ifelse(test$diff == 0, 1, 0)
ggplot(test) + geom_sf(data=states, inherit.aes = FALSE) + geom_tile(aes(x, y, fill=factor(diff)))
