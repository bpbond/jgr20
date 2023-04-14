# Estimates of global Rh and Rs based on ResCOM spreadsheet

library(readr)
library(ggplot2)
theme_set(theme_bw())

dat <- read_csv("ResCom_ Global Respiration estimates and trends - Sheet1.csv", skip = 1)

# Extract "xxxx (yyyy)" into a year column
year_locs <- regexpr("[0-9]{4}", dat$Citation)
years <- sapply(seq_along(year_locs), function(x) {
  if(isTRUE(year_locs[x] > 0)) {
    substr(dat$Citation[x], year_locs[x], year_locs[x] + attr(year_locs, "match.length")[x] - 1)
  } else {
    NA
  }
})
dat$Year <- as.integer(years)

# Handle flux entries that are "x to y"
dat$Flux_mean <- as.numeric(dat$`Mean or median\n(Pg C yr-1)`)
dat$Flux_sd <- as.numeric(dat$`Std\n(PgC yr-1)`)

have_to <- grep("to", dat$`Mean or median\n(Pg C yr-1)`)
split_vals <- strsplit(dat$`Mean or median\n(Pg C yr-1)`[have_to], " to ")
lows <- sapply(split_vals, function(x) as.numeric(x[1]))
highs <- sapply(split_vals, function(x) as.numeric(x[2]))

dat$Flux_mean[have_to] <- (lows + highs) / 2
dat$Flux_sd[have_to] <- (highs - lows) / 2

# Plots

dat <- subset(dat, !is.na(Flux_mean))

rs <- subset(dat, Flux == "Total Rs")

ggplot(rs, aes(Year, Flux_mean)) +  geom_linerange(aes(ymin = Flux_mean - Flux_sd, ymax = Flux_mean + Flux_sd), color= "grey") + geom_point(aes(color = `Model type`))

rh <- subset(dat, Flux == "Belowground Rh")

ggplot(rh, aes(Year, Flux_mean)) +  geom_linerange(aes(ymin = Flux_mean - Flux_sd, ymax = Flux_mean + Flux_sd), color= "grey") + geom_point(aes(color = `Model type`))

