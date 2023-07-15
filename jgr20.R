# Estimates of global Rh and Rs based on ResCOM spreadsheet
# Figure for the JGR 20-year review paper
# Ben Bond-Lamberty July 2023

library(readr)
library(ggplot2)
library(tibble)
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

# Replace missing values with median of present values
impute_median <- function(x) {
  x[is.na(x)] <- median(x[!is.na(x)])
  x
}

# Plots

dat <- dat[!is.na(dat$Flux_mean),]

rs <- dat[dat$Flux %in% c("Belowground Rh", "Total Rs"),]
rs$Flux[rs$Flux != "Total Rs"] <- "Rh"
rs$Flux[rs$Flux == "Total Rs"] <- "Rs"

# We cut off one very early value
MINYEAR <- 1992
print(rs[rs$Year < MINYEAR,])

label_dat <- tibble(Flux = c("Rh", "Rs"), lab = c("(a)~R[H]", "(b)~R[S]"), 
                    # Position where these labels should appear
                    Year = MINYEAR + 2,
                    Flux_mean = c(72, 167))

# Compute the median of last years of estimates
recent_rs <- aggregate(Flux_mean ~ Flux, FUN = median, data = rs[rs$Year > 2003,])
print(recent_rs)

ggplot(rs, aes(Year, Flux_mean)) +  
  geom_linerange(aes(ymin = Flux_mean - Flux_sd, ymax = Flux_mean + Flux_sd), color= "grey") + 
  geom_point(aes(color = `Model type`)) +
  ylab(expression(Flux~(Pg~C~yr^-1))) +
  facet_grid(Flux ~ ., scales = "free_y") +
  coord_cartesian(xlim = c(MINYEAR, 2023)) +
  scale_color_discrete("") +
  geom_text(data = label_dat, aes(label = lab), parse = TRUE, fontface = "bold") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
#  geom_segment(data = recent_rs, x = 2004, xend = 2023, aes(yend = Flux_mean), linetype = 2)
ggsave("resp.png", width = 8, height = 7)

