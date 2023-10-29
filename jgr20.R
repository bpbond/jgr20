# Estimates of global Rh and Rs based on ResCOM spreadsheet
# Figure for the JGR 20-year review paper
# Ben Bond-Lamberty July 2023

library(readr)
library(ggplot2)
theme_set(theme_bw())
library(tibble)
library(dplyr)
library(patchwork)
library(mblm)

dat <- read_csv("ResCom_ Global Respiration estimates and trends - Sheet1.csv", skip = 1)

# Extract "xxxx (yyyy)" into a year column
# First isolate the four-digit numbers in the strings; regexpr returns
# a starting position and, as an attribute, the length of the match (always 4)
year_locs <- regexpr("[0-9]{4}", dat$Citation)
years <- sapply(seq_along(year_locs), function(x) {
  if(isTRUE(year_locs[x] > 0)) {
    substr(dat$Citation[x], year_locs[x], year_locs[x] + 3) # 3 = 4-1
  } else {
    NA_character_
  }
})
dat$Year <- as.integer(years)

# Create easier-to-work-with Flux_mean and Flux_sd columns
dat$Flux_mean <- as.numeric(dat$`Mean or median\n(Pg C yr-1)`)
dat$Flux_sd <- as.numeric(dat$`Std\n(PgC yr-1)`)

# Handle flux entries that are "x to y"
have_to <- grep("to", dat$`Mean or median\n(Pg C yr-1)`)
# For those entries, split the strings and convert each to numeric
split_vals <- strsplit(dat$`Mean or median\n(Pg C yr-1)`[have_to], " to ")
lows <- sapply(split_vals, function(x) as.numeric(x[1]))
highs <- sapply(split_vals, function(x) as.numeric(x[2]))
# The flux mean and sd are calculated based on the provided range
dat$Flux_mean[have_to] <- (lows + highs) / 2 # average
dat$Flux_sd[have_to] <- (highs - lows) / 2   # half the range 

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

p <- ggplot(rs, aes(Year, Flux_mean)) +  
  geom_pointrange(aes(ymin = Flux_mean - Flux_sd, ymax = Flux_mean + Flux_sd, color = `Model type`),
                  position = position_dodge(width = 1)) + 
  ylab(expression(Flux~(Pg~C~yr^-1))) +
  facet_grid(Flux ~ ., scales = "free_y") +
  coord_cartesian(xlim = c(MINYEAR, 2023)) +
  scale_color_discrete("") +
  geom_text(data = label_dat, aes(label = lab), parse = TRUE, fontface = "bold") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

print(p)

#  geom_segment(data = recent_rs, x = 2004, xend = 2023, aes(yend = Flux_mean), linetype = 2)
ggsave("resp.png", width = 8, height = 7)



## plot length of observations needed to observe climate-driven change at a single site
srdb <- read_csv("srdb-data.csv")
srdb %>% 
  select(YearsOfData) %>% 
  filter(!is.na(YearsOfData)) %>% 
  ggplot(aes(YearsOfData)) +
  geom_histogram(col = "black") +
  annotate("text", 1.5, 11100, label = paste("(a)")) +
  ylim(0, 11900)+
  labs(x = "Record length (year)", y = "Count (n)") -> length_dens


# when would expect to see significance? toy-example
# prepare some data first

median_error <- 0.199  # this is Rs measuring error, this value is from COSORE

# set.seed(1234)

trend_emergence <- function(rd, theilsen = FALSE) {
  Year <- seq_len(length(rd))
  trend_p <- rep(NA, length(rd))
  for(i in seq_along(trend_p)) {
    if(i > 2) {
      if(theilsen) {
        df <- tibble(Year = Year[1:i], rd = rd[1:i])
        suppressWarnings(m <- mblm::mblm(rd ~ Year, data = df))  # mblm doesn't like form below
      } else {
        m <- suppressWarnings(lm(rd[1:i] ~ Year[1:i]))
      }
      # Extract 2nd row (Year) and 4th column (Pr>[t] or Pr>|V|)
      trend_p[i] <- summary(m)$coefficients[2, 4]
    }
  }
  trend_p
}

fuzz <- function(x, error) {
  x * rnorm(length(x), mean = 1, sd = error)
}

# Temperature has risen 0.9 C in 40 years, more or less
dTdt <- round(0.9 / 40.0, 3)
q10 <- 2
R0 = 1.0
respdata <- tibble(Year = 1:150,
                   Temp = dTdt * Year,
                   Resp = R0 * q10 ^ (Temp / 10),
                   # This is interannual variability: from SRDB
                   Resp_iav = fuzz(Resp, 0.098),  # this is SRDB Rs_interannual_err
                   Resp_fuzz = fuzz(Resp_iav, median_error))

do_sim <- function(i, respdata, error = 0.0) {
  # This is observational error
  respdata$Resp_fuzz <- fuzz(respdata$Resp_iav, error)
  respdata$trend_p <- trend_emergence(respdata$Resp_fuzz)
  respdata
}

# run the analysis and store the results
results <- list()
library(parallel)
n_sims <- 150
results <- mclapply(seq_len(n_sims), do_sim, respdata, error = median_error)

# summarise results
results %>% 
  bind_rows %>% 
  group_by(Year) %>% 
  summarise(n = n(), 
            Temp = mean(Temp), 
            Resp = mean(Resp),
            Resp_iav_sd = sd(Resp_iav),
            Resp_iav = mean(Resp_iav),
            Resp_fuzz_sd = sd(Resp_fuzz),          
            Resp_fuzz = mean(Resp_fuzz), 
            trend_p_sd = sd(trend_p), 
            trend_p = mean(trend_p)) %>% 
  filter(!is.na(trend_p)) ->
  results_summary

# plot the trend analysis result
p_TheilSen <- ggplot(results_summary, aes(Year, Resp_fuzz, color = trend_p < 0.05)) +
  geom_point() +
  geom_line(aes(y = Resp), color = "grey") +
  geom_ribbon(aes(ymin = Resp_fuzz - Resp_fuzz_sd, 
                  ymax = Resp_fuzz + Resp_fuzz_sd, 
                  fill = trend_p < 0.05), color = NA, alpha = I(0.35)) +
  guides(color = FALSE, fill = FALSE) +
  xlab("Years of observations") +
  ylab("Soil respiration (normalized)")

print(p_TheilSen)

length_dens / p_TheilSen

ggsave("dens_length_obs.png", width = 8, height = 7)

