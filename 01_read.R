library(sf)
library(xml2)
library(zoo)
library(dplyr)
library(ggplot2)

# ---- config -----------------------------------------------------------------
ride_id  <- "5835787118"
bike_kg  <- 10
CdA      <- 0.40  # m^2, drag area
Crr      <- 0.005 # rolling resistance coefficient
rho      <- 1.225 # kg/m^3, air density
eta      <- 0.97  # drivetrain efficiency
g        <- 9.81
bin_s    <- 10    # seconds per bin for the model

gpx <- file.path("data", "activities", paste0(ride_id, ".gpx"))
dir.create("output", showWarnings = FALSE)

# ---- read track points with sf ----------------------------------------------
pts <- st_read(gpx, layer = "track_points", quiet = TRUE)[, c("ele", "time")]

# GDAL's GPX driver doesn't expose the nested <gpxtpx:hr> value, so pull HR
# from the XML directly; trkpt order matches the sf rows.
trkpt <- xml_find_all(read_xml(gpx), "//*[local-name()='trkpt']")
hr <- as.numeric(xml_text(xml_find_first(trkpt, ".//*[local-name()='hr']")))
stopifnot(length(hr) == nrow(pts))
pts$hr <- hr

# ---- rider mass on the ride date (from 00_weight.R) ---------------------------
# Linear interpolation between daily weigh-ins, held flat beyond either end.
weight <- read.csv(file.path("output", "weight_daily.csv"))
weight$date <- as.Date(weight$date)
ride_date <- as.Date(pts$time[1])
rider_kg <- approx(weight$date, weight$kg, xout = ride_date, rule = 2)$y
gap_days <- min(abs(as.numeric(weight$date - ride_date)))
if (gap_days > 30) {
  warning(sprintf("Nearest weigh-in is %d days from the ride", gap_days))
}
cat(sprintf("Ride date %s: rider mass %.1f kg (nearest weigh-in %d days away)\n",
            ride_date, rider_kg, gap_days))

# ---- distance, speed, grade, acceleration -----------------------------------
n <- nrow(pts)
step_m <- c(0, as.numeric(st_distance(pts[-n, ], pts[-1, ], by_element = TRUE)))

ride <- st_drop_geometry(pts) |>
  mutate(
    t    = as.numeric(difftime(time, time[1], units = "secs")),
    dist = cumsum(step_m)
  ) |>
  filter(!duplicated(t))

k <- 5 # half-window (points) for grade
ride <- ride |>
  mutate(
    ele_s = rollmean(ele, 11, fill = "extend"),
    dt    = c(NA, diff(t)),
    v_raw = c(0, diff(dist) / diff(t)),
    # median then mean filter: removes GPS spikes before differencing for acc
    v     = rollmean(rollmedian(v_raw, 5, fill = "extend"), 9, fill = "extend"),
    d_run = lead(dist, k) - lag(dist, k),
    grade = (lead(ele_s, k) - lag(ele_s, k)) / d_run,
    grade = ifelse(is.na(grade) | d_run < 5, 0, pmin(pmax(grade, -0.2), 0.2)),
    acc   = pmin(pmax(c(0, diff(v) / diff(t)), -2), 2)
  )

# ---- estimated power (physics model) ----------------------------------------
m <- rider_kg + bike_kg
ride <- ride |>
  mutate(
    theta = atan(grade),
    P = (m * g * v * (Crr * cos(theta) + sin(theta)) +
           m * acc * v +
           0.5 * rho * CdA * v^3) / eta,
    P = ifelse(v < 1, 0, P)
  )
# P is left unfloored per point: GPS noise in acceleration swings it both ways,
# and clipping each point at 0 would bias the mean upward. Clip after binning.

# ---- bin for modelling -------------------------------------------------------
bins <- ride |>
  mutate(bin = floor(t / bin_s)) |>
  group_by(bin) |>
  summarise(t = mean(t), P = mean(P), hr = mean(hr, na.rm = TRUE),
            v = mean(v), .groups = "drop") |>
  filter(!is.na(hr)) |>
  mutate(P = pmax(P, 0), dt = c(bin_s, diff(t)))

# ---- sanity check vs Strava ---------------------------------------------------
act <- read.csv("data/activities.csv", check.names = FALSE)
strava <- act[act[["Activity ID"]] == as.numeric(ride_id), ]
cat(sprintf("Points: %d (HR missing: %d)  Bins: %d\n", n, sum(is.na(hr)), nrow(bins)))
cat(sprintf("Mean est. power (moving): %.0f W | Strava avg watts: %s W\n",
            mean(bins$P[bins$v >= 1]), strava[["Average Watts"]][1]))
cat(sprintf("Mean HR: %.0f bpm | Strava avg HR: %s bpm\n",
            mean(ride$hr, na.rm = TRUE), strava[["Average Heart Rate"]][1]))

# ---- plot and save -----------------------------------------------------------
long <- bind_rows(
  data.frame(t = bins$t, value = bins$hr, series = "Heart rate (bpm)"),
  data.frame(t = bins$t, value = bins$P,  series = "Estimated power (W)")
)
p <- ggplot(long, aes(t / 60, value)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~series, ncol = 1, scales = "free_y") +
  labs(x = "Elapsed time (min)", y = NULL, title = paste("Ride", ride_id)) +
  theme_minimal()
ggsave(file.path("output", paste0("ride_", ride_id, "_hr_power.png")), p,
       width = 9, height = 5, dpi = 120)

saveRDS(
  list(bins = bins, ride = ride,
       config = list(ride_id = ride_id, rider_kg = rider_kg, bike_kg = bike_kg,
                     CdA = CdA, Crr = Crr, rho = rho, eta = eta, bin_s = bin_s)),
  file.path("output", paste0("ride_", ride_id, ".rds"))
)
