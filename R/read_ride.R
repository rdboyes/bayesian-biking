library(sf)
library(xml2)
library(zoo)
library(dplyr)

ride_defaults <- list(
  CdA   = 0.40,  # m^2, drag area
  Crr   = 0.005, # rolling resistance coefficient
  rho   = 1.225, # kg/m^3, air density
  eta   = 0.97,  # drivetrain efficiency
  g     = 9.81,
  bin_s = 10     # seconds per bin for the model
)

read_weights <- function(path = file.path("output", "weight_daily.csv")) {
  w <- read.csv(path)
  w$date <- as.Date(w$date)
  w
}

# Rider mass on a date: linear interpolation between daily weigh-ins, held flat
# beyond either end. Returns kg and the distance (days) to the nearest weigh-in.
rider_mass <- function(date, weight) {
  list(
    kg  = approx(weight$date, weight$kg, xout = date, rule = 2)$y,
    gap = min(abs(as.numeric(weight$date - date)))
  )
}

# Read one GPX ride, estimate power from GPS, and bin it for the model.
read_ride <- function(gpx, weight, bike_kg = 10, cfg = ride_defaults) {
  # ---- read track points with sf --------------------------------------------
  pts <- st_read(gpx, layer = "track_points", quiet = TRUE)[, c("ele", "time")]

  # GDAL's GPX driver doesn't expose the nested <gpxtpx:hr> value, so pull HR
  # from the XML directly; trkpt order matches the sf rows.
  trkpt <- xml_find_all(read_xml(gpx), "//*[local-name()='trkpt']")
  hr <- as.numeric(xml_text(xml_find_first(trkpt, ".//*[local-name()='hr']")))
  stopifnot(length(hr) == nrow(pts))
  pts$hr <- hr

  ride_date <- as.Date(pts$time[1])
  mass <- rider_mass(ride_date, weight)

  # ---- distance, speed, grade, acceleration ---------------------------------
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

  # ---- estimated power (physics model) --------------------------------------
  m <- mass$kg + bike_kg
  ride <- ride |>
    mutate(
      theta = atan(grade),
      P = (m * cfg$g * v * (cfg$Crr * cos(theta) + sin(theta)) +
             m * acc * v +
             0.5 * cfg$rho * cfg$CdA * v^3) / cfg$eta,
      P = ifelse(v < 1, 0, P)
    )
  # P is left unfloored per point: GPS noise in acceleration swings it both
  # ways, and clipping each point at 0 would bias the mean upward. Clip after
  # binning.

  # ---- bin for modelling ----------------------------------------------------
  bins <- ride |>
    mutate(bin = floor(t / cfg$bin_s)) |>
    group_by(bin) |>
    summarise(t = mean(t), P = mean(P), hr = mean(hr, na.rm = TRUE),
              v = mean(v), .groups = "drop") |>
    filter(!is.na(hr)) |>
    mutate(P = pmax(P, 0), dt = c(cfg$bin_s, diff(t)))

  list(
    bins = bins, ride = ride,
    info = list(date = ride_date, start = pts$time[1], n_points = n,
                hr_missing = sum(is.na(hr)), rider_kg = mass$kg,
                weight_gap_days = mass$gap, bike_kg = bike_kg),
    config = cfg
  )
}
