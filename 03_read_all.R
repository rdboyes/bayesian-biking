source(file.path("R", "read_ride.R"))

# Read every outdoor GPX ride with heart rate, estimate power, and stack the
# binned data for the hierarchical model.

min_bins       <- 60    # drop rides shorter than ~10 minutes of HR data
max_hr_missing <- 0.2   # drop rides missing HR on more than 20% of points

# Physics per bike. The Trek 520 is a touring bike ridden upright on wide tyres,
# so it gets more drag and rolling resistance than the Miele road bike.
# Rides with no gear recorded use the defaults (road bike).
bike_cfg <- list(
  "Trek 520" = modifyList(ride_defaults, list(CdA = 0.55, Crr = 0.008))
)

act <- read.csv("data/activities.csv", check.names = FALSE)
names(act) <- make.unique(names(act))
weight <- read_weights()

cand <- act |>
  filter(`Activity Type` == "Ride",
         grepl("\\.gpx$", Filename),
         !is.na(`Average Heart Rate`)) |>
  transmute(
    ride_id       = as.character(`Activity ID`),
    file          = file.path("data", Filename),
    name          = `Activity Name`,
    bike          = ifelse(`Activity Gear` == "", "Miele", `Activity Gear`),
    bike_kg       = ifelse(is.na(`Bike Weight`), 10, `Bike Weight`),
    strava_watts  = `Average Watts`,
    strava_hr     = `Average Heart Rate`
  )
cat(sprintf("Candidate GPX rides with HR: %d\n", nrow(cand)))

results <- vector("list", nrow(cand))
for (i in seq_len(nrow(cand))) {
  c_i <- cand[i, ]
  results[[i]] <- tryCatch({
    cfg <- if (c_i$bike %in% names(bike_cfg)) bike_cfg[[c_i$bike]] else ride_defaults
    r <- read_ride(c_i$file, weight, bike_kg = c_i$bike_kg, cfg = cfg)
    b <- r$bins |> filter(hr >= 40, hr <= 220)
    list(
      bins = mutate(b, ride_id = c_i$ride_id),
      info = data.frame(
        ride_id = c_i$ride_id, name = c_i$name, bike = c_i$bike,
        date = r$info$date,
        n_bins = nrow(b), hours = max(b$t) / 3600,
        hr_missing_frac = r$info$hr_missing / r$info$n_points,
        rider_kg = r$info$rider_kg, weight_gap_days = r$info$weight_gap_days,
        est_watts = mean(b$P[b$v >= 1]), strava_watts = c_i$strava_watts,
        mean_hr = mean(b$hr), strava_hr = c_i$strava_hr,
        error = NA_character_
      )
    )
  }, error = function(e) {
    list(bins = NULL,
         info = data.frame(ride_id = c_i$ride_id, name = c_i$name,
                           error = conditionMessage(e)))
  })
  if (i %% 20 == 0) cat(sprintf("  read %d / %d\n", i, nrow(cand)))
}

info <- bind_rows(lapply(results, `[[`, "info")) |>
  mutate(keep = is.na(error) & n_bins >= min_bins &
           hr_missing_frac <= max_hr_missing)

cat(sprintf("Read errors: %d | too short: %d | too much HR missing: %d | kept: %d\n",
            sum(!is.na(info$error)),
            sum(is.na(info$error) & info$n_bins < min_bins),
            sum(is.na(info$error) & info$hr_missing_frac > max_hr_missing),
            sum(info$keep)))
if (any(!is.na(info$error))) print(info[!is.na(info$error), c("ride_id", "error")])

keep_ids <- info$ride_id[info$keep]
bins <- bind_rows(lapply(results, `[[`, "bins")) |>
  filter(ride_id %in% keep_ids)

kept <- info |> filter(keep)
cat(sprintf("Kept rides: %s to %s, %.0f hours, %d bins\n",
            min(kept$date), max(kept$date), sum(kept$hours), nrow(bins)))
cat("Est. vs Strava watts by bike:\n")
print(kept |>
  group_by(bike) |>
  summarise(rides = n(),
            ratio_median = median(est_watts / strava_watts, na.rm = TRUE),
            ratio_q25 = quantile(est_watts / strava_watts, 0.25, na.rm = TRUE),
            ratio_q75 = quantile(est_watts / strava_watts, 0.75, na.rm = TRUE)))
cat(sprintf("Rides with nearest weigh-in > 30 days away: %d\n",
            sum(kept$weight_gap_days > 30)))

saveRDS(list(bins = bins, info = info,
             config = list(default = ride_defaults, bike = bike_cfg)),
        file.path("output", "all_rides.rds"))
