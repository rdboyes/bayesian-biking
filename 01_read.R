library(ggplot2)
source(file.path("R", "read_ride.R"))

# Read and inspect a single ride. See 03_read_all.R for every GPX ride.

ride_id <- "5835787118"
bike_kg <- 10

gpx <- file.path("data", "activities", paste0(ride_id, ".gpx"))
dir.create("output", showWarnings = FALSE)

r <- read_ride(gpx, read_weights(), bike_kg = bike_kg)
bins <- r$bins

if (r$info$weight_gap_days > 30) {
  warning(sprintf("Nearest weigh-in is %d days from the ride", r$info$weight_gap_days))
}
cat(sprintf("Ride date %s: rider mass %.1f kg (nearest weigh-in %d days away)\n",
            r$info$date, r$info$rider_kg, r$info$weight_gap_days))

# ---- sanity check vs Strava ---------------------------------------------------
act <- read.csv("data/activities.csv", check.names = FALSE)
strava <- act[act[["Activity ID"]] == as.numeric(ride_id), ]
cat(sprintf("Points: %d (HR missing: %d)  Bins: %d\n",
            r$info$n_points, r$info$hr_missing, nrow(bins)))
cat(sprintf("Mean est. power (moving): %.0f W | Strava avg watts: %s W\n",
            mean(bins$P[bins$v >= 1]), strava[["Average Watts"]][1]))
cat(sprintf("Mean HR: %.0f bpm | Strava avg HR: %s bpm\n",
            mean(r$ride$hr, na.rm = TRUE), strava[["Average Heart Rate"]][1]))

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
  list(bins = bins, ride = r$ride,
       config = c(list(ride_id = ride_id, rider_kg = r$info$rider_kg,
                       bike_kg = bike_kg), r$config)),
  file.path("output", paste0("ride_", ride_id, ".rds"))
)
