library(dplyr)

# Extract body-mass records from the Apple Health export into a daily table.
# export.xml is ~600 MB, so stream it in chunks and keep only matching lines
# rather than parsing the full XML tree.

src <- file.path("data", "health", "apple_health_export", "export.xml")
dir.create("output", showWarnings = FALSE)

con <- file(src, open = "r", encoding = "UTF-8")
hits <- character()
repeat {
  chunk <- readLines(con, n = 200000, warn = FALSE)
  if (length(chunk) == 0) break
  hits <- c(hits, grep('type="HKQuantityTypeIdentifierBodyMass"', chunk,
                       value = TRUE, fixed = TRUE))
}
close(con)

attr_val <- function(x, name) {
  sub(paste0('.*\\b', name, '="([^"]*)".*'), "\\1", x)
}

weights <- data.frame(
  source = attr_val(hits, "sourceName"),
  unit   = attr_val(hits, "unit"),
  start  = attr_val(hits, "startDate"),
  value  = as.numeric(attr_val(hits, "value"))
) |>
  mutate(
    kg   = case_when(unit == "lb" ~ value * 0.45359237,
                     unit == "kg" ~ value,
                     unit == "g"  ~ value / 1000),
    date = as.Date(substr(start, 1, 10))
  )
stopifnot(!anyNA(weights$kg))

# One value per day (median handles same-day duplicates across apps).
daily <- weights |>
  group_by(date) |>
  summarise(kg = median(kg), n = n(), .groups = "drop") |>
  arrange(date)

cat(sprintf("Body-mass records: %d -> %d days (%s to %s), %.1f-%.1f kg\n",
            nrow(weights), nrow(daily), min(daily$date), max(daily$date),
            min(daily$kg), max(daily$kg)))

write.csv(daily, file.path("output", "weight_daily.csv"), row.names = FALSE)
