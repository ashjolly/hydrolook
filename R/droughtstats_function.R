# Function for assembling drought statistics and data
# Directly from Tidyhydat data
# 11Aug2020 Ashlee Jollymore
# =====================



## Calculate percentiles from tidyhydat extracted data
## Limited to current day calculation
calc_percentiles <- function(historical_flow, realtime_data, expected,...) {

  if(class(realtime_data$Date) != "Date") stop("Date column is not in date class")

  df <- historical_flow %>%
    dplyr::filter(yday(Date) == yday(Sys.Date())) %>%
    dplyr::group_by(STATION_NUMBER) %>%
    tidyr::nest() %>%
    dplyr::left_join(realtime_data, by = c("STATION_NUMBER")) %>%
    dplyr::mutate(prctile = map2_dbl(data, Value, ~ecdf(.x$Value)(.y))) %>%
    dplyr::left_join(allstations, by = c("STATION_NUMBER")) %>%
    dplyr::mutate(pct_bin = case_when(
      is.na(prctile) ~ "Not ranked",
      prctile >= 0 & prctile <= 0.01 ~ "Low",
      prctile > 0.01 & prctile <= 0.10 ~ "Much below normal (<10)",
      prctile > 0.10 & prctile <= 0.24 ~ "Below Normal (10-24)",
      prctile > 0.24 & prctile <= 0.75 ~ "Normal (25-75)",
      prctile > 0.75 & prctile <= 0.90 ~ "Above normal (76-90)",
      prctile > 0.90 & prctile < 1 ~ "Much above normal (>90)",
      prctile == 1 ~ "High"
    )) %>%
    dplyr::mutate(pct_bin = factor(pct_bin, levels = expected)) %>%
    dplyr::mutate(prctile = prctile*100)
  #st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
  #         crs= 4326,
  #         agr= "constant") %>%
  #transform_bc_albers(sf_object)

  sp::coordinates(df) <- ~LONGITUDE + LATITUDE
  sp::proj4string(df) <- sp::CRS("+proj=longlat +datum=NAD83")

  spatial = st_as_sf(df, coords = c("LONGITUDE", "LATITUDE"),
                     crs = 4326, agr= "constant")

  spatial_spp <- st_transform(spatial, 4326)
  spatial_spp <- transform_bc_albers(spatial)
}


# Function for calculating the drought statistics
drought_statistics <- function(stations){
 q_stns <- unique(stations) %>%
  hy_stn_data_range() %>%
  filter(DATA_TYPE == "Q") %>%
  filter(RECORD_LENGTH >=20) %>%
  pull(STATION_NUMBER)

 # Query realtime data
 rl_data <- realtime_dd(q_stns)

 ## Find most recent instantaneous discharge value
 rl_data_instant <- rl_data %>%
  dplyr::filter(Parameter == "Flow") %>%
  dplyr::group_by(STATION_NUMBER) %>%
  dplyr::filter(Date == max(Date)) %>%
  dplyr::select(STATION_NUMBER, Date, Value) %>%
  dplyr::mutate(Date = as.Date(Date)) %>%
  dplyr::filter(Date == Sys.Date()) %>% ## drop max values that aren't today
  dplyr::ungroup()

 ## Query historical data
 ## NOTE: Should this be done with rl_data_recent$STATION_NUMBER?
 hist_flow <- hy_daily_flows(q_stns)

 ## Find the average of the last 24 hours
 rl_data_last24 <- rl_data %>%
  dplyr::filter(Parameter == "Flow") %>%
  dplyr::group_by(STATION_NUMBER) %>%
  dplyr::filter(Date >= Sys.time() - 60*60*24) %>% ## all data from last 24 hours
  dplyr::select(STATION_NUMBER, Date, Value) %>%
  dplyr::mutate(Date = Sys.Date()) %>% ## label last twenty four hours as from today
  dplyr::group_by(STATION_NUMBER, Date) %>%
  dplyr::summarise(Value = mean(Value)) %>%
  dplyr::ungroup()

 ## Realtime 7 day average
 ## I took the realtime daily averages
 rl_data_7day_mean <- rl_data %>%
  dplyr::filter(Parameter == "Flow") %>%
  dplyr::mutate(Date = as.Date(Date)) %>%
  dplyr::group_by(STATION_NUMBER, Date) %>%
  dplyr::summarise(DailyMean = mean(Value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STATION_NUMBER) %>%
  dplyr::mutate(Value = zoo::rollapply(DailyMean, 7, align = 'right', fill=NA, mean, na.rm=TRUE, partial = TRUE)) %>%
  dplyr::ungroup()

 # Historic 7 day average
 hist_flow_7day_mean <- hist_flow %>%
  dplyr::group_by(STATION_NUMBER, Date) %>% # some values in the historic record with multime values for one day
  dplyr::summarise(DailyMean = mean(Value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STATION_NUMBER) %>%
  dplyr::mutate(Value = zoo::rollapply(DailyMean, 7, align = 'right', fill=NA, mean, na.rm=TRUE, partial = TRUE)) %>%
  dplyr::ungroup()

 ## Expected NAWW percentile bins
 expected <- c("Not ranked", "Low","Much below normal (<10)",
               "Below Normal (10-24)", "Normal (25-75)",
               "Above normal (76-90)","Much above normal (>90)", "High")

 ## Calculate instantaneous percentiles
 pct_flow_instant <- calc_percentiles(historical_flow = hist_flow,
                                     realtime_data = rl_data_instant,
                                     expected)

 ## Calculate 24 hours percentiles
 pct_flow_last24 <- calc_percentiles(historical_flow = hist_flow,
                                    realtime_data = rl_data_last24,
                                    expected)

 # Calculate 7 day mean percentiles
 pct_flow_7day_mean <- calc_percentiles(hist_flow_7day_mean,
                                       rl_data_7day_mean,
                                       expected)

 num_year_data <- hy_stn_data_range(unique(pct_flow_instant$STATION_NUMBER)) %>%
  filter(DATA_TYPE == "Q") %>%
  select(STATION_NUMBER, RECORD_LENGTH) %>%
  rename(`Record Length` = RECORD_LENGTH)

 ## Grab only the latest flow and merge into one data frame
 pct_flow_instant_tbl_data <- pct_flow_instant %>%
  #st_set_geometry(NULL) %>%
  dplyr::select(STATION_NAME, STATION_NUMBER, Value, prctile) %>%
  dplyr::rename(`%tile-instant` = prctile, `Latest Q` = Value) %>%
  select(-`%tile-instant`) %>% ##remove instant %tile
  st_join(pct_flow_last24 %>%
            #st_set_geometry(NULL) %>%
            dplyr::select(STATION_NUMBER, prctile, pct_bin) %>%
            dplyr::rename(`%tile-last24` = prctile, pct_bin_last24 = pct_bin)) %>%
  st_join(pct_flow_7day_mean %>%
            #st_set_geometry(NULL) %>%
            dplyr::filter(Date == Sys.Date()) %>%
            dplyr::select(STATION_NUMBER, prctile, pct_bin) %>%
            dplyr::rename(`%tile-7day_mean` = prctile, pct_bin_7day = pct_bin)) %>%
  left_join(num_year_data, by = c("STATION_NUMBER")) %>%
  dplyr::arrange(STATION_NUMBER) %>%
  dplyr::select(-STATION_NUMBER.x, -STATION_NUMBER.y, -geometry)

 # Add in regulation status
 regulation <- tidyhydat::hy_stn_regulation(station_number = q_stns) %>%
  dplyr::mutate(regulation = ifelse(REGULATED == "TRUE", "Regulated", "Natural"))

 pct_flow_instant_tbl_data <- full_join(pct_flow_instant_tbl_data, regulation %>% dplyr::select(STATION_NUMBER, regulation))

 ## Add in % median and mean flow
 # Find the mean for this day
 hist_flow_7day_mean_day <- hist_flow_7day_mean %>%
  dplyr::mutate(day_month = paste0(day(Date), "-", month(Date))) %>%
  dplyr::filter(day_month == paste0(day(Sys.Date()), "-", month(Sys.Date()))) %>%
  dplyr::group_by(STATION_NUMBER) %>%
  dplyr::summarize(mean_Q7_forthisdate = round(mean(Value, na.rm = TRUE), digits = 2))

 # current percent of mean
 pct_flow_instant_tbl_data <- full_join(pct_flow_instant_tbl_data, hist_flow_7day_mean_day)

 hist_flow_7day_median_day <- hist_flow_7day_mean %>%
  dplyr::mutate(day_month = paste0(day(Date), "-", month(Date))) %>%
  dplyr::filter(day_month == paste0(day(Sys.Date()), "-", month(Sys.Date()))) %>%
  dplyr::group_by(STATION_NUMBER) %>%
  dplyr::summarize(median_Q7_forthisdate = round(median(Value, na.rm = TRUE), digits = 2))

 # current percent of mean
 pct_flow_instant_tbl_data <- full_join(pct_flow_instant_tbl_data, hist_flow_7day_median_day)

 # Add in watershed area
 description <- hy_stations(station_number = pct_flow_instant_tbl_data$STATION_NUMBER) %>%
  as.list()

 pct_flow_instant_tbl_data$basin_area <- description$DRAINAGE_AREA_GROSS

 # Add in the temperature
 #token = token_ws()
 temp_station <- realtime_ws(station_number = pct_flow_instant_tbl_data$STATION_NUMBER,
                            token = token_ws()) %>%
  dplyr::filter(Code == "TW") %>%
  dplyr::mutate(date_dmy = as.Date(Date)) %>%
  dplyr::filter(date_dmy == as.Date(Sys.Date())) %>% # filter for only today's date
  dplyr::group_by(STATION_NUMBER, date_dmy) %>%
  dplyr::summarise(mean_daily_temp = round(mean(Value, na.rm = TRUE), digits = 2))

 pct_flow_instant_tbl_data <- full_join(pct_flow_instant_tbl_data, temp_station)
 #unique(test_1$Code)

 # Add in value of Q7
 Q7 <- pct_flow_7day_mean %>%
  dplyr::filter(Date == Sys.Date()) %>%
  dplyr::select(STATION_NUMBER, Date, Value, prctile, STATION_NAME) %>%
  dplyr::rename(Q7_value = Value, Q7_prctile = prctile)

 pct_flow_instant_tbl_data <- full_join(pct_flow_instant_tbl_data, data.frame(Q7))

 # Calculate the percent of median and mean Q7
 pct_flow_instant_tbl_data <- pct_flow_instant_tbl_data %>%
  dplyr::mutate(Per_Q7_median = round(Q7_value/median_Q7_forthisdate*100, digits = 0)) %>%
  dplyr::mutate(Per_Q7_mean = round(Q7_value/mean_Q7_forthisdate*100, digits = 0))

 # Add in the Q min 7 day value for today's date
 hist_flow_7day_mean_day <- hist_flow_7day_mean %>%
  dplyr::mutate(day_month = paste0(day(Date), "-", month(Date))) %>%
  dplyr::filter(day_month == paste0(day(Sys.Date()), "-", month(Sys.Date()))) %>%
  dplyr::group_by(STATION_NUMBER) %>%
  summarize(min_Q7 = round(min(Value, na.rm = TRUE), digits = 2))

 pct_flow_instant_tbl_data <- full_join(pct_flow_instant_tbl_data, hist_flow_7day_mean_day)

 # Add in MAD
 long_term_mad <- hy_annual_stats(q_stns) %>%
  filter(Sum_stat == "MEAN") %>%
  spread(Sum_stat, Value) %>%
  group_by(STATION_NUMBER) %>%
  summarise(`MEAN MAD (m^3/s)` = round(mean(MEAN, na.rm = TRUE), digits = 2)) %>%
  left_join(rl_data_instant, by = c("STATION_NUMBER")) %>%
  mutate(`% MAD` = round((Value/`MEAN MAD (m^3/s)`)*100, digits = 2)) %>%
  left_join(allstations[,1:2], by = c("STATION_NUMBER")) %>%
  select(STATION_NAME, STATION_NUMBER, Value, `MEAN MAD (m^3/s)`, `% MAD`) %>%
  rename(`Latest Q` = Value) %>%
  arrange(STATION_NUMBER)

 pct_flow_instant_tbl_data <- full_join(pct_flow_instant_tbl_data, long_term_mad)

 # Assemble table of drought - relavent meta data and statistics
 table_out <- pct_flow_instant_tbl_data %>%
  dplyr::select(-geometry,  -date_dmy, -`%tile-7day_mean`) %>%
  dplyr::arrange(`%tile-last24`) %>%
  dplyr::rename(`Last 24 hour Q (m3/s)` = `Latest Q`,
                `Latest 7 Day Q (m3/s)` = Q7_value,
                `Percentile - Last 24 Hour Q` = `%tile-last24`,
                `ID` = STATION_NUMBER,
                `Station Name` = STATION_NAME,
                `Historic Min 7 Day Q (m3/s)` = min_Q7,
                `Mean Daily Water Temp (degC)` = mean_daily_temp,
                `Record Length` = `Record Length`,
                `Percentile - Q7` = Q7_prctile,
                `Historic Mean Q7 for today` = mean_Q7_forthisdate,
                `Historic Median Q7 for today` = median_Q7_forthisdate,
                `Percent of Daily Median Q7 (%)` = Per_Q7_median,
                `Percent of Daily Mean Q7 (%)` = Per_Q7_mean,
                `Basin Area (km2)` = basin_area,
                `Regulation Status` = regulation) %>%
  dplyr::mutate(`Percent of Daily Mean Q7 (%; Historic Mean Q7 in m3/s)` = paste0(`Percent of Daily Mean Q7 (%)`, " (", `Historic Mean Q7 for today`, ")")) %>%
  dplyr::mutate(`Percent of Daily Median Q7 (%; Historic Median Q7 in m3/s)` = paste0(`Percent of Daily Median Q7 (%)`, " (", `Historic Mean Q7 for today`, ")")) %>%
  dplyr::select(-`Percent of Daily Mean Q7 (%)`, -`Historic Mean Q7 for today`, -`Percent of Daily Median Q7 (%)`, -`Historic Median Q7 for today`) %>%
  dplyr::select(`Station Name`, `ID`, `Record Length`, `Basin Area (km2)`, `Regulation Status`,
                `Last 24 hour Q (m3/s)`, `Percentile - Last 24 Hour Q`, pct_bin_last24,
                `Latest 7 Day Q (m3/s)`, `Percentile - Q7`,  pct_bin_7day,
                `Percent of Daily Mean Q7 (%; Historic Mean Q7 in m3/s)`,
                `Percent of Daily Median Q7 (%; Historic Median Q7 in m3/s)`,
                `Historic Min 7 Day Q (m3/s)`,
                `MEAN MAD (m^3/s)`, `% MAD`, `Mean Daily Water Temp (degC)`, `Date`)
}
