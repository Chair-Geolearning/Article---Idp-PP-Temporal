selected_vars <- c(
  "Temperature", "Humidity", "Windspeed",
  "RG", "VPD", "FWI", "API"
)

reshape_to_long <- function(df, name_region) {
  df %>%
    dplyr::mutate(
      coord_id = as.integer(factor(paste(lon, lat))),
      nuts_id  = name_region,
      date     = Date
    ) %>%
    dplyr::select(lon, lat, date, coord_id, nuts_id, all_of(selected_vars)) %>%
    tidyr::pivot_longer(
      cols = all_of(selected_vars),
      names_to = "variable",
      values_to = "value"
    )
}

compute_extremes <- function(df_long, quantile_map) {
  
  df_long %>%
    left_join(quantile_map, by = "variable") %>%
    group_by(coord_id, variable) %>%
    mutate(
      threshold = quantile(value, unique(prob), na.rm = TRUE),
      is_extreme = dplyr::case_when(
        direction == "high" ~ value > threshold,
        direction == "low"  ~ value < threshold,
        TRUE ~ FALSE
      )
    ) %>%
    ungroup() %>%
    filter(is_extreme) %>%
    split(.$variable)
}


group_single_variable <- function(df_var, temporal_width, distance_threshold) {
  df_var <- df_var %>% arrange(date)
  df_var$event_id <- cumsum(
    c(1, diff(as.numeric(df_var$date)) > temporal_width)
  )
  df_var
}

compute_API <- function(df,
                        precip_var = "Precipitation",
                        date_var   = "Date",
                        coord_vars = c("lon", "lat"),
                        m = 14,
                        k = 0.9) {
  
  df %>%
    arrange(across(all_of(c(coord_vars, date_var)))) %>%
    group_by(across(all_of(coord_vars))) %>%
    mutate(
      API = purrr::map_dbl(
        seq_len(n()),
        ~ {
          idx <- .x
          if (idx <= 1) return(NA_real_)
          
          i_max <- min(m, idx - 1)
          sum(
            k^(0:(i_max - 1)) *
              .data[[precip_var]][(idx - 1):(idx - i_max)],
            na.rm = TRUE
          )
        }
      )
    ) %>%
    ungroup()
}
quantile_map <- tibble::tribble(
  ~variable,        ~prob,  ~direction,
  "FWI",             0.99,  "high",
  "Temperature",     0.99,  "high",
  "Windspeed",       0.99,  "high",
  "RG",              0.99,  "high",
  "VPD",             0.99,  "high",
  "Humidity",        0.01,  "low",
  "API",             0.01,  "low"
)

