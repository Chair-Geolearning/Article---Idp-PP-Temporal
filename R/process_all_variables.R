process_all_variables <- function(years = 1941:2024,
                                  input_dir = "data/processed",
                                  output_dir = "data/extremes",
                                  quantile_precip = 0.90,
                                  quantile_cape = 0.99,
                                  quantile_wind = 0.99) {
  
  cat(sprintf("\n%s\n", strrep("=", 80)))
  cat(sprintf(" CONFIGURATION DE LA DÉTECTION D'EXTRÊMES\n"))
  cat(sprintf("%s\n\n", strrep("=", 80)))
  cat(sprintf("  Période d'analyse: %d-%d (%d ans)\n", 
              min(years), max(years), length(years)))
  cat(sprintf("  Période de référence: TOUTE LA PÉRIODE (1941-2024)\n\n"))
  cat(sprintf("Méthodes par variable:\n"))
  cat(sprintf("  • Précipitations totales: Q%.0f sur jours humides (> 1mm)\n", 
              quantile_precip * 100))
  cat(sprintf("  • Précipitations convectives: Tous les jours positifs (> 0)\n"))
  cat(sprintf("  • CAPE: Q%.0f\n", quantile_cape * 100))
  cat(sprintf("  • Vent: Q%.0f\n", quantile_wind * 100))
  cat(sprintf("\n%s\n\n", strrep("=", 80)))
  
  results <- list()
  
  cat("\n  VARIABLE 1/4: Précipitations totales\n")
  results[["total_precipitation_daily_sum"]] <- process_variable_pipeline(
    var_name = "total_precipitation_daily_sum",
    years = years,
    input_dir = input_dir,
    output_dir = output_dir,
    quantile_threshold = quantile_precip
  )
  gc(verbose = FALSE)
  
  cat("\n  VARIABLE 2/4: Précipitations convectives\n")
  results[["convective_precipitation_daily_sum"]] <- process_variable_pipeline(
    var_name = "convective_precipitation_daily_sum",
    years = years,
    input_dir = input_dir,
    output_dir = output_dir,
    quantile_threshold = 0  
  )
  gc(verbose = FALSE)
  
  cat("\n VARIABLE 3/4: CAPE\n")
  results[["convective_available_potential_energy_daily_maximum"]] <- process_variable_pipeline(
    var_name = "convective_available_potential_energy_daily_maximum",
    years = years,
    input_dir = input_dir,
    output_dir = output_dir,
    quantile_threshold = quantile_cape
  )
  gc(verbose = FALSE)
  
  cat("\n VARIABLE 4/4: Vent (combinaison u + v)\n")
  dt_wind <- process_wind_components(years, input_dir)
  extremes_wind <- compute_extreme_events(
    dt_wind, 
    "wind_speed_daily_maximum",
    quantile_threshold = quantile_wind
  )
  
  output_file <- file.path(output_dir, "extremes_wind_speed_daily_maximum.rds")
  saveRDS(extremes_wind, output_file)
  cat(sprintf("\n✓ Résultats sauvegardés: %s\n", basename(output_file)))
  results[["wind_speed_daily_maximum"]] <- extremes_wind
  gc(verbose = FALSE)
  
  cat(sprintf("\n%s\n", strrep("=", 80)))
  cat(sprintf("✅ TRAITEMENT TERMINÉ POUR TOUTES LES VARIABLES\n"))
  cat(sprintf("%s\n\n", strrep("=", 80)))
  
  for (var in names(results)) {
    cat(sprintf("  • %s: %.0f événements extrêmes\n", 
                var, nrow(results[[var]])))
  }
  
  cat(sprintf("\n%s\n\n", strrep("=", 80)))
  
  return(results)
}
process_variable_pipeline <- function(var_name, years, 
                                      input_dir = "data/processed",
                                      output_dir = "data/extremes",
                                      quantile_threshold = 0.99) {
  
  cat(sprintf("\n%s\n", strrep("=", 80)))
  cat(sprintf("PIPELINE COMPLET: %s\n", var_name))
  cat(sprintf("%s\n\n", strrep("=", 80)))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  df <- load_single_variable(var_name, years, input_dir)
  
  dt <- preprocess_single_variable(df, var_name)
  rm(df)
  gc(verbose = FALSE)
  
  extremes <- compute_extreme_events(
    dt, 
    var_name, 
    quantile_threshold = quantile_threshold
  )
  
  rm(dt)
  gc(verbose = FALSE)
  
  output_file <- file.path(output_dir, paste0("extremes_", var_name, ".rds"))
  saveRDS(extremes, output_file)
  cat(sprintf("\n✓ Résultats sauvegardés: %s\n", basename(output_file)))
  
  return(extremes)
}
compute_extreme_events <- function(dt_var, var_name, quantile_threshold = 0.99) {
  
  cat(sprintf("\n=== Détection des extrêmes pour %s ===\n", var_name))
  
  if (!is.data.table(dt_var)) {
    dt_var <- as.data.table(dt_var)
  }
  
  required_cols <- c("lon", "lat", "value", "date", "coord_id")
  missing_cols <- required_cols[!required_cols %in% names(dt_var)]
  if (length(missing_cols) > 0) {
    stop(paste("Colonnes manquantes:", paste(missing_cols, collapse = ", ")))
  }
  
  setkeyv(dt_var, "coord_id")
  
  # -------------------------------------------------------------------------
  # -------------------------------------------------------------------------
  
  if (var_name == "total_precipitation_daily_sum") {
    # ------------------------
    # ------------------------
    cat("Méthode: Q90 sur jours humides (> 1mm)\n")
    
    dt_wet <- dt_var[value > 0]
    cat(sprintf("  - Jours humides: %.0f / %.0f (%.1f%%)\n", 
                nrow(dt_wet), nrow(dt_var), 100 * nrow(dt_wet) / nrow(dt_var)))
    
    quantiles_dt <- dt_wet[, .(
      threshold = quantile(value, quantile_threshold, na.rm = TRUE),
      n_wet_days = .N
    ), by = coord_id]
    
    setkey(quantiles_dt, coord_id)
    dt_var <- dt_var[quantiles_dt]
    
    extreme_events <- dt_var[value >0 & value >= threshold]
    extreme_events[, excess := value - threshold]
    extreme_events[, method :=  paste0("precip_q", quantile_threshold * 100)]
    
  } else if (var_name == "convective_precipitation_daily_sum") {
    # ------------------------
    # ------------------------
    cat("Méthode: Tous les jours avec précipitations convectives (> 0)\n")
    
    extreme_events <- dt_var[value > 0]
    extreme_events[, excess := value]
    extreme_events[, threshold := 0]
    extreme_events[, method := "convective_positive"]
    
    cat(sprintf("  - Jours avec précip. convectives: %.0f / %.0f (%.1f%%)\n", 
                nrow(extreme_events), nrow(dt_var), 
                100 * nrow(extreme_events) / nrow(dt_var)))
    
  } else if (grepl("cape|convective_available", var_name)) {
    
    cat(sprintf("Méthode: Q%.0f sur CAPE\n", quantile_threshold * 100))
    
    quantiles_dt <- dt_var[, .(
      threshold = quantile(value, quantile_threshold, na.rm = TRUE),
      n_obs = .N
    ), by = coord_id]
    
    setkey(quantiles_dt, coord_id)
    dt_var <- dt_var[quantiles_dt]
    
    extreme_events <- dt_var[value >= threshold]
    extreme_events[, excess := value - threshold]
    extreme_events[, method := paste0("cape_q", quantile_threshold * 100)]
    
  } else if (grepl("wind", var_name)) {
    
    cat(sprintf("Méthode: Q%.0f sur vitesse du vent\n", quantile_threshold * 100))
    
    quantiles_dt <- dt_var[, .(
      threshold = quantile(value, quantile_threshold, na.rm = TRUE),
      n_obs = .N
    ), by = coord_id]
    
    setkey(quantiles_dt, coord_id)
    dt_var <- dt_var[quantiles_dt]
    
    extreme_events <- dt_var[value >= threshold]
    extreme_events[, excess := value - threshold]
    extreme_events[, method := paste0("wind_q", quantile_threshold * 100)]
    
  } else {
    cat(sprintf("Méthode: Q%.0f générique\n", quantile_threshold * 100))
    
    quantiles_dt <- dt_var[, .(
      threshold = quantile(value, quantile_threshold, na.rm = TRUE)
    ), by = coord_id]
    
    setkey(quantiles_dt, coord_id)
    dt_var <- dt_var[quantiles_dt]
    
    extreme_events <- dt_var[value >= threshold]
    extreme_events[, excess := value - threshold]
    extreme_events[, method := paste0("generic_q", quantile_threshold * 100)]
  }
  
  extreme_events[, variable := var_name]
  
  cat(sprintf("✓ %.0f événements extrêmes détectés (%.2f%% des observations)\n", 
              nrow(extreme_events), 
              100 * nrow(extreme_events) / nrow(dt_var)))
  
  return(extreme_events)
}
preprocess_single_variable <- function(df, var_name) {
  cat(sprintf("\n=== Prétraitement de %s ===\n", var_name))
  
  dt <- as.data.table(df)
  
  if ("lyr.1" %in% names(dt)) {
    setnames(dt, "lyr.1", "value")
  }
  
  if (grepl("precipitation", var_name) && !grepl("convective", var_name)) {
    dt[value < 0.001, value := 0]
    cat("✓ Précipitations < 0.001 mises à 0\n")
  }
  
  cat(sprintf("✓ %.0f lignes prétraitées\n", nrow(dt)))
  return(dt)
}
load_single_variable <- function(var_name, years, input_dir = "data/processed") {
  cat(sprintf("\n=== Chargement de %s ===\n", var_name))
  annual_dfs <- list()
  
  for (year in years) {
    file_path <- file.path(input_dir, paste0(var_name, "_", year, ".rds"))
    
    if (file.exists(file_path)) {
      annual_dfs[[as.character(year)]] <- readRDS(file_path)
    } else {
      warning(sprintf("Fichier manquant: %s", basename(file_path)))
    }
  }
  
  if (length(annual_dfs) > 0) {
    result <- do.call(rbind, annual_dfs)
    cat(sprintf("✓ %d années chargées, %.0f lignes\n", 
                length(annual_dfs), nrow(result)))
    return(result)
  } else {
    stop(paste("Aucune donnée trouvée pour", var_name))
  }
}
