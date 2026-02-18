process_netcdf_files <- function(var_names, years = NULL, nuts_level = 2, 
                                 save_intermediate = TRUE, output_dir = "data/processed",
                                 n_cores = NULL, use_parallel = TRUE) {
  

  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
 
  if (save_intermediate && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Configuration:\n"))
  cat(sprintf("  - Variables: %d\n", length(var_names)))
  cat(sprintf("  - Années: %d (%d-%d)\n", length(years), min(years), max(years)))
  cat(sprintf("  - Fichiers totaux: %d\n", length(var_names) * length(years)))
  cat(sprintf("  - Parallélisation: %s (%d cœurs)\n", ifelse(use_parallel, "OUI", "NON"), n_cores))
  cat(sprintf("  - Sauvegarde intermédiaire: %s\n", ifelse(save_intermediate, "OUI", "NON")))
  cat(sprintf("========================================\n\n"))
  
  process_single_file <- function(var_name, year, nuts, nuts_level, output_dir) {
    
    if (grepl("temperature", var_name)) {
      filename <- paste0("data/2m_temperature_daily_mean_", year, ".nc")
    } else {
      filename <- paste0("data/", var_name, "_", year, ".nc")
    }
    
    if (!file.exists(filename)) {
      return(list(
        success = FALSE,
        var_name = var_name,
        year = year,
        message = "Fichier introuvable"
      ))
    }
    
    output_file <- file.path(output_dir, paste0(var_name, "_", year, ".rds"))
    if (file.exists(output_file)) {
      return(list(
        success = TRUE,
        var_name = var_name,
        year = year,
        message = "Déjà traité (fichier existant)",
        rows = 0,
        skipped = TRUE
      ))
    }
    
    tryCatch({
      raw_data <- load_nc_data(filename)
      
      df <- preprocess_data(raw_data, var_name)
      
      points_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
      points_nuts <- st_join(points_sf, nuts, join = st_intersects)
      df$nuts_id <- points_nuts$NUTS_ID
      df$nuts_name <- points_nuts$NAME_LATN
      
      df <- df[complete.cases(df), ]
      
      saveRDS(df, output_file)
      
      return(list(
        success = TRUE,
        var_name = var_name,
        year = year,
        message = "Traité avec succès",
        rows = nrow(df),
        skipped = FALSE
      ))
      
    }, error = function(e) {
      return(list(
        success = FALSE,
        var_name = var_name,
        year = year,
        message = as.character(e$message)
      ))
    })
  }
  
  load_nc_data <- function(file_path) {
    nc_data <- nc_open(file_path)
    on.exit(nc_close(nc_data))
    
    vars <- names(nc_data$var)
    
    main_var <- NULL
    
    if (grepl("precipitation", file_path)) {
      if ("tp" %in% vars) main_var <- "tp"
      else if (any(grepl("precip|pr", vars))) main_var <- vars[grepl("precip|pr", vars)][1]
    } else if (grepl("temperature", file_path)) {
      if ("t2m" %in% vars) main_var <- "t2m"
      else if (any(grepl("t2m|temp|temperature|2m_temperature", vars))) main_var <- vars[grepl("t2m|temp|temperature|2m_temperature", vars)][1]
    } else if (grepl("wind.*u_component", file_path)) {
      if ("u10" %in% vars) main_var <- "u10"
      else if (any(grepl("10m_u|u_component|u10", vars))) main_var <- vars[grepl("10m_u|u_component|u10", vars)][1]
    } else if (grepl("wind.*v_component", file_path)) {
      if ("v10" %in% vars) main_var <- "v10"
      else if (any(grepl("10m_v|v_component|v10", vars))) main_var <- vars[grepl("10m_v|v_component|v10", vars)][1]
    } else if (grepl("cape|convective_available_potential", file_path)) {
      if ("cape" %in% vars) main_var <- "cape"
      else if (any(grepl("convective_available|cape", vars))) main_var <- vars[grepl("convective_available|cape", vars)][1]
    } else if (grepl("convective_precipitation", file_path)) {
      if ("cp" %in% vars) main_var <- "cp"
      else if (any(grepl("convective_precipitation|cp", vars))) main_var <- vars[grepl("convective_precipitation|cp", vars)][1]
    }
    
    if (is.null(main_var) && length(vars) > 0) {
      coord_like <- vars[grepl("^(lon|lat|time|x|y|t)$", vars)]
      main_candidates <- vars[!vars %in% coord_like]
      
      if (length(main_candidates) > 0) {
        main_var <- main_candidates[1]
      } else {
        main_var <- vars[1]
      }
    }
    
    main_data <- NULL
    tryCatch({
      main_data <- ncvar_get(nc_data, main_var)
    }, error = function(e) {
      tryCatch({
        main_data <<- ncvar_get(nc_data, 0)
      }, error = function(e2) {
        for (v in vars) {
          if (!is.null(main_data)) break
          
          tryCatch({
            temp_data <- ncvar_get(nc_data, v)
            if (!is.null(temp_data) && is.array(temp_data) && length(dim(temp_data)) >= 2) {
              main_data <<- temp_data
              main_var <<- v
              break
            }
          }, error = function(e3) {})
        }
      })
    })
    
    if (is.null(main_data)) {
      stop(paste("Impossible d'extraire les données du fichier:", file_path))
    }
    
    find_coord_var <- function(possible_names) {
      match <- vars[vars %in% possible_names]
      if (length(match) > 0) return(match[1]) else return(NA)
    }
    
    lon_var <- find_coord_var(c("lon", "longitude", "x"))
    lat_var <- find_coord_var(c("lat", "latitude", "y"))
    time_var <- find_coord_var(c("time", "valid_time", "t"))
    
    get_coord_safely <- function(var_name) {
      if (is.na(var_name)) return(NULL)
      tryCatch({
        return(ncvar_get(nc_data, var_name))
      }, error = function(e) {
        return(NULL)
      })
    }
    
    lon_data <- get_coord_safely(lon_var)
    lat_data <- get_coord_safely(lat_var)
    time_data <- get_coord_safely(time_var)
    
    data_dims <- dim(main_data)
    
    if (is.null(lon_data) || length(lon_data) != data_dims[1]) {
      lon_data <- seq(-10, 40, length.out = data_dims[1])
    }
    
    if (is.null(lat_data) || length(lat_data) != data_dims[2]) {
      lat_data <- seq(35, 70, length.out = data_dims[2])
    }
    
    dates <- NULL
    if (!is.null(time_data) && !is.na(time_var)) {
      time_units <- tryCatch(
        ncatt_get(nc_data, time_var, "units")$value,
        error = function(e) NULL
      )
      
      if (!is.null(time_units)) {
        if (grepl("hours since", time_units)) {
          origin_date <- sub("hours since ", "", time_units)
          origin_date <- sub(" .*", "", origin_date)
          dates <- as.Date(time_data / 24, origin = origin_date)
        } else if (grepl("days since", time_units)) {
          origin_date <- sub("days since ", "", time_units)
          origin_date <- sub(" .*", "", origin_date)
          dates <- as.Date(time_data, origin = origin_date)
        }
      }
    }
    
    if (is.null(dates)) {
      year_match <- regexpr("\\d{4}", file_path)
      if (year_match > 0) {
        year <- substr(file_path, year_match, year_match + 3)
        start_date <- as.Date(paste0(year, "-01-01"))
      } else {
        start_date <- as.Date("1960-01-01")
      }
      
      n_times <- ifelse(length(data_dims) >= 3, data_dims[3], 1)
      dates <- seq.Date(from = start_date, by = "day", length.out = n_times)
    }
    
    list(
      data = main_data,
      lon = lon_data,
      lat = lat_data,
      dates = dates
    )
  }
  
  preprocess_data <- function(raw_data, var_name) {
    data_dims <- dim(raw_data$data)
    
    if (length(data_dims) >= 3) {
      r_explicit <- rast(
        nrows = length(raw_data$lat),
        ncols = length(raw_data$lon),
        xmin = min(raw_data$lon),
        xmax = max(raw_data$lon),
        ymin = min(raw_data$lat),
        ymax = max(raw_data$lat),
        crs = "+proj=longlat +datum=WGS84"
      )
      
      result_list <- list()
      for (t in 1:length(raw_data$dates)) {
        if (length(data_dims) == 3) {
          layer_data <- raw_data$data[,,t]
        } else if (length(data_dims) == 4) {
          layer_data <- raw_data$data[,,1,t]
        } else {
          layer_data <- raw_data$data[,,rep(1, length(data_dims)-3),t]
        }
        
        values(r_explicit) <- as.vector(layer_data)
        
        cells <- cells(r_explicit)
        coords <- xyFromCell(r_explicit, cells)
        
        df_t <- data.frame(
          lon = coords[, 1],
          lat = coords[, 2],
          value = values(r_explicit),
          date = raw_data$dates[t],
          variable = var_name,
          coord_id = paste0(formatC(round(coords[, 1] * 4) / 4, format = "f", digits = 2), "_", 
                            formatC(round(coords[, 2] * 4) / 4, format = "f", digits = 2))
        )
        
        result_list[[t]] <- df_t
      }
      
      return(do.call(rbind, result_list))
    } else if (length(data_dims) == 2) {
      r_explicit <- rast(
        nrows = length(raw_data$lat),
        ncols = length(raw_data$lon),
        xmin = min(raw_data$lon),
        xmax = max(raw_data$lon),
        ymin = min(raw_data$lat),
        ymax = max(raw_data$lat),
        crs = "+proj=longlat +datum=WGS84"
      )
      
      values(r_explicit) <- as.vector(raw_data$data)
      
      cells <- cells(r_explicit)
      coords <- xyFromCell(r_explicit, cells)
      
      return(data.frame(
        lon = coords[, 1],
        lat = coords[, 2],
        value = values(r_explicit),
        date = raw_data$dates[1],
        variable = var_name,
        coord_id = paste0(formatC(round(coords[, 1] * 4) / 4, format = "f", digits = 2), "_", 
                          formatC(round(coords[, 2] * 4) / 4, format = "f", digits = 2))
      ))
    } else {
      stop(paste("Format de données non supporté: dimensions =", 
                 paste(data_dims, collapse = " × ")))
    }
  }
  
  nuts <- gisco_get_nuts(year = "2021", resolution = "60", nuts_level = nuts_level) %>%
    st_transform(4326)
  
  tasks <- expand.grid(
    var_name = var_names,
    year = years,
    stringsAsFactors = FALSE
  )
  
  if (use_parallel && n_cores > 1) {
    library(parallel)
    library(doParallel)
    
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("process_single_file", "load_nc_data", "preprocess_data", 
                        "nuts", "nuts_level", "output_dir"), envir = environment())
    clusterEvalQ(cl, {
      library(ncdf4)
      library(terra)
      library(sf)
      library(dplyr)
    })
    
    cat(sprintf("Démarrage du traitement parallèle avec %d cœurs...\n\n", n_cores))
    
    results <- parLapply(cl, 1:nrow(tasks), function(i) {
      process_single_file(
        var_name = tasks$var_name[i],
        year = tasks$year[i],
        nuts = nuts,
        nuts_level = nuts_level,
        output_dir = output_dir
      )
    })
    
    stopCluster(cl)
  } else {
    cat("Traitement séquentiel...\n\n")
    results <- lapply(1:nrow(tasks), function(i) {
      result <- process_single_file(
        var_name = tasks$var_name[i],
        year = tasks$year[i],
        nuts = nuts,
        nuts_level = nuts_level,
        output_dir = output_dir
      )
      
      if (result$success && !result$skipped) {
        cat(sprintf("[%d/%d] ✓ %s - %s (%d rows)\n", 
                    i, nrow(tasks), result$var_name, result$year, result$rows))
      } else if (result$success && result$skipped) {
        cat(sprintf("[%d/%d] ⊙ %s - %s (skipped)\n", 
                    i, nrow(tasks), result$var_name, result$year))
      } else {
        cat(sprintf("[%d/%d] ✗ %s - %s: %s\n", 
                    i, nrow(tasks), result$var_name, result$year, result$message))
      }
      
      return(result)
    })
  }
  
  n_success <- sum(sapply(results, function(x) x$success && !x$skipped))
  n_skipped <- sum(sapply(results, function(x) x$success && x$skipped))
  n_failed <- sum(sapply(results, function(x) !x$success))
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("TRAITEMENT TERMINÉ:\n"))
  cat(sprintf("  - Réussis: %d\n", n_success))
  cat(sprintf("  - Déjà traités: %d\n", n_skipped))
  cat(sprintf("  - Échoués: %d\n", n_failed))
  cat(sprintf("  - Fichiers sauvegardés dans: %s\n", output_dir))
  cat(sprintf("========================================\n\n"))
  
  summary_list <- list()
  for (var_name in var_names) {
    var_results <- results[sapply(results, function(x) x$var_name == var_name)]
    years_processed <- sapply(var_results, function(x) if(x$success) x$year else NULL)
    years_processed <- unlist(years_processed[!sapply(years_processed, is.null)])
    
    summary_list[[var_name]] <- data.frame(
      variable = var_name,
      years_processed = paste(sort(years_processed), collapse = ","),
      total_files = length(years_processed),
      output_dir = output_dir
    )
  }
  
  return(summary_list)
}


