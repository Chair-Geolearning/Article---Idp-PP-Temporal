group_extreme_events_by_proximity <- function(extreme_events,
                                              temporal_width = 1,
                                              distance_threshold = 50000,
                                              parallel = TRUE,
                                              num_cores = parallel::detectCores() - 1) {
  
  required_cols <- c("lon", "lat", "date", "coord_id", "nuts_id", "value")
  missing_cols <- setdiff(required_cols, names(extreme_events))
  if (length(missing_cols) > 0) {
    stop("Colonnes manquantes: ", paste(missing_cols, collapse = ", "))
  }
  
  library(data.table)
  library(sf)
  library(igraph)
  if (parallel) {
    library(doParallel)
    library(foreach)
  }
  
  if (!data.table::is.data.table(extreme_events)) {
    extreme_events <- data.table::as.data.table(extreme_events)
  }
  
  data.table::setorder(extreme_events, coord_id, date)
  extreme_events[, `:=`(
    date_diff = as.numeric(date - data.table::shift(date, 1)),
    new_event = FALSE
  ), by = coord_id]
  extreme_events[is.na(date_diff), `:=`(date_diff = 0, new_event = TRUE)]
  extreme_events[date_diff > 1, new_event := TRUE]
  extreme_events[, event_id := paste0(coord_id, "_", cumsum(new_event)), by = coord_id]
  
  events_summary <- extreme_events[, .(
    start_date  = min(date),
    end_date    = max(date),
    duration    = as.numeric(max(date) - min(date) + 1),
    obs_days    = .N,
    total_value = sum(value),
    max_value   = max(value),
    lon         = lon[1],
    lat         = lat[1],
    nuts_id     = nuts_id[1]
  ), by = event_id]
  
  events_sf <- sf::st_as_sf(events_summary, coords = c("lon", "lat"), crs = 4326)
  events_sf <- sf::st_transform(events_sf, 3035)
  
  unique_nuts <- unique(events_sf$nuts_id)
  n_nuts <- length(unique_nuts)
  
  message("\n=== Début du traitement ===")
  message("Nombre total de régions NUTS: ", n_nuts)
  message("Nombre total d'événements (pixel-runs): ", nrow(events_sf))
  message("Seuil de distance: ", distance_threshold, " m")
  message("Fenêtre temporelle (gap entre intervalles): ≤ ", temporal_width, " jours")
  message("Parallélisation: ", ifelse(parallel, "OUI", "NON"))
  if (parallel) message("Nombre de cœurs: ", num_cores)
  message("===========================\n")
  
  start_time <- Sys.time()
  
  process_nuts_region <- function(current_nuts, nuts_index = NULL) {
    if (!is.null(nuts_index)) {
      message(sprintf("[%d/%d] Traitement NUTS: %s", nuts_index, n_nuts, current_nuts))
    }
    
    region_start <- Sys.time()
    nuts_events <- events_sf[events_sf$nuts_id == current_nuts, ]
    n_events <- nrow(nuts_events)
    
    if (n_events == 0) return(NULL)
    
    if (n_events == 1) {
      nuts_events$cluster_id <- paste0(current_nuts, "_1")
      if (!is.null(nuts_index)) {
        message(sprintf("  └─ 1 événement, 1 cluster (%.2f sec)",
                        as.numeric(Sys.time() - region_start, units = "secs")))
      }
      return(nuts_events)
    }
    
    if (!is.null(nuts_index)) {
      message(sprintf("  └─ %d événements à traiter...", n_events))
    }
    
    ord <- order(nuts_events$start_date)
    nuts_events <- nuts_events[ord, ]
    
    start_dates <- nuts_events$start_date
    end_dates   <- nuts_events$end_date
    event_ids   <- nuts_events$event_id
   
    compute_edges <- function(i_vals) {
      edges_list <- vector("list", length(i_vals))
      
      for (idx in seq_along(i_vals)) {
        i <- i_vals[idx]
        
        j_indices <- (i + 1):n_events
        if (length(j_indices) == 0) next
        
        temporal_filter <- start_dates[j_indices] <= (end_dates[i] + temporal_width + 1L)
        if (!any(temporal_filter)) next
        
        j_filtered <- j_indices[temporal_filter]
        
        gap_days <- as.integer(start_dates[j_filtered] - end_dates[i] - 1L)
        gap_days[gap_days < 0L] <- 0L
        
        j_close_time <- j_filtered[gap_days <= temporal_width]
        if (length(j_close_time) == 0) next
        
        distances <- as.numeric(sf::st_distance(nuts_events[i, ], nuts_events[j_close_time, ]))
        close_indices <- j_close_time[distances <= distance_threshold]
        
        if (length(close_indices) > 0) {
          edges_list[[idx]] <- data.frame(
            from = event_ids[i],
            to = event_ids[close_indices],
            stringsAsFactors = FALSE
          )
        }
      }
      
      do.call(rbind, edges_list[!sapply(edges_list, is.null)])
    }
    
    if (parallel && n_events > 100) {
      cl <- parallel::makeCluster(min(num_cores, n_events - 1))
      on.exit(parallel::stopCluster(cl), add = TRUE)
      doParallel::registerDoParallel(cl)
      
      chunk_size <- max(10, ceiling((n_events - 1) / (num_cores * 4)))
      i_vals <- 1:(n_events - 1)
      chunks <- split(i_vals, ceiling(seq_along(i_vals) / chunk_size))
      
      edges <- foreach::foreach(
        chunk = chunks,
        .combine = rbind,
        .packages = c("sf"),
        .export = c("start_dates", "end_dates", "event_ids", "n_events",
                    "nuts_events", "temporal_width", "distance_threshold")
      ) %dopar% {
        edges_list <- vector("list", length(chunk))
        
        for (idx in seq_along(chunk)) {
          i <- chunk[idx]
          
          j_indices <- (i + 1):n_events
          if (length(j_indices) == 0) next
          
          temporal_filter <- start_dates[j_indices] <= (end_dates[i] + temporal_width + 1L)
          if (!any(temporal_filter)) next
          
          j_filtered <- j_indices[temporal_filter]
          
          gap_days <- as.integer(start_dates[j_filtered] - end_dates[i] - 1L)
          gap_days[gap_days < 0L] <- 0L
          
          j_close_time <- j_filtered[gap_days <= temporal_width]
          if (length(j_close_time) == 0) next
          
          distances <- as.numeric(sf::st_distance(nuts_events[i, ], nuts_events[j_close_time, ]))
          close_indices <- j_close_time[distances <= distance_threshold]
          
          if (length(close_indices) > 0) {
            edges_list[[idx]] <- data.frame(
              from = event_ids[i],
              to = event_ids[close_indices],
              stringsAsFactors = FALSE
            )
          }
        }
        
        do.call(rbind, edges_list[!sapply(edges_list, is.null)])
      }
    } else {
      edges <- compute_edges(1:(n_events - 1))
    }
    
    if (is.null(edges) || nrow(edges) == 0) {
      nuts_events$cluster_id <- paste0(current_nuts, "_", seq_len(n_events))
      if (!is.null(nuts_index)) {
        message(sprintf("     Aucune connexion → %d clusters isolés (%.2f sec)",
                        n_events, as.numeric(Sys.time() - region_start, units = "secs")))
      }
      return(nuts_events)
    }
    
    g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = event_ids)
    comps <- igraph::components(g)
    
    n_clusters <- comps$no
    if (!is.null(nuts_index)) {
      message(sprintf("     %d connexions trouvées → %d clusters (%.2f sec)",
                      nrow(edges), n_clusters,
                      as.numeric(Sys.time() - region_start, units = "secs")))
    }
    
    cluster_map <- setNames(
      paste0(current_nuts, "_", comps$membership),
      names(comps$membership)
    )
    
    nuts_events$cluster_id <- cluster_map[nuts_events$event_id]
    
    isolated <- is.na(nuts_events$cluster_id)
    if (any(isolated)) {
      max_cluster <- max(comps$membership, na.rm = TRUE)
      nuts_events$cluster_id[isolated] <- paste0(
        current_nuts, "_",
        max_cluster + seq_len(sum(isolated))
      )
    }
    
    nuts_events
  }
  
  results <- vector("list", n_nuts)
  nuts_completed <- 0
  region_start_time <- Sys.time()
  
  for (i in seq_along(unique_nuts)) {
    results[[i]] <- process_nuts_region(unique_nuts[i], nuts_index = i)
    
    nuts_completed <- nuts_completed + 1
    elapsed <- as.numeric(Sys.time() - region_start_time, units = "secs")
    avg_time <- elapsed / nuts_completed
    remaining <- n_nuts - nuts_completed
    eta <- avg_time * remaining
    pct <- round(100 * nuts_completed / n_nuts, 1)
    
    if (nuts_completed %% max(1, floor(n_nuts / 20)) == 0 || nuts_completed == n_nuts) {
      message(sprintf("\n[%d/%d] %.1f%% complété | Temps écoulé: %.1fs | ETA: ~%.1fs (%d régions restantes)",
                      nuts_completed, n_nuts, pct, elapsed, eta, remaining))
    }
  }
  
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) stop("Aucun résultat trouvé")
  
  message("\n--- Assemblage des résultats ---")
  final_result <- do.call(rbind, results)
  final_result_df <- sf::st_drop_geometry(final_result)
  
  message("Calcul des statistiques par cluster...")
  final_dt <- data.table::as.data.table(final_result_df)
  
  cluster_summary <- final_dt[, .(
    event_count = .N,
    start_date = min(start_date),
    end_date   = max(end_date),
    duration   = as.numeric(max(end_date) - min(start_date) + 1),
    obs_days   = sum(obs_days),
    total_value = sum(total_value),
    max_value   = max(max_value),
    spatial_extent = NA_real_
  ), by = cluster_id]
  
  total_time <- as.numeric(Sys.time() - start_time, units = "secs")
  
  message("\n=== Traitement terminé ===")
  message("Temps total: ", sprintf("%.1f sec (%.1f min)", total_time, total_time/60))
  message("Événements traités: ", nrow(final_result_df))
  message("Clusters identifiés: ", nrow(cluster_summary))
  message("Temps moyen par région NUTS: ", sprintf("%.2f sec", total_time/n_nuts))
  message("==========================\n")
  
  list(
    events = final_result_df,
    clusters = as.data.frame(cluster_summary)
  )
}
