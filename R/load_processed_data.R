load_processed_data <- function(var_names, years, input_dir = "data/processed",
                                sample_dates = NULL, 
                                nuts_filter = NULL) { 
  result_list <- list()
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("CHARGEMENT DES DONNÉES TRAITÉES\n"))
  cat(sprintf("========================================\n\n"))
  
  for (var_name in var_names) {
    cat(sprintf("Chargement de %s...\n", var_name))
    annual_dfs <- list()
    
    for (year in years) {
      file_path <- file.path(input_dir, paste0(var_name, "_", year, ".rds"))
      
      if (file.exists(file_path)) {
        df <- readRDS(file_path)
        
        if (!is.null(sample_dates)) {
          unique_dates <- unique(df$date)
          sampled_dates <- unique_dates[seq(1, length(unique_dates), by = sample_dates)]
          df <- df[df$date %in% sampled_dates, ]
        }
        
        if (!is.null(nuts_filter)) {
          df <- df[substr(df$nuts_id, 1, 2) %in% nuts_filter, ]
        }
        
        annual_dfs[[as.character(year)]] <- df
        rm(df)
        gc(verbose = FALSE)
      }
    }
    
    if (length(annual_dfs) > 0) {
      result_list[[var_name]] <- do.call(rbind, annual_dfs)
      cat(sprintf("  ✓ %d années, %.0f lignes\n", 
                  length(annual_dfs), nrow(result_list[[var_name]])))
      rm(annual_dfs)
      gc(verbose = FALSE)
    }
  }
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("CHARGEMENT TERMINÉ - %d variables\n", length(result_list)))
  cat(sprintf("========================================\n\n"))
  
  return(result_list)
}
