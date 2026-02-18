extract_all_nuts <- function(variables_list) {
  # Extraire NUTS du premier élément de la liste
  nuts_list <- unique(variables_list[[1]]$events$nuts_id)
  nuts_list <- nuts_list[!is.na(nuts_list) & nuts_list != ""]
  nuts_sorted <- sort(nuts_list)
  
  cat("Total NUTS found:", length(nuts_sorted), "\n")
  cat("First 10:", paste(head(nuts_sorted, 10), collapse = ", "), "\n")
  
  return(nuts_sorted)
}
