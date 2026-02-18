plot_idp_maps_pw_multi <- function(
    inter_df_all,
    intra_df,
    out_dir = "output/testidp_precip_wind_cape",
    nuts_refs = c("FRL0", "FI1D", "DE21"),
    alpha = 0.05,
    nuts_level = 2,
    nuts_year = 2016,
    nuts_resolution = "20",
    width = 2400, height = 1800, res = 200,
    p_floor = 1e-6,
    log10_cap = 6
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(sf)
    library(ggplot2)
    library(giscoR)
    library(stringr)
    library(viridis)
    library(tibble)
  })
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  safe <- function(x) {
    x |>
      as.character() |>
      str_replace_all("[^A-Za-z0-9_\\-]+", "_") |>
      str_replace_all("_+", "_") |>
      str_replace_all("^_|_$", "")
  }
  
  n_more_extreme <- function(p, n_sim) {
    ifelse(
      is.na(p) | is.na(n_sim),
      NA_real_,
      pmax(0, pmin(n_sim, round((1 - p) * (n_sim + 1))))
    )
  }
  
  neglog10p <- function(p) {
    out <- -log10(pmax(p, p_floor))
    pmin(out, log10_cap)
  }
  
  nuts_raw <- giscoR::gisco_get_nuts(
    year = nuts_year,
    nuts_level = nuts_level,
    resolution = nuts_resolution
  ) |> sf::st_as_sf()
  
  id_candidates <- c("nuts_id", "NUTS_ID", "id", "ID")
  id_col <- id_candidates[id_candidates %in% names(nuts_raw)][1]
  if (is.na(id_col)) {
    stop("Could not find a NUTS id column in gisco_get_nuts(). Found: ",
         paste(names(nuts_raw), collapse = ", "))
  }
  
  nuts_sf <- nuts_raw |>
    mutate(nuts_id = .data[[id_col]]) |>
    select(nuts_id, geometry)
  
  out_files <- c()
  
  for (ref in nuts_refs) {
    
    inter_df_ref <- inter_df_all |> filter(nuts_ref == ref)
    if (nrow(inter_df_ref) == 0) next
    
    inter_vars <- sort(unique(inter_df_ref$variable))
    
    for (vv in inter_vars) {
      
      df_v <- inter_df_ref |>
        filter(variable == vv) |>
        transmute(
          nuts_id = nuts_target,
          p_value = p_value,
          n_sim   = n_sim,
          dep = ifelse(is.na(p_value), NA_character_,
                       ifelse(p_value < alpha, "Dependence", "Independence")),
          n_out = n_more_extreme(p_value, n_sim),
          neglogp = neglog10p(p_value)
        )
      
      df_ref <- tibble(
        nuts_id = ref,
        p_value = NA_real_,
        n_sim   = unique(inter_df_ref$n_sim)[1],
        dep     = "Reference",
        n_out   = NA_real_,
        neglogp = NA_real_
      )
      
      df_plot <- bind_rows(df_v, df_ref)
      
      sf_v <- nuts_sf |>
        filter(nuts_id %in% df_plot$nuts_id) |>
        left_join(df_plot, by = "nuts_id")
      
      p_bin <- ggplot(sf_v) +
        geom_sf(aes(fill = dep), colour = "grey35", linewidth = 0.15) +
        scale_fill_manual(
          values = c(
            "Independence" = "#4CAF50",
            "Dependence"   = "#D32F2F",
            "Reference"    = "#FB8C00"
          ),
          na.value = "grey85",
          drop = FALSE
        ) +
        labs(
          title = paste0("INTER (", ref, " vs all) | ", vv),
          subtitle = paste0("Translation | periodic yearly | r = 3.5 | two-sided | alpha = ", alpha),
          fill = NULL
        ) +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(),
              plot.title = element_text(face = "bold"))
      
      f_bin <- file.path(out_dir, paste0("MAP_INTER_", safe(ref), "_", safe(vv), "_binary.png"))
      ggsave(f_bin, p_bin, width = width / res, height = height / res, dpi = res)
      out_files <- c(out_files, f_bin)
      
      p_out <- ggplot(sf_v) +
        geom_sf(aes(fill = n_out), colour = "grey35", linewidth = 0.15) +
        geom_sf(data = sf_v |> filter(nuts_id == ref),
                fill = NA, colour = "#FB8C00", linewidth = 0.9) +
        scale_fill_viridis_c(na.value = "grey85") +
        labs(
          title = paste0("INTER (", ref, " vs all) | ", vv),
          subtitle = "Proxy extremeness: round((1 - p) * (n_sim + 1))",
          fill = "n_more_extreme"
        ) +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(),
              plot.title = element_text(face = "bold"))
      
      f_out <- file.path(out_dir, paste0("MAP_INTER_", safe(ref), "_", safe(vv), "_n_more_extreme.png"))
      ggsave(f_out, p_out, width = width / res, height = height / res, dpi = res)
      out_files <- c(out_files, f_out)
      
      p_log <- ggplot(sf_v) +
        geom_sf(aes(fill = neglogp), colour = "grey35", linewidth = 0.15) +
        geom_sf(data = sf_v |> filter(nuts_id == ref),
                fill = NA, colour = "#FB8C00", linewidth = 0.9) +
        scale_fill_viridis_c(na.value = "grey85", limits = c(0, log10_cap)) +
        labs(
          title = paste0("INTER (", ref, " vs all) | ", vv),
          subtitle = paste0("-log10(p), floored at ", p_floor, ", capped at ", log10_cap),
          fill = expression(-log[10](p))
        ) +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(),
              plot.title = element_text(face = "bold"))
      
      f_log <- file.path(out_dir, paste0("MAP_INTER_", safe(ref), "_", safe(vv), "_neglog10p.png"))
      ggsave(f_log, p_log, width = width / res, height = height / res, dpi = res)
      out_files <- c(out_files, f_log)
    }
  }
  intra_df2 <- intra_df
  
  if (!"nuts_code" %in% names(intra_df2) && "nuts_id" %in% names(intra_df2)) {
    intra_df2 <- intra_df2 |> mutate(nuts_code = nuts_id)
  }
  
  if (!"direction" %in% names(intra_df2)) {
    if (all(c("var1","var2") %in% names(intra_df2))) {
      intra_df2 <- intra_df2 |> mutate(direction = paste0(var1, "_to_", var2))
    } else {
      stop("INTRA df must contain either `direction` or (`var1`,`var2`).")
    }
  }
  
  if (!"alternative" %in% names(intra_df2)) {
    intra_df2 <- intra_df2 |> mutate(alternative = "greater")
  }
  
  intra_keys <- intra_df2 |>
    distinct(direction, alternative) |>
    arrange(direction, alternative)
  
  for (i in seq_len(nrow(intra_keys))) {
    dir_lab <- intra_keys$direction[i]
    alt_lab <- intra_keys$alternative[i]
    
    df_v <- intra_df2 |>
      filter(direction == dir_lab, alternative == alt_lab) |>
      transmute(
        nuts_id = nuts_code,
        p_value = p_value,
        n_sim   = n_sim,
        dep = ifelse(is.na(p_value), NA_character_,
                     ifelse(p_value < alpha, "Dependence", "Independence")),
        n_out = n_more_extreme(p_value, n_sim),
        neglogp = neglog10p(p_value)
      )
    
    if (nrow(df_v) == 0) next
    
    sf_v <- nuts_sf |>
      filter(nuts_id %in% df_v$nuts_id) |>
      left_join(df_v, by = "nuts_id")
    
    title_base <- paste0("INTRA directional | ", gsub("_to_", " \u2192 ", dir_lab),
                         " | ", alt_lab)
    
    # (1) Binary
    p_bin <- ggplot(sf_v) +
      geom_sf(aes(fill = dep), colour = "grey35", linewidth = 0.15) +
      scale_fill_manual(
        values = c(
          "Independence" = "#4CAF50",
          "Dependence"   = "#D32F2F"
        ),
        na.value = "grey85",
        drop = FALSE
      ) +
      labs(
        title = title_base,
        subtitle = paste0("Translation | periodic yearly | r = 3.5 | alpha = ", alpha),
        fill = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold"))
    
    f_bin <- file.path(out_dir, paste0("MAP_INTRA_DIR_", safe(dir_lab), "_", safe(alt_lab), "_binary.png"))
    ggsave(f_bin, p_bin, width = width / res, height = height / res, dpi = res)
    out_files <- c(out_files, f_bin)
    
    p_out <- ggplot(sf_v) +
      geom_sf(aes(fill = n_out), colour = "grey35", linewidth = 0.15) +
      scale_fill_viridis_c(na.value = "grey85") +
      labs(
        title = title_base,
        subtitle = "Proxy extremeness: round((1 - p) * (n_sim + 1))",
        fill = "n_more_extreme"
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold"))
    
    f_out <- file.path(out_dir, paste0("MAP_INTRA_DIR_", safe(dir_lab), "_", safe(alt_lab), "_n_more_extreme.png"))
    ggsave(f_out, p_out, width = width / res, height = height / res, dpi = res)
    out_files <- c(out_files, f_out)
    
    p_log <- ggplot(sf_v) +
      geom_sf(aes(fill = neglogp), colour = "grey35", linewidth = 0.15) +
      scale_fill_viridis_c(na.value = "grey85", limits = c(0, log10_cap)) +
      labs(
        title = title_base,
        subtitle = paste0("-log10(p), floored at ", p_floor, ", capped at ", log10_cap),
        fill = expression(-log[10](p))
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold"))
    
    f_log <- file.path(out_dir, paste0("MAP_INTRA_DIR_", safe(dir_lab), "_", safe(alt_lab), "_neglog10p.png"))
    ggsave(f_log, p_log, width = width / res, height = height / res, dpi = res)
    out_files <- c(out_files, f_log)
  }
  
  unique(out_files)
}
