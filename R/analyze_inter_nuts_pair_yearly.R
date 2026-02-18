analyze_inter_nuts_pair_yearly <- function(
    df_ref,
    df_target,
    var_name,
    nuts_ref,
    nuts_target,
    sim_method = "translation",
    alternative = "two.sided",
    r_value = 3,
    n_sim = 999,
    kde_bw = 0.08,
    homogeneous = FALSE,
    periodic = TRUE
) {
  
  if (nuts_ref == nuts_target) return(NULL)
  
  df_ref <- df_ref[!is.na(df_ref$start_date), ]
  df_target <- df_target[!is.na(df_target$start_date), ]
  
  df_ref <- ensure_jour(df_ref)
  df_target <- ensure_jour(df_target)
  
  df_ref_var <- df_ref[df_ref$Variable == var_name, ]
  df_target_var <- df_target[df_target$Variable == var_name, ]
  
  if (nrow(df_ref_var) < 5 || nrow(df_target_var) < 5) return(NULL)
  
  t_min <- max(min(df_ref$start_date), min(df_target$start_date))
  t_max <- min(max(df_ref$start_date), max(df_target$start_date))
  if (!is.finite(t_min) || !is.finite(t_max) || t_min >= t_max) return(NULL)
  
  df_ref_f <- df_ref[df_ref$start_date >= t_min & df_ref$start_date <= t_max, ]
  df_tar_f <- df_target[df_target$start_date >= t_min & df_target$start_date <= t_max, ]
  
  df_ref_f <- df_ref_f[df_ref_f$Variable == var_name, ]
  df_tar_f <- df_tar_f[df_tar_f$Variable == var_name, ]
  
  if (nrow(df_ref_f) < 5 || nrow(df_tar_f) < 5) return(NULL)
  
  res <- tryCatch({
    calculate_bivariate_K_GET_yearly(
      df1 = df_ref_f,
      df2 = df_tar_f,
      var1_name = var_name,
      var2_name = var_name,
      max_r = r_value,
      time_col = "Jour",
      T_total = 365,
      n_sim = n_sim,
      directional = FALSE,
      sim_method = sim_method,
      periodic = periodic,
      homogeneous = homogeneous,
      alternative = alternative,
      kde_bw = kde_bw,
      save_plot = FALSE
    )
  }, error = function(e) NULL)
  
  if (is.null(res)) return(NULL)
  
  pv <- res$p_value
  pv_num <- if (length(pv) == 2) pv[2] else pv
  
  data.frame(
    nuts_ref = nuts_ref,
    nuts_target = nuts_target,
    variable = var_name,
    type = "inter",
    sim_method = sim_method,
    alternative = alternative,
    r_value = r_value,
    n_sim = n_sim,
    periodic = periodic,
    homogeneous = homogeneous,
    kde_bw = kde_bw,
    p_value = pv_num,
    stringsAsFactors = FALSE
  )
}


analyze_intra_nuts_yearly <- function(
    df_nuts,
    var_names,
    nuts_code,
    sim_method = "translation",
    alternative = "two.sided",
    r_value = 3,
    n_sim = 999,
    kde_bw = 0.08,
    homogeneous = FALSE,
    periodic = TRUE
) {
  
  df_nuts <- df_nuts[!is.na(df_nuts$start_date), ]
  df_nuts <- ensure_jour(df_nuts)
  
  df_nuts <- df_nuts[df_nuts$Variable %in% var_names, ]
  if (nrow(df_nuts) == 0) return(NULL)
  
  pairs <- combn(var_names, 2, simplify = FALSE)
  
  out <- lapply(pairs, function(pr) {
    v1 <- pr[1]; v2 <- pr[2]
    d1 <- df_nuts[df_nuts$Variable == v1, ]
    d2 <- df_nuts[df_nuts$Variable == v2, ]
    
    if (nrow(d1) < 5 || nrow(d2) < 5) return(NULL)
    
    res <- tryCatch({
      calculate_bivariate_K_GET_yearly(
        df1 = df_nuts,  
        df2 = df_nuts,
        var1_name = v1,
        var2_name = v2,
        max_r = r_value,
        time_col = "Jour",
        T_total = 365,
        n_sim = n_sim,
        directional = FALSE,
        sim_method = sim_method,
        periodic = periodic,
        homogeneous = homogeneous,
        alternative = alternative,
        kde_bw = kde_bw,
        save_plot = FALSE
      )
    }, error = function(e) NULL)
    
    if (is.null(res)) return(NULL)
    
    pv <- res$p_value
    pv_num <- if (length(pv) == 2) pv[2] else pv
    
    data.frame(
      nuts_code = nuts_code,
      var1 = v1,
      var2 = v2,
      type = "intra",
      sim_method = sim_method,
      alternative = alternative,
      r_value = r_value,
      n_sim = n_sim,
      periodic = periodic,
      homogeneous = homogeneous,
      kde_bw = kde_bw,
      p_value = pv_num,
      stringsAsFactors = FALSE
    )
  })
  
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}
