calculate_bivariate_K_GET_yearly <- function(
    df1, df2,
    var1_name, var2_name,
    max_r,
    time_col = "Jour",
    T_total = 365,
    n_sim = 999,
    directional = FALSE,
    sim_method = c("translation", "permutation", "nhpp"),
    homogeneous = TRUE,
    periodic = TRUE,
    alternative = c("two.sided", "less", "greater"),
    seed = NULL,
    kde_bw = "nrd0",
    kde_n = 365,
    lambda_floor = 1e-6,
    label_x = NULL,
    label_y = NULL,
    center = FALSE,
    plot_title = NULL,
    save_plot = FALSE,
    save_path = NULL,
    width = 900,
    height = 650,
    res = 120
) {
  sim_method <- match.arg(sim_method)
  alternative <- match.arg(alternative)
  if (!is.null(seed)) set.seed(seed)
  
  if (!requireNamespace("GET", quietly = TRUE)) stop("Package 'GET' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package 'lubridate' is required.")
  
  has_stringr <- requireNamespace("stringr", quietly = TRUE)
  
  .wrap <- function(x, width = 70) {
    if (has_stringr) stringr::str_wrap(x, width = width) else paste(strwrap(x, width = width), collapse = "\n")
  }
  
  .check_cols <- function(df, needed) {
    miss <- setdiff(needed, names(df))
    if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))
  }
  
  .dcirc <- function(u, v) {
    d <- abs(outer(u, v, "-"))
    pmin(d, 1 - d)
  }
  
  .dforward <- function(u, v) {
    (outer(v, u, "-") %% 1) |> t()
  }
  
  .estimate_lambda_doy_kde <- function(u, years, bw = kde_bw, n = kde_n, floor = lambda_floor) {
    u <- u %% 1
    u_ext <- c(u - 1, u, u + 1)
    
    dens <- stats::density(u_ext, bw = bw, from = 0, to = 1, n = n)
    fhat <- pmax(dens$y, 1e-12)
    
    nyears <- length(unique(years))
    if (nyears < 1) nyears <- 1
    nbar <- length(u) / nyears              # mean events/year
    
    lambda_grid <- nbar * fhat              # events/year on [0,1)
    lambda_fun <- function(u_query) {
      uq <- (u_query %% 1)
      stats::approx(dens$x, pmax(lambda_grid, floor), xout = uq, rule = 2)$y
    }
    
    list(x = dens$x, lambda = lambda_grid, lambda_fun = lambda_fun, nbar = nbar)
  }
  
  .build_rescaling_maps <- function(lambda_est) {
    x <- lambda_est$x
    lam <- pmax(lambda_est$lambda, lambda_floor)
    
    dx <- diff(x)
    mid <- (lam[-1] + lam[-length(lam)]) / 2
    cum <- c(0, cumsum(mid * dx))
    total <- cum[length(cum)]
    if (total <= 0) stop("Invalid lambda integral (<=0).")
    
    F <- cum / total
    eps <- 1e-12
    for (i in 2:length(F)) if (F[i] <= F[i - 1]) F[i] <- F[i - 1] + eps
    F <- F / max(F)
    
    u_to_z <- function(u) {
      uu <- u %% 1
      stats::approx(x, F, xout = uu, rule = 2)$y
    }
    z_to_u <- function(z) {
      zz <- z %% 1
      stats::approx(F, x, xout = zz, rule = 2)$y %% 1
    }
    
    list(u_to_z = u_to_z, z_to_u = z_to_u)
  }
  
  .simulate_nhpp_by_year <- function(years_common, lambda_est) {
    w <- pmax(lambda_est$lambda, lambda_floor)
    cdf <- cumsum(w) / sum(w)
    xgrid <- lambda_est$x
    
    y_year_sim <- integer(0)
    y_u_sim <- numeric(0)
    
    for (yr in years_common) {
      nyr <- stats::rpois(1, lambda_est$nbar)
      if (nyr <= 0) next
      u01 <- runif(nyr)
      idx <- findInterval(u01, cdf) + 1
      idx[idx < 1] <- 1
      idx[idx > length(xgrid)] <- length(xgrid)
      uu <- (xgrid[idx] + runif(nyr, -0.5 / kde_n, 0.5 / kde_n)) %% 1
      y_u_sim <- c(y_u_sim, uu)
      y_year_sim <- c(y_year_sim, rep(yr, nyr))
    }
    list(year = y_year_sim, u = y_u_sim)
  }
  
 
  .K_total <- function(x_year, x_u, y_year, y_u, r_frac,
                       directional = FALSE,
                       periodic = TRUE,
                       homogeneous = TRUE,
                       lambda_x_fun = NULL,
                       lambda_y_fun = NULL,
                       floor = lambda_floor) {
    
    y0 <- min(c(x_year, y_year), na.rm = TRUE)
    tx <- (x_year - y0) + (x_u %% 1)
    ty <- (y_year - y0) + (y_u %% 1)
    
    nx <- length(tx); ny <- length(ty)
    if (nx == 0 || ny == 0) return(rep(NA_real_, length(r_frac)))
    
    L <- (max(c(x_year, y_year), na.rm = TRUE) - y0) + 1
    if (L <= 0) return(rep(NA_real_, length(r_frac)))
    
    dcirc_L <- function(a, b, L) {
      d <- abs(outer(a, b, "-"))
      pmin(d, L - d)
    }
    dforward_L <- function(a, b, L) {
      (outer(b, a, "-") %% L) |> t()
    }
    
    out <- numeric(length(r_frac))
    
    for (rr in seq_along(r_frac)) {
      r <- r_frac[rr]
      
      if (periodic) {
        if (!directional) {
          ind <- (dcirc_L(tx, ty, L) <= r)
        } else {
          ind <- (dforward_L(tx, ty, L) <= r)
        }
      } else {
        d <- abs(outer(tx, ty, "-"))
        if (!directional) {
          ind <- (d <= r)
        } else {
          ind <- (outer(ty, tx, "-") >= 0) & (d <= r)
        }
      }
      
      if (homogeneous) {
      
        out[rr] <- (L / (nx * ny)) * sum(ind)
      } else {
        ux <- x_u %% 1
        uy <- y_u %% 1
        lamx <- pmax(lambda_x_fun(ux), floor)
        lamy <- pmax(lambda_y_fun(uy), floor)
        W <- 1 / (outer(lamx, lamy, "*"))
        out[rr] <- (1 / L) * sum(W * ind)
      }
    }
    
    out
  }
  
  needed <- c("Variable", "start_date", time_col)
  .check_cols(df1, needed)
  .check_cols(df2, needed)
  
  d1 <- df1 |>
    dplyr::filter(.data$Variable == var1_name) |>
    dplyr::mutate(
      year = lubridate::year(.data$start_date),
      doy = as.numeric(.data[[time_col]]),
      u = ((doy - 1) / T_total) %% 1
    )
  
  d2 <- df2 |>
    dplyr::filter(.data$Variable == var2_name) |>
    dplyr::mutate(
      year = lubridate::year(.data$start_date),
      doy = as.numeric(.data[[time_col]]),
      u = ((doy - 1) / T_total) %% 1
    )
  
  if (nrow(d1) == 0) stop("No events found in df1 for var1_name = ", var1_name)
  if (nrow(d2) == 0) stop("No events found in df2 for var2_name = ", var2_name)
  
  x_year <- d1$year
  x_u <- d1$u
  y_year <- d2$year
  y_u <- d2$u
  
  years_common <- sort(intersect(unique(x_year), unique(y_year)))
  if (length(years_common) == 0) stop("No common years between X and Y.")
  
  r_days <- seq_len(max_r)
  r_frac <- r_days / T_total
  
  if (!homogeneous) {
    lambda_x_est <- .estimate_lambda_doy_kde(x_u, x_year)
    lambda_y_est <- .estimate_lambda_doy_kde(y_u, y_year)
    lambda_x_fun <- lambda_x_est$lambda_fun
    lambda_y_fun <- lambda_y_est$lambda_fun
    y_maps <- .build_rescaling_maps(lambda_y_est) 
  } else {
    lambda_x_est <- NULL; lambda_y_est <- NULL
    lambda_x_fun <- NULL; lambda_y_fun <- NULL
    y_maps <- NULL
  }
  
  K_obs <- .K_total(
    x_year, x_u, y_year, y_u, r_frac,
    directional = directional,
    periodic = periodic,
    homogeneous = homogeneous,
    lambda_x_fun = lambda_x_fun,
    lambda_y_fun = lambda_y_fun
  )
  
  K_sim <- matrix(NA_real_, nrow = length(r_frac), ncol = n_sim)
  n_sim_events <- integer(n_sim)
  
  for (i in seq_len(n_sim)) {
    if (sim_method == "translation") {
      y_year_sim <- y_year
      y_u_sim <- y_u
      shifts <- stats::runif(length(years_common), 0, 1)
      names(shifts) <- as.character(years_common)
      
      for (yr in years_common) {
        idx <- which(y_year == yr)
        if (length(idx) == 0) next
        if (homogeneous) {
          y_u_sim[idx] <- (y_u_sim[idx] + shifts[as.character(yr)]) %% 1
        } else {
          z <- y_maps$u_to_z(y_u[idx])
          z2 <- (z + shifts[as.character(yr)]) %% 1
          y_u_sim[idx] <- y_maps$z_to_u(z2)
        }
      }
      
    } else if (sim_method == "permutation") {
      y_year_sim <- y_year
      y_u_sim <- y_u
      for (yr in years_common) {
        idx <- which(y_year == yr)
        if (length(idx) == 0) next
        y_u_sim[idx] <- sample(y_u_sim[idx], replace = FALSE)
      }
      
    } else if (sim_method == "nhpp") {
      if (homogeneous) {
        nbar <- length(y_u) / length(unique(y_year))
        y_year_sim <- integer(0)
        y_u_sim <- numeric(0)
        for (yr in years_common) {
          nyr <- stats::rpois(1, nbar)
          if (nyr <= 0) next
          y_u_sim <- c(y_u_sim, runif(nyr))
          y_year_sim <- c(y_year_sim, rep(yr, nyr))
        }
      } else {
        sim <- .simulate_nhpp_by_year(years_common, lambda_y_est)
        y_year_sim <- sim$year
        y_u_sim <- sim$u
      }
    } else {
      stop("Unknown sim_method: ", sim_method)
    }
    
    n_sim_events[i] <- length(y_u_sim)
    
    K_sim[, i] <- .K_total(
      x_year, x_u, y_year_sim, y_u_sim, r_frac,
      directional = directional,
      periodic = periodic,
      homogeneous = homogeneous,
      lambda_x_fun = lambda_x_fun,
      lambda_y_fun = lambda_y_fun
    )
  }
  
  curve_r <- r_days
  
  K_theo <- if (directional) r_frac else 2 * r_frac
  
  if (isTRUE(center)) {
    obs_for_test <- K_obs - K_theo
    sim_for_test <- sweep(K_sim, 1, K_theo, "-")
    center_line <- rep(0, length(r_days))
  } else {
    obs_for_test <- K_obs
    sim_for_test <- K_sim
    center_line <- K_theo
  }
  
  m <- length(obs_for_test)
  
  if (is.null(dim(sim_for_test))) stop("sim_for_test must be a matrix.")
  if (nrow(sim_for_test) != m && ncol(sim_for_test) == m) {
    sim_for_test <- t(sim_for_test)
  }
  if (nrow(sim_for_test) != m) {
    stop(sprintf("sim_for_test has wrong dims: %d x %d (expected %d x n_sim)",
                 nrow(sim_for_test), ncol(sim_for_test), m))
  }
  
  curve_set <- GET::create_curve_set(list(
    r = r_days,
    obs = obs_for_test,
    sim_m = sim_for_test
  ))
  
  global_test <- GET::global_envelope_test(curve_set, type = "rank", alternative = alternative)
  p_value <- attr(global_test, "p")
  
  get_lohi <- function(gt) {
    if (!is.null(gt$lo) && !is.null(gt$hi)) return(list(lo = gt$lo, hi = gt$hi))
    if (!is.null(gt$global_envelope) &&
        !is.null(gt$global_envelope$lo) && !is.null(gt$global_envelope$hi)) {
      return(list(lo = gt$global_envelope$lo, hi = gt$global_envelope$hi))
    }
    if (!is.null(gt$envelope) &&
        !is.null(gt$envelope$lo) && !is.null(gt$envelope$hi)) {
      return(list(lo = gt$envelope$lo, hi = gt$envelope$hi))
    }
    return(NULL)
  }
  
  eh <- get_lohi(global_test)
  
  if (is.null(eh)) {
    mean_sim <- rowMeans(sim_for_test, na.rm = TRUE)
    sd_sim   <- apply(sim_for_test, 1, stats::sd, na.rm = TRUE)
    warning("Could not extract GET lo/hi; using mean +/- 1.96*sd as fallback.")
    z <- stats::qnorm(0.975)
    lo <- mean_sim - z * sd_sim
    hi <- mean_sim + z * sd_sim
  } else {
    lo <- eh$lo
    hi <- eh$hi
  }
  
  plot_data <- data.frame(
    r   = r_days,
    obs = obs_for_test,
    lo  = lo,
    hi  = hi
  )
  

  p <- NULL
  intensity_plot <- NULL
  
  hom_short <- if (homogeneous) "H" else "NH"
  dir_short <- if (directional) "X→Y" else "SYM"
  sim_short <- switch(sim_method, "translation" = "trans", "permutation" = "perm", "nhpp" = "nhpp")
  
  pv <- p_value
  p_txt <- if (length(pv) == 2) sprintf("[%.3f, %.3f]", pv[1], pv[2]) else
    ifelse(is.na(pv), "NA", ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv)))
  
  short_subtitle <- .wrap(sprintf("%s | sim=%s | %s | p=%s", hom_short, sim_short, dir_short, p_txt), width = 70)
  short_caption <- .wrap(sprintf("r≤%gd | n_sim=%d", max_r, n_sim), width = 70)
  
  if (is.null(plot_title)) {
    plot_title <- sprintf("%s vs %s", var1_name, var2_name)
    if (!isTRUE(center)) plot_title <- paste0(plot_title, " (K)")
    if (isTRUE(center))  plot_title <- paste0(plot_title, " (K - ref)")
  }
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    subtitle_color <- switch(
      alternative,
      "two.sided" = "#1976D2",
      "greater"   = "#2E7D32",
      "less"      = "#C62828"
    )
    
    ylab <- if (center) {
      if (directional) "K(r) - r" else "K(r) - 2r"
    } else {
      "K(r)"
    }
    
    xlab <- if (directional) {
      sprintf("Temporal distance r (days): %s precedes %s (all years; continuous time)", var1_name, var2_name)
    } else {
      "Temporal distance r (days) (all years; continuous time)"
    }
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = r)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "lightblue", alpha = 0.5) +
      ggplot2::geom_line(ggplot2::aes(y = lo), color = "blue", linewidth = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = hi), color = "blue", linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(y = obs), color = "red", size = 2.2) +
      ggplot2::labs(
        title = plot_title,
        subtitle = short_subtitle,
        caption = short_caption,
        x = xlab,
        y = ylab
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 11, face = "bold", color = subtitle_color),
        plot.caption  = ggplot2::element_text(size = 9, hjust = 0, lineheight = 1.15),
        axis.title    = ggplot2::element_text(size = 11),
        axis.text     = ggplot2::element_text(size = 10),
        plot.margin   = ggplot2::margin(10, 22, 10, 10)
      )
    
    if (!homogeneous) {
      lab_x <- if (!is.null(label_x)) label_x else var1_name
      lab_y <- if (!is.null(label_y)) label_y else var2_name
      
      grid_doy <- seq(1, T_total, length.out = kde_n)
      grid_u <- (grid_doy - 1) / T_total
      
      lx <- lambda_x_fun(grid_u)
      ly <- lambda_y_fun(grid_u)
      
      df_intensity <- data.frame(
        DOY = rep(grid_doy, 2),
        lambda = c(lx, ly),
        series = rep(c(lab_x, lab_y), each = length(grid_doy))
      )
      
      intensity_plot <- ggplot2::ggplot(df_intensity,
                                        ggplot2::aes(x = DOY, y = lambda, color = series)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::labs(
          title = "Seasonal intensity (KDE)",
          subtitle = .wrap(sprintf("bw=%s | periodic DOY | %s vs %s",
                                   if (is.character(kde_bw)) kde_bw else sprintf("%.3f", kde_bw),
                                   lab_x, lab_y), 70),
          x = "Day of Year (DOY)",
          y = expression(hat(lambda)(DOY)~~"(events/year)"),
          color = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 12, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 10),
          legend.position = "bottom",
          plot.margin = ggplot2::margin(10, 22, 10, 10)
        )
    }
    
    if (isTRUE(save_plot) || !is.null(save_path)) {
      if (is.null(save_path)) {
        save_path <- sprintf("K_%s_%s_%s.png", var1_name, var2_name, sim_method)
      }
      
      if (!requireNamespace("Cairo", quietly = TRUE)) {
        warning("Package 'Cairo' not available: using png() fallback.")
        grDevices::png(filename = save_path, width = width, height = height, res = res)
        print(p); grDevices::dev.off()
      } else {
        Cairo::CairoPNG(filename = save_path, width = width, height = height, res = res)
        print(p); grDevices::dev.off()
      }
      
      if (!is.null(intensity_plot)) {
        int_path <- sub("\\.png$", "_intensity.png", save_path)
        if (requireNamespace("Cairo", quietly = TRUE)) {
          Cairo::CairoPNG(filename = int_path, width = width, height = height, res = res)
          print(intensity_plot); grDevices::dev.off()
        } else {
          grDevices::png(filename = int_path, width = width, height = height, res = res)
          print(intensity_plot); grDevices::dev.off()
        }
      }
    }
  } else {
    warning("ggplot2 not installed: skipping plot creation.")
  }
  
  return(list(
    K_obs = K_obs,
    K_sim = K_sim,
    center = center,
    obs_used = obs_for_test,
    sim_used = sim_for_test,
    ref_used = center_line,
    n_sim_events = n_sim_events,
    r_days = r_days,
    r_frac = r_frac,
    K_theo = K_theo,
    curve_set = curve_set,
    global_test = global_test,
    p_value = p_value,
    alternative = alternative,
    var1_name = var1_name,
    var2_name = var2_name,
    directional = directional,
    periodic = periodic,
    sim_method = sim_method,
    homogeneous = homogeneous,
    plot = p,
    intensity_plot = intensity_plot,
    lambda_x = lambda_x_est,
    lambda_y = lambda_y_est
  ))
}
