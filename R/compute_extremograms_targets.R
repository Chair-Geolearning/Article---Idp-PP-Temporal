compute_extremograms_targets <- function(
    all_extremes98,
    nuts_ids = c("FRL0", "FI1D"),
    variables = c(
      "total_precipitation_daily_sum",
      "convective_available_potential_energy_daily_maximum",
      "wind_speed_daily_maximum"
    ),
    lag_max = 30,
    n_pixels_max = Inf,
    bootstrap_B = 300,
    signif_B = 300,
    alpha = 0.05,
    signif_method = c("circular_shift", "permute"),
    seed = 42,
    start_date = NULL,
    end_date = NULL,
    out_dir = file.path("output", "extremogramme"),
    make_plots = TRUE,
    recur_K = 7,
    make_fboxplot = TRUE,
    fbplot_max_curves = 400,     
    fbplot_min_finite = 8,       
    fbplot_nbasis = 12           
  stopifnot(is.list(all_extremes98))
  signif_method <- match.arg(signif_method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  make_I_from_dates <- function(dates, start_date, end_date) {
    start_date <- as.Date(start_date)
    end_date   <- as.Date(end_date)
    Tn <- as.integer(end_date - start_date) + 1L
    I <- integer(Tn)
    if (length(dates) == 0) return(I)
    idx <- as.integer(as.Date(dates) - start_date) + 1L
    idx <- idx[idx >= 1L & idx <= Tn]
    if (length(idx) > 0) I[unique(idx)] <- 1L
    I
  }
  
  chi_binary <- function(I, lag_max = 30) {
    Tn <- length(I)
    out <- rep(NA_real_, lag_max)
    for (h in 1:lag_max) {
      if (Tn - h <= 0) next
      I0 <- I[1:(Tn - h)]
      Ih <- I[(1 + h):Tn]
      den <- sum(I0, na.rm = TRUE)
      if (!is.finite(den) || den == 0) out[h] <- NA_real_
      else out[h] <- sum(I0 * Ih, na.rm = TRUE) / den
    }
    out
  }
  
  bootstrap_mean <- function(mat, B = 300, seed = 1, probs = c(0.05, 0.95)) {
    set.seed(seed)
    n <- nrow(mat); L <- ncol(mat)
    if (n == 0) return(list(mean=rep(NA_real_, L), lo=rep(NA_real_, L), hi=rep(NA_real_, L)))
    boot <- matrix(NA_real_, nrow = B, ncol = L)
    for (b in 1:B) {
      idx <- sample.int(n, size = n, replace = TRUE)
      boot[b, ] <- colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
    }
    list(
      mean = colMeans(mat, na.rm = TRUE),
      lo   = apply(boot, 2, stats::quantile, probs = probs[1], na.rm = TRUE),
      hi   = apply(boot, 2, stats::quantile, probs = probs[2], na.rm = TRUE)
    )
  }
  
  nullize_I <- function(I, method = c("circular_shift", "permute")) {
    method <- match.arg(method)
    if (all(is.na(I))) return(I)
    Tn <- length(I)
    if (Tn <= 1) return(I)
    if (method == "permute") sample(I, size = Tn, replace = FALSE)
    else {
      k <- sample.int(Tn, 1)
      c(I[(k+1):Tn], I[1:k])
    }
  }
  
  signif_envelope <- function(I_list, lag_max, B, alpha, method, seed) {
    set.seed(seed)
    n <- length(I_list)
    if (n == 0 || B <= 0) return(list(lo=rep(NA_real_, lag_max), hi=rep(NA_real_, lag_max)))
    sims <- matrix(NA_real_, nrow = B, ncol = lag_max)
    for (b in 1:B) {
      chi_mat_b <- matrix(NA_real_, nrow = n, ncol = lag_max)
      for (i in 1:n) {
        Istar <- nullize_I(I_list[[i]], method = method)
        chi_mat_b[i, ] <- chi_binary(Istar, lag_max = lag_max)
      }
      sims[b, ] <- colMeans(chi_mat_b, na.rm = TRUE)
    }
    list(
      lo = apply(sims, 2, stats::quantile, probs = alpha/2, na.rm = TRUE),
      hi = apply(sims, 2, stats::quantile, probs = 1 - alpha/2, na.rm = TRUE)
    )
  }
  
  condprob_any_1_to_k <- function(chi_h) {
    1 - prod(1 - pmin(pmax(chi_h, 0), 1))
  }
  
  plot_chi_curve <- function(x, y, lo_ci, hi_ci, p_base,
                             lo_sig = NULL, hi_sig = NULL,
                             main, out_png, alpha=0.05) {
    Cairo::CairoPNG(out_png, width = 1600, height = 1000, res = 150)
    graphics::par(mar=c(5,5,4,2)+0.1)
    
    ylim <- range(c(y, lo_ci, hi_ci, p_base, lo_sig, hi_sig), na.rm = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(0, 1)
    
    graphics::plot(x, y, type="l", lwd=3,
                   xlab="Lag (jours)",
                   ylab="chi(h) = P(I[t+h]=1 | I[t]=1)",
                   main=main, ylim=ylim)
    
    graphics::polygon(c(x, rev(x)), c(lo_ci, rev(hi_ci)),
                      border=NA, col=grDevices::adjustcolor("grey70", alpha.f=0.6))
    graphics::lines(x, y, lwd=3)
    graphics::abline(h=p_base, lty=2)
    
    if (!is.null(lo_sig) && !is.null(hi_sig)) {
      graphics::lines(x, lo_sig, lty=3)
      graphics::lines(x, hi_sig, lty=3)
      sig_idx <- which(is.finite(y) & is.finite(hi_sig) & y > hi_sig)
      if (length(sig_idx) > 0) graphics::points(x[sig_idx], y[sig_idx], pch=16, cex=1.2)
      graphics::mtext(sprintf("Enveloppe H0 (%.0f%%) + points si χ(h) > borne sup.", 100*(1-alpha)),
                      side=3, line=0.2, cex=0.9)
    }
    
    grDevices::dev.off()
  }
  
  impute_curve <- function(z) {
    x <- seq_along(z)
    ok <- is.finite(z)
    if (sum(ok) == 0) return(rep(NA_real_, length(z)))
    if (sum(ok) == 1) return(rep(z[ok][1], length(z)))
    z_imp <- z
    z_imp[!ok] <- stats::approx(x[ok], z[ok], xout = x[!ok], rule = 2)$y
    z_imp
  }
  
  plot_chi_fda_fbplot <- function(chi_mat_pixels, main, out_png,
                                  max_curves, min_finite, nbasis, seed_local) {
    Cairo::CairoPNG(out_png, width = 1600, height = 1000, res = 150)
    graphics::par(mar=c(5,5,4,2)+0.1)
    
    if (!requireNamespace("fda", quietly = TRUE)) {
      graphics::plot.new(); graphics::title(main = main)
      graphics::text(0.5, 0.5, "Package 'fda' non installé -> fbplot indisponible.")
      grDevices::dev.off()
      return(invisible(NULL))
    }
    
    finite_counts <- rowSums(is.finite(chi_mat_pixels))
    keep <- which(finite_counts >= min_finite)
    chi <- if (length(keep) > 0) chi_mat_pixels[keep, , drop = FALSE] else chi_mat_pixels[0, , drop = FALSE]
    
    if (nrow(chi) < 3) {
      graphics::plot.new(); graphics::title(main = main)
      graphics::text(0.5, 0.5, "Pas assez de pixels valides (>=3) après filtrage.")
      grDevices::dev.off()
      return(invisible(NULL))
    }
    
    if (is.finite(max_curves) && nrow(chi) > max_curves) {
      set.seed(seed_local)
      chi <- chi[sample.int(nrow(chi), max_curves), , drop = FALSE]
    }
    
    chi_imp <- t(apply(chi, 1, impute_curve))
    keep2 <- which(rowSums(is.finite(chi_imp)) == ncol(chi_imp))
    chi_imp <- if (length(keep2) > 0) chi_imp[keep2, , drop=FALSE] else chi_imp[0, , drop=FALSE]
    
    if (nrow(chi_imp) < 3) {
      graphics::plot.new(); graphics::title(main = main)
      graphics::text(0.5, 0.5, "Pas assez de courbes complètes après imputation.")
      grDevices::dev.off()
      return(invisible(NULL))
    }
    
    chi_imp <- pmin(pmax(chi_imp, 0), 1)
    
    xgrid <- 1:ncol(chi_imp)
    Y <- t(chi_imp) 
    nbasis <- max(3, min(as.integer(nbasis), length(xgrid)))
    
    tryCatch({
      basis <- fda::create.bspline.basis(rangeval = range(xgrid), nbasis = nbasis)
      sm <- fda::smooth.basis(argvals = xgrid, y = Y, fdParobj = basis)
      fit_fd <- sm$fd
      
      if (is.null(fit_fd$coefs) || length(fit_fd$coefs) == 0) stop("fit_fd vide (coefs longueur 0)")
      
      fda::fbplot(
        fit  = fit_fd,
        main = main,
        xlab = "Lag (jours)",
        ylab = "chi(h) par pixel"
      )
    }, error = function(e) {
      graphics::plot.new(); graphics::title(main = main)
      graphics::text(0.5, 0.5, paste("Erreur fda::fbplot:", conditionMessage(e)))
    })
    
    grDevices::dev.off()
  }
  
  results <- list()
  summary_table <- data.table::data.table()
  recur_table <- data.table::data.table()
  
  for (v in variables) {
    if (!v %in% names(all_extremes98)) stop("Variable absente: ", v)
    dtv <- all_extremes98[[v]]
    data.table::setDT(dtv)
    
    for (n in nuts_ids) {
      dt <- dtv[nuts_id == n]
      if (nrow(dt) == 0) next
      
      sd <- if (!is.null(start_date)) as.Date(start_date) else min(dt$date)
      ed <- if (!is.null(end_date))   as.Date(end_date)   else max(dt$date)
      
      coord_all <- unique(dt$coord_id)
      set.seed(seed)
      coord_use <- if (length(coord_all) > n_pixels_max) sample(coord_all, n_pixels_max) else coord_all
      
      chi_mat <- matrix(NA_real_, nrow = length(coord_use), ncol = lag_max)
      p_vec   <- rep(NA_real_, length(coord_use))
      I_list  <- vector("list", length(coord_use))
      
      for (i in seq_along(coord_use)) {
        cid <- coord_use[i]
        dti <- dt[coord_id == cid]
        I <- make_I_from_dates(dti$date, sd, ed)
        I_list[[i]] <- I
        
        p <- mean(I)
        p_vec[i] <- p
        if (!is.finite(p) || p <= 0) next
        
        chi_mat[i, ] <- chi_binary(I, lag_max = lag_max)
      }
      
      keep <- which(rowSums(is.finite(chi_mat)) > 0)
      chi2   <- if (length(keep) > 0) chi_mat[keep, , drop = FALSE] else chi_mat[0, , drop = FALSE]
      p_keep <- p_vec[keep]
      I_keep <- I_list[keep]
      p_mean <- mean(p_keep, na.rm = TRUE)
      
      boot_chi <- if (bootstrap_B > 0) bootstrap_mean(chi2, B = bootstrap_B, seed = seed,
                                                      probs = c(0.05, 0.95)) else
                                                        list(mean = colMeans(chi2, na.rm = TRUE),
                                                             lo = rep(NA_real_, lag_max),
                                                             hi = rep(NA_real_, lag_max))
      
      sig_env <- if (signif_B > 0) signif_envelope(I_keep, lag_max, signif_B, alpha,
                                                   method = signif_method, seed = seed + 1) else
                                                     list(lo = rep(NA_real_, lag_max), hi = rep(NA_real_, lag_max))
      
      K <- min(recur_K, lag_max)
      any_probs <- sapply(1:K, function(k) condprob_any_1_to_k(boot_chi$mean[1:k]))
      
      safe_v <- gsub("[^A-Za-z0-9]+", "_", v)
      out_chi <- file.path(out_dir, sprintf("chi_%s_%s.png", n, safe_v))
      out_fbx <- file.path(out_dir, sprintf("fbplot_fda_%s_%s.png", n, safe_v))
      
      if (isTRUE(make_plots)) {
        plot_chi_curve(
          x = 1:lag_max,
          y = boot_chi$mean,
          lo_ci = boot_chi$lo,
          hi_ci = boot_chi$hi,
          p_base = p_mean,
          lo_sig = sig_env$lo,
          hi_sig = sig_env$hi,
          main = sprintf("Extremogramme χ(h) | %s | %s", n, v),
          out_png = out_chi,
          alpha = alpha
        )
      }
      
      if (isTRUE(make_fboxplot)) {
        plot_chi_fda_fbplot(
          chi_mat_pixels = chi2,
          main = sprintf("Functional boxplot (fda::fbplot) des χ(h) | %s | %s", n, v),
          out_png = out_fbx,
          max_curves = fbplot_max_curves,
          min_finite = fbplot_min_finite,
          nbasis = fbplot_nbasis,
          seed_local = seed + 123
        )
      }
      
      key <- paste(v, n, sep="__")
      results[[key]] <- list(
        nuts_id = n,
        variable = v,
        start_date = sd,
        end_date = ed,
        n_pixels_used = nrow(chi2),
        lag_max = lag_max,
        p_mean = p_mean,
        chi_mean = boot_chi$mean,
        chi_lo = boot_chi$lo,
        chi_hi = boot_chi$hi,
        signif_lo = sig_env$lo,
        signif_hi = sig_env$hi,
        condprob_any_1_to_k = any_probs,
        plot_chi = if (isTRUE(make_plots)) out_chi else NA_character_,
        plot_fda_fbplot = if (isTRUE(make_fboxplot)) out_fbx else NA_character_
      )
      
      summary_table <- data.table::rbindlist(
        list(summary_table,
             data.table::data.table(
               nuts_id = n,
               variable = v,
               start_date = sd,
               end_date = ed,
               n_pixels_used = nrow(chi2),
               p_mean = p_mean,
               chi_h1 = boot_chi$mean[1],
               chi_h2 = if (lag_max >= 2) boot_chi$mean[2] else NA_real_,
               any_within_1 = any_probs[1],
               any_within_2 = if (K >= 2) any_probs[2] else NA_real_,
               any_within_3 = if (K >= 3) any_probs[3] else NA_real_,
               plot_chi = if (isTRUE(make_plots)) out_chi else NA_character_,
               plot_fda_fbplot = if (isTRUE(make_fboxplot)) out_fbx else NA_character_
             )),
        use.names=TRUE, fill=TRUE
      )
      
      recur_table <- data.table::rbindlist(
        list(recur_table,
             data.table::data.table(
               nuts_id = n,
               variable = v,
               k = 1:length(any_probs),
               any_within_k = as.numeric(any_probs)
             )),
        use.names=TRUE, fill=TRUE
      )
    }
  }
  
  list(results = results, summary = summary_table, recur_summary = recur_table)
}
