plot_event_timeline_beeswarm <- function(variables_list, var_names, nuts_target, year, output_file) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggbeeswarm)
    library(Cairo)
    library(scico) 
  })
  
  stopifnot(length(variables_list) == length(var_names))
  
  get_events_df <- function(var_obj) {
    if (is.null(var_obj)) return(NULL)
    
    if (is.list(var_obj) && "events" %in% names(var_obj) && is.data.frame(var_obj$events)) {
      ev <- var_obj$events
    } else if (is.data.frame(var_obj)) {
      ev <- var_obj
    } else if (is.list(var_obj) && length(var_obj) == 1 && is.data.frame(var_obj[[1]])) {
      ev <- var_obj[[1]]
    } else {
      return(NULL)
    }
    
    if (!("start_date" %in% names(ev)) && ("start" %in% names(ev))) ev <- ev %>% rename(start_date = start)
    if (!("end_date" %in% names(ev)) && ("end" %in% names(ev))) ev <- ev %>% rename(end_date = end)
    
    if (!all(c("nuts_id", "start_date", "end_date") %in% names(ev))) return(NULL)
    
    ev %>%
      mutate(
        start_date = as.Date(.data$start_date),
        end_date   = as.Date(.data$end_date)
      )
  }
  
  process_events <- function(event_df) {
    if (is.null(event_df) || nrow(event_df) == 0) {
      return(tibble(
        cluster_id  = character(),
        nuts_id     = character(),
        start_date  = as.Date(character()),
        end_date    = as.Date(character()),
        duration    = numeric()
      ))
    }
    
    if (!("cluster_id" %in% names(event_df))) {
      event_df <- event_df %>% mutate(cluster_id = paste0("event_", dplyr::row_number()))
    }
    
    event_df %>%
      group_by(cluster_id, nuts_id) %>%
      summarise(
        start_date = min(start_date, na.rm = TRUE),
        end_date   = max(end_date, na.rm = TRUE),
        duration   = as.numeric(end_date - start_date) + 1,
        .groups    = "drop"
      )
  }
  
  all_events <- vector("list", length(variables_list))
  for (i in seq_along(variables_list)) {
    label <- var_names[[i]]
    
    ev_raw <- get_events_df(variables_list[[i]])
    ev <- process_events(ev_raw) %>%
      filter(.data$nuts_id == nuts_target) %>%
      mutate(Variable = label)
    
    all_events[[i]] <- ev
  }
  
  events_df <- bind_rows(all_events) %>%
    mutate(
      Jour      = as.integer(format(.data$start_date, "%j")),
      Magnitude = .data$duration,
      Variable  = factor(.data$Variable, levels = rev(var_names))
    )
  
  message("Nombre d'événements pour ", nuts_target, " : ", nrow(events_df))
  
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  Cairo::CairoPNG(filename = output_file, width = 1200, height = 800, res = 110, bg = "white")
  
  pal <- scico::scico(n = length(var_names), palette = "batlow")
  names(pal) <- rev(var_names)
  
  pl <- ggplot(events_df, aes(x = Jour, y = Variable)) +
    theme_bw(base_size = 20, base_family = "sans") +
    scale_x_continuous(
      limits = c(0, 366),
      breaks = c(15,45,75,105,135,165,195,225,255,285,315,345),
      labels = month.abb
    ) +
    geom_quasirandom(
      aes(size = Magnitude, color = Variable),
      shape = 1,
      alpha = 0.6,
      stroke = 1.2,       
      groupOnX = FALSE,
      dodge.width = 0.5,
      varwidth = TRUE,
      method = "quasirandom"
    ) +
    scale_size_continuous(
      range = c(0.5, 25),  
      name = "Duration (days)"
    ) +
    scale_color_manual(
      values = pal,
      guide = "none"
    ) +
    labs(
      x = "Month",
      y = "Type of events"
    ) +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  print(pl)
  dev.off()
  
  invisible(events_df)
}
plot_event_timeline_beeswarm_year <- function(variables_list, var_names,
                                              nuts_target, year,
                                              output_file,
                                              width = 1200, height = 800, res = 110,
                                              size_range = c(1, 40)) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggbeeswarm)
    library(lubridate)
    library(Cairo)
    
    library(scico)
  })
  
  stopifnot(length(variables_list) == length(var_names))
  if (!is.numeric(year) || length(year) != 1) {
    stop("'year' must be a single numeric value, e.g. year = 2019")
  }
  
  get_events_df <- function(var_obj) {
    if (is.null(var_obj)) return(NULL)
    
    if (is.list(var_obj) && "events" %in% names(var_obj) && is.data.frame(var_obj$events)) {
      ev <- var_obj$events
    } else if (is.data.frame(var_obj)) {
      ev <- var_obj
    } else if (is.list(var_obj) && length(var_obj) == 1 && is.data.frame(var_obj[[1]])) {
      ev <- var_obj[[1]]
    } else {
      return(NULL)
    }
    
    if (!("start_date" %in% names(ev)) && ("start" %in% names(ev))) ev <- ev %>% rename(start_date = start)
    if (!("end_date" %in% names(ev)) && ("end" %in% names(ev))) ev <- ev %>% rename(end_date = end)
    
    if (!all(c("nuts_id", "start_date", "end_date") %in% names(ev))) return(NULL)
    
    ev %>%
      mutate(
        start_date = as.Date(.data$start_date),
        end_date   = as.Date(.data$end_date)
      )
  }
  
  process_events <- function(event_df) {
    if (is.null(event_df) || nrow(event_df) == 0) {
      return(tibble(
        cluster_id = character(),
        nuts_id = character(),
        start_date = as.Date(character()),
        end_date = as.Date(character()),
        duration = numeric()
      ))
    }
    
    if (!("cluster_id" %in% names(event_df))) {
      event_df <- event_df %>% mutate(cluster_id = paste0("event_", dplyr::row_number()))
    }
    
    event_df %>%
      group_by(cluster_id, nuts_id) %>%
      summarise(
        start_date = min(start_date, na.rm = TRUE),
        end_date   = max(end_date, na.rm = TRUE),
        duration   = as.numeric(end_date - start_date) + 1,
        .groups = "drop"
      )
  }
  
  all_events <- vector("list", length(variables_list))
  
  for (i in seq_along(variables_list)) {
    label <- var_names[[i]]
    ev_raw  <- get_events_df(variables_list[[i]])
    ev_proc <- process_events(ev_raw)
    
    all_events[[i]] <- ev_proc %>%
      filter(.data$nuts_id == nuts_target) %>%
      filter(lubridate::year(.data$start_date) == year) %>%
      mutate(Variable = label)
  }
  
  events_df <- bind_rows(all_events) %>%
    mutate(
      Jour = as.integer(format(.data$start_date, "%j")),
      Magnitude = .data$duration,
      Variable = factor(.data$Variable, levels = rev(var_names))
    )
  
  message("Events for ", nuts_target, " in ", year, ": ", nrow(events_df))
  
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  Cairo::CairoPNG(filename = output_file, width = width, height = height, res = res, bg = "white")
  
  pal <- scico::scico(n = length(var_names), palette = "batlow")
  names(pal) <- rev(var_names)
  
  pl <- ggplot(events_df, aes(x = Jour, y = Variable)) +
    theme_bw(base_size = 20) +
    scale_x_continuous(
      limits = c(0, 366),
      breaks = c(15,45,75,105,135,165,195,225,255,285,315,345),
      labels = month.abb
    ) +
    geom_quasirandom(
      aes(size = Magnitude, colour = Variable),
      shape = 1,
      alpha = 1,
      stroke = 1.2,
      groupOnX = FALSE,
      dodge.width = 0.5,
      varwidth = FALSE,
      method = "quasirandom"
    ) +
    scale_size_continuous(
      range = c(2, 15),
      breaks = c(10, 20, 30),
      name = "Duration (days)"
    ) +
    scale_colour_manual(
      values = pal,
      guide = "none"   
    ) +
    labs(
      x = "Month",
      y = "Event type"
    ) +
    theme(
      legend.position = "right",  
      axis.text.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  print(pl)
  dev.off()
  
  invisible(events_df)
}

plot_event_timeline_beeswarm_year <- function(variables_list, var_names,
                                              nuts_target, year,
                                              output_file,
                                              width = 1200, height = 800, res = 110,
                                              size_range = c(1, 40)) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggbeeswarm)
    library(lubridate)
    library(Cairo)
   
    library(scico)
  })
  
  stopifnot(length(variables_list) == length(var_names))
  if (!is.numeric(year) || length(year) != 1) {
    stop("'year' must be a single numeric value, e.g. year = 2019")
  }
  
  get_events_df <- function(var_obj) {
    if (is.null(var_obj)) return(NULL)
    
    if (is.list(var_obj) && "events" %in% names(var_obj) && is.data.frame(var_obj$events)) {
      ev <- var_obj$events
    } else if (is.data.frame(var_obj)) {
      ev <- var_obj
    } else if (is.list(var_obj) && length(var_obj) == 1 && is.data.frame(var_obj[[1]])) {
      ev <- var_obj[[1]]
    } else {
      return(NULL)
    }
    
    if (!("start_date" %in% names(ev)) && ("start" %in% names(ev))) ev <- ev %>% rename(start_date = start)
    if (!("end_date" %in% names(ev)) && ("end" %in% names(ev))) ev <- ev %>% rename(end_date = end)
    
    if (!all(c("nuts_id", "start_date", "end_date") %in% names(ev))) return(NULL)
    
    ev %>%
      mutate(
        start_date = as.Date(.data$start_date),
        end_date   = as.Date(.data$end_date)
      )
  }
  
  process_events <- function(event_df) {
    if (is.null(event_df) || nrow(event_df) == 0) {
      return(tibble(
        cluster_id = character(),
        nuts_id = character(),
        start_date = as.Date(character()),
        end_date = as.Date(character()),
        duration = numeric()
      ))
    }
    
    if (!("cluster_id" %in% names(event_df))) {
      event_df <- event_df %>% mutate(cluster_id = paste0("event_", dplyr::row_number()))
    }
    
    event_df %>%
      group_by(cluster_id, nuts_id) %>%
      summarise(
        start_date = min(start_date, na.rm = TRUE),
        end_date   = max(end_date, na.rm = TRUE),
        duration   = as.numeric(end_date - start_date) + 1,
        .groups = "drop"
      )
  }
  
  all_events <- vector("list", length(variables_list))
  
  for (i in seq_along(variables_list)) {
    label <- var_names[[i]]
    ev_raw  <- get_events_df(variables_list[[i]])
    ev_proc <- process_events(ev_raw)
    
    all_events[[i]] <- ev_proc %>%
      filter(.data$nuts_id == nuts_target) %>%
      filter(lubridate::year(.data$start_date) == year) %>%
      mutate(Variable = label)
  }
  
  events_df <- bind_rows(all_events) %>%
    mutate(
      Jour = as.integer(format(.data$start_date, "%j")),
      Magnitude = .data$duration,
      Variable = factor(.data$Variable, levels = rev(var_names))
    )
  
  message("Events for ", nuts_target, " in ", year, ": ", nrow(events_df))
  
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  Cairo::CairoPNG(filename = output_file, width = width, height = height, res = res, bg = "white")
  
  pal <- scico::scico(n = length(var_names), palette = "batlow")
  names(pal) <- rev(var_names)
  
  pl <- ggplot(events_df, aes(x = Jour, y = Variable)) +
    theme_bw(base_size = 20) +
    scale_x_continuous(
      limits = c(0, 366),
      breaks = c(15,45,75,105,135,165,195,225,255,285,315,345),
      labels = month.abb
    ) +
    geom_quasirandom(
      aes(size = Magnitude, colour = Variable),
      shape = 1,
      alpha = 1,
      stroke = 1.2,
      groupOnX = FALSE,
      dodge.width = 0.5,
      varwidth = FALSE,
      method = "quasirandom"
    ) +
    scale_size_continuous(
      range = c(2, 15),
      breaks = c(10, 20, 30),
      name = "Duration (days)"
    ) +
    scale_colour_manual(
      values = pal,
      guide = "none"  
    ) +
    labs(
      x = "Month",
      y = "Event type"
    ) +
    theme(
      legend.position = "right",  
      axis.text.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  print(pl)
  dev.off()
  
  invisible(events_df)
}
plot_marked_cumulative_trajectories_yearly <- function(
    df,
    nuts_target,
    variable_target,
    year_col = NULL,
    date_col = "start_date",
    doy_col = NULL,
    mark_col = "Magnitude",          
    agg_same_day = TRUE,
    highlight_year = 2022,
    top_k_extremes = 5,
    x_breaks = c(15,45,75,105,135,165,195,225,255,285,315,345),
    x_labels = month.abb,
    output_file = NULL,
    width = 1400, height = 900, res = 110,
    title = NULL,
    subtitle = NULL,
    add_year_zero = TRUE,
    
    show_steps = TRUE,
    step_direction = "hv",         
    step_size_other = 0.35,
    step_size_extreme = 0.55,
    step_size_highlight = 1.10,
    step_alpha_other = 0.18,
    step_alpha_extreme = 0.45,
    step_alpha_highlight = 0.95,
    
    show_points = TRUE,
    point_size_other = 0.9,
    point_size_extreme = 1.3,
    point_size_highlight = 2.2,
    point_alpha_other = 0.20,
    point_alpha_extreme = 0.55,
    point_alpha_highlight = 1.00,
    
    shape_other = 16,              
    shape_extreme = 16,
    shape_highlight = 23,          
    highlight_fill = "#56B4E9",       
    highlight_outline = "black",
    highlight_stroke = 0.7,
    
    show_rug = TRUE,
    rug_alpha_other = 0.10,
    rug_alpha_extreme = 0.25,
    rug_alpha_highlight = 0.55
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(lubridate)
    library(Cairo)
    library(tibble)
  })
  
  stopifnot(is.data.frame(df))
  if (!all(c("nuts_id", "Variable", date_col) %in% names(df))) {
    stop("df must contain columns: nuts_id, Variable, and ", date_col)
  }
  
  d <- df %>%
    filter(.data$nuts_id == nuts_target,
           as.character(.data$Variable) == variable_target) %>%
    mutate(
      .date = as.Date(.data[[date_col]]),
      .year = if (is.null(year_col)) lubridate::year(.date) else as.integer(.data[[year_col]]),
      .doy  = if (is.null(doy_col)) as.integer(format(.date, "%j")) else as.integer(.data[[doy_col]])
    )
  
  if (nrow(d) == 0) stop("No events after filtering.")
  
  if (is.null(mark_col)) {
    d <- d %>% mutate(.mark = 1)
    y_lab <- "Cumulative number of events"
  } else {
    if (!mark_col %in% names(d)) stop("mark_col not found in df: ", mark_col)
    d <- d %>% mutate(.mark = suppressWarnings(as.numeric(.data[[mark_col]])))
    y_lab <- paste0("Cumulative sum of duration")
  }
  
  if (isTRUE(agg_same_day)) {
    d <- d %>%
      group_by(.year, .doy) %>%
      summarise(.mark = sum(.mark, na.rm = TRUE), .groups = "drop")
  } else {
    d <- d %>% select(.year, .doy, .mark)
  }
  
  traj <- d %>%
    group_by(.year) %>%
    arrange(.doy, .by_group = TRUE) %>%
    mutate(.cum = cumsum(.mark)) %>%
    ungroup()
  
  if (isTRUE(add_year_zero)) {
    traj <- traj %>%
      group_by(.year) %>%
      group_modify(~{
        first_doy <- suppressWarnings(min(.x$.doy, na.rm = TRUE))
        if (is.finite(first_doy) && first_doy > 1) {
          bind_rows(
            tibble(.year = .x$.year[1], .doy = 1L, .mark = 0, .cum = 0),
            .x
          )
        } else .x
      }) %>%
      arrange(.year, .doy) %>%
      ungroup()
  }
  
  final_by_year <- traj %>%
    group_by(.year) %>%
    summarise(final = max(.cum, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(final))
  
  extreme_years <- head(final_by_year$.year, top_k_extremes)
  
  traj <- traj %>%
    mutate(.group = case_when(
      .year == highlight_year ~ "highlight",
      .year %in% extreme_years ~ "extreme",
      TRUE ~ "other"
    ))
  
  p <- ggplot(traj, aes(x = .doy, y = .cum)) +
    theme_bw(base_size = 18, base_family = "sans") +
    scale_x_continuous(limits = c(1, 365), breaks = x_breaks, labels = x_labels) +
    labs(x = "Month (day of year)", y = y_lab) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  
  if (isTRUE(show_steps)) {
    p <- p +
      geom_step(
        data = dplyr::filter(traj, .group == "other"),
        aes(group = .year),
        direction = step_direction,
        linewidth = step_size_other,
        alpha = step_alpha_other,
        colour = "grey45"
      ) +
      geom_step(
        data = dplyr::filter(traj, .group == "extreme"),
        aes(group = .year),
        direction = step_direction,
        linewidth = step_size_extreme,
        alpha = step_alpha_extreme,
        colour = "grey15"
      ) +
      geom_step(
        data = dplyr::filter(traj, .group == "highlight"),
        aes(group = .year),
        direction = step_direction,
        linewidth = step_size_highlight,
        alpha = step_alpha_highlight,
        colour = highlight_outline
      )
  }
  
  if (isTRUE(show_points)) {
    p <- p +
      geom_point(
        data = dplyr::filter(traj, .group == "other" & .mark > 0),
        shape = shape_other,
        colour = "grey45",
        size = point_size_other,
        alpha = point_alpha_other
      ) +
      geom_point(
        data = dplyr::filter(traj, .group == "extreme" & .mark > 0),
        shape = shape_extreme,
        colour = "grey10",
        size = point_size_extreme,
        alpha = point_alpha_extreme
      ) +
      geom_point(
        data = dplyr::filter(traj, .group == "highlight" & .mark > 0),
        shape = shape_highlight,
        fill = highlight_fill,
        colour = highlight_outline,
        stroke = highlight_stroke,
        size = point_size_highlight,
        alpha = point_alpha_highlight
      )
  }
  
  if (isTRUE(show_rug)) {
    p <- p +
      geom_rug(
        data = dplyr::filter(traj, .group == "other" & .mark > 0),
        aes(x = .doy), inherit.aes = FALSE,
        sides = "b", alpha = rug_alpha_other, colour = "grey45"
      ) +
      geom_rug(
        data = dplyr::filter(traj, .group == "extreme" & .mark > 0),
        aes(x = .doy), inherit.aes = FALSE,
        sides = "b", alpha = rug_alpha_extreme, colour = "grey10"
      ) +
      geom_rug(
        data = dplyr::filter(traj, .group == "highlight" & .mark > 0),
        aes(x = .doy), inherit.aes = FALSE,
        sides = "b", alpha = rug_alpha_highlight, colour = highlight_outline
      )
  }
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    Cairo::CairoPNG(filename = output_file,
                    width = width, height = height, res = res, bg = "white")
    print(p)
    dev.off()
  }
  
  invisible(list(plot = p, data = traj, final_by_year = final_by_year))
}
