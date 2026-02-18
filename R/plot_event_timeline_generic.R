plot_event_timeline_generic <- function(
    variables_list, var_names, nuts_target, year,
    output_file = NULL,
    make_plot = FALSE
) {
  library(dplyr)
  library(data.table)
  library(ggplot2)
  
  get_events_df <- function(var_obj) {
    if (is.null(var_obj)) return(NULL)
    if (is.list(var_obj) && "events" %in% names(var_obj)) {
      ev <- var_obj$events
    } else if (is.data.frame(var_obj)) {
      ev <- var_obj
    } else if (is.list(var_obj) && length(var_obj) == 1 && is.data.frame(var_obj[[1]])) {
      ev <- var_obj[[1]]
    } else {
      ev <- NULL
    }
    
    if (!is.null(ev)) {
      if (!("start_date" %in% names(ev)) && ("start" %in% names(ev))) ev <- ev %>% rename(start_date = start)
      if (!("end_date" %in% names(ev)) && ("end" %in% names(ev))) ev <- ev %>% rename(end_date = end)
      if (!all(c("nuts_id","start_date","end_date") %in% names(ev))) return(NULL)
      if (!("cluster_id" %in% names(ev))) ev$cluster_id <- NA_character_
    }
    ev
  }
  
  process_events <- function(event_df) {
    if (is.null(event_df) || nrow(event_df) == 0) {
      return(data.frame(
        cluster_id = character(),
        nuts_id = character(),
        start_date = as.Date(character()),
        end_date = as.Date(character()),
        duration = numeric()
      ))
    }
    
    event_df %>%
      group_by(cluster_id, nuts_id) %>%
      summarise(
        start_date = min(as.Date(start_date), na.rm = TRUE),
        end_date   = max(as.Date(end_date), na.rm = TRUE),
        duration   = as.numeric(max(as.Date(end_date), na.rm = TRUE) -
                                  min(as.Date(start_date), na.rm = TRUE)) + 1,
        .groups = "drop"
      )
  }
  
  all_events <- vector("list", length(variables_list))
  for (i in seq_along(variables_list)) {
    ev_raw <- get_events_df(variables_list[[i]])
    label  <- var_names[i]
    
    ev <- process_events(ev_raw) %>%
      filter(nuts_id == nuts_target) %>%
      mutate(Variable = label)
    
    all_events[[i]] <- ev
  }
  
  events_df <- bind_rows(all_events) %>%
    mutate(
      Jour = as.numeric(format(as.Date(start_date), "%j")),
      Magnitude = duration,
      Variable = factor(Variable, levels = rev(var_names))
    )
  
  message("Nombre d'événements pour ", nuts_target, " : ", nrow(events_df))
  
  if (isTRUE(make_plot) || !is.null(output_file)) {
    if (!is.null(output_file)) {
      dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
      Cairo::CairoPNG(filename = output_file, width = 1200, height = 800, bg = "white")
      on.exit(dev.off(), add = TRUE)
    }
    
    pl <- ggplot(events_df, aes(x = Jour, y = Variable)) +
      theme_bw(base_size = 20) +
      scale_x_continuous(
        limits = c(0, 365),
        breaks = c(15,45,75,105,135,165,195,225,255,285,315,345),
        labels = month.abb
      ) +
      geom_point(aes(size = Magnitude, color = Variable), shape = 1, alpha = 1, stroke = 2) +
      scale_size_continuous(range = c(1, 40)) +
      labs(
        title = paste("Events with Duration -", nuts_target, year),
        x = "Month", y = "Type of events", size = "Duration (days)"
      ) +
      theme(
        legend.position = "right",
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      ) +
      scale_color_discrete(guide = guide_legend(reverse = FALSE, override.aes = list(size = 15)))
    
    print(pl)
  }
  
  return(events_df)
}
plot_event_timeline_beeswarm_year_precip <- function(variables_list, var_names,
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
  pal["Severe precipitation"] <- "#E69F00"  
  
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
