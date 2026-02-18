#################################################
library(targets)
library(future)
library(future.callr)

# plan(multisession, workers = 30)

tar_option_set(
  packages = c(
    "dplyr", "purrr", "sf", "foreach",
    "ggplot2", "IndTestPP", "Cairo", "giscoR", "tarchetypes", "data.table"
  ),
  memory = "transient",
  format = "qs",
  seed = 1995
)

#################################################################################################
#########                               2. Source script                                #########
#################################################################################################
purrr::walk(list.files("R", full.names = TRUE), source)
#################################################################################################
#########                               3. Targets plan                                 #########
################################################################################################

list(
   tar_target(
    process_summary,
    process_netcdf_files(
      var_names = c(
        "total_precipitation_daily_sum",
        "10m_u_component_of_wind_daily_maximum",
        "10m_v_component_of_wind_daily_maximum",
        "convective_available_potential_energy_daily_maximum",
        "convective_precipitation_daily_sum",
        "2m_temperature_daily_mean"
      ),
      years = 1941:2024,
      nuts_level = 2,
      save_intermediate = TRUE,
      output_dir = "data/processed",
      n_cores = 20,
      use_parallel = TRUE
    )
  ),

  tar_target(
    raw_data_19412024,
    load_processed_data(
      var_names = c(
        "total_precipitation_daily_sum",
        "10m_u_component_of_wind_daily_maximum",
        "10m_v_component_of_wind_daily_maximum",
        "convective_available_potential_energy_daily_maximum",
        "convective_precipitation_daily_sum",
        "2m_temperature_daily_mean"
      ),
      years = 1941:2024,
      input_dir = "data/processed"
    )
  ),


  tar_target(
    all_extremes999,
    process_all_variables(
      years = 1941:2024,
      input_dir = "data/processed",
      output_dir = "data/extremes",
      quantile_precip = 0.995, 
      quantile_cape = 0.995,  
      quantile_wind = 0.995    
    )
  ),
  # tar_target(
  #   extremograms_FRL0_FI1D_annual,
  #   compute_extremograms_targets(
  #     all_extremes98 = all_extremes999,
  #     nuts_ids = c("FRL0", "FI1D"),
  #     variables = c("total_precipitation_daily_sum", "wind_speed_daily_maximum"),
  #     out_dir = file.path("output", "extremogramme", "annual"),
  #     make_plots = TRUE
  #   )
  # ),

  # tar_target(
  #   extremograms_FRL0_FI1D_cape_JJA,
  #   compute_extremograms_targets(
  #     all_extremes98 = all_extremes999,
  #     nuts_ids = c("FRL0", "FI1D"),
  #     variables = c("convective_available_potential_energy_daily_maximum"),
  #     out_dir = file.path("output", "extremogramme", "CAPE_JJA"),
  #     make_plots = TRUE
  #   )
  # ),
  tar_target(
    grouped_extreme_events_precip_19412024_q999,
    group_extreme_events_by_proximity(
      extreme_events = all_extremes999$total_precipitation_daily_sum,
      temporal_width = 2,
      distance_threshold = 50000,  # 50 km
      parallel = TRUE,
      num_cores = 40
    ),
    deployment = "main"
  ),

  tar_target(
    grouped_extreme_events_cape_19412024,
    group_extreme_events_by_proximity(
      extreme_events = all_extremes999$convective_available_potential_energy_daily_maximum,
      temporal_width = 2,
      distance_threshold = 50000,  # 50 km
      parallel = TRUE,
      num_cores = 40
    ),
    deployment = "main"
  ),

  tar_target(
    grouped_extreme_events_wind_speed_19412024,
    group_extreme_events_by_proximity(
      extreme_events = all_extremes999$wind_speed_daily_maximum,
      temporal_width = 2,
      distance_threshold = 50000,  # 50 km
      parallel = TRUE,
      num_cores = 40
    ),
    deployment = "main"
  ),
  tar_target(
    all_nuts_codes,
    extract_all_nuts(list(
      grouped_extreme_events_precip_19412024_q999,
      grouped_extreme_events_wind_speed_19412024,
      grouped_extreme_events_cape_19412024
    )),
    deployment = "main"
  ),
  tar_target(
    plot_timeline_all_nuts,
    {
      plot_event_timeline_generic(
        variables_list = list(
          grouped_extreme_events_precip_19412024_q999,
          grouped_extreme_events_wind_speed_19412024,
          grouped_extreme_events_cape_19412024
        ),
        var_names = c("Severe precipitation", "Severe wind", "CAPE"),
        nuts_target = all_nuts_codes,
        year = "1941-2024",
        output_file = NULL   # <- pas d'écriture
      )
    },
    pattern = map(all_nuts_codes),
    iteration = "list"
  ),

  tar_target(
    plot_timeline19412024_FRL0,
      plot_event_timeline_generic(
        variables_list = list(
          grouped_extreme_events_precip_19412024_q999,
          grouped_extreme_events_wind_speed_19412024,
          grouped_extreme_events_cape_19412024
        ),
        var_names = c("Severe precipitation", "Severe wind", "CAPE"),
        nuts_target = "FRL0",
        year = "1941-2024",
        output_file = "output/graphe/timeline_19412024_FRL0.png",
        make_plot = F
  )),
  tar_target(
    plot_timeline2024_FRL0_bees,
    plot_event_timeline_beeswarm_year_precip(
      variables_list = list(
        grouped_extreme_events_precip_19412024_q999,
        grouped_extreme_events_wind_speed_19412024,
        grouped_extreme_events_cape_19412024
      ),
      var_names = c("Severe precipitation", "Severe wind", "CAPE"),
      nuts_target = "FRL0",
      year = 2024,
      output_file = "output/graphe/timeline_2024_FRL0.png"
    )
  ),

    tar_target(
      var_names_pw,
      c("Severe precipitation", "Severe wind", "CAPE"),
      deployment = "main"
    ),

    tar_target(
      nuts_ref_codes_pw,
      c("FRL0", "FI1D", "DE21"),
      deployment = "main"
    ),

    tar_target(
      r_value_pw,
      3.5,
      deployment = "main"
    ),

    tar_target(
      out_dir_pw,
      "output/testidp_precip_wind_cape",
      deployment = "main"
    ),

    tar_target(
      alpha_pw,
      0.05,
      deployment = "main"
    ),

    # =============================================================================
    # INTER jobs: (ref in {FRL0,FI1D,DE21}) vs each nuts_id, for each variable
    # =============================================================================
    tar_target(
      inter_jobs_pw,
      {
        tidyr::expand_grid(
          nuts_ref    = nuts_ref_codes_pw,
          nuts_target = all_nuts_codes,
          var_name    = var_names_pw
        ) |>
          dplyr::filter(nuts_target != nuts_ref)
      }
    ),

    tar_target(
      inter_results_branches_pw,
      {
        dir.create(out_dir_pw, recursive = TRUE, showWarnings = FALSE)

        job <- inter_jobs_pw

        df_ref <- get_nuts_df(plot_timeline_all_nuts, job$nuts_ref)
        df_tar <- get_nuts_df(plot_timeline_all_nuts, job$nuts_target)

        if (is.null(df_ref) || nrow(df_ref) == 0) return(NULL)
        if (is.null(df_tar) || nrow(df_tar) == 0) return(NULL)

        analyze_inter_nuts_pair_yearly(
          df_ref = df_ref,
          df_target = df_tar,
          var_name = job$var_name,
          nuts_ref = job$nuts_ref,
          nuts_target = job$nuts_target,
          sim_method = "translation",
          alternative = "two.sided",
          r_value = r_value_pw,
          n_sim = 999,
          kde_bw = 0.08,
          homogeneous = FALSE,
          periodic = TRUE
        )
      },
      pattern = map(inter_jobs_pw),
      iteration = "list"
    ),

    tar_target(
      inter_results_pw,
      {
        lst <- inter_results_branches_pw[!vapply(inter_results_branches_pw, is.null, logical(1))]
        if (length(lst) == 0) return(NULL)
        dplyr::bind_rows(lst)
      }
    ),

    tar_target(
      inter_results_pw_csv,
      {
        dir.create(out_dir_pw, recursive = TRUE, showWarnings = FALSE)
        path <- file.path(out_dir_pw, "inter_multiRef_vs_all.csv")
        readr::write_csv(inter_results_pw, path)
        path
      },
      format = "file"
    ),

    # =============================================================================
    # INTRA directional jobs (greater only)
    # =============================================================================
    tar_target(
      intra_dir_pairs_pw,
      tibble::tibble(
        var1 = c("CAPE", "CAPE", "Severe precipitation", "Severe wind"),
        var2 = c("Severe precipitation", "Severe wind", "Severe wind", "Severe precipitation")
      )
    ),

    tar_target(
      intra_jobs_pw,
      {
        tidyr::expand_grid(
          nuts_code = all_nuts_codes,
          pair_id   = seq_len(nrow(intra_dir_pairs_pw))
        ) |>
          dplyr::left_join(
            dplyr::mutate(intra_dir_pairs_pw, pair_id = seq_len(dplyr::n())),
            by = "pair_id"
          ) |>
          dplyr::select(nuts_code, var1, var2)
      }
    ),

    tar_target(
      intra_results_branches_pw,
      {
        dir.create(out_dir_pw, recursive = TRUE, showWarnings = FALSE)

        job <- intra_jobs_pw

        df_n <- get_nuts_df(plot_timeline_all_nuts, job$nuts_code)
        if (is.null(df_n) || nrow(df_n) == 0) return(NULL)

        # Directionnel: var1 -> var2, alternative = greater
        analyze_intra_nuts_directional_yearly(
          df_nuts = df_n,
          var1_name = job$var1,
          var2_name = job$var2,
          nuts_code = job$nuts_code,
          sim_method = "translation",
          alternative = "greater",
          r_value = r_value_pw,
          n_sim = 999,
          kde_bw = 0.08,
          homogeneous = FALSE,
          periodic = TRUE
        )
      },
      pattern = map(intra_jobs_pw),
      iteration = "list"
    ),

    tar_target(
      intra_results_pw,
      {
        lst <- intra_results_branches_pw[!vapply(intra_results_branches_pw, is.null, logical(1))]
        if (length(lst) == 0) return(NULL)
        dplyr::bind_rows(lst)
      }
    ),

    tar_target(
      intra_results_pw_csv,
      {
        dir.create(out_dir_pw, recursive = TRUE, showWarnings = FALSE)
        path <- file.path(out_dir_pw, "intra_directional_greater.csv")
        readr::write_csv(intra_results_pw, path)
        path
      },
      format = "file"
    ),

  tar_target(
    inter_refs_pw,
    c("FRL0", "FI1D", "DE21"),
    deployment = "main"
  ),

  tar_target(
    pw_maps_files,
    {
      inter_df_all <- readr::read_csv(inter_results_pw_csv, show_col_types = FALSE)
      intra_df     <- readr::read_csv(intra_results_pw_csv, show_col_types = FALSE)

      plot_idp_maps_pw_multi(
        inter_df_all = inter_df_all,
        intra_df = intra_df,
        out_dir  = out_dir_pw,
        nuts_refs = inter_refs_pw,
        alpha    = alpha_pw
      )
    },
    format = "file"
  )

  ,

tar_target(
  var_precip_only,
  "Severe precipitation",
  deployment = "main"
),


tar_target(
  nuts_ref_code,
  "FRL0",
  deployment = "main"
),

tar_target(
  inter_jobs_precip_less,
  tidyr::expand_grid(
    nuts_target = all_nuts_codes
  ) |>
    dplyr::filter(nuts_target != nuts_ref_code)
),

tar_target(
  inter_results_precip_less_branches,
  {
    dir.create(out_dir_pw, recursive = TRUE, showWarnings = FALSE)
    
    job <- inter_jobs_precip_less
    df_ref <- get_nuts_df(plot_timeline_all_nuts, nuts_ref_code)
    df_tar <- get_nuts_df(plot_timeline_all_nuts, job$nuts_target)
    
    if (is.null(df_ref) || nrow(df_ref) == 0) return(NULL)
    if (is.null(df_tar) || nrow(df_tar) == 0) return(NULL)
    
    analyze_inter_nuts_pair_yearly(
      df_ref      = df_ref,
      df_target   = df_tar,
      var_name    = var_precip_only,
      nuts_ref    = nuts_ref_code,
      nuts_target = job$nuts_target,
      sim_method  = "translation",
      alternative = "less",   
      r_value     = r_value_pw,
      n_sim       = 999,
      kde_bw      = 0.08,
      homogeneous = FALSE,
      periodic    = TRUE
    )
  },
  pattern = map(inter_jobs_precip_less),
  iteration = "list"
),

tar_target(
  inter_results_precip_less,
  {
    lst <- inter_results_precip_less_branches[!vapply(inter_results_precip_less_branches, is.null, logical(1))]
    if (length(lst) == 0) return(NULL)
    dplyr::bind_rows(lst)
  }
),

tar_target(
  inter_results_precip_less_csv,
  {
    dir.create(out_dir_pw, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(out_dir_pw, "inter_FRL0_vs_all_precip_less.csv")
    readr::write_csv(inter_results_precip_less, path)
    path
  },
  format = "file"
)
)

#FWI ==============================================================================


FWI_COVARIATES <- c("Temperature", "Humidity", "Windspeed", "Radiation", "API", "VPD")
#
# # ======================================================================

list(

  # ---------------------------
  # ---------------------------

  tar_target(
    load_aquitaine,
    { load("~/postdoc_temporal_cee/data/df_ERA5-Aquitaine.RData"); df_ERA5 }
  ),
  tar_target(df_api_aquitaine, compute_API(load_aquitaine, m = 14, k = 0.9)),
  tar_target(df_long_aquitaine, reshape_to_long(df_api_aquitaine, name_region = "Aquitaine")),
  tar_target(extreme_events_aquitaine, compute_extremes(df_long_aquitaine, quantile_map)),

  tar_target(
    load_catalonia,
    { load("~/postdoc_temporal_cee/data/df_ERA5-Catalonia.RData"); df_ERA5 }
  ),
  tar_target(df_api_catalonia, compute_API(load_catalonia, m = 14, k = 0.9)),
  tar_target(df_long_catalonia, reshape_to_long(df_api_catalonia, name_region = "Catalonia")),
  tar_target(extreme_events_catalonia, compute_extremes(df_long_catalonia, quantile_map)),

  tar_target(
    load_portugal,
    { load("~/postdoc_temporal_cee/data/df_ERA5-Portugal.RData"); df_ERA5 }
  ),
  tar_target(df_api_portugal, compute_API(load_portugal, m = 14, k = 0.9)),
  tar_target(df_long_portugal, reshape_to_long(df_api_portugal, name_region = "Portugal")),
  tar_target(extreme_events_portugal, compute_extremes(df_long_portugal, quantile_map)),

    # --- Aquitaine
  tar_target(grouped_temperature_AQUITAINE,
             group_extreme_events_by_proximity(extreme_events_aquitaine$Temperature,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_humidity_AQUITAINE,
             group_extreme_events_by_proximity(extreme_events_aquitaine$Humidity,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_windspeed_AQUITAINE,
             group_extreme_events_by_proximity(extreme_events_aquitaine$Windspeed,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_rg_AQUITAINE,
             group_extreme_events_by_proximity(extreme_events_aquitaine$RG,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_api_AQUITAINE,
             group_extreme_events_by_proximity(extreme_events_aquitaine$API,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_vpd_AQUITAINE,
             group_extreme_events_by_proximity(extreme_events_aquitaine$VPD,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_fwi_AQUITAINE,
             group_extreme_events_by_proximity(extreme_events_aquitaine$FWI,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),

  # --- Catalonia
  tar_target(grouped_temperature_CATALONIA,
             group_extreme_events_by_proximity(extreme_events_catalonia$Temperature,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_humidity_CATALONIA,
             group_extreme_events_by_proximity(extreme_events_catalonia$Humidity,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_windspeed_CATALONIA,
             group_extreme_events_by_proximity(extreme_events_catalonia$Windspeed,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_rg_CATALONIA,
             group_extreme_events_by_proximity(extreme_events_catalonia$RG,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_api_CATALONIA,
             group_extreme_events_by_proximity(extreme_events_catalonia$API,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_vpd_CATALONIA,
             group_extreme_events_by_proximity(extreme_events_catalonia$VPD,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_fwi_CATALONIA,
             group_extreme_events_by_proximity(extreme_events_catalonia$FWI,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),

  # --- Portugal
  tar_target(grouped_temperature_PORTUGAL,
             group_extreme_events_by_proximity(extreme_events_portugal$Temperature,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_humidity_PORTUGAL,
             group_extreme_events_by_proximity(extreme_events_portugal$Humidity,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_windspeed_PORTUGAL,
             group_extreme_events_by_proximity(extreme_events_portugal$Windspeed,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_rg_PORTUGAL,
             group_extreme_events_by_proximity(extreme_events_portugal$RG,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_api_PORTUGAL,
             group_extreme_events_by_proximity(extreme_events_portugal$API,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_vpd_PORTUGAL,
             group_extreme_events_by_proximity(extreme_events_portugal$VPD,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),
  tar_target(grouped_fwi_PORTUGAL,
             group_extreme_events_by_proximity(extreme_events_portugal$FWI,
                                               temporal_width = 7, distance_threshold = 50000, parallel = TRUE, num_cores = 40)),

  # common var list and names
  tar_target(
    var_names_common,
    c("Temperature","Humidity","Windspeed","Radiation","API","VPD","FWI")
  ),

  # Aquitaine timelines
  tar_target(
    plot_timeline_AQUITAINE_beeswarm,
    plot_event_timeline_beeswarm(
      variables_list = list(
        grouped_temperature_AQUITAINE,
        grouped_humidity_AQUITAINE,
        grouped_windspeed_AQUITAINE,
        grouped_rg_AQUITAINE,
        grouped_api_AQUITAINE,
        grouped_vpd_AQUITAINE,
        grouped_fwi_AQUITAINE
      ),
      var_names = var_names_common,
      nuts_target = "Aquitaine",
      year = "2000-2022",
      output_file = "output/graphe/timeline_2000_2022_AQUITAINE.png"
    )
  ),
  tar_target(
    plot_timeline_AQUITAINE_beeswarm_2022,
    plot_event_timeline_beeswarm_year(
      variables_list = list(
        grouped_temperature_AQUITAINE,
        grouped_humidity_AQUITAINE,
        grouped_windspeed_AQUITAINE,
        grouped_rg_AQUITAINE,
        grouped_api_AQUITAINE,
        grouped_vpd_AQUITAINE,
        grouped_fwi_AQUITAINE
      ),
      var_names = var_names_common,
      nuts_target = "Aquitaine",
      year = 2022,
      output_file = "output/graphe/timeline_2022_AQUITAINE.png"
    )
  ),
  tar_target(
    plot_cumu_2022_aqui,
    plot_marked_cumulative_trajectories_yearly(
      df = plot_timeline_AQUITAINE_beeswarm,
      nuts_target = "Aquitaine",
      variable_target = "FWI",
      mark_col = "Magnitude",
      highlight_year = 2022,
      output_file = "output/graphe/fwi_cumtraj_Aquitaine_2022.png"
    )
  ),

  # Catalonia timelines
  tar_target(
    plot_timeline_CATALONIA_beeswarm,
    plot_event_timeline_beeswarm(
      variables_list = list(
        grouped_temperature_CATALONIA,
        grouped_humidity_CATALONIA,
        grouped_windspeed_CATALONIA,
        grouped_rg_CATALONIA,
        grouped_api_CATALONIA,
        grouped_vpd_CATALONIA,
        grouped_fwi_CATALONIA
      ),
      var_names = var_names_common,
      nuts_target = "Catalonia",
      year = "2000-2022",
      output_file = "output/graphe/timeline_2000_2022_CATALONIA.png"
    )),
  tar_target(
    plot_timeline_CATALONIA_beeswarm_2022,
    plot_event_timeline_beeswarm_year(
      variables_list = list(
        grouped_temperature_CATALONIA,
        grouped_humidity_CATALONIA,
        grouped_windspeed_CATALONIA,
        grouped_rg_CATALONIA,
        grouped_api_CATALONIA,
        grouped_vpd_CATALONIA,
        grouped_fwi_CATALONIA
      ),
      var_names = var_names_common,
      nuts_target = "Catalonia",
      year = 2022,
      output_file = "output/graphe/timeline_2022_CATALONIA.png"
    )
  ),
  tar_target(
    plot_cumu_2022_cat,
    plot_marked_cumulative_trajectories_yearly(
      df = plot_timeline_CATALONIA_beeswarm,
      nuts_target = "Catalonia",
      variable_target = "FWI",
      mark_col = "Magnitude",
      highlight_year = 2022,
      output_file = "output/graphe/fwi_cumtraj_Catalonia_2022.png"
    )
  ),

  # Portugal timelines
  tar_target(
    plot_timeline_PORTUGAL_beeswarm,
    plot_event_timeline_beeswarm(
      variables_list = list(
        grouped_temperature_PORTUGAL,
        grouped_humidity_PORTUGAL,
        grouped_windspeed_PORTUGAL,
        grouped_rg_PORTUGAL,
        grouped_api_PORTUGAL,
        grouped_vpd_PORTUGAL,
        grouped_fwi_PORTUGAL
      ),
      var_names = var_names_common,
      nuts_target = "Portugal",
      year = "2000-2022",
      output_file = "output/graphe/timeline_2000_2022_PORTUGAL.png"
    )
  ),
  tar_target(
    plot_timeline_PORTUGAL_beeswarm_2022,
    plot_event_timeline_beeswarm_year(
      variables_list = list(
        grouped_temperature_PORTUGAL,
        grouped_humidity_PORTUGAL,
        grouped_windspeed_PORTUGAL,
        grouped_rg_PORTUGAL,
        grouped_api_PORTUGAL,
        grouped_vpd_PORTUGAL,
        grouped_fwi_PORTUGAL
      ),
      var_names = var_names_common,
      nuts_target = "Portugal",
      year = 2022,
      output_file = "output/graphe/timeline_2022_PORTUGAL.png"
    )
  ),
  tar_target(
    plot_cumu_2022_pt,
    plot_marked_cumulative_trajectories_yearly(
      df = plot_timeline_PORTUGAL_beeswarm,
      nuts_target = "Portugal",
      variable_target = "FWI",
      mark_col = "Magnitude",
      highlight_year = 2022,
      output_file = "output/graphe/fwi_cumtraj_Portugal_2022.png"
    )
  ),

  
  tar_target(
    fwi_covariates_vec,
    c("Temperature","Humidity","Windspeed","Radiation","API","VPD")
  ),

  tar_target(
    fwi_jobs,
    {
      covs <- fwi_covariates_vec

      jobs_sym <- expand.grid(
        covariate = covs,
        directional = FALSE,
        alternative = "two.sided",
        stringsAsFactors = FALSE
      )

      jobs_dir <- expand.grid(
        covariate = covs,
        directional = TRUE,
        alternative = c("greater","less"),
        direction = c("FWI_to_cov","cov_to_FWI"),
        stringsAsFactors = FALSE
      ) |>
        dplyr::mutate(covariate = paste0(covariate, "__", direction)) |>
        dplyr::select(-direction)

      dplyr::bind_rows(jobs_sym, jobs_dir)
    }
  ),

  # --- Aquitaine branches
  tar_target(
    fwi_results_AQUITAINE_branches,
    {
      job <- fwi_jobs
      run_one_FWI_job(
        df = plot_timeline_AQUITAINE_beeswarm,
        region = "AQUITAINE",
        covariate = job$covariate,
        max_r = 10.5,
        n_sim = 999,
        kde_bw = 0.08,
        sim_method = "translation",
        homogeneous = FALSE,
        periodic = TRUE,
        center = TRUE,
        directional = job$directional,
        alternative = job$alternative,
        out_dir = "output/graphe_k/FWI/AQUITAINE"
      )
    },
    pattern = map(fwi_jobs)
  ),
  tar_target(
    fwi_results_AQUITAINE_all,
    dplyr::bind_rows(fwi_results_AQUITAINE_branches)
  ),
  tar_target(
    fwi_results_AQUITAINE_csv,
    {
      dir.create("output/graphe_k/FWI/AQUITAINE", recursive = TRUE, showWarnings = FALSE)
      path <- "output/graphe_k/FWI/AQUITAINE/summary_FWI_AQUITAINE.csv"
      readr::write_csv(fwi_results_AQUITAINE_all, path)
      path
    },
    format = "file"
  ),

  # --- Catalonia branches
  tar_target(
    fwi_results_CATALONIA_branches,
    {
      job <- fwi_jobs
      run_one_FWI_job(
        df = plot_timeline_CATALONIA_beeswarm,
        region = "CATALONIA",
        covariate = job$covariate,
        max_r = 10.5,
        n_sim = 999,
        kde_bw = 0.08,
        sim_method = "translation",
        homogeneous = FALSE,
        periodic = TRUE,
        center = TRUE,
        directional = job$directional,
        alternative = job$alternative,
        out_dir = "output/graphe_k/FWI/CATALONIA"
      )
    },
    pattern = map(fwi_jobs)
  ),
  tar_target(
    fwi_results_CATALONIA_all,
    dplyr::bind_rows(fwi_results_CATALONIA_branches)
  ),
  tar_target(
    fwi_results_CATALONIA_csv,
    {
      dir.create("output/graphe_k/FWI/CATALONIA", recursive = TRUE, showWarnings = FALSE)
      path <- "output/graphe_k/FWI/CATALONIA/summary_FWI_CATALONIA.csv"
      readr::write_csv(fwi_results_CATALONIA_all, path)
      path
    },
    format = "file"
  ),

  # --- Portugal branches
  tar_target(
    fwi_results_PORTUGAL_branches,
    {
      job <- fwi_jobs
      run_one_FWI_job(
        df = plot_timeline_PORTUGAL_beeswarm,
        region = "PORTUGAL",
        covariate = job$covariate,
        max_r = 10.5,
        n_sim = 999,
        kde_bw = 0.08,
        sim_method = "translation",
        homogeneous = FALSE,
        periodic = TRUE,
        center = TRUE,
        directional = job$directional,
        alternative = job$alternative,
        out_dir = "output/graphe_k/FWI/PORTUGAL"
      )
    },
    pattern = map(fwi_jobs)
  ),
  tar_target(
    fwi_results_PORTUGAL_all,
    dplyr::bind_rows(fwi_results_PORTUGAL_branches)
  ),
  tar_target(
    fwi_results_PORTUGAL_csv,
    {
      dir.create("output/graphe_k/FWI/PORTUGAL", recursive = TRUE, showWarnings = FALSE)
      path <- "output/graphe_k/FWI/PORTUGAL/summary_FWI_PORTUGAL.csv"
      readr::write_csv(fwi_results_PORTUGAL_all, path)
      path
    },
    format = "file"
  ),

  tar_target(
    fwi_region_pairs,
    tibble::tibble(
      r1 = c("AQUITAINE","AQUITAINE","PORTUGAL"),
      r2 = c("CATALONIA","PORTUGAL","CATALONIA")
    )
  ),

  tar_target(
    fwi_pairwise_branches,
    {
      suppressPackageStartupMessages({
        library(dplyr)
        library(stringr)
        library(readr)
      })

      pair <- fwi_region_pairs

    
      df1 <- switch(pair$r1,
                    "AQUITAINE" = plot_timeline_AQUITAINE_beeswarm,
                    "CATALONIA" = plot_timeline_CATALONIA_beeswarm,
                    "PORTUGAL"  = plot_timeline_PORTUGAL_beeswarm
      )

      df2 <- switch(pair$r2,
                    "AQUITAINE" = plot_timeline_AQUITAINE_beeswarm,
                    "CATALONIA" = plot_timeline_CATALONIA_beeswarm,
                    "PORTUGAL"  = plot_timeline_PORTUGAL_beeswarm
      )

      out_dir <- "output/graphe_k/FWI/BETWEEN"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      safe <- function(x) {
        x |>
          as.character() |>
          str_replace_all("[^A-Za-z0-9_\\-]+", "_") |>
          str_replace_all("_+", "_") |>
          str_replace_all("^_|_$", "")
      }

      tag <- paste0("FWI_between__", pair$r1, "__vs__", pair$r2)
      plot_path <- file.path(out_dir, paste0(safe(tag), ".png"))

      # --- run K test (sym only) ---
      res <- calculate_bivariate_K_GET_yearly(
        df1 = df1,
        df2 = df2,
        var1_name = "FWI",
        var2_name = "FWI",
        max_r = 10.5,
        time_col = "Jour",
        T_total = 365,
        n_sim = 999,
        kde_bw = 0.08,
        periodic = TRUE,
        homogeneous = FALSE,
        directional = FALSE,
        alternative = "two.sided",
        sim_method = "translation",
        center = TRUE,
        save_plot = TRUE,
        save_path = plot_path,
        plot_title = paste0("FWI ↔ FWI | ", pair$r1, " vs ", pair$r2)
      )

      pv <- res$p_value
      pv_num <- if (length(pv) == 2) pv[2] else pv

      tibble::tibble(
        test = "FWI_between_regions",
        region1 = pair$r1,
        region2 = pair$r2,
        var1 = "FWI",
        var2 = "FWI",
        directional = FALSE,
        alternative = "two.sided",
        sim_method = "translation",
        homogeneous = FALSE,
        periodic = TRUE,
        center = TRUE,
        max_r = 10.5,
        n_sim = 999,
        kde_bw = 0.08,
        p_value = pv_num,
        plot_path = plot_path
      )
    },
    pattern = map(fwi_region_pairs)
  ),

  tar_target(
    fwi_pairwise_all,
    dplyr::bind_rows(fwi_pairwise_branches) |>
      dplyr::arrange(region1, region2)
  ),

  tar_target(
    fwi_pairwise_csv,
    {
      dir.create("output/graphe_k/FWI/BETWEEN", recursive = TRUE, showWarnings = FALSE)
      path <- "output/graphe_k/FWI/BETWEEN/summary_FWI_pairwise_regions.csv"
      readr::write_csv(fwi_pairwise_all, path)
      path
    },
    format = "file"
  ),

tar_target(
  fwi_jobs_dir_cov_to_fwi,
  {
    covs <- c("Temperature","Humidity","Windspeed","Radiation","API","VPD")

    tidyr::expand_grid(
      covariate   = covs,
      alternative = c("two.sided","greater","less"),
      sim_method  = c("translation","nhpp")
    )
  }
),
tar_target(
  fwi_dir_cov_to_fwi_AQUITAINE_branches,
  {
    job <- fwi_jobs_dir_cov_to_fwi
    run_one_FWI_job_dir_cov_to_FWI(
      df = plot_timeline_AQUITAINE_beeswarm,
      region = "AQUITAINE",
      covariate = job$covariate,
      max_r = 10.5,
      n_sim = 999,
      kde_bw = 0.08,
      sim_method = job$sim_method,
      alternative = job$alternative,
      out_dir = "output/graphe_k/FWI/AQUITAINE_DIR_COV_TO_FWI"
    )
  },
  pattern = map(fwi_jobs_dir_cov_to_fwi),
  iteration = "list"
),

tar_target(
  fwi_dir_cov_to_fwi_AQUITAINE_all,
  {
    lst <- fwi_dir_cov_to_fwi_AQUITAINE_branches
    lst <- lst[!vapply(lst, is.null, logical(1))]
    if (length(lst) == 0) return(NULL)
    dplyr::bind_rows(lst) |>
      dplyr::arrange(covariate, sim_method, alternative)
  }
),

tar_target(
  fwi_dir_cov_to_fwi_AQUITAINE_csv,
  {
    dir.create("output/graphe_k/FWI/AQUITAINE_DIR_COV_TO_FWI", recursive = TRUE, showWarnings = FALSE)
    path <- "output/graphe_k/FWI/AQUITAINE_DIR_COV_TO_FWI/summary_dir_cov_to_fwi_AQUITAINE.csv"
    readr::write_csv(fwi_dir_cov_to_fwi_AQUITAINE_all, path)
    path
  },
  format = "file"
),
tar_target(
  fwi_dir_cov_to_fwi_PORTUGAL_branches,
  {
    job <- fwi_jobs_dir_cov_to_fwi
    run_one_FWI_job_dir_cov_to_FWI(
      df = plot_timeline_PORTUGAL_beeswarm,
      region = "PORTUGAL",
      covariate = job$covariate,
      max_r = 10.5,
      n_sim = 999,
      kde_bw = 0.08,
      sim_method = job$sim_method,
      alternative = job$alternative,
      out_dir = "output/graphe_k/FWI/PORTUGAL_DIR_COV_TO_FWI"
    )
  },
  pattern = map(fwi_jobs_dir_cov_to_fwi),
  iteration = "list"
),

tar_target(
  fwi_dir_cov_to_fwi_PORTUGAL_all,
  {
    lst <- fwi_dir_cov_to_fwi_PORTUGAL_branches
    lst <- lst[!vapply(lst, is.null, logical(1))]
    if (length(lst) == 0) return(NULL)
    dplyr::bind_rows(lst) |>
      dplyr::arrange(covariate, sim_method, alternative)
  }
),

tar_target(
  fwi_dir_cov_to_fwi_PORTUGAL_csv,
  {
    dir.create("output/graphe_k/FWI/PORTUGAL_DIR_COV_TO_FWI", recursive = TRUE, showWarnings = FALSE)
    path <- "output/graphe_k/FWI/PORTUGAL_DIR_COV_TO_FWI/summary_dir_cov_to_fwi_PORTUGAL.csv"
    readr::write_csv(fwi_dir_cov_to_fwi_PORTUGAL_all, path)
    path
  },
  format = "file"
),
tar_target(
  fwi_dir_cov_to_fwi_CATALONIA_branches,
  {
    job <- fwi_jobs_dir_cov_to_fwi

    run_one_FWI_job_dir_cov_to_FWI(
      df = plot_timeline_CATALONIA_beeswarm,
      region = "CATALONIA",
      covariate = job$covariate,
      max_r = 10.5,
      n_sim = 999,
      kde_bw = 0.08,
      sim_method = job$sim_method,
      alternative = job$alternative,
      out_dir = "output/graphe_k/FWI/CATALONIA_DIR_COV_TO_FWI"
    )
  },
  pattern = map(fwi_jobs_dir_cov_to_fwi),  # <-- OK, car fwi_jobs_dir_cov_to_fwi est un tar_target
  iteration = "list"
),
tar_target(
  fwi_dir_cov_to_fwi_CATALONIA_all,
  {
    lst <- fwi_dir_cov_to_fwi_CATALONIA_branches
    lst <- lst[!vapply(lst, is.null, logical(1))]
    if (length(lst) == 0) return(NULL)
    dplyr::bind_rows(lst) |>
      dplyr::arrange(covariate, sim_method, alternative)
  }
),

tar_target(
  fwi_dir_cov_to_fwi_CATALONIA_csv,
  {
    dir.create("output/graphe_k/FWI/CATALONIA_DIR_COV_TO_FWI", recursive = TRUE, showWarnings = FALSE)
    path <- "output/graphe_k/FWI/CATALONIA_DIR_COV_TO_FWI/summary_dir_cov_to_fwi_CATALONIA.csv"
    readr::write_csv(fwi_dir_cov_to_fwi_CATALONIA_all, path)
    path
  },
  format = "file"
),
tar_target(
  fwi_dir_cov_to_fwi_all_regions,
  dplyr::bind_rows(
    fwi_dir_cov_to_fwi_AQUITAINE_all,
    fwi_dir_cov_to_fwi_CATALONIA_all,
    fwi_dir_cov_to_fwi_PORTUGAL_all
  )
),

tar_target(
  fwi_dir_cov_to_fwi_all_regions_csv,
  {
    dir.create("output/graphe_k/FWI/DIR_COV_TO_FWI_ALL", recursive = TRUE, showWarnings = FALSE)
    path <- "output/graphe_k/FWI/DIR_COV_TO_FWI_ALL/summary_dir_cov_to_fwi_all_regions.csv"
    readr::write_csv(fwi_dir_cov_to_fwi_all_regions, path)
    path
  },
  format = "file"
)
)
