#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  check_required_packages <- function(pkgs) {
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop(
        "Missing required packages: ",
        paste(missing_pkgs, collapse = ", "),
        ". Please install them before running this script."
      )
    }
  }

  required_pkgs <- c("terra", "dplyr", "purrr", "tibble", "readr", "tidyr", "yaml")
  check_required_packages(required_pkgs)

  library(terra)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(yaml)
})

normalize_positive <- function(x, eps = 1e-12) {
  x[is.na(x)] <- 0
  x[x < 0] <- 0
  s <- sum(x)
  if (s <= eps) rep(1 / length(x), length(x)) else x / s
}

calc_critic_weights <- function(df_vars) {
  var_names <- colnames(df_vars)
  p <- ncol(df_vars)

  if (p == 1) {
    return(tibble(
      Variable = var_names,
      CRITIC_score = 1,
      Weight_CRITIC = 1
    ))
  }

  disp <- apply(df_vars, 2, IQR, na.rm = TRUE)
  corr_abs <- abs(stats::cor(df_vars, method = "spearman", use = "pairwise.complete.obs"))
  diag(corr_abs) <- 1

  c_score <- disp * rowSums(1 - corr_abs, na.rm = TRUE)
  c_score[is.na(c_score)] <- 0
  w <- normalize_positive(c_score)

  tibble(
    Variable = var_names,
    CRITIC_score = as.numeric(c_score),
    Weight_CRITIC = as.numeric(w)
  )
}

calc_entropy_weights <- function(df_vars, eps = 1e-6) {
  var_names <- colnames(df_vars)
  p <- ncol(df_vars)

  if (p == 1) {
    return(tibble(
      Variable = var_names,
      Entropy = 0,
      Divergence = 1,
      Weight_Entropy = 1
    ))
  }

  x <- as.matrix(df_vars) + eps
  p_mat <- sweep(x, 2, colSums(x, na.rm = TRUE), FUN = "/")
  n <- nrow(p_mat)

  e <- -colSums(p_mat * log(p_mat), na.rm = TRUE) / log(n)
  d <- 1 - e
  d[is.na(d)] <- 0
  w <- normalize_positive(d)

  tibble(
    Variable = var_names,
    Entropy = as.numeric(e),
    Divergence = as.numeric(d),
    Weight_Entropy = as.numeric(w)
  )
}

mer_fuse_weights <- function(candidate_list,
                             candidate_weights = NULL,
                             lambda = 0.01,
                             eps = 0.001,
                             maxit = 2000) {
  p <- length(candidate_list[[1]])
  if (p * eps >= 1) stop("MER constraint violated: p * eps must be < 1.")

  if (p == 1) {
    return(list(
      weight = 1,
      objective = 0,
      entropy = 0,
      convergence = 0,
      method = "single-variable"
    ))
  }

  cand_mat <- do.call(cbind, lapply(candidate_list, normalize_positive))
  rownames(cand_mat) <- names(candidate_list[[1]])

  if (is.null(candidate_weights)) {
    candidate_weights <- rep(1 / ncol(cand_mat), ncol(cand_mat))
  } else {
    candidate_weights <- normalize_positive(candidate_weights)
  }

  theta_to_weight <- function(theta) {
    q <- exp(theta - max(theta))
    q <- q / sum(q)
    eps + (1 - p * eps) * q
  }

  objective_fun <- function(theta) {
    w <- theta_to_weight(theta)
    residual_term <- 0
    for (m in seq_len(ncol(cand_mat))) {
      residual_term <- residual_term + candidate_weights[m] * sum((w - cand_mat[, m])^2)
    }
    entropy_term <- -sum(w * log(pmax(w, 1e-15)))
    residual_term - lambda * entropy_term
  }

  avg_cand <- normalize_positive(rowMeans(cand_mat))
  theta_init <- log(pmax(avg_cand, 1e-8))

  opt <- tryCatch(
    stats::optim(
      par = theta_init,
      fn = objective_fun,
      method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-10)
    ),
    error = function(e) NULL
  )

  if (is.null(opt)) {
    q <- normalize_positive(avg_cand)
    w <- eps + (1 - p * eps) * q
    w <- normalize_positive(w)
    return(list(
      weight = as.numeric(w),
      objective = NA_real_,
      entropy = -sum(w * log(pmax(w, 1e-15))),
      convergence = 999,
      method = "fallback-average"
    ))
  }

  w_final <- theta_to_weight(opt$par)
  list(
    weight = as.numeric(w_final),
    objective = opt$value,
    entropy = -sum(w_final * log(pmax(w_final, 1e-15))),
    convergence = opt$convergence,
    method = "BFGS"
  )
}

calc_fixed_from_annual_series <- function(weight_wide_df,
                                          vars,
                                          lambda = 0.01,
                                          eps = 0.001,
                                          use_stability = TRUE) {
  annual_mat <- weight_wide_df[, vars, drop = FALSE]

  critic_long <- calc_critic_weights(annual_mat)
  entropy_long <- calc_entropy_weights(annual_mat)

  mean_v <- colMeans(annual_mat, na.rm = TRUE)
  sd_v <- apply(annual_mat, 2, sd, na.rm = TRUE)
  cv_v <- sd_v / pmax(mean_v, 1e-8)
  stability_raw <- exp(-cv_v)
  stability_w <- normalize_positive(stability_raw)

  vC <- critic_long$Weight_CRITIC[match(vars, critic_long$Variable)]
  names(vC) <- vars

  vE <- entropy_long$Weight_Entropy[match(vars, entropy_long$Variable)]
  names(vE) <- vars

  candidates <- list(LongTermCRITIC = vC, LongTermEntropy = vE)

  if (use_stability) {
    vS <- stability_w
    names(vS) <- vars
    candidates$Stability <- vS
    cand_weights <- rep(1 / 3, 3)
  } else {
    vS <- rep(NA_real_, length(vars))
    names(vS) <- vars
    cand_weights <- rep(1 / 2, 2)
  }

  mer_fixed <- mer_fuse_weights(
    candidate_list = candidates,
    candidate_weights = cand_weights,
    lambda = lambda,
    eps = eps
  )

  tibble(
    Variable = vars,
    LongTermWeight_CRITIC = vC,
    LongTermWeight_Entropy = vE,
    StabilityWeight = vS,
    FixedWeight_MER = mer_fixed$weight,
    FixedMER_objective = mer_fixed$objective,
    FixedMER_entropy = mer_fixed$entropy,
    FixedMER_convergence = mer_fixed$convergence,
    FixedMER_method = mer_fixed$method,
    MeanSeriesWeight = mean_v[vars],
    SDSeriesWeight = sd_v[vars],
    CVSeriesWeight = cv_v[vars]
  )
}

build_file_path <- function(var_name, year, variables_cfg) {
  pattern <- variables_cfg[[var_name]]$pattern
  file.path(variables_cfg[[var_name]]$dir, gsub("\\{year\\}", as.character(year), pattern))
}

align_to_reference <- function(r, ref, method = "bilinear") {
  if (!terra::same.crs(r, ref)) {
    r <- terra::project(r, ref, method = method)
  }
  if (!terra::compareGeom(ref, r, stopOnError = FALSE)) {
    r <- terra::resample(r, ref, method = method)
  }
  if (!terra::compareGeom(ref, r, stopOnError = FALSE)) {
    stop("Raster geometry mismatch remains after alignment.")
  }
  r
}

read_year_stack <- function(year, variables_cfg, variable_names) {
  files <- purrr::map_chr(variable_names, ~ build_file_path(.x, year, variables_cfg))
  rasters <- purrr::map(files, terra::rast)

  ref <- rasters[[1]]
  aligned <- vector("list", length(rasters))
  aligned[[1]] <- ref

  if (length(rasters) > 1) {
    for (k in 2:length(rasters)) {
      aligned[[k]] <- align_to_reference(rasters[[k]], ref, method = "bilinear")
    }
  }

  s <- do.call(c, aligned)
  names(s) <- variable_names
  s
}

assign_block_id <- function(df_xy, ref_raster, block_size_lonlat, block_size_projected) {
  if (terra::is.lonlat(ref_raster)) {
    bx <- floor((df_xy$x + 180) / block_size_lonlat)
    by <- floor((df_xy$y + 90) / block_size_lonlat)
  } else {
    bx <- floor(df_xy$x / block_size_projected)
    by <- floor(df_xy$y / block_size_projected)
  }
  paste0("B", bx, "_", by)
}

sample_year_values <- function(year,
                               variables_cfg,
                               variable_names,
                               sample_n,
                               seed,
                               block_size_lonlat,
                               block_size_projected) {
  set.seed(seed + year)
  s <- read_year_stack(year, variables_cfg, variable_names)
  ref <- s[[1]]

  vals <- terra::spatSample(
    s,
    size = sample_n,
    method = "random",
    na.rm = TRUE,
    values = TRUE,
    xy = TRUE,
    as.data.frame = TRUE
  )

  vals <- vals %>%
    dplyr::rename(x = x, y = y) %>%
    dplyr::select(x, y, all_of(variable_names)) %>%
    dplyr::filter(stats::complete.cases(.))

  vals$block_id <- assign_block_id(
    vals[, c("x", "y")],
    ref,
    block_size_lonlat = block_size_lonlat,
    block_size_projected = block_size_projected
  )
  vals
}

calc_weights_one_year <- function(df_vars,
                                  year,
                                  groups_fixed,
                                  mer_lambda_within = 0.01,
                                  mer_lambda_between = 0.01,
                                  mer_eps = 0.001,
                                  mer_candidate_weights = c(0.5, 0.5)) {
  within_rows <- list()
  group_score_df <- tibble(.rows = nrow(df_vars))

  for (g in names(groups_fixed)) {
    vars_g <- groups_fixed[[g]]
    dat_g <- df_vars[, vars_g, drop = FALSE]

    critic_g <- calc_critic_weights(dat_g)
    entropy_g <- calc_entropy_weights(dat_g)

    wC <- critic_g$Weight_CRITIC
    names(wC) <- critic_g$Variable

    wE <- entropy_g$Weight_Entropy
    names(wE) <- entropy_g$Variable

    mer_g <- mer_fuse_weights(
      candidate_list = list(CRITIC = wC, Entropy = wE),
      candidate_weights = mer_candidate_weights,
      lambda = mer_lambda_within,
      eps = mer_eps
    )

    wM <- mer_g$weight
    names(wM) <- vars_g

    res_g <- tibble(
      Year = year,
      Group = g,
      Variable = vars_g,
      CRITIC_score = critic_g$CRITIC_score[match(vars_g, critic_g$Variable)],
      Weight_CRITIC = critic_g$Weight_CRITIC[match(vars_g, critic_g$Variable)],
      Entropy = entropy_g$Entropy[match(vars_g, entropy_g$Variable)],
      Divergence = entropy_g$Divergence[match(vars_g, entropy_g$Variable)],
      Weight_Entropy = entropy_g$Weight_Entropy[match(vars_g, entropy_g$Variable)],
      Weight_MER = as.numeric(wM)
    )

    within_rows[[g]] <- res_g
    group_score_df[[g]] <- as.numeric(as.matrix(dat_g) %*% wM)
  }

  group_critic <- calc_critic_weights(group_score_df)
  group_entropy <- calc_entropy_weights(group_score_df)

  WC <- group_critic$Weight_CRITIC
  names(WC) <- group_critic$Variable
  WE <- group_entropy$Weight_Entropy
  names(WE) <- group_entropy$Variable

  mer_between <- mer_fuse_weights(
    candidate_list = list(CRITIC = WC, Entropy = WE),
    candidate_weights = mer_candidate_weights,
    lambda = mer_lambda_between,
    eps = mer_eps
  )

  WM <- mer_between$weight
  names(WM) <- names(groups_fixed)

  between_tbl <- tibble(
    Year = year,
    Group = names(groups_fixed),
    CRITIC_score = group_critic$CRITIC_score[match(names(groups_fixed), group_critic$Variable)],
    Weight_CRITIC = group_critic$Weight_CRITIC[match(names(groups_fixed), group_critic$Variable)],
    Entropy = group_entropy$Entropy[match(names(groups_fixed), group_entropy$Variable)],
    Divergence = group_entropy$Divergence[match(names(groups_fixed), group_entropy$Variable)],
    Weight_Entropy = group_entropy$Weight_Entropy[match(names(groups_fixed), group_entropy$Variable)],
    Weight_MER = as.numeric(WM)
  )

  within_tbl <- bind_rows(within_rows)

  annual_indicator_tbl <- within_tbl %>%
    left_join(
      between_tbl %>% select(Group, BetweenGroupWeight_MER = Weight_MER),
      by = "Group"
    ) %>%
    rename(WithinGroupWeight_MER = Weight_MER) %>%
    mutate(
      AnnualIndicatorWeight = WithinGroupWeight_MER * BetweenGroupWeight_MER
    ) %>%
    select(
      Year, Group, Variable,
      WithinGroupWeight_MER, BetweenGroupWeight_MER, AnnualIndicatorWeight
    )

  list(
    within = within_tbl,
    between = between_tbl,
    annual_indicator = annual_indicator_tbl
  )
}

calc_final_fixed_weights_hier <- function(within_tbl,
                                          between_tbl,
                                          groups_fixed,
                                          fixed_lambda_between = 0.01,
                                          fixed_lambda_group = NULL,
                                          mer_eps = 0.001) {
  between_wide <- between_tbl %>%
    select(Year, Group, Weight_MER) %>%
    pivot_wider(names_from = Group, values_from = Weight_MER) %>%
    arrange(Year)

  fixed_between <- calc_fixed_from_annual_series(
    weight_wide_df = between_wide,
    vars = names(groups_fixed),
    lambda = fixed_lambda_between,
    eps = mer_eps,
    use_stability = TRUE
  ) %>%
    rename(Group = Variable, FixedBetweenWeight_HierMER = FixedWeight_MER)

  fixed_within_list <- list()

  for (g in names(groups_fixed)) {
    vars_g <- groups_fixed[[g]]
    lambda_g <- fixed_lambda_group[[g]]

    within_wide_g <- within_tbl %>%
      filter(Group == g) %>%
      select(Year, Variable, Weight_MER) %>%
      pivot_wider(names_from = Variable, values_from = Weight_MER) %>%
      arrange(Year)

    fixed_within_g <- calc_fixed_from_annual_series(
      weight_wide_df = within_wide_g,
      vars = vars_g,
      lambda = lambda_g,
      eps = mer_eps,
      use_stability = TRUE
    ) %>%
      mutate(Group = g, FixedWithinLambda = lambda_g) %>%
      rename(FixedWithinWeight_HierMER = FixedWeight_MER)

    fixed_within_list[[g]] <- fixed_within_g
  }

  fixed_within <- bind_rows(fixed_within_list)

  final_fixed <- fixed_within %>%
    left_join(
      fixed_between %>% select(Group, FixedBetweenWeight_HierMER),
      by = "Group"
    ) %>%
    mutate(
      FinalFixedWeight_HierMER = FixedWithinWeight_HierMER * FixedBetweenWeight_HierMER
    ) %>%
    select(
      Group, Variable,
      FixedBetweenWeight_HierMER,
      FixedWithinWeight_HierMER,
      FinalFixedWeight_HierMER,
      LongTermWeight_CRITIC,
      LongTermWeight_Entropy,
      StabilityWeight,
      MeanSeriesWeight, SDSeriesWeight, CVSeriesWeight,
      FixedMER_objective, FixedMER_entropy, FixedMER_convergence, FixedMER_method,
      FixedWithinLambda
    ) %>%
    arrange(desc(FinalFixedWeight_HierMER)) %>%
    mutate(Rank = dplyr::row_number()) %>%
    select(Rank, everything())

  list(
    fixed_between = fixed_between,
    fixed_within = fixed_within,
    final_fixed = final_fixed
  )
}

bootstrap_resample_one_year <- function(df_year, max_points = 60000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  unique_blocks <- unique(df_year$block_id)
  sampled_blocks <- sample(unique_blocks, size = length(unique_blocks), replace = TRUE)

  resampled <- purrr::map_dfr(seq_along(sampled_blocks), function(i) {
    b <- sampled_blocks[i]
    df_year %>%
      dplyr::filter(block_id == b) %>%
      dplyr::mutate(.boot_block_rep = i)
  })

  if (nrow(resampled) > max_points) {
    idx <- sample(seq_len(nrow(resampled)), size = max_points, replace = FALSE)
    resampled <- resampled[idx, , drop = FALSE]
  }
  resampled
}

run_hdci_weighting <- function(config_path = "config.yml") {
  cfg <- yaml::read_yaml(config_path)

  years <- seq(cfg$years$start, cfg$years$end)
  variable_names <- names(cfg$variables)
  groups_fixed <- cfg$groups

  out_root <- cfg$output_dir
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  path_df <- tidyr::crossing(
    Year = years,
    Variable = variable_names
  ) %>%
    mutate(
      File = purrr::map2_chr(Variable, Year, ~ build_file_path(.x, .y, cfg$variables)),
      Exists = file.exists(File)
    )

  readr::write_csv(path_df, file.path(out_root, "all_file_paths_and_existence.csv"))

  complete_years <- path_df %>%
    group_by(Year) %>%
    summarise(AllFilesExist = all(Exists), .groups = "drop") %>%
    filter(AllFilesExist) %>%
    pull(Year)

  if (length(complete_years) < 3) {
    stop("Fewer than 3 complete years are available. Check input files and filename patterns.")
  }

  readr::write_csv(
    tibble(Year = complete_years),
    file.path(out_root, "available_years.csv")
  )

  year_samples <- list()
  within_all <- list()
  between_all <- list()
  annual_indicator_all <- list()

  for (yr in complete_years) {
    message("Processing year: ", yr)

    df_year <- sample_year_values(
      year = yr,
      variables_cfg = cfg$variables,
      variable_names = variable_names,
      sample_n = cfg$sampling$sample_n_per_year,
      seed = cfg$sampling$seed,
      block_size_lonlat = cfg$bootstrap$block_size_lonlat,
      block_size_projected = cfg$bootstrap$block_size_projected
    )
    year_samples[[as.character(yr)]] <- df_year

    df_vars <- df_year %>% select(all_of(variable_names))

    res_year <- calc_weights_one_year(
      df_vars = df_vars,
      year = yr,
      groups_fixed = groups_fixed,
      mer_lambda_within = cfg$mer$within_lambda,
      mer_lambda_between = cfg$mer$between_lambda,
      mer_eps = cfg$mer$eps,
      mer_candidate_weights = unlist(cfg$mer$candidate_weights)
    )

    within_all[[as.character(yr)]] <- res_year$within
    between_all[[as.character(yr)]] <- res_year$between
    annual_indicator_all[[as.character(yr)]] <- res_year$annual_indicator
  }

  within_tbl <- bind_rows(within_all)
  between_tbl <- bind_rows(between_all)
  annual_indicator_tbl <- bind_rows(annual_indicator_all)

  readr::write_csv(within_tbl, file.path(out_root, "yearly_within_group_weights.csv"))
  readr::write_csv(between_tbl, file.path(out_root, "yearly_between_group_weights.csv"))
  readr::write_csv(annual_indicator_tbl, file.path(out_root, "yearly_indicator_weights.csv"))
  saveRDS(year_samples, file.path(out_root, "year_samples.rds"))

  fixed_res <- calc_final_fixed_weights_hier(
    within_tbl = within_tbl,
    between_tbl = between_tbl,
    groups_fixed = groups_fixed,
    fixed_lambda_between = cfg$fixed$between_lambda,
    fixed_lambda_group = cfg$fixed$group_lambda,
    mer_eps = cfg$mer$eps
  )

  readr::write_csv(fixed_res$fixed_between, file.path(out_root, "fixed_between_group_weights.csv"))
  readr::write_csv(fixed_res$fixed_within, file.path(out_root, "fixed_within_group_weights.csv"))
  readr::write_csv(fixed_res$final_fixed, file.path(out_root, "final_fixed_weights.csv"))

  if (isTRUE(cfg$bootstrap$run)) {
    message("Running spatial block bootstrap ...")
    boot_list <- vector("list", cfg$bootstrap$n_reps)

    for (b in seq_len(cfg$bootstrap$n_reps)) {
      if (b %% 10 == 0) message("Bootstrap replicate: ", b, " / ", cfg$bootstrap$n_reps)

      within_boot_all <- list()
      between_boot_all <- list()

      for (yr in complete_years) {
        df_year <- year_samples[[as.character(yr)]]
        boot_df <- bootstrap_resample_one_year(
          df_year = df_year,
          max_points = cfg$bootstrap$max_points_per_year,
          seed = cfg$sampling$seed + b + yr
        )

        df_vars_boot <- boot_df %>% select(all_of(variable_names))

        res_boot <- calc_weights_one_year(
          df_vars = df_vars_boot,
          year = yr,
          groups_fixed = groups_fixed,
          mer_lambda_within = cfg$mer$within_lambda,
          mer_lambda_between = cfg$mer$between_lambda,
          mer_eps = cfg$mer$eps,
          mer_candidate_weights = unlist(cfg$mer$candidate_weights)
        )

        within_boot_all[[as.character(yr)]] <- res_boot$within
        between_boot_all[[as.character(yr)]] <- res_boot$between
      }

      within_boot_tbl <- bind_rows(within_boot_all)
      between_boot_tbl <- bind_rows(between_boot_all)

      fixed_boot <- calc_final_fixed_weights_hier(
        within_tbl = within_boot_tbl,
        between_tbl = between_boot_tbl,
        groups_fixed = groups_fixed,
        fixed_lambda_between = cfg$fixed$between_lambda,
        fixed_lambda_group = cfg$fixed$group_lambda,
        mer_eps = cfg$mer$eps
      )$final_fixed %>%
        select(Variable, FinalFixedWeight_HierMER) %>%
        mutate(Bootstrap = b)

      boot_list[[b]] <- fixed_boot
    }

    bootstrap_tbl <- bind_rows(boot_list)
    bootstrap_summary_tbl <- bootstrap_tbl %>%
      group_by(Variable) %>%
      summarise(
        median_bootstrap_weight = median(FinalFixedWeight_HierMER, na.rm = TRUE),
        mean_bootstrap_weight = mean(FinalFixedWeight_HierMER, na.rm = TRUE),
        sd_bootstrap_weight = sd(FinalFixedWeight_HierMER, na.rm = TRUE),
        ci_lower_95 = quantile(FinalFixedWeight_HierMER, 0.025, na.rm = TRUE),
        ci_upper_95 = quantile(FinalFixedWeight_HierMER, 0.975, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(median_bootstrap_weight))

    readr::write_csv(bootstrap_tbl, file.path(out_root, "bootstrap_fixed_weights_all_replicates.csv"))
    readr::write_csv(bootstrap_summary_tbl, file.path(out_root, "bootstrap_fixed_weights_summary.csv"))
  }

  invisible(
    list(
      yearly_within = within_tbl,
      yearly_between = between_tbl,
      yearly_indicator = annual_indicator_tbl,
      final_fixed = fixed_res$final_fixed
    )
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0) {
  cfg_path <- if (length(args) >= 1) args[1] else "config.yml"
  run_hdci_weighting(cfg_path)
}
