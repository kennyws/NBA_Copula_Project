# NBA player boxscore joint-distribution comparison:
# Three forecasters of joint threshold events
#   1. Gaussian copula
#   2. Empirical joint rate
#   3. Independence (product of empirical marginal rates)
# Each run uses a separate two-season window and samples one starter per team.
# The 2020-21 NBA season (hoopR season label 2021) is intentionally excluded.


library(gsl)
library(copula)
library(hoopR)
library(dplyr)
library(rstanarm)
library(tidyr)
library(purrr)
library(ggplot2)
library(Matrix)
library(MASS)

set.seed(1)

# --------------------------------------------------
# RUN FIVE PERIOD-SPECIFIC SAMPLES
# --------------------------------------------------

run_periods <- list(
  Run_1_2025_2026 = c(2025, 2026),
  Run_2_2023_2024 = c(2023, 2024),
  Run_3_2018_2019 = c(2018, 2019),
  Run_4_2016_2017 = c(2016, 2017),
  Run_5_2014_2015 = c(2014, 2015),
  Run_6_2012_2013 = c(2012, 2013),
  Run_7_2010_2011 = c(2010, 2011)
)

seasons_to_load <- sort(unique(unlist(run_periods)))

boxscores <- load_nba_player_box(seasons = seasons_to_load)

# Filtering for regular season games, where players played non-zero minutes
boxscores <- boxscores %>%
  filter(
    season_type == 2,
    team_id %in% 1:30,
    !is.na(minutes),
    minutes > 0
  )

full_data <- boxscores %>%
  dplyr::select(
    season, game_date, game_id, athlete_id, athlete_display_name,
    team_id, starter, points, assists, rebounds
  )

# --------------------------------------------------
# QUANTILE GRID
# --------------------------------------------------

grid_quantiles <- c(
  Q25 = 0.25,
  Q50 = 0.50,
  Q75 = 0.75
)

lv <- names(grid_quantiles)

# --------------------------------------------------
# LOG LOSS HELPER - avoiding 0s and 1s
# --------------------------------------------------

ll <- function(p, y) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  -(y * log(p) + (1 - y) * log(1 - p))
}

# --------------------------------------------------
# PER-PLAYER WALK-FORWARD PIPELINE
# --------------------------------------------------

run_player_pipeline <- function(player_id, player_name, team_id_val, full_data) {
  
  data <- full_data %>%
    filter(athlete_id == player_id) %>%
    arrange(game_date)
  
  if (nrow(data) < 20) {
    stop(
      "Not enough games (", nrow(data),
      ") for a stable walk-forward fit"
    )
  }
  
  player_results <- map_dfr(21:nrow(data), function(t) {
    
    train <- data[1:(t - 1), ]
    actual <- data[t, ]
    
    # Three-dimensional training data
    train_stats <- train %>%
      dplyr::select(points, assists, rebounds)
    
    u <- pobs(as.matrix(train_stats))
    
    # ---- Multivariate Gaussian copula ----
    
    gauss_result <- tryCatch({
      
      z <- qnorm(u)
      
      if (any(!is.finite(z))) {
        stop("non-finite z-scores (likely a constant column)")
      }
      
      corr_matrix <- cor(z)
      
      if (any(!is.finite(corr_matrix))) {
        stop("non-finite correlation matrix")
      }
      
      corr_matrix <- as.matrix(
        Matrix::nearPD(corr_matrix, corr = TRUE)$mat
      )
      
      sim_gauss <- MASS::mvrnorm(
        1000,
        mu = rep(0, ncol(corr_matrix)),
        Sigma = corr_matrix
      )
      
      sim_gauss <- pnorm(sim_gauss)
      colnames(sim_gauss) <- colnames(train_stats)
      
      list(
        pts = quantile(train$points,   sim_gauss[, "points"],   type = 7),
        ast = quantile(train$assists,  sim_gauss[, "assists"],  type = 7),
        reb = quantile(train$rebounds, sim_gauss[, "rebounds"], type = 7),
        ok = TRUE
      )
      
    }, error = function(e) {
      
      message(
        "    Gaussian fit failed at t=", t,
        " for ", player_name, ": ",
        conditionMessage(e)
      )
      
      list(ok = FALSE)
    })
    
    # ---- Joint threshold grid: all (PTS, REB, AST) quantile combos ----
    
    grid <- expand.grid(
      pts_level = lv,
      reb_level = lv,
      ast_level = lv,
      stringsAsFactors = FALSE
    )
    
    # Thresholds computed once per game, then looked up by level name
    pts_ths <- quantile(train$points,   grid_quantiles, type = 7, names = FALSE)
    reb_ths <- quantile(train$rebounds, grid_quantiles, type = 7, names = FALSE)
    ast_ths <- quantile(train$assists,  grid_quantiles, type = 7, names = FALSE)
    names(pts_ths) <- names(reb_ths) <- names(ast_ths) <- lv
    
    grid <- grid %>%
      mutate(
        pts_th = unname(pts_ths[pts_level]),
        reb_th = unname(reb_ths[reb_level]),
        ast_th = unname(ast_ths[ast_level]),
        game_id = actual$game_id,
        game_date = actual$game_date,
        n_training_games = t - 1
      )
    
    # 1. Gaussian copula probability for each grid row
    grid$p_gaussian <- if (gauss_result$ok) {
      mapply(
        function(p, r, a) {
          mean(gauss_result$pts >= p &
                 gauss_result$reb >= r &
                 gauss_result$ast >= a)
        },
        grid$pts_th, grid$reb_th, grid$ast_th
      )
    } else {
      1e-6
    }
    
    # 2. Empirical joint probability for each grid row
    grid$p_empirical <- mapply(
      function(p, r, a) {
        mean(train$points >= p &
               train$rebounds >= r &
               train$assists >= a)
      },
      grid$pts_th, grid$reb_th, grid$ast_th
    )
    
    # 3. Independence: product of the empirical marginal exceedance rates
    #    P(PTS >= p) * P(REB >= r) * P(AST >= a)
    grid$p_independent <- mapply(
      function(p, r, a) {
        mean(train$points   >= p) *
          mean(train$rebounds >= r) *
          mean(train$assists  >= a)
      },
      grid$pts_th, grid$reb_th, grid$ast_th
    )
    
    # Actual joint outcome for each grid row
    grid$outcome <- as.integer(
      actual$points   >= grid$pts_th &
        actual$rebounds >= grid$reb_th &
        actual$assists  >= grid$ast_th
    )
    
    grid
  })
  
  player_results %>%
    mutate(
      athlete_id = player_id,
      athlete_display_name = player_name,
      team_id = team_id_val
    )
}

# --------------------------------------------------
# RUN PLAYERS WITH A FALLBACK OPTION
# --------------------------------------------------

run_team_with_fallback <- function(team_candidates, period_data, n_players = 1) {
  
  team_id_val <- team_candidates$team_id[1]
  
  candidates <- team_candidates %>%
    distinct(athlete_id, athlete_display_name)
  
  results <- list()
  
  for (i in seq_len(nrow(candidates))) {
    
    if (length(results) >= n_players) break
    
    player_id <- candidates$athlete_id[i]
    player_name <- candidates$athlete_display_name[i]
    
    message("Running: ", player_name, " (team ", team_id_val, ")")
    
    result <- tryCatch({
      run_player_pipeline(
        player_id = player_id,
        player_name = player_name,
        team_id_val = team_id_val,
        full_data = period_data
      )
    }, error = function(e) {
      message("  Skipped ", player_name, " - ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(result)) {
      results[[length(results) + 1]] <- result
    }
  }
  
  if (length(results) < n_players) {
    message(
      "  Only ", length(results), " of ", n_players,
      " players succeeded for team_id ", team_id_val
    )
  }
  
  if (length(results) == 0) return(NULL)
  
  bind_rows(results)
}

# --------------------------------------------------
# RUN ONE TWO-SEASON PERIOD
# Randomly sample two games per team, then take one successful starter.
# --------------------------------------------------

make_candidates <- function(period_data) {
  
  sampled_games <- period_data %>%
    distinct(team_id, game_id) %>%
    group_by(team_id) %>%
    slice_sample(prop = 1) %>%
    slice_head(n = 2) %>%
    ungroup()
  
  sampled_games %>%
    left_join(
      period_data %>%
        filter(starter == TRUE) %>%
        distinct(
          team_id, game_id, athlete_id,
          athlete_display_name
        ),
      by = c("team_id", "game_id")
    ) %>%
    filter(!is.na(athlete_id)) %>%
    distinct(team_id, athlete_id, athlete_display_name) %>%
    group_by(team_id) %>%
    slice_sample(prop = 1) %>%
    ungroup() %>%
    group_by(team_id) %>%
    group_split()
}

run_one_period <- function(period_name, period_seasons) {
  
  message("\n========== ", period_name, " ==========")
  
  period_data <- full_data %>%
    filter(season %in% period_seasons) %>%
    arrange(game_date)
  
  if (nrow(period_data) == 0) {
    stop("No data found for ", period_name)
  }
  
  candidates_by_team <- make_candidates(period_data)
  
  period_results <- map_dfr(
    candidates_by_team,
    run_team_with_fallback,
    period_data = period_data,
    n_players = 1,
    .progress = TRUE
  )
  
  if (nrow(period_results) == 0) {
    stop("No players completed successfully for ", period_name)
  }
  
  period_results %>%
    mutate(
      run_id = period_name,
      seasons = paste(period_seasons, collapse = "-")
    )
}

all_results <- imap_dfr(
  run_periods,
  ~ run_one_period(.y, .x)
)

# --------------------------------------------------
# SANITY CHECKS
# --------------------------------------------------

# Number of teams and players successfully sampled in each period.
run_counts <- all_results %>%
  distinct(run_id, team_id, athlete_id, athlete_display_name) %>%
  group_by(run_id) %>%
  summarise(
    n_teams = n_distinct(team_id),
    n_players = n(),
    .groups = "drop"
  )

print(run_counts, n = Inf)

# ==================================================
# Validation: Gaussian copula vs empirical joint vs independence
#
# - All runs pooled; every observation is one player-game forecast.
# - Observations are not averaged within player first.
# - Uncertainty comes from a cluster bootstrap over player-runs (one
#   player sampled in one run).
#
# Differences are always (loss of the first) - (loss of the second),
# so a POSITIVE value means the second method is better:
#   diff_emp_gauss = empirical   - Gaussian    (positive = Gaussian wins)
#   diff_ind_gauss = independent - Gaussian    (positive = Gaussian wins)
#   diff_ind_emp   = independent - empirical   (positive = empirical wins)
# ==================================================

set.seed(1)
B <- 10000   # bootstrap resamples

diag_results <- all_results %>%
  mutate(
    across(c(pts_level, reb_level, ast_level), ~ factor(.x, levels = lv)),
    across(c(p_gaussian, p_empirical, p_independent),
           ~ pmin(pmax(.x, 1e-6), 1 - 1e-6)),
    log_loss_gaussian     = ll(p_gaussian, outcome),
    log_loss_empirical    = ll(p_empirical, outcome),
    log_loss_independent  = ll(p_independent, outcome),
    diff_emp_gauss = log_loss_empirical   - log_loss_gaussian,
    diff_ind_gauss = log_loss_independent - log_loss_gaussian,
    diff_ind_emp   = log_loss_independent - log_loss_empirical,
    combo = paste0("PTS ", pts_level, " | REB ", reb_level, " | AST ", ast_level),
    # cluster = one player in one run
    cluster_id = paste(run_id, athlete_id, sep = "_")
  )

n_players <- n_distinct(diag_results$cluster_id)
n_runs    <- n_distinct(diag_results$run_id)

combo_vars <- c("pts_level", "reb_level", "ast_level", "combo")

# --------------------------------------------------
# Diagnostic: how often does each estimator say "never" or "always"?
# --------------------------------------------------

diag_results %>%
  group_by(across(all_of(combo_vars))) %>%
  summarise(
    n_obs = n(),
    observed_event_rate        = mean(outcome),
    mean_p_gaussian            = mean(p_gaussian),
    mean_p_empirical           = mean(p_empirical),
    mean_p_independent         = mean(p_independent),
    share_empirical_zero       = mean(p_empirical <= 1e-6),
    share_empirical_one        = mean(p_empirical >= 1 - 1e-6),
    share_independent_zero     = mean(p_independent <= 1e-6),
    share_independent_one      = mean(p_independent >= 1 - 1e-6),
    .groups = "drop"
  ) %>%
  print(n = Inf, width = Inf)

# --------------------------------------------------
# TEST 1: paired t-test per combination, for each of the 3 comparisons
# (Holm across the combinations within each comparison)
# --------------------------------------------------

comparisons <- tibble::tribble(
  ~comparison,        ~loss_a,                 ~loss_b,               ~diff_col,
  "Empirical vs Gaussian",   "log_loss_empirical",   "log_loss_gaussian",  "diff_emp_gauss",
  "Independent vs Gaussian", "log_loss_independent", "log_loss_gaussian",  "diff_ind_gauss",
  "Independent vs Empirical","log_loss_independent", "log_loss_empirical", "diff_ind_emp"
)

tests_t <- pmap_dfr(comparisons, function(comparison, loss_a, loss_b, diff_col) {
  diag_results %>%
    group_by(across(all_of(combo_vars))) %>%
    group_modify(~ {
      tt <- t.test(.x[[loss_a]], .x[[loss_b]], paired = TRUE)
      tibble(
        n_obs = nrow(.x),
        mean_diff = mean(.x[[diff_col]]),
        t_stat = unname(tt$statistic),
        p_t = tt$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(comparison = comparison,
           p_t_holm = p.adjust(p_t, method = "holm"))
})

# --------------------------------------------------
# TEST 2: cluster bootstrap per combination
# Resample player-runs (with all of their games), recompute the pooled mean.
# Gives a 95% CI and a two-sided p-value that respect within-player dependence.
# --------------------------------------------------

cluster_boot <- function(df, col) {
  per_cluster <- tibble(cluster_id = df$cluster_id, x = df[[col]]) %>%
    group_by(cluster_id) %>%
    summarise(s = sum(x), n = n(), .groups = "drop")
  
  P   <- nrow(per_cluster)
  idx <- matrix(sample(P, P * B, replace = TRUE), nrow = B)
  means <- rowSums(matrix(per_cluster$s[idx], nrow = B)) /
    rowSums(matrix(per_cluster$n[idx], nrow = B))
  
  tibble(
    estimate = sum(per_cluster$s) / sum(per_cluster$n),
    lo = unname(quantile(means, 0.025)),
    hi = unname(quantile(means, 0.975)),
    p_boot = max(1 / B, min(1, 2 * min(mean(means <= 0), mean(means >= 0))))
  )
}

boot_measures <- c(
  "log_loss_gaussian", "log_loss_empirical", "log_loss_independent",
  "diff_emp_gauss", "diff_ind_gauss", "diff_ind_emp"
)

boot_tbl <- map_dfr(boot_measures, function(col) {
  diag_results %>%
    group_by(across(all_of(combo_vars))) %>%
    group_modify(~ cluster_boot(.x, col)) %>%
    ungroup() %>%
    mutate(measure = col)
})

diff_measures <- c("diff_emp_gauss", "diff_ind_gauss", "diff_ind_emp")

tests_boot <- boot_tbl %>%
  filter(measure %in% diff_measures) %>%
  group_by(measure) %>%
  mutate(p_boot_holm = p.adjust(p_boot, method = "holm")) %>%
  ungroup() %>%
  left_join(comparisons %>% dplyr::select(comparison, diff_col),
            by = c("measure" = "diff_col")) %>%
  dplyr::select(
    comparison, all_of(combo_vars),
    mean_diff_boot = estimate, ci_lo = lo, ci_hi = hi,
    p_boot, p_boot_holm
  )

cat("\n--- Pairwise comparisons, pooled across runs, by quantile combination ---\n")
cat("(positive mean_diff = the SECOND method named is better)\n")
tests_t %>%
  left_join(tests_boot, by = c("comparison", combo_vars)) %>%
  arrange(comparison, desc(mean_diff)) %>%
  print(n = Inf, width = Inf)

# --------------------------------------------------
# Summary: how many of the combinations does each method win?
# (a "win" requires the bootstrap CI to exclude 0)
# --------------------------------------------------

cat("\n--- Number of combinations with a bootstrap CI excluding 0 ---\n")
tests_boot %>%
  group_by(comparison) %>%
  summarise(
    first_method_better  = sum(ci_hi < 0),
    second_method_better = sum(ci_lo > 0),
    no_clear_difference  = sum(ci_lo <= 0 & ci_hi >= 0),
    .groups = "drop"
  ) %>%
  print(width = Inf)

# --------------------------------------------------
# PLOT 1: mean log loss by method, all combinations
# (95% cluster bootstrap CI)
# --------------------------------------------------

plot_loss <- boot_tbl %>%
  filter(measure %in% c("log_loss_gaussian", "log_loss_empirical",
                        "log_loss_independent")) %>%
  mutate(model = recode(measure,
                        log_loss_gaussian    = "Gaussian copula",
                        log_loss_empirical   = "Empirical joint rate",
                        log_loss_independent = "Independence"))

ggplot(plot_loss, aes(reorder(combo, estimate), estimate, color = model)) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.7), size = 0.3) +
  coord_flip() +
  labs(
    title = "Log loss of joint threshold forecasts by quantile combination",
    subtitle = paste0("Pooled across ", n_runs, " runs and ", n_players,
                      " player-runs; bars are 95% cluster bootstrap CIs"),
    x = NULL, y = "Mean log loss (lower is better)", color = NULL
  ) +
  theme_minimal()

# --------------------------------------------------
# PLOT 2: pairwise heatmaps (facet = assist level)
# --------------------------------------------------

make_heatmap <- function(measure_name, title, subtitle, legend_name) {
  
  plot_diff <- boot_tbl %>%
    filter(measure == measure_name) %>%
    mutate(sig = ifelse(lo > 0 | hi < 0, "*", ""))
  
  ggplot(plot_diff, aes(pts_level, reb_level, fill = estimate)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(sprintf("%.4f", estimate), sig)), size = 3) +
    facet_wrap(~ paste("AST", ast_level)) +
    scale_fill_gradient2(
      low = "firebrick", mid = "white", high = "steelblue",
      midpoint = 0, name = legend_name
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Points quantile level", y = "Rebounds quantile level"
    ) +
    theme_minimal()
}

# (a) Gaussian copula vs empirical joint rate
make_heatmap(
  "diff_emp_gauss",
  "Log loss improvement of the Gaussian copula over the empirical joint rate",
  "Blue = copula better; * = 95% cluster bootstrap CI excludes 0",
  "Empirical \u2212\nGaussian"
)

# (b) Gaussian copula vs independence
make_heatmap(
  "diff_ind_gauss",
  "Log loss improvement of the Gaussian copula over independence",
  "Blue = copula better; * = 95% cluster bootstrap CI excludes 0",
  "Independent \u2212\nGaussian"
)

# (c) Empirical joint rate vs independence
make_heatmap(
  "diff_ind_emp",
  "Log loss improvement of the empirical joint rate over independence",
  "Blue = empirical joint better; * = 95% cluster bootstrap CI excludes 0",
  "Independent \u2212\nEmpirical"
)

# --------------------------------------------------
# RELATIVE LOG LOSS IMPROVEMENT OF THE GAUSSIAN COPULA
#   rel_improvement = 1 - (Gaussian log loss / comparison log loss)
#   > 0 : Gaussian better (e.g. 0.05 = 5% lower log loss)
#   < 0 : Gaussian worse
# Computed per quantile combination as a ratio of pooled mean log losses.
# The cluster bootstrap resamples player-runs and recomputes the whole
# ratio each time, so the CI reflects uncertainty in the ratio itself.
# --------------------------------------------------

cluster_boot_rel <- function(df, num_col, den_col) {
  per_cluster <- df %>%
    group_by(cluster_id) %>%
    summarise(
      num = sum(.data[[num_col]]),
      den = sum(.data[[den_col]]),
      .groups = "drop"
    )
  
  P   <- nrow(per_cluster)
  idx <- matrix(sample(P, P * B, replace = TRUE), nrow = B)
  
  # n cancels in the ratio of means, so summed losses are enough
  stats <- 1 - rowSums(matrix(per_cluster$num[idx], nrow = B)) /
    rowSums(matrix(per_cluster$den[idx], nrow = B))
  
  tibble(
    estimate = 1 - sum(per_cluster$num) / sum(per_cluster$den),
    lo = unname(quantile(stats, 0.025)),
    hi = unname(quantile(stats, 0.975))
  )
}

rel_comparisons <- tibble::tribble(
  ~comparison,              ~den_col,
  "Independence",           "log_loss_independent",
  "Empirical joint rate",   "log_loss_empirical"
)

rel_tbl <- pmap_dfr(rel_comparisons, function(comparison, den_col) {
  diag_results %>%
    group_by(across(all_of(combo_vars))) %>%
    group_modify(~ cluster_boot_rel(.x, "log_loss_gaussian", den_col)) %>%
    ungroup() %>%
    mutate(comparison = comparison)
})

make_rel_plot <- function(comp) {
  
  d <- rel_tbl %>%
    filter(comparison == comp) %>%
    mutate(sig = lo > 0 | hi < 0)
  
  ggplot(d, aes(reorder(combo, estimate), estimate, color = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3) +
    coord_flip() +
    scale_color_manual(
      values = c(`TRUE` = "steelblue", `FALSE` = "gray60"),
      labels = c(`TRUE` = "CI excludes 0", `FALSE` = "CI includes 0"),
      name = NULL
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = paste0("Relative log loss improvement of the Gaussian copula over ",
                     tolower(comp)),
      subtitle = paste0("1 \u2212 (Gaussian log loss / ", tolower(comp),
                        " log loss)"),
      x = NULL,
      y = "Relative improvement (positive = Gaussian better)"
    ) +
    theme_minimal()
}

# Plot A: Gaussian copula vs independence
make_rel_plot("Independence")

# Plot B: Gaussian copula vs empirical joint rate
make_rel_plot("Empirical joint rate")

# --------------------------------------------------
# Summary: Average log loss for selected combinations
# --------------------------------------------------

reference_combos <- c(
  "PTS Q25 | REB Q25 | AST Q25",
  "PTS Q50 | REB Q50 | AST Q50",
  "PTS Q75 | REB Q75 | AST Q75"
)

reference_table <- boot_tbl %>%
  filter(
    combo %in% reference_combos,
    measure %in% c(
      "log_loss_gaussian",
      "log_loss_empirical",
      "log_loss_independent"
    )
  ) %>%
  mutate(
    model = recode(
      measure,
      log_loss_gaussian = "Gaussian",
      log_loss_empirical = "Empirical",
      log_loss_independent = "Independence"
    )
  ) %>%
  dplyr::select(combo, model, estimate) %>%
  tidyr::pivot_wider(
    names_from = model,
    values_from = estimate
  ) %>%
  arrange(match(combo, reference_combos))

print(reference_table, width = Inf)

# ==================================================
# Calibration curves (all three models)
# Equal-count bins per model; CIs from a cluster bootstrap over player-runs.
# ==================================================

n_bins <- 10
B_cal  <- 2000
set.seed(1)

calib_long <- diag_results %>%
  dplyr::select(cluster_id, combo, outcome,
                p_gaussian, p_empirical, p_independent) %>%
  pivot_longer(c(p_gaussian, p_empirical, p_independent),
               names_to = "model", values_to = "p") %>%
  mutate(model = recode(model,
                        p_gaussian    = "Gaussian copula",
                        p_empirical   = "Empirical joint rate",
                        p_independent = "Independence"))

calib_one <- function(d, n_bins = 10, B = 2000) {
  # Equal-count bins; unique() handles models with many tied forecasts (e.g. zeros)
  breaks <- unique(quantile(d$p, probs = seq(0, 1, length.out = n_bins + 1)))
  d$bin  <- cut(d$p, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  per <- d %>%
    group_by(cluster_id, bin) %>%
    summarise(y = sum(outcome), n = n(), ps = sum(p), .groups = "drop")
  
  clusters <- unique(per$cluster_id)
  bins     <- sort(unique(per$bin))
  Y <- N <- PS <- matrix(0, length(clusters), length(bins))
  ij <- cbind(match(per$cluster_id, clusters), match(per$bin, bins))
  Y[ij] <- per$y; N[ij] <- per$n; PS[ij] <- per$ps
  
  P <- length(clusters)
  boot_obs <- replicate(B, {
    idx <- sample.int(P, P, replace = TRUE)
    colSums(Y[idx, , drop = FALSE]) / colSums(N[idx, , drop = FALSE])
  })
  if (is.null(dim(boot_obs))) boot_obs <- matrix(boot_obs, nrow = 1)
  
  tibble(
    bin  = bins,
    n    = colSums(N),
    pred = colSums(PS) / colSums(N),   # mean forecast in the bin
    obs  = colSums(Y) / colSums(N),    # observed event rate in the bin
    lo   = apply(boot_obs, 1, quantile, probs = 0.025, na.rm = TRUE),
    hi   = apply(boot_obs, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

make_calib <- function(data, title) {
  tbl <- data %>%
    group_by(model) %>%
    group_modify(~ calib_one(.x, n_bins, B_cal)) %>%
    ungroup()
  
  print(
    tbl %>%
      group_by(model) %>%
      summarise(ECE = sum(n * abs(obs - pred)) / sum(n), .groups = "drop")
  )  # expected calibration error (lower is better)
  
  ggplot(tbl, aes(pred, obs, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_line() +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3) +
    facet_wrap(~ model) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    guides(color = "none") +
    labs(
      title = title,
      subtitle = paste0(n_bins, " equal-count bins; bars are 95% cluster bootstrap CIs"),
      x = "Mean predicted probability", y = "Observed frequency"
    ) +
    theme_minimal()
}

# (a) All combinations pooled
make_calib(calib_long, "Calibration of joint threshold forecasts (all combinations)")

# (b) The high-threshold tail, where the models differ most
make_calib(
  calib_long %>% filter(combo == "PTS Q75 | REB Q75 | AST Q75"),
  "Calibration: PTS Q75 | REB Q75 | AST Q75"
)

# ==================================================
# Simple correlation matrix: PTS, REB, AST
# Uses the raw stat lines of every evaluated player-game.
# ==================================================

eval_games <- all_results %>%
  distinct(run_id, athlete_id, game_id) %>%
  inner_join(
    full_data %>% dplyr::select(athlete_id, game_id, points, rebounds, assists),
    by = c("athlete_id", "game_id")
  ) %>%
  mutate(cluster_id = paste(run_id, athlete_id, sep = "_")) %>%
  rename(PTS = points, REB = rebounds, AST = assists)

stat_cols <- c("PTS", "REB", "AST")

# Within-player: correlation computed for each player-run, then averaged.

within <- eval_games %>%
  group_by(cluster_id) %>%
  filter(n() >= 20) %>%
  summarise(
    PTS_REB = suppressWarnings(cor(PTS, REB, method = "spearman")),
    PTS_AST = suppressWarnings(cor(PTS, AST, method = "spearman")),
    REB_AST = suppressWarnings(cor(REB, AST, method = "spearman")),
    .groups = "drop"
  )

pair_summary <- within %>%
  pivot_longer(-cluster_id, names_to = "pair", values_to = "rho") %>%
  group_by(pair) %>%
  summarise(
    mean = mean(rho, na.rm = TRUE),
    p10  = quantile(rho, 0.10, na.rm = TRUE),
    p90  = quantile(rho, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  separate(pair, into = c("a", "b"), sep = "_")

# Fill both triangles symmetrically, plus a diagonal of 1s
corr_long <- bind_rows(
  pair_summary,
  pair_summary %>% rename(a = b, b = a),
  tibble(a = stat_cols, b = stat_cols, mean = 1, p10 = NA_real_, p90 = NA_real_)
) %>%
  mutate(
    a = factor(a, levels = stat_cols),
    b = factor(b, levels = rev(stat_cols)),
    label = ifelse(
      is.na(p10),
      "1.00",
      sprintf("%.2f\n(%.2f, %.2f)", mean, p10, p90)
    )
  )

ggplot(corr_long, aes(a, b, fill = mean)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(
    low = "firebrick", mid = "white", high = "steelblue",
    midpoint = 0, limits = c(-1, 1), name = "Mean\nSpearman"
  ) +
  coord_equal() +
  labs(
    title = "Within-player correlation of box score statistics",
    subtitle = paste0("Mean Spearman correlation across all players; (10th, 90th percentile) across players"),
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# Correlation Box Plots

library(tidyverse)

player_cor <- eval_games %>%
  group_by(cluster_id) %>%
  filter(n() >= 20) %>%
  summarise(
    `PTS vs REB` = suppressWarnings(cor(PTS, REB, method = "spearman")),
    `PTS vs AST` = suppressWarnings(cor(PTS, AST, method = "spearman")),
    `REB vs AST` = suppressWarnings(cor(REB, AST, method = "spearman")),
    .groups = "drop"
  ) %>%
  pivot_longer(-cluster_id, names_to = "pair", values_to = "rho")

ggplot(player_cor, aes(x = pair, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.shape = NA, fill = "steelblue", alpha = 0.25,
               color = "steelblue4", width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1, color = "grey30") +
  labs(
    title = "Within-player correlation of box score statistics",
    subtitle = "Each point is one player's Spearman correlation of the pairwise relationship",
    x = NULL, y = "Spearman correlation"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())


