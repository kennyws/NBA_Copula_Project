# 3D plots for the Gaussian copula analysis
#
# Run this AFTER the main script (it uses full_data, all_results, run_periods).
# Extra packages: install.packages(c("plotly", "ks"))

library(plotly)
library(ks)
library(dplyr)
library(purrr)
library(tibble)
library(copula)
library(Matrix)

set.seed(1)

# --------------------------------------------------
# SHARED HELPERS
# --------------------------------------------------

stat_cols <- c("points", "rebounds", "assists")

# All games for one player (name pattern, regex) in the given seasons
get_player_data <- function(name_pattern, seasons) {
  full_data %>%
    filter(grepl(name_pattern, athlete_display_name, ignore.case = TRUE),
           season %in% seasons) %>%
    arrange(game_date)
}

# Draws from the fitted Gaussian copula, mapped back to the player's marginals
simulate_gaussian_copula <- function(df, n = 1e5) {
  X <- as.matrix(df[, stat_cols])
  z <- qnorm(pobs(X))
  R <- as.matrix(Matrix::nearPD(cor(z), corr = TRUE)$mat)
  s <- pnorm(MASS::mvrnorm(n, rep(0, 3), R))
  tibble(
    points   = quantile(df$points,   s[, 1], type = 7, names = FALSE),
    rebounds = quantile(df$rebounds, s[, 2], type = 7, names = FALSE),
    assists  = quantile(df$assists,  s[, 3], type = 7, names = FALSE)
  )
}

# Draws with the same marginals but independent columns
simulate_independent <- function(df, n = 1e5) {
  tibble(
    points   = sample(df$points,   n, replace = TRUE),
    rebounds = sample(df$rebounds, n, replace = TRUE),
    assists  = sample(df$assists,  n, replace = TRUE)
  )
}

# ==================================================
# PART 1: JOINT-PROBABILITY RATIO SURFACE
# ==================================================

surface_matrices <- function(df,
                             ast_q = 0.75,
                             q_seq = seq(0.05, 0.95, by = 0.05),
                             n_sim = 1e5,
                             baseline = c("independent", "empirical")) {
  baseline <- match.arg(baseline)
  
  sim <- simulate_gaussian_copula(df, n_sim)
  
  pt <- quantile(df$points,   q_seq, type = 7, names = FALSE)
  rt <- quantile(df$rebounds, q_seq, type = 7, names = FALSE)
  at <- quantile(df$assists,  ast_q, type = 7, names = FALSE)
  
  sim_ast_ok <- sim$assists >= at
  df_ast_ok  <- df$assists  >= at
  
  # rows = rebounds level (y axis), columns = points level (x axis)
  k   <- length(q_seq)
  p_g <- p_b <- matrix(NA_real_, k, k)
  
  for (i in seq_len(k)) {
    sim_pts_ok <- sim$points >= pt[i]
    df_pts_ok  <- df$points  >= pt[i]
    for (j in seq_len(k)) {
      p_g[j, i] <- mean(sim_pts_ok & sim$rebounds >= rt[j] & sim_ast_ok)
      p_b[j, i] <- if (baseline == "independent") {
        mean(df_pts_ok) * mean(df$rebounds >= rt[j]) * mean(df_ast_ok)
      } else {
        mean(df_pts_ok & df$rebounds >= rt[j] & df_ast_ok)
      }
    }
  }
  
  # Floors keep the ratio finite: a simulated zero is worth 1/n_sim,
  # an empirical zero is worth half an observation.
  p_g <- pmax(p_g, 1 / n_sim)
  p_b <- pmax(p_b, 0.5 / nrow(df))
  
  list(q_seq = q_seq, ast_q = ast_q, baseline = baseline,
       p_gauss = p_g, p_base = p_b, ratio = p_g / p_b)
}

plot_ratio_surface <- function(m, title = NULL) {
  
  base_label <- if (m$baseline == "independent") "independence" else "empirical joint rate"
  
  if (is.null(title)) {
    title <- paste0("Gaussian copula vs ", base_label,
                    ": joint probability ratio (assists \u2265 Q",
                    round(100 * m$ast_q), ")")
  }
  
  plot_ly() %>%
    add_surface(
      x = m$q_seq, y = m$q_seq, z = m$ratio,
      colorscale = "Viridis",
      colorbar = list(title = "Gaussian /<br>baseline")
    ) %>%
    add_surface(
      x = m$q_seq, y = m$q_seq,
      z = matrix(1, length(m$q_seq), length(m$q_seq)),
      showscale = FALSE, opacity = 0.25,
      colorscale = list(c(0, "gray"), c(1, "gray"))
    ) %>%
    plotly::layout(
      title = list(
        text = title,
        y = 0.98,
        pad = list(b = 20)
      ),
      margin = list(t = 100, b = 50, l = 50, r = 50),
      scene = list(
        aspectratio = list(x = 1, y = 1, z = 1),
        xaxis = list(
          title = "Points quantile level",
          titlefont = list(size = 12),
          tickfont = list(size = 10)
        ),
        yaxis = list(
          title = "Rebounds quantile level",
          titlefont = list(size = 12),
          tickfont = list(size = 10)
        ),
        zaxis = list(title = "P(Gaussian) / P(baseline)")
      )
    )
}

# ---- Pooled surface across every player-run in all_results ----
# Geometric mean of the per-player ratios (average of the log ratios).
# Uses a coarser grid and fewer simulations to keep runtime reasonable.
pooled_surface <- function(ast_q = ast_q,
                           baseline = "independent",
                           q_seq = seq(0.1, 0.9, by = 0.1),
                           n_sim = 2e4) {
  
  keys <- all_results %>% distinct(run_id, athlete_id)
  
  log_ratios <- map2(keys$run_id, keys$athlete_id, function(rid, aid) {
    d <- full_data %>%
      filter(athlete_id == aid, season %in% run_periods[[rid]]) %>%
      arrange(game_date)
    tryCatch(
      log(surface_matrices(d, ast_q, q_seq, n_sim, baseline)$ratio),
      error = function(e) NULL
    )
  }) %>%
    compact()
  
  list(q_seq = q_seq, ast_q = ast_q, baseline = baseline,
       ratio = exp(Reduce(`+`, log_ratios) / length(log_ratios)),
       n_player_runs = length(log_ratios))
}

pooled <- pooled_surface(ast_q = 0.75, baseline = "independent")
plot_ratio_surface(
  pooled,
  title = paste0("Pooled over ", pooled$n_player_runs,
                 " players: Gaussian copula / independence (assists \u2265 Q75)")
)

pooled <- pooled_surface(ast_q = 0.25, baseline = "independent")
plot_ratio_surface(
  pooled,
  title = paste0("Pooled over ", pooled$n_player_runs,
                 " players: Gaussian copula / independence (assists \u2265 Q25)")
)
