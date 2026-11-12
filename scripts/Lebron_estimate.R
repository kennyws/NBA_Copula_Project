library(hoopR)
library(dplyr)
library(MASS)
library(copula)

set.seed(1)

# hoopR season labels: 2026 = 2025-26, 2027 = 2026-27
lebron <- load_nba_player_box(seasons = c(2026, 2027)) %>%
  filter(
    athlete_display_name == "LeBron James",
    season_type == 2,
    !is.na(minutes), minutes > 0
  ) %>%
  arrange(game_date) %>%
  dplyr::select(game_date, points, rebounds, assists)

n <- nrow(lebron)
cat("Games in sample:", n, "\n")

# Thresholds: P(PTS >= 27, REB >= 7, AST >= 7)
th <- c(points = 27, rebounds = 7, assists = 7)

# ---- Gaussian copula (same recipe as your pipeline) ----
stats <- as.matrix(lebron[, c("points", "rebounds", "assists")])

z    <- qnorm(pobs(stats))
corr <- as.matrix(Matrix::nearPD(cor(z), corr = TRUE)$mat)

n_sim <- 1e6   # large, since this joint event is fairly rare
u_sim <- pnorm(MASS::mvrnorm(n_sim, mu = rep(0, 3), Sigma = corr))

sim_pts <- quantile(lebron$points,   u_sim[, 1], type = 1, names = FALSE)
sim_reb <- quantile(lebron$rebounds, u_sim[, 2], type = 1, names = FALSE)
sim_ast <- quantile(lebron$assists,  u_sim[, 3], type = 1, names = FALSE)

p_gaussian <- mean(sim_pts >= th["points"] &
                     sim_reb >= th["rebounds"] &
                     sim_ast >= th["assists"])

# ---- Benchmarks for context ----
p_empirical <- mean(lebron$points   >= th["points"] &
                      lebron$rebounds >= th["rebounds"] &
                      lebron$assists  >= th["assists"])

p_independent <- mean(lebron$points   >= th["points"]) *
  mean(lebron$rebounds >= th["rebounds"]) *
  mean(lebron$assists  >= th["assists"])

round(c(gaussian_copula = p_gaussian,
        empirical       = p_empirical,
        independence    = p_independent), 4)