library(hoopR)
library(dplyr)
library(rstanarm)
library(VineCopula)
library(tidyr)
library(purrr)
library(ggplot2)

set.seed(1)  # reproducible player sampling

# --------------------------------------------------
# LOAD + PREP FULL LEAGUE DATA (once, not per player)
# --------------------------------------------------

boxscores <- load_nba_player_box(seasons = c(2026))

boxscores <- boxscores %>%
  filter(season_type == 2) %>%
  filter(team_id %in% c(1:30)) %>%
  filter(!is.na(minutes), minutes > 0)

# Team Level Possessions

team_game_poss <- boxscores %>%
  group_by(season, game_id, team_id, opponent_team_id) %>%
  summarise(
    fga  = sum(field_goals_attempted, na.rm = TRUE),
    fta  = sum(free_throws_attempted, na.rm = TRUE),
    oreb = sum(offensive_rebounds, na.rm = TRUE),
    tov  = sum(turnovers, na.rm = TRUE),
    game_date = first(game_date),
    .groups = "drop"
  ) %>%
  mutate(
    possessions = fga + 0.44 * fta + tov - oreb
  ) %>%
  arrange(season, team_id, game_date)

# Team Rolling Possessions

team_game_poss <- team_game_poss %>%
  group_by(season, team_id) %>%
  mutate(
    poss_per_min = possessions / 48,
    rolling_poss_per_min = lag(cummean(poss_per_min))
  ) %>%
  ungroup()

# Joining and Adding Player Estimated Possessions

boxscores <- boxscores %>%
  left_join(
    team_game_poss %>% dplyr::select(season, game_id, team_id, team_possessions = possessions,
                                     team_rolling_poss_per_min = rolling_poss_per_min),
    by = c("season", "game_id", "team_id")
  ) %>%
  mutate(
    player_possessions = (team_possessions / 48) * minutes
  )

full_data <- boxscores %>%
  dplyr::select(game_date, game_id, player_possessions, athlete_id, athlete_display_name,
                team_id, starter, points, assists, rebounds, three_point_field_goals_made)

# --------------------------------------------------
# BUILD SHUFFLED CANDIDATE LIST PER TEAM
# (not just one sampled starter - a ranked fallback list per team)
# --------------------------------------------------

candidate_pool <- full_data %>%
  # First randomly select one game for each team
  distinct(team_id, game_id) %>%
  group_by(team_id) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # Get the players from those selected games
  left_join(
    full_data %>%
      filter(starter == TRUE) %>%
      distinct(team_id, game_id, athlete_id, athlete_display_name),
    by = c("team_id", "game_id")
  )

candidates_by_team <- candidate_pool %>%
  group_by(team_id) %>%
  group_split()

# --------------------------------------------------
# PER-PLAYER PIPELINE (unchanged core logic)
# --------------------------------------------------

grid_quantiles <- c(Q25 = 0.25, Q50 = 0.50, Q75 = 0.75)

run_player_pipeline <- function(player_name, team_id_val, full_data) {
  
  data <- full_data %>%
    filter(athlete_display_name == player_name) %>%
    arrange(game_date)
  
  if (nrow(data) < 20) {
    stop("Not enough games (", nrow(data), ") for a stable walk-forward fit")
  }
  
  stat_cols <- data %>%
    dplyr::select(points, assists, rebounds,
                  three_point_field_goals_made, player_possessions)
  
  vars <- colnames(stat_cols)
  poss_idx <- which(vars == "player_possessions")
  order <- c(poss_idx, setdiff(seq_along(vars), poss_idx))
  
  d <- length(vars)
  npairs <- d * (d - 1) / 2
  
  skeleton <- C2RVine(
    order,
    family = rep(0, npairs),
    par = rep(0, npairs),
    par2 = rep(0, npairs)
  )
  
  player_results <- map_dfr(11:nrow(data), function(t) {
    
    train <- data[1:(t - 1), ]
    actual <- data[t, ]
    
    u <- pobs(as.matrix(
      train %>%
        dplyr::select(points, assists, rebounds,
                      three_point_field_goals_made,
                      player_possessions)
    ))
    
    vine_fit <- RVineCopSelect(
      u,
      familyset = NA,
      Matrix = skeleton$Matrix,
      indeptest = TRUE
    )
    
    sim <- RVineSim(10000, vine_fit)
    
    pts_sim <- quantile(train$points, sim[, 2], type = 7)
    ast_sim <- quantile(train$assists, sim[, 3], type = 7)
    reb_sim <- quantile(train$rebounds, sim[, 1], type = 7)
    
    grid <- expand.grid(
      pts_level = names(grid_quantiles),
      reb_level = names(grid_quantiles),
      ast_level = names(grid_quantiles),
      stringsAsFactors = FALSE
    )
    
    grid %>%
      rowwise() %>%
      mutate(
        pts_th = quantile(train$points, grid_quantiles[pts_level], type = 7, names = FALSE),
        reb_th = quantile(train$rebounds, grid_quantiles[reb_level], type = 7, names = FALSE),
        ast_th = quantile(train$assists, grid_quantiles[ast_level], type = 7, names = FALSE),
        
        game_id = actual$game_id,
        game_date = actual$game_date,
        n_training_games = t - 1,
        
        p_vine = mean(
          pts_sim >= pts_th &
            reb_sim >= reb_th &
            ast_sim >= ast_th
        ),
        
        p_empirical = mean(
          train$points >= pts_th &
            train$rebounds >= reb_th &
            train$assists >= ast_th
        ),
        
        outcome = as.integer(
          actual$points >= pts_th &
            actual$rebounds >= reb_th &
            actual$assists >= ast_th
        )
      ) %>%
      ungroup()
  })
  
  player_results %>%
    mutate(
      athlete_display_name = player_name,
      team_id = team_id_val
    )
}

# --------------------------------------------------
# RUN PIPELINE PER TEAM, FALLING BACK THROUGH CANDIDATES
# UNTIL ONE SUCCEEDS (or the team's list is exhausted)
# --------------------------------------------------

run_team_with_fallback <- function(team_candidates) {
  
  team_id_val <- team_candidates$team_id[1]
  names_to_try <- team_candidates$athlete_display_name
  
  for (nm in names_to_try) {
    message("Running: ", nm)
    
    result <- tryCatch({
      run_player_pipeline(nm, team_id_val, full_data)
    }, error = function(e) {
      message("  Skipped ", nm, " - ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(result)) {
      return(result)  # success - stop trying more candidates for this team
    }
  }
  
  message("  ALL candidates failed for team_id ", team_id_val, " - team excluded")
  NULL
}

all_results <- map_dfr(candidates_by_team, run_team_with_fallback, .progress = TRUE)

# Sanity check: how many teams actually made it in
n_distinct(all_results$team_id)

# --------------------------------------------------
# LOG LOSS
# --------------------------------------------------

all_results <- all_results %>%
  mutate(
    p_vine = pmin(pmax(p_vine, 1e-6), 1 - 1e-6),
    p_empirical = pmin(pmax(p_empirical, 1e-6), 1 - 1e-6),
    
    log_loss_vine = -(
      outcome * log(p_vine) + (1 - outcome) * log(1 - p_vine)
    ),
    log_loss_empirical = -(
      outcome * log(p_empirical) + (1 - outcome) * log(1 - p_empirical)
    ),
    
    log_loss_diff = log_loss_empirical - log_loss_vine  # positive = vine wins
  ) %>%
  mutate(
    pts_level = factor(pts_level, levels = c("Q25", "Q50", "Q75")),
    reb_level = factor(reb_level, levels = c("Q25", "Q50", "Q75")),
    ast_level = factor(ast_level, levels = c("Q25", "Q50", "Q75"))
  )

# --------------------------------------------------
# PER-PLAYER SUMMARY (one row per successfully sampled starter)
# --------------------------------------------------

player_summary <- all_results %>%
  group_by(team_id, athlete_display_name) %>%
  summarise(
    mean_ll_vine = mean(log_loss_vine),
    mean_ll_empirical = mean(log_loss_empirical),
    mean_diff = mean(log_loss_diff),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_diff)

print(player_summary, n = 30)

# --------------------------------------------------
# VARIANCE + MEAN-DIFFERENCE ANALYSIS
# --------------------------------------------------

paired_t <- t.test(player_summary$mean_ll_empirical,
                   player_summary$mean_ll_vine,
                   paired = TRUE)
print(paired_t)

wilcox_test <- wilcox.test(player_summary$mean_ll_empirical,
                           player_summary$mean_ll_vine,
                           paired = TRUE)
print(wilcox_test)

var_test <- var.test(player_summary$mean_ll_vine,
                     player_summary$mean_ll_empirical)
print(var_test)

player_summary %>%
  summarise(
    across(c(mean_ll_vine, mean_ll_empirical, mean_diff),
           list(mean = mean, sd = sd, var = var),
           .names = "{.col}_{.fn}")
  )

# --------------------------------------------------
# VISUALIZATION 1 — distribution of per-player mean log loss by model
# --------------------------------------------------

player_summary_long <- player_summary %>%
  pivot_longer(cols = c(mean_ll_vine, mean_ll_empirical),
               names_to = "model", values_to = "mean_log_loss") %>%
  mutate(model = recode(model,
                        mean_ll_vine = "Vine",
                        mean_ll_empirical = "Empirical"))

ggplot(player_summary_long, aes(x = model, y = mean_log_loss, fill = model)) +
  geom_boxplot(width = 0.4, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.6, size = 2) +
  labs(
    title = "Per-starter mean log loss: Vine vs Empirical",
    subtitle = paste0("n = ", nrow(player_summary), " sampled starters (one per team)"),
    x = NULL, y = "Mean log loss (across threshold grid)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# --------------------------------------------------
# VISUALIZATION 2 — paired difference per player (dumbbell-style)
# --------------------------------------------------

ggplot(player_summary, aes(x = reorder(athlete_display_name, mean_diff), y = mean_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = mean_diff > 0), size = 3) +
  geom_segment(aes(xend = athlete_display_name, y = 0, yend = mean_diff,
                   color = mean_diff > 0)) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "Vine wins", "FALSE" = "Empirical wins"),
                     name = NULL) +
  coord_flip() +
  labs(
    title = "Log loss difference by starter (Empirical - Vine)",
    subtitle = "Positive = vine model outperforms empirical baseline for that player",
    x = NULL, y = "Mean log loss difference"
  ) +
  theme_minimal()

# --------------------------------------------------
# VISUALIZATION 3 — full distribution of raw grid-cell log losses
# --------------------------------------------------

all_results_long <- all_results %>%
  pivot_longer(cols = c(log_loss_vine, log_loss_empirical),
               names_to = "model", values_to = "log_loss") %>%
  mutate(model = recode(model,
                        log_loss_vine = "Vine",
                        log_loss_empirical = "Empirical"))

ggplot(all_results_long, aes(x = log_loss, fill = model)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of log loss across all players x grid cells",
    x = "Log loss", y = "Density"
  ) +
  theme_minimal()
