# =========================================================================
# Title: Efficiency Index and PCA 
#
# Description: This code shows each player's efficiency that is obtained by
# a modified formula to balance the efficiency between offensive and 
# defensive positions by using PCA.
# It also shows all players' values(efficiency / salary)
# =========================================================================

# The functions below change scientific notation (e) to numbers.
options(scipen = 5)
options(digits = 10)

# Warn.conflicts option helps us remove any warning while learning library.
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)

# Set my directory
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")

# Read the database
cleandata <- read.csv("data/cleandata/roster-salary-stats.csv",
                      row.names = 1,
                      stringsAsFactors = FALSE)

# Efficiency
# Caclulate and put negative sign(MFG, MFT, TO). Add them into each subset
# subset by position
center <- cleandata %>%
  filter(Position == 'C') %>%
  mutate(MFT = Free.Throws - Free.Throw.Attempts) %>%
  mutate(MFG = Field.Goals - Field.Goal.Attempts) %>%
  mutate(TO = -1 * Turnovers)

power_forward <- cleandata %>%
  filter(Position == 'PF') %>%
  mutate(MFT = Free.Throws - Free.Throw.Attempts) %>%
  mutate(MFG = Field.Goals - Field.Goal.Attempts) %>%
  mutate(TO = -1 * Turnovers)

small_forward <- cleandata %>%
  filter(Position == 'SF') %>%
  mutate(MFT = Free.Throws - Free.Throw.Attempts) %>%
  mutate(MFG = Field.Goals - Field.Goal.Attempts) %>%
  mutate(TO = -1 * Turnovers)

shooting_guard <- cleandata %>%
  filter(Position == 'SG') %>%
  mutate(MFT = Free.Throws - Free.Throw.Attempts) %>%
  mutate(MFG = Field.Goals - Field.Goal.Attempts) %>%
  mutate(TO = -1 * Turnovers)

point_guard <- cleandata %>%
  filter(Position == 'PG') %>%
  mutate(MFT = Free.Throws - Free.Throw.Attempts) %>%
  mutate(MFG = Field.Goals - Field.Goal.Attempts) %>%
  mutate(TO = -1 * Turnovers)

# Statistics for efficiency
# Order of the datatable ('PTS', 'TRB', 'AST', 'STL', 'BLK',
# 'MFG', 'MFT', 'TO')
stats <- c('Points', 'Total.Rebounds', 'Assists', 'Steals',
           'Blocks', 'MFG', 'MFT', 'TO')

# X part of the EFF equation.
X_center <- as.matrix(center[, stats] / center$Games)
print(round(cor(X_center), 2), print.gap = 2)

X_power_forward <- as.matrix(power_forward[, stats] / power_forward$Games)
print(round(cor(X_power_forward), 2), print.gap = 2)

X_small_forward <- as.matrix(small_forward[, stats] / small_forward$Games)
print(round(cor(X_small_forward), 2), print.gap = 2)

X_shooting_guard <- as.matrix(shooting_guard[, stats] / shooting_guard$Games)
print(round(cor(X_shooting_guard), 2), print.gap = 2)

X_point_guard <- as.matrix(point_guard[, stats] / point_guard$Games)
print(round(cor(X_point_guard), 2), print.gap = 2)


# The following codes are weight part of the equation
# balacing the offense oriented formula by using prc
center_pca <- prcomp(X_center, center = TRUE, scale. = TRUE)
center_weights <- center_pca$rotation[, 1]
if (sum(sign(center_weights)) < 0) {
  center_weights = -1 * center_weights
}

power_forward_pca <- prcomp(X_power_forward, center = TRUE, scale. = TRUE)
power_forward_weights <- power_forward_pca$rotation[, 1]
if (sum(sign(power_forward_weights)) < 0) {
  power_forward_weights = -1 * power_forward_weights
}

small_forward_pca <- prcomp(X_small_forward, center = TRUE, scale. = TRUE)
small_forward_weights <- small_forward_pca$rotation[, 1]
if (sum(sign(small_forward_weights)) < 0) {
  small_forward_weights = -1 * small_forward_weights
}

shooting_guard_pca <- prcomp(X_shooting_guard, center = TRUE, scale. = TRUE)
shooting_guard_weights <- shooting_guard_pca$rotation[, 1]
if (sum(sign(shooting_guard_weights)) < 0) {
  shooting_guard_weights = -1 * shooting_guard_weights
}

point_guard_pca <- prcomp(X_point_guard, center = TRUE, scale. = TRUE)
point_guard_weights <- point_guard_pca$rotation[, 1]
if (sum(sign(point_guard_weights)) < 0) {
  point_guard_weights = -1 * point_guard_weights
}

# Get standard deviations
center_sigmas <- apply(X_center, 2, sd)
power_forward_sigmas <- apply(X_power_forward, 2, sd)
small_forward_sigmas <- apply(X_small_forward, 2, sd)
shooting_guard_sigmas <- apply(X_shooting_guard, 2, sd)
point_guard_sigmas <- apply(X_point_guard, 2, sd)

# Modified efficiency and integrate it to each subset
center_eff <- X_center %*% (center_weights / center_sigmas)
center$EFF <- center_eff

power_forward_eff <- X_power_forward %*% (power_forward_weights /
                                            power_forward_sigmas)
power_forward$EFF <- power_forward_eff

small_forward_eff <- X_small_forward %*% (small_forward_weights /
                                            small_forward_sigmas)
small_forward$EFF <- small_forward_eff

shooting_guard_eff <- X_shooting_guard %*% (shooting_guard_weights /
                                              shooting_guard_sigmas)
shooting_guard$EFF <- shooting_guard_eff

point_guard_eff <- X_point_guard %*% (point_guard_weights /
                                        point_guard_sigmas)
point_guard$EFF <- point_guard_eff


# Essence = player's name + total points + total rebounds ~~ salary
# change negative sign to positive
center_essence <- center[, c("Player", "Points", "Total.Rebounds", "Assists",
                             "Steals", "Blocks", "MFG", "MFT", "Turnovers",
                             "Games", "EFF", "Salary", "Position")]
center_essence$MFG <- with(center_essence, -1 * MFG)
center_essence$MFT <- with(center_essence, -1 * MFT)


power_forward_essence <- power_forward[c("Player", "Points", "Total.Rebounds",
                                         "Assists", "Steals","Blocks", "MFG",
                                         "MFT", "Turnovers", "Games",
                                         "EFF", "Salary", "Position")]
power_forward_essence$MFG <- with(power_forward_essence, -1 * MFG)
power_forward_essence$MFT <- with(power_forward_essence, -1 * MFT)

small_forward_essence <- small_forward[c("Player", "Points", "Total.Rebounds",
                                         "Assists", "Steals","Blocks", "MFG",
                                         "MFT", "Turnovers", "Games",
                                         "EFF", "Salary", "Position")]
small_forward_essence$MFG <- with(small_forward_essence, -1 * MFG)
small_forward_essence$MFT <- with(small_forward_essence, -1 * MFT)

shooting_guard_essence <- shooting_guard[c("Player", "Points", 
                                           "Total.Rebounds", "Assists",
                                           "Steals","Blocks", "MFG", "MFT",
                                           "Turnovers", "Games",
                                           "EFF", "Salary", "Position")]
shooting_guard_essence$MFG <- with(shooting_guard_essence, -1 * MFG)
shooting_guard_essence$MFT <- with(shooting_guard_essence, -1 * MFT)

point_guard_essence <- point_guard[c("Player", "Points", "Total.Rebounds",
                                     "Assists", "Steals","Blocks", "MFG",
                                     "MFT", "Turnovers", "Games",
                                     "EFF", "Salary", "Position")]
point_guard_essence$MFG <- with(point_guard_essence, -1 * MFG)
point_guard_essence$MFT <- with(point_guard_essence, -1 * MFT)


# Integrate information of each position into one datatable
all_eff <- rbind(center_essence, power_forward_essence, small_forward_essence,
                 shooting_guard_essence, point_guard_essence)

names(all_eff) <- c("Player's name", "Total Points", "Total Rebounds",
                    "Assists", "Steals", "Blocks","Missed Field goals",
                    "Missed Free Throws", "Turnovers",
                    "Games Played", "Efficiency Index", "Salary", "Position")

# Extract all_eff to csv file
write.csv(all_eff, "data/cleandata/eff-stats-salary.csv")


# Value of a player
center$value <- with(center, EFF / Salary)
power_forward$value <- with(power_forward, EFF / Salary)
small_forward$value <- with(small_forward, EFF / Salary)
shooting_guard$value <- with(shooting_guard, EFF / Salary)
point_guard$value <- with(point_guard, EFF / Salary)

# Merge all the players
value_all_players <- rbind(center, power_forward, small_forward,
                           shooting_guard, point_guard)

# Extract best 20 and worst 20 value
best_20_value <- value_all_players[with(value_all_players, order(-value)),]
best_20_value <- best_20_value[1:20,]
best_20_value$identify <- c("most valuable")
worst_20_value <- value_all_players[with(value_all_players, order(value)),]
worst_20_value <- worst_20_value[1:20,]
worst_20_value$identify <- c("worst value")


# Integrate information of each position to one datatable
total_40_value <- rbind(best_20_value, worst_20_value)

# Select player, position, value, identify columns
total_40 <- data.frame(total_40_value$Player, total_40_value$Position,
                       total_40_value$value, total_40_value$identify)

# Rename columns(unify a title style)
names(total_40) <- c("Player's name", "Position", "Value", "Identify")


# Extract all_eff to csv file
sink(file = "data/cleandata/best-worst-value-players.txt")
total_40
sink()

# Codes below are for getting the correlation between salaries
# and each skills.
sal = all_eff$"Salary"
pt1 = all_eff$"Total Points"
cor_point <- cor(sal, pt1)
cor_point

pt2 = all_eff$"Total Rebounds"
cor_rebounds <- cor(sal, pt2)
cor_rebounds

pt3 = all_eff$"Assists"
cor_assists <- cor(sal, pt3)
cor_assists

pt4 = all_eff$"Steals"
cor_steals <- cor(sal, pt4)
cor_steals

pt5 = all_eff$"Blocks"
cor_blocks <- cor(sal, pt5)
cor_blocks

pt6 = all_eff$"Missed Field goals"
cor_missed_field_goals <- cor(sal, pt6)
cor_missed_field_goals

pt7 = all_eff$"Missed Free Throws"
cor_missed_free_throws <- cor(sal, pt7)
cor_missed_free_throws

pt8 = all_eff$"Turnovers"
cor_turnovers <- cor(sal, pt8)
cor_turnovers

pt9 = all_eff$"Games Played"
cor_games <- cor(sal, pt9)
cor_games

# Make it to print out the name of the max and min of correlations
y <- c("cor_point", "cor_rebounds", "cor_assists", "cor_steals", "cor_blocks", 
       "cor_missed_field_goals", "cor_missed_free_throws", "cor_turnovers",
       "cor_games")

# Make it to see what skills are the most and least correlated
x <- c(cor_point, cor_rebounds, cor_assists, cor_steals, cor_blocks, 
       cor_missed_field_goals, cor_missed_free_throws, cor_turnovers,
       cor_games)

# Try to find the best correlated skills with the salaries
print(y[which.max(x)])

# Try to find the least correlated skills with the salaries
print(y[which.min(x)])

# Make regression for best correlated skill: points
z <- ggplot (all_eff, aes(x = pt1,y = sal)) +
  geom_point() + theme_bw()
z
plot (sal ~ pt1, all_eff, main = "regression salary and points",
      xlab = "points", ylab = "salary")
lin = lm (sal ~ pt1, all_eff)
lin
abline(lin, col = "red")

# Make regression for least correlated skill: blocks
z2 <- ggplot (all_eff, aes(x = pt5, y = sal)) +
  geom_point() + theme_bw()
z2
plot (sal ~ pt5, all_eff, main = "regression salary and blocks",
      xlab = "blocks", ylab = "salary")
lin2 = lm (sal ~ pt5, all_eff)
lin2
abline(lin2, col = "red")

# Position_essence contains positive values of all skills and
# other information.
sal_center = center_essence$"Salary"

# Calculate the correlation between salary and each skill 
# First, Do center.
pt1_center = center_essence$"Points"
cor_point_center <- cor(sal_center, pt1_center)

pt2_center = center_essence$"Total.Rebounds"
cor_rebounds_center <- cor(sal_center, pt2_center)

pt3_center = center_essence$"Assists"
cor_assists_center <- cor(sal_center, pt3_center)

pt4_center = center_essence$"Steals"
cor_steals_center <- cor(sal_center, pt4_center)

pt5_center = center_essence$"Blocks"
cor_blocks_center <- cor(sal_center, pt5_center)

pt6_center = center_essence$"MFG"
cor_missed_field_goals_center <- cor(sal_center, pt6_center)

pt7_center = center_essence$"MFT"
cor_missed_free_throws_center <-cor(sal_center, pt7_center)

pt8_center = center_essence$"Turnovers"
cor_turnovers_center <- cor(sal_center, pt8_center)

pt9_center = center_essence$"Games"
cor_games_center <- cor(sal_center, pt9_center)

# Do Power forward
sal_power_forward = power_forward_essence$"Salary"

pt1_power_forward = power_forward_essence$"Points"
cor_point_power_forward <- cor(sal_power_forward, pt1_power_forward)

pt2_power_forward = power_forward_essence$"Total.Rebounds"
cor_rebounds_power_forward <- cor(sal_power_forward, pt2_power_forward)

pt3_power_forward = power_forward_essence$"Assists"
cor_assists_power_forward <- cor(sal_power_forward, pt3_power_forward)

pt4_power_forward = power_forward_essence$"Steals"
cor_steals_power_forward <- cor(sal_power_forward, pt4_power_forward)

pt5_power_forward = power_forward_essence$"Blocks"
cor_blocks_power_forward <- cor(sal_power_forward, pt5_power_forward)

pt6_power_forward = power_forward_essence$"MFG"
cor_missed_field_goals_power_forward <- cor(sal_power_forward,
                                            pt6_power_forward)

pt7_power_forward = power_forward_essence$"MFT"
cor_missed_free_throws_power_forward <-cor(sal_power_forward,
                                           pt7_power_forward)

pt8_power_forward = power_forward_essence$"Turnovers"
cor_turnovers_power_forward <- cor(sal_power_forward,
                                   pt8_power_forward)

pt9_power_forward = power_forward_essence$"Games"
cor_games_power_forward <- cor(sal_power_forward,
                               pt9_power_forward)

# Do Small_forward
sal_small_forward = small_forward_essence$"Salary"

pt1_small_forward = small_forward_essence$"Points"
cor_point_small_forward <- cor(sal_small_forward, pt1_small_forward)

pt2_small_forward = small_forward_essence$"Total.Rebounds"
cor_rebounds_small_forward <- cor(sal_small_forward, pt2_small_forward)

pt3_small_forward = small_forward_essence$"Assists"
cor_assists_small_forward <- cor(sal_small_forward, pt3_small_forward)

pt4_small_forward = small_forward_essence$"Steals"
cor_steals_small_forward <- cor(sal_small_forward, pt4_small_forward)

pt5_small_forward = small_forward_essence$"Blocks"
cor_blocks_small_forward <- cor(sal_small_forward, pt5_small_forward)

pt6_small_forward = small_forward_essence$"MFG"
cor_missed_field_goals_small_forward <- cor(sal_small_forward,
                                            pt6_small_forward)

pt7_small_forward = small_forward_essence$"MFT"
cor_missed_free_throws_small_forward <-cor(sal_small_forward,
                                           pt7_small_forward)

pt8_small_forward = small_forward_essence$"Turnovers"
cor_turnovers_small_forward <- cor(sal_small_forward, pt8_small_forward)

pt9_small_forward = small_forward_essence$"Games"
cor_games_small_forward <- cor(sal_small_forward, pt9_small_forward)

# Do Shooting_guard
sal_shooting_guard = shooting_guard_essence$"Salary"

pt1_shooting_guard = shooting_guard_essence$"Points"
cor_point_shooting_guard <- cor(sal_shooting_guard, pt1_shooting_guard)

pt2_shooting_guard = shooting_guard_essence$"Total.Rebounds"
cor_rebounds_shooting_guard <- cor(sal_shooting_guard, pt2_shooting_guard)

pt3_shooting_guard = shooting_guard_essence$"Assists"
cor_assists_shooting_guard <- cor(sal_shooting_guard, pt3_shooting_guard)

pt4_shooting_guard = shooting_guard_essence$"Steals"
cor_steals_shooting_guard <- cor(sal_shooting_guard, pt4_shooting_guard)

pt5_shooting_guard = shooting_guard_essence$"Blocks"
cor_blocks_shooting_guard <- cor(sal_shooting_guard, pt5_shooting_guard)

pt6_shooting_guard = shooting_guard_essence$"MFG"
cor_missed_field_goals_shooting_guard <- cor(sal_shooting_guard,
                                             pt6_shooting_guard)

pt7_shooting_guard = shooting_guard_essence$"MFT"
cor_missed_free_throws_shooting_guard <-cor(sal_shooting_guard,
                                            pt7_shooting_guard)

pt8_shooting_guard = shooting_guard_essence$"Turnovers"
cor_turnovers_shooting_guard <- cor(sal_shooting_guard,
                                    pt8_shooting_guard)

pt9_shooting_guard = shooting_guard_essence$"Games"
cor_games_shooting_guard <- cor(sal_shooting_guard,
                                pt9_shooting_guard)

# Do Point_guard
sal_point_guard = point_guard_essence$"Salary"

pt1_point_guard = point_guard_essence$"Points"
cor_point_point_guard <- cor(sal_point_guard, pt1_point_guard)

pt2_point_guard = point_guard_essence$"Total.Rebounds"
cor_rebounds_point_guard <- cor(sal_point_guard, pt2_point_guard)

pt3_point_guard = point_guard_essence$"Assists"
cor_assists_point_guard <- cor(sal_point_guard, pt3_point_guard)

pt4_point_guard = point_guard_essence$"Steals"
cor_steals_point_guard <- cor(sal_point_guard, pt4_point_guard)

pt5_point_guard = point_guard_essence$"Blocks"
cor_blocks_point_guard <- cor(sal_point_guard, pt5_point_guard)

pt6_point_guard = point_guard_essence$"MFG"
cor_missed_field_goals_point_guard <- cor(sal_point_guard, pt6_point_guard)

pt7_point_guard = point_guard_essence$"MFT"
cor_missed_free_throws_point_guard <-cor(sal_point_guard, pt7_point_guard)

pt8_point_guard = point_guard_essence$"Turnovers"
cor_turnovers_point_guard <- cor(sal_point_guard, pt8_point_guard)

pt9_point_guard = point_guard_essence$"Games"
cor_games_point_guard <- cor(sal_point_guard, pt9_point_guard)

# Print the correlation of all positions
cor_point_center
cor_rebounds_center
cor_assists_center
cor_steals_center
cor_blocks_center
cor_missed_field_goals_center
cor_missed_free_throws_center
cor_turnovers_center
cor_games_center

cor_point_power_forward
cor_rebounds_power_forward
cor_assists_power_forward
cor_steals_power_forward
cor_blocks_power_forward
cor_missed_field_goals_power_forward
cor_missed_free_throws_power_forward
cor_turnovers_power_forward
cor_games_power_forward

cor_point_small_forward
cor_rebounds_small_forward
cor_assists_small_forward
cor_steals_small_forward
cor_blocks_small_forward
cor_missed_field_goals_small_forward
cor_missed_free_throws_small_forward
cor_turnovers_small_forward
cor_games_small_forward

cor_point_shooting_guard
cor_rebounds_shooting_guard
cor_assists_shooting_guard
cor_steals_shooting_guard
cor_blocks_shooting_guard
cor_missed_field_goals_shooting_guard
cor_missed_free_throws_shooting_guard
cor_turnovers_shooting_guard
cor_games_shooting_guard

cor_point_point_guard
cor_rebounds_point_guard
cor_assists_point_guard
cor_steals_point_guard
cor_blocks_point_guard
cor_missed_field_goals_point_guard
cor_missed_free_throws_point_guard
cor_turnovers_point_guard
cor_games_point_guard

# Make it to print out the name of the max and min of correlations
y_center <- c("cor_point_center", "cor_rebounds_center", "cor_assists_center",
              "cor_steals_center", "cor_blocks_center",
              "cor_missed_field_goals_center",
              "cor_missed_free_throws_center", "cor_turnovers_center",
              "cor_games_center")

y_power_forward <- c("cor_point_power_forward", "cor_rebounds_power_forward",
                     "cor_assists_power_forward", "cor_steals_power_forward",
                     "cor_blocks_power_forward",
                     "cor_missed_field_goals_power_forward",
                     "cor_missed_free_throws_power_forward",
                     "cor_turnovers_power_forward", "cor_games_power_forward")

y_small_forward <- c("cor_point_small_forward", "cor_rebounds_small_forward",
                     "cor_assists_small_forward", 
                     "cor_steals_small_forward",
                     "cor_blocks_small_forward", 
                     "cor_missed_field_goals_small_forward",
                     "cor_missed_free_throws_small_forward",
                     "cor_turnovers_small_forward", "cor_games_small_forward")

y_shooting_guard <- c("cor_point_shooting_guard", 
                      "cor_rebounds_shooting_guard",
                      "cor_assists_shooting_guard", 
                      "cor_steals_shooting_guard",
                      "cor_blocks_shooting_guard", 
                      "cor_missed_field_goals_shooting_guard",
                      "cor_missed_free_throws_shooting_guard",
                      "cor_turnovers_shooting_guard",
                      "cor_games_shooting_guard")

y_point_guard <- c("cor_point_point_guard", "cor_rebounds_point_guard",
                   "cor_assists_point_guard", "cor_steals_point_guard",
                   "cor_blocks_point_guard", 
                   "cor_missed_field_goals_point_guard",
                   "cor_missed_free_throws_point_guard",
                   "cor_turnovers_point_guard", "cor_games_point_guard")


# Make it to see what skills are the most and least correlated
x_center <- c(cor_point_center, cor_rebounds_center, cor_assists_center,
              cor_steals_center, cor_blocks_center,
              cor_missed_field_goals_center, cor_missed_free_throws_center,
              cor_turnovers_center, cor_games_center)

x_power_forward <- c(cor_point_power_forward, cor_rebounds_power_forward,
                     cor_assists_power_forward, cor_steals_power_forward,
                     cor_blocks_power_forward,
                     cor_missed_field_goals_power_forward,
                     cor_missed_free_throws_power_forward,
                     cor_turnovers_power_forward, cor_games_power_forward)

x_small_forward <- c(cor_point_small_forward, cor_rebounds_small_forward,
                     cor_assists_small_forward, cor_steals_small_forward,
                     cor_blocks_small_forward,
                     cor_missed_field_goals_small_forward,
                     cor_missed_free_throws_small_forward,
                     cor_turnovers_small_forward, cor_games_small_forward)

x_shooting_guard <- c(cor_point_shooting_guard, cor_rebounds_shooting_guard,
                      cor_assists_shooting_guard, cor_steals_shooting_guard,
                      cor_blocks_shooting_guard,
                      cor_missed_field_goals_shooting_guard,
                      cor_missed_free_throws_shooting_guard,
                      cor_turnovers_shooting_guard, cor_games_shooting_guard)

x_point_guard <- c(cor_point_point_guard, cor_rebounds_point_guard,
                   cor_assists_point_guard,
                   cor_steals_point_guard, cor_blocks_point_guard,
                   cor_missed_field_goals_point_guard,
                   cor_missed_free_throws_point_guard,
                   cor_turnovers_point_guard, cor_games_point_guard)


# Try to find the best and least correlated skills with the salaries 
# for center
print(y_center[which.max(x_center)])
print(y_center[which.min(x_center)])

# Try to find the best and least correlated skills with the salaries
# for point forward
print(y_power_forward[which.max(x_power_forward)])
print(y_power_forward[which.min(x_power_forward)])

# Try to find the best and least correlated skills with the salaries
# for small forward
print(y_small_forward[which.max(x_small_forward)])
print(y_small_forward[which.min(x_small_forward)])

# Try to find the best and least correlated skills with the salaries
# for shooting guard
print(y_shooting_guard[which.max(x_shooting_guard)])
print(y_shooting_guard[which.min(x_shooting_guard)])

# Try to find the best and least correlated skills with the
# salaries for point guard
print(y_point_guard[which.max(x_point_guard)])
print(y_point_guard[which.min(x_point_guard)])

# Make regression for worthy players
worthy_eff <- value_all_players$EFF
worthy_salary <- value_all_players$Salary

ggplot_worthy <- ggplot (value_all_players,
                         aes(x = worthy_salary, y = worthy_eff)) +
  geom_point() + theme_bw()

plot (worthy_eff ~ worthy_salary, value_all_players,
      main = "Regression Salary and Eff", xlab = "Salary",
      ylab = "Efficiency")
lin = lm(worthy_eff ~ worthy_salary, value_all_players)
lin
abline(lin, col = "red")

# Ordering players by the value of players
valued_players_descending <- value_all_players[with(value_all_players,
                                                    order(-value)),]
under_valued_30_players <- valued_players_descending[1:30,]

# Extracting neccessary columns from the data table
under_valued_players <- data.frame(under_valued_30_players$Player,
                                   under_valued_30_players$EFF,
                                   under_valued_30_players$Salary,
                                   under_valued_30_players$value)

# Ordering players by the value of players
valued_players_ascending <- value_all_players[with(value_all_players,
                                                   order(value)),]
over_valued_30_players <- valued_players_ascending[1:30,]

# Extracting neccessary columns from the data table
over_valued_players <- data.frame(over_valued_30_players$Player,
                                  over_valued_30_players$EFF,
                                  over_valued_30_players$Salary,
                                  over_valued_30_players$value)

# Printing undervalued and overvalued players
under_valued_players
over_valued_players