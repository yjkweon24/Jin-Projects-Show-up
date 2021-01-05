#================================================================
# Title: Creating Data Dictionary for the data we used in this project
#
# Description: This script generate data dictionary for each data and generate
# it as csv files. The data dictionary contains exact variable name as in 
# the data file, a longer explanation of what the variable means, the measure-
# ment units, and etc. 
#================================================================

library(readr)
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")


name <-c(0)
description <- c(0)
explanataion <- c(0)

# make a data dictionary for eff-stats-salary.csv

cd1 <- data.frame(
  name = c("Total Points", "Total Rebounds", "Assists", "Steals", "Blocks",
           "Missed Field Goals", "Missed Free Throws", 
          "Turnovers", "Games Played", "Efficiency Index", "Salary",
          "Position"),
  description = c("total points number of each player",
                  "total rebounds number for each player", 
                  "number of assists by each player",
                  "number of steals by each player",
                  "number of blocks by each player",
                  "number of missed field goals by each player",
                  "number of missed free throw by each player", 
                  "number of turnovers by each player",
                  "games played by each player", 
                  "efficiency index for each player" , 
                  "salary for each player", 
                  "position for each player"),
  explanataion = c(" how many total points that each player gets",
                   "a statistic awarded to a player who retrieves
                   the ball after a missed field goal or free throw",
                   "attributed to a player who passes the ball to a
                   teammate in a way that leads to a score by field goal",
                   "a steal occurs when a defensive player legally
                   causes a turnover by his positive, aggressive action", 
                   "a block occurs when a defensive player legally 
                   deflects a field goal attempt from an offensive player",
                   "number of missed basket scored on any shot or
                   tap other than a free throw", 
                   "missed unopposed attempts to score points from
                   a restricted area on the court the free throw line",
                   "a turnover occurs when a team loses possession of 
                   the ball to the opposing team before a player takes
                   a shot at his team's basket", 
                   "total number of games that each player has played before",
                   "efficiency index we calculated by Eff Model",
                   "salary earned by each player", "position for each player")
)

write.csv(cd1, file = "data/cleandata/dictionary-cleaned-eff-stats-salary.csv")

# make a data dictionary for roster-salary-stats.csv
cd2 <- data.frame(
  name = c("Field Goal Percentage", "Three Point Field Goal Percentage", 
           "Two Point Field Goal Percentage", 
           "Effective Field Goal Percentage", "Free Throws Percentage"), 
  description = c("the percentage of field goal",
                  "the percentage of three point field goal",
                  "the percentage of two point field goal",
                  "adjusts field goal percentage to 
                  account for the fact that three-point field goals
                  count for three points while field goals only
                  count for two points", "the percentage of free throws"), 
  explanation = c("Field goal divided by Field Goal Attempt", 
                  "Three Point Goal divided by Three Point Goal Attempt",
                  "Two Point Field Goal divided by
                  Two Point Field Goal Attempt",
                  "Field goal + 0.5 times Three Point Field
                  Goals and then divided by Field Goal Attempt",
                  "Free Throws divided by Free Throws Attempt")
)

write.csv(cd2, file = "data/cleandata/dictionary-cleaned-roster-salary-stats.csv")


# make a data dictionary for team-salaries.csv
cd3 <- data.frame(
  name = c("first_quartile_salary", "median_salary", "third_quartile_salary", 
           "average_salary", "interquartile_range"),
  description = c("the first quartile of salary", "the median of salary",
                  "the third quartile of salary", 
                  "the average of salary",
                  "the interquartile range of salary"),
  explanation = c("the 25% of salaries for players in each team",
                  "the 50% of salaries for players in each team",
                  "the 75% of salaries for players in each team", 
                  "for each team, add salary of each players and 
                  divided by number of players in that team",
                  "the third quartile - the first quartile for the salaries
                  in each team")
)
  
write.csv(cd3, file = "data/cleandata/dictionary-cleaned-team-salaries.csv")
  
# make a data dictionary for roster data
cd4 <- data.frame(
  name = c("roster.No.", "roster.Player", "roster.Pos", "roster.Ht",
           "roster.Wt", "roster.Birth.Date", "roster.Country", "roster.Exp",
           "roster.College"),
  description = c("roster number of each player in this team", "name", 
                  "position", "height", 
                  "the weight of each player in this team", "birth date", 
                  "country", 
                  "years of experience", 
                  "college name"),
  explanation = c("the number where this player is at in the original
                  roster table", 
                  "the name of each player in this team", 
                  "the specific position for each player in this team, 
                  there are five types of positions, which are C, PF, SF,
                  SG, and PG", 
                  "the height of each player in the team in ft",
                  "the weight of each player in the team in pounds", 
                  "the birth date of each player in the team", 
                  "the country that each player comes from in this team", 
                  "the experience of playing basketball games for each player
                  in this team, unit is year",
                  "the name of college that each player in this team was
                  graduated from")
)
write.csv(cd4, file = "data/cleandata/dictionary-raw-roster.csv")
  
# make a data dictionary for salary-data
cd5 <- data.frame(
  name = c("salaries.Rk", "salaries.Player", "salaries.Salary"),
  description = c("the ranking of salary for each player in this team",
                  "the name of player corresponding to this salary
                  in this team",
                 "the salary of this player in this team"),
  explanation = c("this number shows the ranking for this specific player
                  in this team",
                 "this is the name of the player who earns this salary in
                 this team", 
                 "this is the actual salary the player earns in this team,
                 unit is U.S. dollars")
)
write.csv(cd5, file = "data/cleandata/dictionary-raw-salary.csv")

# make a data dictionary for stat-data
cd6 <- data.frame(
  name = c("totals.Rk", "totals.Player", "totals.Age", "totals.G",
           "totals.GS", "totals.MP", "totals.FG", 
           "totals.FGA", "totals.3P", "totals.3PA", "totals.3P.",
           "totals.2P", "totals.2PA", "totals.2P.",
           "totals.eFG.", "totals.FT", "totals.FTA", "totals.FT.",
           "totals.ORB", "totals.DRB", "totals.TRB", 
           "totals.AST", "totals.STL", "totals.BLK", "totals.TOV",
           "totals.PF", "totals.PTS"), 
  description = c("ranking", "name", "age", "games", "games started",
                  "minutes played", 
                  "field goals", "field goals attempted",
                  "three-point field goals",
                  "three-point field goal attempted", 
                  "three-point field goals percentage", 
                  "two-point field goals", 
                  "two-point field goals attempted", 
                  "two-point field goal percentage",
                  "effective field goals percentage", "free throws",
                  "free throws attempted", 
                  "free throws percentage", "offensive rebound",
                  "defensive rebound", "total rebound", 
                  "assists", "steals", "blocks", "turnovers",
                  "personal fouls", "points"),
  explanation = c("the ranking of the player in this team", 
                  "the name of the player in this team", 
                  "the age of the player in this team", 
                  "total number of games played by this player",
                  "total number of games started of this player",  
                  "the total minutes played by this player", 
                  "total number of field goals of this player", 
                  "total number of field goals attempted of this player", 
                  "total number of three-point field goal of this player", 
                  "total number of three-point field goal attempt of
                  this player",
                  "three-point field goals percentage of this player,
                  calculated from dividing three-point field goals by 
                  three-point field goals attempted", 
                  "total number of two-point field goal of this player", 
                  "total number of two-point field goal attempted
                  of this player", 
                  "two-point field goal percentage,
                  calculated from dividing two-point field goals by
                  two-point field goals attempted", 
                  "effective field goals percentage of this player,
                  formula is (FG + 0.5 * 3P) / FGA", 
                  "total number of free throws by this player", 
                  "total number of free throw attempted by this player", 
                  "free throw percentage of this player, calculated from 
                  dividing number of free throws by number of free throws
                  attempted", 
                  "total number of offensive rebound of this player", 
                  "totalnumber of defensive rebound of this player", 
                  "total number of total rebound of this player", 
                  "total number of assists to turnover of this player", 
                  "total number of steals made by this player", 
                  "total number of blocks made by this player", 
                  "total number of turnovers made by this player", 
                  "total number of personal fouls of this player", 
                  "total number of points that this player got")
)
write.csv(cd6, file = "data/cleandata/dictionary-raw-stat.csv")

  
