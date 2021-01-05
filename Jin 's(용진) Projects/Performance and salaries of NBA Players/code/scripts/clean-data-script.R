# =========================================================================
# Title: Cleaning the Raw Data Files
# 
# Description:
# This script contains R code to help clean raw data files and save it as a 
# csv file to the corresponding subdirectory in the 'cleandata/' folder
#
# Imported Functions : 'readfiles', 'orderfiles', 'jointable', and 'jointable2'
# are the imported functions from 'data-cleaning-functions.R' script. Each
# function's description is explained in 'data-cleaning-functions.R' script.
# =========================================================================

# The following functions help change the scientific notation (e) to numbers
options(scipen = 5)
options(digits = 10)

# Warn.conflicts help us erase any error while use library functions.
library(dplyr, warn.conflicts = FALSE)
library(XML, warn.conflicts = FALSE)
library(readr)
library(stringr)
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")
source("code/scripts/download-data-script2.R")
source("code/functions/data-cleaning-functions.R")

# Set the vectors for the first row to have appropriate names.
a <- c(0)
b <- c(0)
cleaned <- c(0)

for(i in 1:30) {
  # Save each team's roster, totals, and salaries data into each variables
  roster_raw = readfiles('data/rawdata/roster-data/roster-', team_names[i])
  totals_raw = readfiles('data/rawdata/stat-data/stats-', team_names[i])
  salaries_raw = readfiles('data/rawdata/salary-data/salaries-', team_names[i])
  
  # Order each team's data by name
  order_roster = orderfiles(roster_raw, 3)
  order_salaries = orderfiles(salaries_raw, 3)
  order_totals = orderfiles(totals_raw, 3)
  
  # variable 'a' get the result of combined roster file's data with salary
  # file's data by the name of the players. Only select and combine the 
  # player's information who are listed on both files.
  a <- jointable(order_roster, order_salaries)
  
  # variable 'b' get the result of combined the data in 'a' with total file's
  # data by the name of the players. Only select and combine the 
  # player's information who are listed on both files.
  b <- jointable2(a, order_totals)
  
  # Add team's name at the last columns
  b$Team <- team_names[i]
  
  # If it examined the first team, then save data in variable 'b' to cleaned.
  # If it's not examined the first team, then combine the data in b with
  # existed cleaned data.
  if(i == 1) {
    cleaned <- b
  }
  else {
    cleaned <- rbind(cleaned, b)
  }
}

# Remove unnecessary columns
cleaned[c(1, 11, 14)] <- list(NULL)

# Change column's name into more descriptive name
colnames(cleaned)[1:38] <- c("Roster Number", "Player", "Position", "Height",
                             "Weight", "Birth Date", "Country", "Experience",
                             "College", "Salary Rank", "Salary", "Total Rank",
                             "Age", "Games" ,"Games Started", "Minutes Played",
                             "Field Goals", "Field Goal Attempts",
                             "Field Goal Percentage", 
                             "Three-Point Field Goals", 
                             "Three-Point Field Goal Attempts",
                             "Three-Point Field Goal Percentage",
                             "Two-Point Field Goals",
                             "Two-Point Field Goal Attempts",
                             "Two-Point Field Goal Percentage",
                             "Effective Field Goal Percentage", "Free Throws",
                             "Free Throw Attempts", "Free Throw Percentage",
                             "Offensive Rebounds", "Defensive Rebounds", 
                             "Total Rebounds", "Assists", "Steals", "Blocks",
                             "Turnovers", "Personal Fouls", "Points")

# Modify Salary column's data into numeric vectors and remove dollar sign
cleaned$Salary = as.numeric(gsub("\\$|,", "", cleaned$Salary))

# Modify Height column's data into numeric vectors and convert into inch
cleaned$Height = as.numeric(str_split_fixed(cleaned$Height,
                                            "-", n = 2)[, 1]) * 12 +
  as.numeric(str_split_fixed(cleaned$Height, "-", n = 2)[, 2])

# Remove dulicated players on the list
cleaned = cleaned[!duplicated(cleaned$Player), ]

# Export cleaned dataframe as a csv file
write.csv(cleaned, "data/cleandata/roster-salary-stats.csv")


