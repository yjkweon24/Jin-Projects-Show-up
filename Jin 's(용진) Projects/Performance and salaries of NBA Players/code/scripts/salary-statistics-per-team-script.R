# =========================================================================
# Title: Salary Statistics per team
# 
# Description: This script analyze salaries aggregated by teams and obtain
# a table 'team-salaries.csv' which contains salary's statistics for each
# team such as total payroll, minimum salary, maximum salary, first quartile
# salary, median salary, third quartile salary, average salary, interquartile
# range, standard deviation.
# =========================================================================

# Change the scientific notation (e) to the number
options(scipen = 5)
options(digits = 10)

setwd("E:/Cal/2016 Fall/Stat 133/Final/project")

data <- read.csv("data/cleandata/roster-salary-stats.csv")

# Set vectors for making a data frame
total_payroll <- c(0)
minimum_salary <- c(0)
maximum_salary <- c(0)
first_quartile_salary <- c(0)
median_salary <- c(0)
third_quartile_salary <- c(0)
average_salary <- c(0)
interquartile_range <- c(0)
standard_deviation <- c(0)

# make team dataframe to save statistics
team <- data.frame(total_payroll, minimum_salary, maximum_salary,
                   first_quartile_salary, median_salary, third_quartile_salary,
                   average_salary, interquartile_range, standard_deviation)


# Assign 1 for the newly made variables
head = 1
tail = 1
count = 1

# analyze each team's salary statistics and save it into dataframe 'team'
for(i in 1:(length(data[,40]) - 1)) {
  
  # if team's name is changed, analyze the data from 'head' row to 'i' row
  # or if i is 470 which is the second last row of 'data', then let i = 471
  # and analyze data from 'head' row to 471 row. 
  if(data[i,40] != data[i + 1,40] | i == 470) {
    if(i == 470){
      i = i + 1
    }
    team[count,1] = sum(data[head:i, 12], na.rm = TRUE)
    team[count,2] = min(data[head:i, 12], na.rm = TRUE)
    team[count,3] = max(data[head:i, 12], na.rm = TRUE)
    team[count,4] = quantile(data[head:i, 12], c(0.25), na.rm = TRUE)
    team[count,5] = quantile(data[head:i, 12], c(0.5), na.rm = TRUE)
    team[count,6] = quantile(data[head:i, 12], c(0.75), na.rm = TRUE)
    team[count,7] = mean(data[head:i, 12], na.rm = TRUE)
    team[count,8] = IQR(data[head:i, 12], na.rm = TRUE)
    team[count,9] = sd(data[head:i, 12], na.rm = TRUE)
    head = i+1
    count = count + 1
  }
}


# Set rownames equals to each team's name
rownames(team) <- c("CLE", "TOR", "MIA", "ATL", "BOS", "CHO", "IND", "DET",
                    "CHI", "WAS", "ORL", "MIL", "NYK", "BRK", "PHI", "GSW",
                    "SAS", "OKC", "LAC", "POR", "DAL", "MEM", "HOU", "UTA",
                    "SAC", "DEN", "NOP", "MIN", "PHO", "LAL")

# Export it as a csv file 
write.csv(team, "data/cleandata/team-salaries.csv")
