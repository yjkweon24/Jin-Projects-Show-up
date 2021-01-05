# =========================================================================
# Title: Scrape raw html tables
#
# Description:
# This script contains R code to help you scrape the tables 
# 'Roster', 'Totals', and 'Salaries', for a specific NBA team.
# Each table is read as a data.frame, which is then exported as a csv file
# to the corresponding subdirectory in the 'rawdata/' folder
#
# Imported Functions : Imported Functions in this script is from 'download-
# data-fruntions.R' script which is 'exportcsv.'
#
# The function 'exportcsv' takes the argument of 'id' which is passed to grep
# function to recognize start line of roster table, 'filepath' to save file,
# 'basketref' for the link of NBA website for scraping, 'team_hrefs' to create
# right URL for each team, and 'team_names' to save csv file as each team's
# name. This function obtain the information of roster, salary, and stat table
# and return it as a csv file.
# =========================================================================

library(XML, warn.conflicts = FALSE)
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")
source("code/functions/download-data-functions.R")
source("code/scripts/download-data-script2.R")


# Call the function in download-data-functions.R to export each 
# team's statistic as a csv file.
exportcsv(id = 'id="roster"', filepath = 'roster-data/roster-',  
          basketref, team_hrefs, team_names)

exportcsv(id = 'id="totals"', filepath = 'stat-data/stats-',
          basketref, team_hrefs, team_names)

exportcsv(id = 'id="salaries"', filepath = 'salary-data/salaries-',
          basketref, team_hrefs, team_names)

