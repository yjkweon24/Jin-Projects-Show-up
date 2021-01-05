# =========================================================================
# Title: Scrape raw html tables function
#
# Description:
# The function 'exportcsv' takes the argument of 'id' which is passed to grep
# function to recognize start line of roster table, 'filepath' to save file,
# 'basketref' for the link of NBA website for scraping, 'team_hrefs' to create
# right URL for each team, and 'team_names' to save csv file as each team's
# name. This function obtain the information of roster, salary, and stat table
# and return it as a csv file.
# Each table is read as a data.frame, which is then exported as a csv file
# to the corresponding subdirectory in the 'rawdata/' folder
# =========================================================================

library(XML, warn.conflicts = FALSE)
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")


  
# The function 'exportcsv' takes the argument of 'id' which is passed to grep
# function to recognize start line of roster table, 'filepath' to save file,
# 'basketref' for the link of NBA website for scraping, 'team_hrefs' to create
# right URL for each team, and 'team_names' to save csv file as each team's
# name. This function obtain the information of roster, salary, and stat table
# and return it as a csv file.

exportcsv <- function(id, filepath, basketref, team_hrefs, team_names) {
  for(i in 1:30) {
    # Read html document (as a character vector) for a given team
    # (first team is "CLE")
    url_team <- paste0(basketref, team_hrefs[i])
    html_doc <- readLines(con = url_team)
    # initial line position of input file(such as roster) html table
    begin_roster <- grep(id, html_doc)
    # find the line where the html ends
    line_counter <- begin_roster
    while (!grepl("</table>", html_doc[line_counter])) {
      line_counter <- line_counter + 1
    }
    # read input argument file's(such as roster) table as data.frame
    # and export it as csv
    roster <- readHTMLTable(html_doc[begin_roster:line_counter])
    write.csv(roster, file = paste0('data/rawdata/', filepath,
                                    team_names[i], '.csv'))
  }
}