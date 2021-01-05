# =========================================================================
# Title: Scrape raw html tables
#
# Description:
# This script contains R code to help you scrape the tables 
# 'Roster', 'Totals', and 'Salaries', for a specific NBA team by their URL.
# This script will be soucred by 'download-data-script.R' to provide information
#  of each team.
# =========================================================================


library(XML, warn.conflicts = FALSE)
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")


# The desired html tables are in URL's having this form:
# "http://www.basketball-reference.com/teams/CLE/2016.html"
# "http://www.basketball-reference.com/teams/TOR/2016.html"
# "http://www.basketball-reference.com/teams/MIA/2016.html"
#
# The first step is to extract the part of the url associate with each
# team, that is: /teams/CLE/2016.html, /teams/TOR/2016.html, ...
#
# To do that, we'll scrape the the page:
# "http://www.basketball-reference.com/leagues/NBA_2016.html"
# in order to extract the 'href' attributes of the anchor tags:
# /teams/CLE/2016.html
# /teams/TOR/2016.html
# /teams/MIA/2016.html
#
# These attributes are extracted as a character vector
# (these will be used to parse each team's page)

# Extract_team_info 
# Base url
basketref <- 'http://www.basketball-reference.com'

# Parse 'http://www.basketball-reference.com/leagues/NBA_2016.html'
url <- paste0(basketref, '/leagues/NBA_2016.html')
doc <- htmlParse(url)

# Identify nodes with anchor tags for each team and
# extract the href attribute from the anchor tags
team_rows <- getNodeSet(doc, "//th[@scope='row']/a")
team_hrefs <- xmlSApply(team_rows, xmlAttrs)

# Just in case, here's the character vector with the team abbreviations
team_names <- substr(team_hrefs, 8, 10)
