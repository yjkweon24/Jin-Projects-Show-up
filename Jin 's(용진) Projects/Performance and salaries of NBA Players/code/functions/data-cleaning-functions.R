#================================================
# Title : Data cleaning Function Script
#
# Description : This script is a collection of functions which is sourced by
# 'clean-data-script' to help data cleaning of NBA players' records
#================================================


# Function 'readfiles' takes the arguments of path of the file and each
# NBA team's name to read the each team's raw data file. It returns the file
# which has the name as filepath+team_names.csv.
readfiles <- function(filepath, team_names) {
  return(read.csv(file = paste0(filepath, team_names, '.csv'),
                  stringsAsFactors = FALSE, header = TRUE ))
}

# function 'orderfiles' takes the arguments of file and index which is the
# column to order the file. It returns the ordered file by the arugment index. 
orderfiles <- function(rawfile, index) {
  return(rawfile[order(rawfile[, index]),])
}

# function 'jointable' takes the arugments of two files and returns the
# joinned table which contains intersection of both files by NBA players' name.
jointable <- function(file1, file2) {
  a <- inner_join(file1, file2, 
                  by = c("roster.Player" = "salaries.Player"))
  return (a)
}

# function 'jointable2' takes the arugments of two files and returns the
# joinned table which contains intersection of both files by NBA players' name.
jointable2 <- function(file1, file2) {
  a <- inner_join(file1, file2, 
                  by = c("roster.Player" = "totals.Player"))
  return (a)
}
