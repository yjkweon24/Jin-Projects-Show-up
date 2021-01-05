# =========================================================================
# Title: eda-script
#
# Description:
# For quantitative variables (e.g. salary, games played, free throws, etc) 
# graph histograms and boxplots
# For qualitative variables (e.g. position, team) graph barcharts
# of their frequencies
# =========================================================================

# To globally change the behavior of console session in favor of ten digits
options(digits = 12)
options(scipen = 5)

# Warn.conflicts option helps us remove any warning while running library.
library(readr)
library(ggplot2, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(grid, warn.conflicts = FALSE)

# For getting a frequency for qualitative variables.
library(MASS, warn.conflicts = FALSE) 

setwd("E:/Cal/2016 Fall/Stat 133/Final/project")

data <- read.csv("data/cleandata/roster-salary-stats.csv")

# Set vectors to make a data frame
minimum <- c(0)
median <- c(0)
ten_quartiles <- c(0)
twentyfive_quartiles <- c(0)
seventyfive_quartiles <- c(0)
ninty_quartiles <- c(0)
maximum <- c(0)
get_mean <- c(0)
Inter_quartiles <- c(0)
std_dev <- c(0)
ranges <- c(0)
difference <- c(0)
sums <- c(0)

# Create quantitative and qualitative vector columns
quant <- c(2, 5, 6, 9, 11:39)
quali <- c(3, 4, 7, 8, 10, 40)

# change Experience column's data 'R' into 0 and other data into
# numeric vectors

data[, quant[4]] = as.character(data[, quant[4]])

for(j in 1:length(data[, quant[4]])) {
  if(data[j, quant[4]]== 'R') {
    data[j, quant[4]] = 0
  }
}
data[, quant[4]] = as.numeric(data[, quant[4]])


# Get minimum for each column
for (i in 1:33) {
  minimum[i] = min(data[, quant[i]], na.rm = TRUE)
}


# Get maximum for each column
for (i in 1:33) {
  maximum[i] = max(data[, quant[i]], na.rm = TRUE)
}


# Get median for each column
for (i in 1:33) {
  median[i] = quantile(data[, quant[i]], c(0.5), na.rm = TRUE)
}

# Get mean for each column
for (i in 1:33) {
  get_mean[i] = mean(data[, quant[i]], na.rm = TRUE)
}

# Get quartiles for 10th
for (i in 1:33) {
  ten_quartiles[i] = quantile(data[, quant[i]], c(0.10), na.rm = TRUE)
}

# Get quartiles for 25th
for (i in 1:33) {
  twentyfive_quartiles[i] = quantile(data[, quant[i]], c(0.25), na.rm = TRUE)
}

# Get quartiles for 75th
for (i in 1:33) {
  seventyfive_quartiles[i] = quantile(data[, quant[i]], c(0.75), na.rm = TRUE)
}

# Get quartiles for 90th
for (i in 1:33) {
  ninty_quartiles[i] = quantile(data[, quant[i]], c(0.90), na.rm = TRUE)
}

# Get inter quartiles 
for (i in 1:33) {
  Inter_quartiles[i] = IQR(data[, quant[i]], na.rm = TRUE)
}

# Get standard deviation
for (i in 1:33) {
  std_dev[i] = sd(data[, quant[i]], na.rm = TRUE)
}


# Get ranges

for(i in 1:33) {
  assign(paste("range", quant[i], sep = ""),range(data[, quant[i]],
                                                na.rm = TRUE))
}


ranges <- c(range2, range5, range6, range9, range11, range12, range13,
            range14, range15, range16, range17, range18, range19, range20,
            range21, range22, range23, range24, range25, range26, range27,
            range28, range29, range30, range31, range32, range33, range34,
            range35, range36, range37, range38, range39)

ranges = as.integer(ranges)

# Get difference of ranges
for (i in 1:33) {
  difference[i] = maximum[i] - minimum[i]
}

# Get sum of each column
for (i in 1:33) {
  sums[i] = sum(as.numeric(data[, quant[i]]), na.rm = TRUE)
}

# Get summary of the data set at the end
summary(data)

# Get frequency distributions for qualitiative variables
for(i in 2:6) {
  assign(paste0("frequencies", quali[i]), table(data[,quali[i]]))
}

frequencies <- c(frequencies4, frequencies7, frequencies8, frequencies10, 
                 frequencies40)

# Get send the statistics to eda-output.txt
sink(file = "data/cleandata/eda-output.txt")

cat("STATISTICS FOR EACH VARAIABLE : Qualitative and Quantative variables
    in the order of columns")
cat("\n") # to give blank space between the sentences.
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Quantitative Data\n")
cat("Analysis: 1. Minimum  2. Maximum  3. Median  4. Mean  
    5. Quartiles (10th, 25th, 75th, 90th)  6. IQR  7. Standard Deviation
    8. Difference  9. Sum")
cat("\n")
cat("\n")
cat("Order: Roster number | Height | Weight | Experience | Salary Rank |
    Salary | Total Rank | Age | Games | Games Started | Minutes Played |
    Field Goals | Field Goal Attempts | Field Goal Percentage |
    Three-Point Field Goals | Three-Point Field Goal Attempts | 
    Three-Point Field Goal Percentage | Two-Point Field Goals | 
    Two-Point Field Goal Attempts | Two-Point Field Goal Percentage |
    Effective Field Goal Percentage | Free Throws | Free Throw Attempts |
    Free Throw Percentage | Offensive Rebounds | Defensive Rebounds |
    Total Rebounds | Assists | Steals | Blocks | Turnovers |
    Personal Fouls | Points")
cat("\n")
cat("\n")
cat("\n")

cat("Qualitative Data\n")
cat("Analysis: 1. Player  2. Position  3. Birthdate  4. Country 
    5. College  6. Team")
cat("\n")
cat("\n")
cat("Order: Position | Birthdate | Country | College | Team")
cat("\n")
cat("\n")

# Get structures of the data
cat("Structure\n")
cat("\n")
str(data)
cat("\n")
cat("\n")

# Get minimum for each column
cat("Minimum\n")
cat("\n")
minimum
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Maximum\n")
cat("\n")
maximum
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Median\n")
cat("\n")
median
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Ten quartiles\n")
cat("\n")
ten_quartiles
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Twentyfive quartiles\n")
cat("\n")
twentyfive_quartiles
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Seventyfive quartiles\n")
cat("\n")
seventyfive_quartiles
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Nintyfive quartiles\n")
cat("\n")
ninty_quartiles
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Mean")
cat("\n")
cat("\n")
get_mean 
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Interquartiles\n")
cat("\n")
Inter_quartiles 
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Standard Deviation")
cat("\n")
cat("\n")
std_dev 
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Ranges\n")
cat("\n")
ranges 
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Difference\n")
cat("\n")
difference 
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Sum\n")
cat("\n")
sums
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Frequency for qualitative variables\n")
cat("\n")
frequencies
cat("\n")
cat("\n")
cat("\n")
cat("\n")

cat("Summary\n")
cat("\n")
summary(data)
cat("\n")
cat("\n")
cat("\n")
cat("\n")

sink()

# Ggplot Begins here!
# I do not make graphs for players' names, birthdata, salary rank,
# and total rank since they do not contribute for exploratory data analysis.


# Change the direcotry
setwd("images")


# Create quantitative and qualitative vector columns
quant2 <- c(2, 5, 6, 9, 14:39)
quali2 <- c(4, 8, 10, 40)

# Get four bar charts for four qualitative variables
for (i in 1:4) {
  g <- ggplot(data, aes(x = data[, quali2[i]], fill = data[, quali2[i]])) +
    geom_bar(width = 0.5) + theme_grey() + theme(legend.position = "none") +
    coord_flip() +xlab(colnames(data)[quali2[i]])
  
  ggsave(sprintf("Bar-chart.%s.png", i), width = 10, height = 6)
}

# Draw separately, since binwidth is much higher than others.
salary = data$Salary
graphsalary <- ggplot(data, aes(x = salary)) +
  geom_histogram(binwidth = 500000) +
  theme_grey() + theme(legend.position = "none")

ggsave(graphsalary, file = "Histogram.0.png", width = 10, height = 6)

# I run ggplots based on the sizes of binwidths (for visual purpose)
# Combine with vectors where binwidth sizes are similar
a <- c(1:3, 12:13, 19:20, 22:30)

# Make histograms where binwidths are 20
for (i in a) {
  g <- ggplot(data, aes(x = data[, quant2[i]], fill = data[, quant2[i]])) +
    geom_histogram(binwidth = 20) + theme_grey() +
    theme(legend.position = "none") + 
    xlab(colnames(data)[quant2[i]])
  
  ggsave(sprintf("Histogram.%s.png", i), width = 10, height = 6)
}

# Make histograms where binwidths are 2
for (i in 4:7) {
  g <- ggplot(data, aes(x = data[, quant2[i]], fill = data[, quant2[i]])) +
    geom_histogram(binwidth = 2) + theme_grey() +
    theme(legend.position ="none") + xlab(colnames(data)[quant2[i]])
  ggsave(sprintf("Histogram.%s.png", i), width = 10, height = 6)
}

# Combine with vectors where binwidth sizes are similar
b <- c(8:10, 15:16)

# Make histograms where binwidths are 40
for (i in b) {
  g <- ggplot(data, aes(x = data[, quant2[i]], fill = data[, quant2[i]])) +
    geom_histogram(binwidth = 40) + theme_grey() +
    theme(legend.position="none") + xlab(colnames(data)[quant2[i]])
  
  ggsave(sprintf("Histogram.%s.png", i), width = 10, height = 6)
}

# Combine with vectors where binwidth sizes are similar
d <- c(11, 14, 17:18, 21)

# Make histograms where binwidths are 0.100
for (i in d) {
  
  # Remove Na values in each columns
  data_subset <- subset(data, !is.na(data[, quant2[i]]))
  
  g <- ggplot(data_subset, 
              aes(x = data_subset[, quant2[i]],
                  fill = data_subset[, quant2[i]])) +
    geom_histogram(binwidth = 0.100) + theme_grey() + 
    theme(legend.position = "none") + xlab(colnames(data)[quant2[i]])
  
  ggsave(sprintf("Histogram.%s.png", i), width = 10, height = 6)
}

# Need to make box plots below:
for (i in 1:30) {
  data_subset <- subset(data, !is.na(data[, quant2[i]]))
  
  g <- ggplot(data_subset, aes(x = data_subset[, quant2[i]],
                               y = data_subset[, quant2[i]], 
                               fill = data_subset[, quant2[i]])) +
    geom_boxplot(aes(group = cut_width(data_subset[,quant2[i]], 
                                       10000000))) +
    theme_grey() + theme(legend.position = "none") +
    xlab(colnames(data)[quant2[i]]) + ylab(colnames(data)[quant2[i]])
  
  ggsave(sprintf("Boxplots.%s.png", i), width = 10, height = 6)
}

g <- ggplot(data, aes(x = salary, y = salary, fill = salary)) +
  geom_boxplot(aes(group = cut_width(data_subset[, quant2[i]],
                                     100000000))) +
  theme_grey() + theme(legend.position = "none") 

ggsave(g, file = "Boxplots.0.png", width = 10, height = 6)