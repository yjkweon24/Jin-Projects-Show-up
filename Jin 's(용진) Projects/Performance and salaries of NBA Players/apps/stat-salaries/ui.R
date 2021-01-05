# =========================================================================
# Title: Visulization of the Relationship Between All the Player Statistics
# 
# Description: This shiny app shows scatterplot in which the x-axis corresponds
# to one statistics, and the y-axis corresponds to one statistics(salary
# by default) by two input widgets. It includes another widget to indicate
# wheter the dots should be colored by position or not. Lastly, the app 
# displays the correlation coefficient between two chosen variables in the text.
# =========================================================================

# Change the scientific notation (e) to the numbers
options(scipen = 5)
options(digits = 10)
library(shiny)

# Set the working directory
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")


# Read csv file
salary <- read.csv("data/cleandata/eff-stats-salary.csv")

# Define UI for application that draws a scatterplot
shinyUI(fluidPage(
  
# Application title
  titlePanel("Relationship Between All the Player Statistics"),
  
# Generate a row with a sidebar
  sidebarLayout(  
    
# Define the sidebar with three inputs
    sidebarPanel(
      selectInput(inputId = "xaxis", label = "Select X-axis:", 
                  choices = colnames(salary[3:13])),
      helpText("Select variable of data from eff-stats-salary.csv
                which contains salary statistics for each team"),
      
      selectInput(inputId = "yaxis", label = "Select Y-axis:",
                  choices = colnames(salary[3:13]),
                  selected = colnames(salary[13])),
      helpText("Select variable of data from eff-stats-salary.csv
                which contains salary statistics for each team"),
      
      selectInput(inputId = "color", label = "Select Color by Position:",
                  choices = list("Color by Position" = "color", 
                                 "No Color" = "no")),
      helpText("Select 'Color by Position' if you want to color each player's
               by color")
      ),
    
    # Create a spot for the scatterplot
    mainPanel(
      plotOutput("salaryPlot", height = "700px"),  
      verbatimTextOutput("cor"),
      helpText("Correlation Coefficient Between two chosen variables")
      )
    )
  )
)
