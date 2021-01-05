# =========================================================================
# Title: Visualization of Different Statistics of Salary
# 
# Description: This shiny app shows the horizontal bar-chart using one widget
# to select the desired statistics, and another widget to indicate whether the
# bar should be displayed in decreasing order or in ascending.
# =========================================================================

# Change the scientific notation (e) to the numbers
options(scipen = 5)
options(digits = 10)

library(shiny)
library(ggplot2)

# Set the working directory
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")


# Read csv file
salary <- read.csv("data/cleandata/team-salaries.csv")

# Define UI for application that draws a bar-chart
shinyUI(fluidPage(
  
# Application title
  titlePanel("Team Salaries Summary by Different Standard"),
  
# Generate a row with a sidebar
  sidebarLayout(  
    
# Define the sidebar with two inputs
    sidebarPanel(
      selectInput(inputId = "standard", label = "Standard:", 
                  choices = colnames(salary[2:10])),
      helpText("Data from team-salaries.csv which contains salary statistics
               for each team"),
      
      selectInput(inputId = "order", label = "Select the type of order:",
                  choices = list("Ascending" = "a" , "Descending" = "d",
                                 "None" = "na")),
      
      helpText("If you want to order the data for each team, select these
               option")
      
      ),
    
# Create a spot for the bar-chart
    mainPanel(
      plotOutput("salaryPlot", height = "700px")  
      )
    )
  )
)

