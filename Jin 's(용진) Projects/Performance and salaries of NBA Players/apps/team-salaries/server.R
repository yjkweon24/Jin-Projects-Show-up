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

# Define server logic required to draw a bar-chart
shinyServer(function(input, output) {
  
# Fill in the spot we created for a plot
  output$salaryPlot <- renderPlot({
    
# if input$order is Ascending, sort data 'salary' by ascending order
# if input$order is Descending, sort data 'salary' by descending order
# otherwise, don't sort
    if(input$order == "a") {
            salary = salary[order(salary[, input$standard]),]
          }
    else if(input$order == "d"){
      salary = salary[order(salary[, input$standard], decreasing = TRUE),]
      }
    
# Render a bar-chart
    barplot(salary[, input$standard],
            main = input$standard,
            xlab = "Statistics of Selected Standard",
            ylab = "Teams", col=c(rainbow(30)),
            xlim = c(-1, max(salary[, input$standard])), horiz = TRUE,
            axes = TRUE, names.arg= as.character(salary[1:30, 1]), cex.names=1,
            las = 1)
  })
})
