# =========================================================================
# Title: Visulization of the Relationship Between All the Player Statistics
# 
# Description: This shiny app shows scatterplot in which the x-axis corresponds
# to one statistics, and the y-axis corresponds to one statistics(salary
# by default) by two input widgets. It includes another widget to indicate
# wheter the dots should be colored by position or not. Lastly, the app 
# displays the correlation coefficient between two chosen variables in
# the text.
# =========================================================================

# Change scientific notation (e) to numbers. 
options(scipen = 5)
options(digits = 10)

library(shiny) 
setwd("E:/Cal/2016 Fall/Stat 133/Final/project")
salary <- read.csv("data/cleandata/eff-stats-salary.csv")

# Define server logic required to draw a scatterplot
shinyServer(function(input, output) {
  
# Fill in the spot we created for a Correlattion coefficient
output$cor <- renderText({
     paste("Correlation Coefficient Between two chosen variable is",
           cor(salary[,input$xaxis], salary[, input$yaxis]))
   })
  
# Fill in the spot we created for a plot
output$salaryPlot <- renderPlot({
    if(input$color == "color") {
       color = salary$Position
     }
    else {
      color = "black"
    }
# Render a scatterplot
    plot(salary[, c(input$xaxis, input$yaxis)], col = color,
            main = paste(input$xaxis, " vs ", input$yaxis),
            xlab= input$xaxis,
            ylab= input$yaxis)
    legend(max(salary[,input$xaxis]) / 90, max(salary[, input$yaxis]),
           unique(salary$Position),
           col = 1:length(color),
           pch = 1)
  })
})
