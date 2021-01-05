## Stat133 Final Project Fall 2016

#### Authors: Jerome Kim, Zimei Yuan, Jiyoon Jeong, Jin Kweon

### **Description: In this project we analyze data about basketball players from the NBA League
in the 2015-2016 season. The central topic has to do with the salary and performance of 
NBA players. We focus on the question that *"In the 2015-2016 season, how do the skills
of a player relate to his salary?"* **

#### Important instructions: 

```
   1. Cleaning and scraping data 2. EDA analysis 3. Using shiny app for "team-salaries" file  
   4. EFf model with PCA(Principle Components Analysis) 5. Using shiny app for "stat-salaries" file 
   6. Calculate value of players 7. Conclusion
   
```
-----

### The directory-files structure
```
  project 
    README.md
    code/
      functions/
      scripts/
        download-data-script.R
        clean-data-script.R
        eda-script.R
        compute-efficiency-script.R
    data/
      rawdata/
        roster-data/
        salary-data/
        stat-data/
      cleandata/
        roster-salary-stats.csv
        eff-stats-salary.csv
        eda-output.txt
    images/
    apps/
        team-salaries/
          ui.R
          server.R
        stat-salaries/
          ui.R
          server.R
      report/
        report.Rmd
        report.pdf
      slides/
        slides.Rmd
        slides.html
  
```
