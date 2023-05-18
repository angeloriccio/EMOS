
library(INLA)
library(fields)
library(abind)
library(lattice)
library(MASS)

##--- Source the function for a specific day
source(paste("./","run_single_instance.R",sep=""))
 
## For which day do you want the prediction?
i_day <- 4

## where input data are?
data.path <- "./"

##  what date (dd-mmm-yyyy format)
Dates.lista <- "04-Jan-2022"

## what pollutant?
poll <- "PM10"

sink(file=paste(data.path, "console", poll, ".txt",sep=""))
run_single_instance(Dates.lista, i_day, poll)
sink()
