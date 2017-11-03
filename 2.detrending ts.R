#STRATAA time series analysis

# times series decomposition and detrending

rm(list=ls())
#load packages
library(lubridate)
library(dplyr)
library(plyr)
library(zoo)
library(xts)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(reshape2)
library(ggfortify)
library(data.table)
library(tidyquant)


#load climate
#===========================================
climate <- read.csv("/Users/dthindwa/Documents/PUBLICATION/Strataa/Time series/Data/Datasets in-use/climate.csv")

#creating mean values for rainfall
climate$rainfall <- (climate$chil_r + climate$chic_r)/2 

#creating mean values for temperature
climate$min_mean_t <- (climate$chil_mint + climate$chil_maxt)/2
climate$max_mean_t <- (climate$chic_mint + climate$chic_maxt)/2
climate$temperature <- (climate$min_mean_t + climate$max_mean_t)/2

#creating mean values for humidity
climate$humidity <- (climate$chil_h + climate$chic_h)/2

#converting to date using lubridate package
climate$climate_date <- dmy(climate$date)

#final climate dataset to use
climate <- subset(climate, select = c(climate_date, rainfall, temperature, humidity))


#load typhoid/iNTS cases
#=======================
case <- read.csv("/Users/dthindwa/Documents/PUBLICATION/Strataa/Time series/Data/Datasets in-use/case.csv")
case$case_date <- dmy(case$date_s)
case$sex <-case$gender
case$age <- case$age_yrs
case$organism_type <- case$org_type
case$organism <- case$org
case$case_count <- c(1)

#Final case dataset to use
case <- subset(case, select = c(case_date, age, sex, organism_type, organism, case_count))

#Time series manipulation for typhi & iNTS
#=========================================
#select cases of typhi only, or iNTS only
case_typhi <- subset(case, organism == "typhi", select = c(case_date, age, sex, organism_type, organism, case_count))
case_iNTS <- subset(case, organism == "iNTS", select = c(case_date, age, sex, organism_type, organism, case_count))

#select counts of typhi or iNTS.
case_typhi.x <- subset(case_typhi, select = c(case_date, case_count))
case_iNTS.x <- subset(case_iNTS, select = c(case_date, case_count))

#convert from data frame to xts for typhi or iNTS
case_typhi.x = as.xts(case_typhi.x[,-1,drop = FALSE], order.by = as.Date(case_typhi.x[,1]))
case_iNTS.x = as.xts(case_iNTS.x[,-1,drop = FALSE], order.by = as.Date(case_iNTS.x[,1]))

#generate monthly cases for typhi or iNTS
case_typhi.x.w <- apply.monthly(case_typhi.x, FUN = sum)
case_iNTS.x.w <- apply.monthly(case_iNTS.x, FUN = sum)


#Time series manipulation for rainfall, humidity, temperature
#============================================================
#select rainfall,  only
climate_rain <- subset(climate, select = c(climate_date, rainfall))
climate_humi <- subset(climate, select = c(climate_date, humidity))
climate_temp <- subset(climate, select = c(climate_date, temperature))

#convert from data frame to xts
climate_rain.x = as.xts(climate_rain[,-1,drop = FALSE], order.by = as.Date(climate_rain[,1]))
climate_humi.x = as.xts(climate_humi[,-1,drop = FALSE], order.by = as.Date(climate_humi[,1]))
climate_temp.x = as.xts(climate_temp[,-1,drop = FALSE], order.by = as.Date(climate_temp[,1]))

#generate monthly rainfall, humidity, temperature values
climate_rain.x.w <- apply.monthly(climate_rain.x, FUN = mean)
climate_humi.x.w <- apply.monthly(climate_humi.x, FUN = mean)
climate_temp.x.w <- apply.monthly(climate_temp.x, FUN = mean)

#generate & keep only time series-format datasets
typhi <-case_typhi.x.w 
iNTS <-case_iNTS.x.w 
rain <-climate_rain.x.w
humi <-climate_humi.x.w
temp <-climate_temp.x.w
rm(list = ls()[grep("^climate", ls())]) #deletes any dataset with word 'climate'
rm(list = ls()[grep("^case", ls())]) #deletes any dataset with word 'case'

#=====================================================================================
#using decompose function to detrend time series. Only if you have evenly (monthly) spaced observations in time series

#decompose iNTS
ts_iNTS = ts(iNTS$case_count, frequency = 12, start = c(2000,1), end = c(2015,12))
decompose_iNTS = decompose(ts_iNTS, "multiplicative")
par(mar=c(2,2,2,2)) #creates bigger plot space
plot(as.ts(decompose_iNTS$seasonal))
plot(as.ts(decompose_iNTS$trend))
plot(as.ts(decompose_iNTS$random))
plot(decompose_iNTS)
#recompose iNTS
recompose_iNTS = decompose_iNTS$seasonal*decompose_iNTS$trend*decompose_iNTS$random
par(mfrow = c(2,1), mar = c(2,2,2,2))
plot(ts_iNTS)
plot(recompose_iNTS)

#decompose rains
ts_rain = ts(rain$rainfall, frequency = 12, start = c(2000,1), end = c(2015,12))
decompose_rain = decompose(ts_rain, "additive")
par(mar=c(2,2,2,2)) #creates bigger plot space
plot(as.ts(decompose_rain$seasonal))
plot(as.ts(decompose_rain$trend))
plot(as.ts(decompose_rain$random))
plot(decompose_rain)
#recompose rains
recompose_rain = decompose_rain$seasonal+decompose_rain$trend+decompose_rain$random
par(mfrow = c(2,1), mar = c(2,2,2,2))
plot(ts_rain)
plot(recompose_rain)

#decompose temperature
ts_temp = ts(temp$temperature, frequency = 12, start = c(2000,1), end = c(2015,12))
decompose_temp = decompose(ts_temp, "additive")
par(mar=c(2,2,2,2)) #creates bigger plot space
plot(as.ts(decompose_temp$seasonal))
plot(as.ts(decompose_temp$trend))
plot(as.ts(decompose_temp$random))
plot(decompose_temp)
#recompose temperature
recompose_temp = decompose_temp$seasonal+decompose_temp$trend+decompose_temp$random
par(mfrow = c(2,1), mar = c(2,2,2,2))
plot(ts_temp)
plot(recompose_temp)

#decompose humidity
ts_humi = ts(humi$humidity, frequency = 12, start = c(2000,1), end = c(2015,12))
decompose_humi = decompose(ts_humi, "additive")
par(mar=c(2,2,2,2)) #creates bigger plot space
plot(as.ts(decompose_humi$seasonal))
plot(as.ts(decompose_humi$trend))
plot(as.ts(decompose_humi$random))
plot(decompose_humi)
#recompose humidity
recompose_humi = decompose_humi$seasonal+decompose_humi$trend+decompose_humi$random
par(mfrow = c(2,1), mar = c(2,2,2,2))
plot(ts_humi)
plot(recompose_humi)


#plot trends of iNTS & rainfall on same graph
par(mfrow = c(1,1), mar = c(5,4,4,6)) #set plot margins enough for captions
plot(decompose_iNTS$trend, axes=TRUE, xlab="", ylab="", main="Non-typhoid Salmonella & Rainfall trends", col="red3")
axis(2, col.axis="red")
mtext("Number of non-typhoid cases", side=2, line=2, col="red3")
box()
par(new = TRUE)
plot(decompose_rain$trend, axes=FALSE, xlab="", ylab="", col="blue")
axis(4, col="blue", col.axis="blue", las=1)
mtext("Rainfall [mm]", side=4, line=4, col="blue")
mtext("Time (year)", side=1, col="black", line=2.5)


#plot trends of iNTS & temperature on same graph
par(mfrow = c(1,1), mar = c(5,4,4,6)) #set plot margins enough for captions
plot(decompose_iNTS$trend, axes=TRUE, xlab="", ylab="", main="Non-typhoid Salmonella & Temperature trends", col="red3")
axis(2, col.axis="red")
mtext("Number of non-typhoid cases", side=2, line=2, col="red3")
box()
par(new = TRUE)
plot(decompose_temp$trend, axes=FALSE, xlab="", ylab="", col="blue")
axis(4, col="blue", col.axis="blue", las=1)
mtext("Temperature [°C]", side=4, line=4, col="blue")
mtext("Time (year)", side=1, col="black", line=2.5)


#plot trends of iNTS & humidity on same graph
par(mfrow = c(1,1), mar = c(5,4,4,6)) #set plot margins enough for captions
plot(decompose_iNTS$trend, axes=TRUE, xlab="", ylab="", main="Non-typhoid Salmonella & Humidity trends", col="red3")
axis(2, col.axis="red")
mtext("Number of non-typhoid cases", side=2, line=2, col="red3")
box()
par(new = TRUE)
plot(decompose_humi$trend, axes=FALSE, xlab="", ylab="", col="blue")
axis(4, col="blue", col.axis="blue", las=1)
mtext("Humidity [%]", side=4, line=4, col="blue")
mtext("Time (year)", side=1, col="black", line=2.5)

