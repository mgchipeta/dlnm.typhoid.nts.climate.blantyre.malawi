#STRATAA time series analysis

#Required packages
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

#Load climatic dataset and data manipulation
#===========================================
climate <- read.csv("/Users/dthindwa/Documents/Rprogramming/Strataa_timeseries/climate.csv")

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
climate <- subset(climate, select = c(climate_date, rainfall, temperature, humidity ))

#Load typhoid/iNTS dataset and data manipulation
#===============================================
case <- read.csv("/Users/dthindwa/Documents/Rprogramming/Strataa_timeseries/case.csv")
case$case_date <- dmy(case$date_s)
case$sex <-case$gender
case$age <- case$age_yrs
case$organism_type <- case$org_type
case$organism <- case$org
case$case_count <- c(1)

#Final case dataset to use
case <- subset(case, select = c(case_date, age, sex, organism_type, organism, case_count))

#Time series manipulation for typhi
#==================================
#select cases of typhi only
case_typhi <- subset(case, organism == "typhi", select = c(case_date, age, sex, organism_type, organism, case_count))

#select counts of typhi.
case_typhi.x <- subset(case_typhi, select = c(case_date, case_count))

#convert from data frame to xts
case_typhi.x = as.xts(case_typhi.x[,-1,drop = FALSE], order.by = as.Date(case_typhi.x[,1]))

#generate daily, quartery and yearly cases
case_typhi.x.d <- apply.daily(case_typhi.x, FUN = sum)
case_typhi.x.q <- apply.quarterly(case_typhi.x, FUN = sum)
case_typhi.x.y <- apply.yearly(case_typhi.x, FUN = sum)

#Time series manipulation for iNTS
#=================================
#select cases of iNTS only
case_iNTS <- subset(case, organism == "iNTS", select = c(case_date, age, sex, organism_type, organism, case_count))

# select counts of iNTS.
case_iNTS.x <- subset(case_iNTS, select = c(case_date, case_count))

#convert from data frame to xts
case_iNTS.x = as.xts(case_iNTS.x[,-1,drop = FALSE], order.by = as.Date(case_iNTS.x[,1]))

#generate daily, quartery and yearly cases
case_iNTS.x.d <- apply.daily(case_iNTS.x, FUN = sum)
case_iNTS.x.q <- apply.quarterly(case_iNTS.x, FUN = sum)
case_iNTS.x.y <- apply.yearly(case_iNTS.x, FUN = sum)

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

#generate daily, quartery and yearly rainfall, humidity, temperature values
climate_rain.x.d <- apply.daily(climate_rain.x, FUN = mean)
climate_rain.x.q <- apply.quarterly(climate_rain.x, FUN = mean)
climate_rain.x.y <- apply.yearly(climate_rain.x, FUN = mean)

climate_humi.x.d <- apply.daily(climate_humi.x, FUN = mean)
climate_humi.x.q <- apply.quarterly(climate_humi.x, FUN = mean)
climate_humi.x.y <- apply.yearly(climate_humi.x, FUN = mean)

climate_temp.x.d <- apply.daily(climate_temp.x, FUN = mean)
climate_temp.x.q <- apply.quarterly(climate_temp.x, FUN = mean)
climate_temp.x.y <- apply.yearly(climate_temp.x, FUN = mean)

#Converting xts objects back to data frames for ggplotting
#========================================================
case_typhi.x.df <- case_typhi.x.d %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

case_typhi.x.qf <- case_typhi.x.q %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

case_typhi.x.yf <- case_typhi.x.y %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

case_iNTS.x.df <- case_iNTS.x.d %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

case_iNTS.x.qf <- case_iNTS.x.q %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

case_iNTS.x.yf <- case_iNTS.x.y %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_rain.x.df <- climate_rain.x.d %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_rain.x.qf <- climate_rain.x.q %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_rain.x.yf <- climate_rain.x.y %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_humi.x.df <- climate_humi.x.d %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_humi.x.qf <- climate_humi.x.q %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_humi.x.yf <- climate_humi.x.y %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_temp.x.df <- climate_temp.x.d %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_temp.x.qf <- climate_temp.x.q %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

climate_temp.x.yf <- climate_temp.x.y %>% as_tibble(preserve_row_names = TRUE) %>% 
mutate(date = ymd(row.names)) %>% select(-row.names) %>% select(date, everything())

#Plotting using ggplot2
#======================
#daily, quarterly iNTS cases + rainfall
ggplot() + geom_bar(data = case_iNTS.x.df, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_rain.x.df, aes(x=date, y=rainfall/10), colour = "blue", alpha = 0.6) + scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Rainfall [mm]")) + ggtitle("Daily iNTS Cases & Daily Raifall (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_iNTS.x.qf, aes(x=date, y=case_count), stat = "identity", colour="red", size = 1.2) + geom_line(data = climate_rain.x.qf, aes(x=date, y=rainfall*30), colour = "blue", size = 1.2) + scale_y_continuous(sec.axis = sec_axis(~./30, name = "Rainfall [mm]")) + ggtitle("Quarterly iNTS Cases & Quarterly Average Raifall (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

#daily, quarterly iNTS cases + temperature
ggplot() + geom_bar(data = case_iNTS.x.df, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_temp.x.df, aes(x=date, y=temperature/3), colour = "blue", alpha = 0.2) + scale_y_continuous(sec.axis = sec_axis(~.*3, name = "Temperature [°C]")) + ggtitle("Daily iNTS Cases & Daily Temperature (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_iNTS.x.qf, aes(x=date, y=case_count), stat = "identity", colour="red", size = 1.2) + geom_line(data = climate_temp.x.qf, aes(x=date, y=temperature*10), colour = "blue", size = 1.2) + scale_y_continuous(sec.axis = sec_axis(~./10, name = "Temperature [°C]")) + ggtitle("Quarterly iNTS Cases & Quarterly Average Temperature (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

#daily, quarterly iNTS cases + humidity
ggplot() + geom_bar(data = case_iNTS.x.df, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_humi.x.df, aes(x=date, y=humidity/5), colour = "blue", alpha = 0.3) + scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Humidity [%]")) + ggtitle("Daily iNTS Cases & Daily Humidity (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_iNTS.x.qf, aes(x=date, y=case_count), stat = "identity", colour="red", size = 1.2) + geom_line(data = climate_humi.x.qf, aes(x=date, y=humidity*3), colour = "blue", size = 1.2) + scale_y_continuous(sec.axis = sec_axis(~./3, name = "Humidity [%]")) + ggtitle("Quarterly iNTS Cases & Quarterly Average Humidity (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

#daily, quarterly typhoid cases + rainfall
ggplot() + geom_bar(data = case_typhi.x.df, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_rain.x.df, aes(x=date, y=rainfall/10), colour = "blue", alpha = 0.6) + scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Rainfall [mm]")) + ggtitle("Daily Typhoid Cases & Daily Raifall (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_typhi.x.qf, aes(x=date, y=case_count), stat = "identity", colour="red", size = 1.2) + geom_line(data = climate_rain.x.qf, aes(x=date, y=rainfall*5), colour = "blue", size = 1.2) + scale_y_continuous(sec.axis = sec_axis(~./5, name = "Rainfall [mm]")) + ggtitle("Quarterly Typhoid Cases & Quarterly Average Raifall (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

#daily, quarterly typhoid cases + temperature
ggplot() + geom_bar(data = case_typhi.x.df, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_line(data = climate_temp.x.df, aes(x=date, y=temperature/15), colour = "blue", alpha = 0.5) + scale_y_continuous(sec.axis = sec_axis(~.*15, name = "Temperature [°C]")) + ggtitle("Daily Typhoid Cases & Daily Temperature (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_typhi.x.qf, aes(x=date, y=case_count), stat = "identity", colour="red", size = 1.2) + geom_line(data = climate_temp.x.qf, aes(x=date, y=temperature*10), colour = "blue", size = 1.2) + scale_y_continuous(sec.axis = sec_axis(~./10, name = "Temperature [°C]")) + ggtitle("Quarterly Typhoid Cases & Quarterly Average Temperature (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

#daily, quarterly tyhpoid cases + humidity
ggplot() + geom_bar(data = case_typhi.x.df, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_humi.x.df, aes(x=date, y=humidity/15), colour = "blue", alpha = 0.3) + scale_y_continuous(sec.axis = sec_axis(~.*15, name = "Humidity [%]")) + ggtitle("Daily Typhoid Cases & Daily Humidity (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_typhi.x.qf, aes(x=date, y=case_count), stat = "identity", colour="red", size = 1.2) + geom_line(data = climate_humi.x.qf, aes(x=date, y=humidity*3), colour = "blue", size = 1.2) + scale_y_continuous(sec.axis = sec_axis(~./3, name = "Humidity [%]")) + ggtitle("Quarterly Typhoid Cases & Quarterly Average Humidity (2000-2015)") + ylab('Cases') + xlab('Time') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

#Zoom-in daily plots for 2000-03, 2013-14
#create datasets for 2002-03, 2013-2014
case_iNTS.x.df.zoom1 <- subset(case_iNTS.x.df, year(date)>=2002 & year(date)<=2003)
case_iNTS.x.df.zoom2 <- subset(case_iNTS.x.df, year(date)>=2013 & year(date)<=2014)

case_typhi.x.df.zoom1 <- subset(case_typhi.x.df, year(date)>=2002 & year(date)<=2003)
case_typhi.x.df.zoom2 <- subset(case_typhi.x.df, year(date)>=2013 & year(date)<=2014)

climate_humi.x.df.zoom1 <- subset(climate_humi.x.df, year(date)>=2002 & year(date)<=2003)
climate_humi.x.df.zoom2 <- subset(climate_humi.x.df, year(date)>=2013 & year(date)<=2014)

climate_rain.x.df.zoom1 <- subset(climate_rain.x.df, year(date)>=2002 & year(date)<=2003)
climate_rain.x.df.zoom2 <- subset(climate_rain.x.df, year(date)>=2013 & year(date)<=2014)

climate_temp.x.df.zoom1 <- subset(climate_temp.x.df, year(date)>=2002 & year(date)<=2003)
climate_temp.x.df.zoom2 <- subset(climate_temp.x.df, year(date)>=2013 & year(date)<=2014)

#ggplotplot iNTS + climate zoomed-in 
ggplot() + geom_line(data = case_iNTS.x.df.zoom1, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_rain.x.df.zoom1, aes(x=date, y=rainfall/10), colour = "blue", alpha = 0.6) + scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Rainfall [mm]")) + ggtitle("Daily iNTS Cases & Daily Raifall (2002-2003)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_iNTS.x.df.zoom2, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_rain.x.df.zoom2, aes(x=date, y=rainfall/10), colour = "blue", alpha = 0.6) + scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Rainfall [mm]")) + ggtitle("Daily iNTS Cases & Daily Raifall (2013-2014)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

ggplot() + geom_line(data = case_iNTS.x.df.zoom1, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_humi.x.df.zoom1, aes(x=date, y=humidity/10), colour = "blue", alpha = 0.3) + scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Humidity [%]")) + ggtitle("Daily iNTS Cases & Daily Humidity (2002-2003)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_iNTS.x.df.zoom2, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_humi.x.df.zoom2, aes(x=date, y=humidity/30), colour = "blue", alpha = 0.3) + scale_y_continuous(sec.axis = sec_axis(~.*30, name = "Humidity [%]")) + ggtitle("Daily iNTS Cases & Daily Humidity (2013-2014)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

ggplot() + geom_line(data = case_iNTS.x.df.zoom1, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_temp.x.df.zoom1, aes(x=date, y=temperature/3), colour = "blue", alpha = 0.2) + scale_y_continuous(sec.axis = sec_axis(~.*3, name = "Temperature [°C]")) + ggtitle("Daily iNTS Cases & Daily Temperature (2002-2003)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_iNTS.x.df.zoom2, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_temp.x.df.zoom2, aes(x=date, y=temperature/10), colour = "blue", alpha = 0.2) + scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Temperature [°C]")) + ggtitle("Daily iNTS Cases & Daily Temperature (2013-2014)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

#ggplotplot typhi + climate zoomed-in
ggplot() + geom_bar(data = case_typhi.x.df.zoom1, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_rain.x.df.zoom1, aes(x=date, y=rainfall/10), colour = "blue", alpha = 0.6) + scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Rainfall [mm]")) + ggtitle("Daily Typhi Cases & Daily Raifall (2002-2003)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_typhi.x.df.zoom2, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_rain.x.df.zoom2, aes(x=date, y=rainfall/10), colour = "blue", alpha = 0.6) + scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Rainfall [mm]")) + ggtitle("Daily Typhi Cases & Daily Raifall (2013-2014)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

ggplot() + geom_bar(data = case_typhi.x.df.zoom1, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_humi.x.df.zoom1, aes(x=date, y=humidity/40), colour = "blue", alpha = 0.3) + scale_y_continuous(sec.axis = sec_axis(~.*40, name = "Humidity [%]")) + ggtitle("Daily Typhi Cases & Daily Humidity (2002-2003)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_typhi.x.df.zoom2, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_humi.x.df.zoom2, aes(x=date, y=humidity/15), colour = "blue", alpha = 0.3) + scale_y_continuous(sec.axis = sec_axis(~.*15, name = "Humidity [%]")) + ggtitle("Daily Typhi Cases & Daily Humidity (2013-2014)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

ggplot() + geom_bar(data = case_typhi.x.df.zoom1, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_temp.x.df.zoom1, aes(x=date, y=temperature/20), colour = "blue", alpha = 0.2) + scale_y_continuous(sec.axis = sec_axis(~.*20, name = "Temperature [°C]")) + ggtitle("Daily Typhi Cases & Daily Temperature (2002-2003)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))
ggplot() + geom_line(data = case_typhi.x.df.zoom2, aes(x=date, y=case_count), stat = "identity", colour="red") + geom_point(data = climate_temp.x.df.zoom2, aes(x=date, y=temperature/5), colour = "blue", alpha = 0.2) + scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Temperature [°C]")) + ggtitle("Daily Typhi Cases & Daily Temperature (2013-2014)") + ylab('Cases') + xlab('Time/day') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue", face = "bold"), axis.title.y = element_text(color="red", face = "bold"))

#END

# Age distribution plots
library(ggplot2)
library(dplyr)
library(tidyr)

# Density/distributions of typhi/iNTS by sex, age or both  
case$sex[case$sex == ""] <- NA
case$age[case$age == ""] <- NA

caseJoy <- case
caseJoy$date <- ymd(caseJoy$case_date)
caseJoy$case_year <- year(caseJoy$date)
caseJoy$agecat[caseJoy$age >=0 & caseJoy$age <=4] <- " 0-4"
caseJoy$agecat[caseJoy$age >4 & caseJoy$age <=9] <- " 5-9"
caseJoy$agecat[caseJoy$age >9 & caseJoy$age <=14] <- "10-14"
caseJoy$agecat[caseJoy$age >14 & caseJoy$age <=19] <- "15-19"
caseJoy$agecat[caseJoy$age >19 & caseJoy$age <=24] <- "20-24"
caseJoy$agecat[caseJoy$age >24 & caseJoy$age <=29] <- "25-29"
caseJoy$agecat[caseJoy$age >29 & caseJoy$age <=34] <- "30-34"
caseJoy$agecat[caseJoy$age >34 & caseJoy$age <=39] <- "35-39"
caseJoy$agecat[caseJoy$age >39 & caseJoy$age <=44] <- "40-44"
caseJoy$agecat[caseJoy$age >44 & caseJoy$age <=49] <- "45-49"
caseJoy$agecat[caseJoy$age >49 & caseJoy$age <=54] <- "50-54"
caseJoy$agecat[caseJoy$age >54 & caseJoy$age <=59] <- "55-59"
caseJoy$agecat[caseJoy$age >59 & caseJoy$age <=64] <- "60-64"
caseJoy$agecat[caseJoy$age >64 & caseJoy$age <=69] <- "65-69"
caseJoy$agecat[caseJoy$age >69 & caseJoy$age <=74] <- "70-74"
caseJoy$agecat[caseJoy$age >74] <- "75 +"

# Age frequency polygons stratified by year
ggplot(subset(caseJoy, organism=="iNTS" & !is.na(age)),aes(x=age)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 15)) + geom_freqpoly(aes(y=..ncount..),binwidth=5) + facet_wrap(~year,nrow=4,dir="v") + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age distribution of iNTS cases by year") + theme(plot.title = element_text(hjust = 0.5))
ggplot(subset(caseJoy, organism=="typhi" & !is.na(age)),aes(x=age)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 15)) + geom_freqpoly(aes(y=..ncount..),binwidth=5) + facet_wrap(~year,nrow=4,dir="v") + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age distribution of typhi cases by year") + theme(plot.title = element_text(hjust = 0.5))
ggplot(subset(caseJoy, !is.na(age)),aes(x=age, fill=organism, color=organism)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 15)) + geom_freqpoly(aes(y=..ncount..),binwidth=5) + facet_wrap(~year,nrow=4,dir="v") + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age distribution of iNTS & typhi cases by year") + theme(plot.title = element_text(hjust = 0.5))

#Age frequency polygons unstratified by year
ggplot(subset(caseJoy, organism=="iNTS" & !is.na(age)),aes(x=age)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 5)) + geom_freqpoly(aes(y=..ncount..),binwidth=5) + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age distribution of iNTS cases") + theme(plot.title = element_text(hjust = 0.5))
ggplot(subset(caseJoy, organism=="typhi" & !is.na(age)),aes(x=age)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 5)) + geom_freqpoly(aes(y=..ncount..),binwidth=5) + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age distribution of typhi cases") + theme(plot.title = element_text(hjust = 0.5))
ggplot(subset(caseJoy, !is.na(age)),aes(x=age, fill=organism, color=organism)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 5)) + geom_freqpoly(aes(y=..ncount..),binwidth=5) + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age distribution iNTS & typhi cases") + theme(plot.title = element_text(hjust = 0.5))
ggplot(subset(caseJoy, organism=="iNTS" & !is.na(age)),aes(x=age, fill=organism_type, color=organism_type)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 5)) + geom_freqpoly(aes(y=..ncount..),binwidth=5) + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age distribution iNTS species") + theme(plot.title = element_text(hjust = 0.5))

#Age/sex frequency polygons unstratified by year
ggplot(subset(caseJoy, organism=="iNTS" & !is.na(age) & sex != "Unknown"),aes(x=age, fill=sex, color=sex)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 5)) + geom_freqpoly(aes(y=..ncount..),binwidth=1) + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age/sex distribution iNTS cases 1998-2016") + theme(plot.title = element_text(hjust = 0.5))
ggplot(subset(caseJoy, organism=="typhi" & !is.na(age) & sex != "Unknown"),aes(x=age, fill=sex, color=sex)) + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 3)) + geom_freqpoly(aes(y=..ncount..),binwidth=1) + labs(x="Age (years)",y="Proportion of cases") + ggtitle("Age distribution typhi cases 1998-2016") + theme(plot.title = element_text(hjust = 0.5))
