#02/08/2018
#by Deus
#STRATAA distributed lag nonlinear analysis & modelling

#===============================================LOAD PACKAGES=====================================================================
#load/attach all required packages for this analysis
library(tidyverse)
library(lubridate)
library(xts)
library(ggthemes, PerformanceAnalytics)
library(reshape2)
library(rugarch)
library(timetk, parallel)
library(timeSeries)
library(tseries)
library(data.table)
library(ggplot2)
library(dlnm)
library(broom)
library(caret)
library(gridExtra)
library(splines)
library(splines2)
library(pspline)
library(cowplot)
library(mgcv)
library(spi)
library(chron)
library(gridExtra)
library(grid)
library(pscl)
library(here)

#===============================================LOAD CASE DATA=====================================================================

#load typhoid and iNTS cases dataset.
case <- read.csv(here("Time.Series", "data", "case.csv"))
case$case_date <- dmy(case$date_s)
case$sex <-case$gender
case$age <- case$age_yrs
case$organism_type <- case$org_type
case$organism <- case$org
case$case_count <- c(1)
demog <- subset(case, select = c(case_date, age, sex, organism_type, organism, case_count))

#create separate datasets for typhi and iNTS cases.
case_typhi <- subset(case, organism == "typhi", select = c(case_date, age, sex, organism_type, organism, case_count))
case_iNTS <- subset(case, organism == "iNTS", select = c(case_date, age, sex, organism_type, organism, case_count))

#wrangle so that consecutive dates appear in the dataset of typhi. assign 0 to case_count when dates have no case.
case_typhi.x <- subset(case_typhi, select = c(case_date, case_count))
case_typhi.x <-aggregate(case_typhi.x$case_count, by=list(case_typhi.x$case_date), FUN=sum, na.rm=TRUE)
names(case_typhi.x) <- c("date", "case_count")
alldates <- data.table(date=seq.Date(min(case$case_date), max(case$case_date), by="day"))
case_typhi.x.w <- merge(case_typhi.x, alldates, by="date", all=TRUE)
case_typhi.x.w[is.na(case_typhi.x.w)] <- 0
dlnm.typhi <- case_typhi.x.w 

#wrangle so that consecutive dates appear in the dataset of iNTS. assign 0 to case_count when dates have no case.
case_iNTS.x <- subset(case_iNTS, select = c(case_date, case_count))
case_iNTS.x <-aggregate(case_iNTS.x$case_count, by=list(case_iNTS.x$case_date), FUN=sum, na.rm=TRUE)
names(case_iNTS.x) <- c("date", "case_count")
alldates <- data.table(date=seq.Date(min(case$case_date), max(case$case_date), by="day"))
case_iNTS.x.w <- merge(case_iNTS.x, alldates, by="date", all=TRUE)
case_iNTS.x.w[is.na(case_iNTS.x.w)] <- 0
dlnm.nts <- case_iNTS.x.w 

#convert typhi and iNTS data frames to xts objects for use in time series plotting.
case_typhi.x = as.xts(case_typhi.x.w[,-1,drop = FALSE], order.by = as.Date(case_typhi.x.w[,1]))
case_iNTS.x = as.xts(case_iNTS.x.w[,-1,drop = FALSE], order.by = as.Date(case_iNTS.x.w[,1]))

#weekly aggregates for typhi and iNTS cases.
wk.typhi <- apply.weekly(case_typhi.x, FUN = sum)
wk.iNTS <- apply.weekly(case_iNTS.x, FUN = sum)

#===============================================LOAD CLIMATE DATA=====================================================================

#load climate dataset.
climate <- read.csv(here("Time.Series", "data", "climate.csv"))

#calculate daily average values for temperature and rainfall.
climate$rainfall <- (climate$chil_r + climate$chic_r)/2 
climate$temperature <- ((climate$chil_mint + climate$chil_maxt)/2 + (climate$chic_mint + climate$chic_maxt)/2)/2

climate$climate_date <- dmy(climate$date)
climate <- subset(climate, select = c(climate_date, rainfall, temperature))
dlnm.climate <- climate

#create separate datasets for temperature, humidity and rainfall. these have initial daily values, no need to wrangle.
climate_rain <- subset(climate, select = c(climate_date, rainfall))
climate_temp <- subset(climate, select = c(climate_date, temperature))

#convert temperature, humidity and rainfall data frames to xts objects for use in time series plotting.
climate_rain.x = as.xts(climate_rain[,-1,drop = FALSE], order.by = as.Date(climate_rain[,1]))
climate_temp.x = as.xts(climate_temp[,-1,drop = FALSE], order.by = as.Date(climate_temp[,1]))

#weekly aggregates for rainfall and temperature
wk.rain <- apply.weekly(climate_rain.x, FUN = mean)
wk.temp <- apply.weekly(climate_temp.x, FUN = mean)

#===============================================DELETE ALL UNNECESSARY DATASETS=====================================================================

rm(list = ls()[grep("^climate", ls())]) 
rm(list = ls()[grep("^case", ls())])
rm(list = ls()[grep("^alldates", ls())]) 

#===============================================CREATE TIBBLES FOR WEEKLY GGPLOT2 PLOTTING====================================================================

wk.typhi <-tk_tbl(wk.typhi, preserve_index = TRUE, rename_index = "date") 
wk.typhi$incid <- wk.typhi$case_count*100000/688490
wk.iNTS <-tk_tbl(wk.iNTS, preserve_index = TRUE, rename_index = "date") 
wk.iNTS$incid <- wk.iNTS$case_count*100000/688490
wk.rain <-tk_tbl(wk.rain, preserve_index = TRUE, rename_index = "date") 
wk.temp <-tk_tbl(wk.temp, preserve_index = TRUE, rename_index = "date") 

#ggplot 1: iNTS v rainfall and temperature.
wk.p1 <- ggplot() + geom_line(data = wk.iNTS, aes(x=date, y=incid), stat = "identity", colour="tan4", size = 0.6) + geom_line(data = wk.rain, aes(x=date, y=rainfall/4), colour = "green4", alpha = 0.5, size = 0.6) + scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Pmean (mm)")) + ggtitle("(A)") + theme_classic() + ylab('NTS incidence rate (/100k pop)') + xlab('Year') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "green4", size = 11), axis.title.y = element_text(color="tan4", size = 11)) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 
wk.p2 <- ggplot() + geom_line(data = wk.iNTS, aes(x=date, y=incid), stat = "identity", colour="tan4", size = 0.6) + geom_line(data = wk.temp, aes(x=date, y=temperature/7), colour = "blue4", alpha = 0.5, size = 0.6) + scale_y_continuous(sec.axis = sec_axis(~.*7, name = "Tmean (°C)")) + ggtitle("(B)") + theme_classic() + ylab('NTS incidence rate (/100k pop)') + xlab('Year') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue4", size = 11), axis.title.y = element_text(color="tan4", size = 11)) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

#ggplot 2: typhi v rainfall and temperature.
wk.p3 <- ggplot() + geom_line(data = wk.typhi, aes(x=date, y=incid), stat = "identity", colour="red4", size = 0.6) + geom_line(data = wk.rain, aes(x=date, y=rainfall*0.3), colour = "green4", alpha = 0.5, size = 0.6) + scale_y_continuous(sec.axis = sec_axis(~./0.3, name = "Pmean (mm)")) + ggtitle("(C)") + theme_classic() + ylab('Typhoid incidence rate (/100k pop)') + xlab('Year') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "green4", size = 11), axis.title.y = element_text(color="red4", size = 11)) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 
wk.p4 <- ggplot() + geom_line(data = wk.typhi, aes(x=date, y=incid), stat = "identity", colour="red4", size = 0.6) + geom_line(data = wk.temp, aes(x=date, y=temperature*0.09), colour = "blue4", alpha = 0.5, size = 0.6) + scale_y_continuous(sec.axis = sec_axis(~./0.09, name = "Tmean (°C)")) + ggtitle("(D)") + theme_classic() + ylab('Typhoid incidence rate (/100k pop)') + xlab('Year') + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.title.y.right = element_text(color = "blue4", size = 11), axis.title.y = element_text(color="red4", size = 11)) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

#ggplot 3: combine all weekly time series plots.
plot_grid(wk.p1, wk.p2, wk.p3, wk.p4)

#ggplot 4: Boxplots for temperature, rainfall, NTS and typhoid
j <- seq(from=1, to=53, by=4)
k1 <-seq(from=16, to=40, by=4)
k2 <-seq(from=0, to=20, by=4)
k3 <-seq(from=0, to=7, by=1)
k4 <-seq(from=0, to=1.5, by=0.5)

pox1 <-ggplot(wk.temp, aes(x=week(date), y=temperature))  + 
  geom_boxplot(aes( group=week(date)), outlier.shape = NA, color="blue4", fill="blue4", alpha=0.2) + 
  labs(title="(A)",x="Week", y = "Tmean (°C)") + 
  scale_x_discrete(limits = j) + 
  scale_y_discrete(limits = k1) +
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) +
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

pox2 <- ggplot(wk.rain, aes(x=week(date), y=rainfall))  + 
  geom_boxplot(aes( group=week(date)), outlier.shape = NA, color="green4", fill="green4", alpha=0.2) + 
  labs(title="(B)",x="Week", y = "Pmean (mm/week)") + 
  scale_x_discrete(limits = j) + 
  scale_y_discrete(limits = k2) + 
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

pox3 <- ggplot(wk.iNTS, aes(x=week(date), y=incid))  + 
  geom_boxplot(aes( group=week(date)), outlier.shape = NA, color="tan4", fill="tan4", alpha=0.2) + 
  labs(title="(C)",x="Week", y = "NTS incidence rate (/100k pop)") + 
  scale_x_discrete(limits = j) + 
  scale_y_discrete(limits = k3) + 
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

pox4 <- ggplot(wk.typhi, aes(x=week(date), y=incid))  + 
  geom_boxplot(aes( group=week(date)), outlier.shape = NA, color="red4", fill="red", alpha=0.2) + 
  labs(title="(D)",x="Week", y = "Typhoid incidence rate (/100k pop)") + 
  scale_x_discrete(limits = j) + 
  scale_y_discrete(limits = k4) + 
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

#ggplot 5: combine all weekly boxplots.
plot_grid(pox1, pox2, pox3, pox4)

#===============================================DEMOGRAPHICS OF STUDY POPULATION GGPLOT2====================================================================

#distributions of typhi and iNTS cases by sex, age or both.
demog$sex[demog$sex == ""] <- NA
demog$age[demog$age == ""] <- NA
demog$date <- ymd(demog$case_date)
demog$year <- year(demog$date)

#categorize and label age into 5 year bands.
cutx <- function(x, lower = 0, upper, by = 5, sep = "-", above.char = "+") {
  labs <- c(paste(seq(lower, upper - by, by = by),
                  seq(lower + by - 1, upper - 1, by = by), 
                  sep = sep),
            paste(upper, above.char, sep = ""))
  
  cut(floor(x), breaks = c(seq(lower, upper, by = by), Inf),
      right = FALSE, labels = labs)
}
demog$agecat <- cutx(demog$age, upper = 90)

#age and sex frequency distribution
agesex.p1 <- ggplot(subset(demog, organism=="iNTS" & !is.na(age) & sex != "Unknown"), aes(x=age, fill=sex, color=sex)) + theme_classic() + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 5)) + geom_freqpoly(aes(y=..ncount..),binwidth=1, size = 1) + labs(x="Age (years)",y="Proportion of NTS cases") + ggtitle("(A)") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + theme(legend.title = element_blank(), legend.justification=c(0.5,0), legend.position = c(0.8, 0.5), legend.text = element_text(size = 16))
agesex.p2 <- ggplot(subset(demog, organism=="typhi" & !is.na(age) & sex != "Unknown"), aes(x=age, fill=sex, color=sex)) + theme_classic() + scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 5)) + geom_freqpoly(aes(y=..ncount..),binwidth=1, size = 1) + labs(x="Age (years)",y="Proportion of typhoid cases") + ggtitle("(B)") + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none') + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11))
plot_grid(agesex.p1, agesex.p2)

#============================================UNIVARIATE NATURAL SPLINE MODEL SELECTION===============================================================================
#prepare weekly dlnm dataset
wk.dlnm <- bind_cols(wk.typhi, wk.iNTS, wk.rain, wk.temp, id=NULL)
wk.dlnm$date1 <- wk.dlnm$date2 <- wk.dlnm$date3 <- NULL
wk.dlnm <- rename(wk.dlnm, typhi=case_count, NTS=case_count1, incidtyp = incid, incidNTS=incid1)
wk.dlnm$time <- seq.int(from = 1, to=836, by=1)
wk.dlnm$year <- year(wk.dlnm$date)
wk.dlnm$month <- month(wk.dlnm$date)
wk.dlnm$week <- week(wk.dlnm$date)

#Selecting best univariate model that fits weekly mean rainfall
wk.umod.rain1 <- gam(wk.dlnm$rainfall ~ ns(wk.dlnm$week, df=2))
wk.umod.rain2 <- gam(wk.dlnm$rainfall ~ ns(wk.dlnm$week, df=3))
wk.umod.rain3 <- gam(wk.dlnm$rainfall ~ ns(wk.dlnm$week, df=4))
wk.umod.rain4 <- gam(wk.dlnm$rainfall ~ ns(wk.dlnm$week, df=5))
wk.umod.rain5 <- gam(wk.dlnm$rainfall ~ ns(wk.dlnm$week, df=6))
wk.umod.rain6 <- gam(wk.dlnm$rainfall ~ ns(wk.dlnm$week, df=7))
wk.umod.rain7 <- gam(wk.dlnm$rainfall ~ ns(wk.dlnm$week, df=8))
wk.umod.rain8 <- gam(wk.dlnm$rainfall ~ ns(wk.dlnm$week, df=9))
AIC(wk.umod.rain1, wk.umod.rain2, wk.umod.rain3, wk.umod.rain4, wk.umod.rain5, wk.umod.rain6, wk.umod.rain7, wk.umod.rain8)
BIC(wk.umod.rain1, wk.umod.rain2, wk.umod.rain3, wk.umod.rain4, wk.umod.rain5, wk.umod.rain6, wk.umod.rain7, wk.umod.rain8)

#Selecting best univariate model that fits weekly mean temperature
wk.umod.temp1 <- gam(wk.dlnm$temperature ~ ns(wk.dlnm$week, df=2))
wk.umod.temp2 <- gam(wk.dlnm$temperature ~ ns(wk.dlnm$week, df=3))
wk.umod.temp3 <- gam(wk.dlnm$temperature ~ ns(wk.dlnm$week, df=4))
wk.umod.temp4 <- gam(wk.dlnm$temperature ~ ns(wk.dlnm$week, df=5))
wk.umod.temp5 <- gam(wk.dlnm$temperature ~ ns(wk.dlnm$week, df=6))
wk.umod.temp6 <- gam(wk.dlnm$temperature ~ ns(wk.dlnm$week, df=7))
wk.umod.temp7 <- gam(wk.dlnm$temperature ~ ns(wk.dlnm$week, df=8))
wk.umod.temp8 <- gam(wk.dlnm$temperature ~ ns(wk.dlnm$week, df=9))
AIC(wk.umod.temp1, wk.umod.temp2, wk.umod.temp3, wk.umod.temp4, wk.umod.temp5, wk.umod.temp6, wk.umod.temp7, wk.umod.temp8)
BIC(wk.umod.temp1, wk.umod.temp2, wk.umod.temp3, wk.umod.temp4, wk.umod.temp5, wk.umod.temp6, wk.umod.temp7, wk.umod.temp8)

#select optimal univariate natural spline function to describe weekly rainfall
wk.rain.spline.p1 <- ggplot(wk.dlnm, aes(week,rainfall)) + 
  geom_point(color="gray50") + 
  stat_smooth(method = gam, formula = y ~ ns(x, df=2), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 2")) +  
  stat_smooth(method = gam, formula = y ~ ns(x, df=3), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 3")) +   
  stat_smooth(method = gam, formula = y ~ ns(x, df=4), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 4")) +   
  stat_smooth(method = gam, formula = y ~ ns(x, df=5), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 5")) +   
  stat_smooth(method = gam, formula = y ~ ns(x, df=6), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 6")) +   
  stat_smooth(method = gam, formula = y ~ ns(x, df=7), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 7")) +   
  stat_smooth(method = gam, formula = y ~ ns(x, df=8), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 8")) +   
  stat_smooth(method = gam, formula = y ~ ns(x, df=9), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 9")) +   
  scale_colour_manual(name="", values=c("df = 2"="black", "df = 3"="green", "df = 4"="blue", "df = 5"="red", "df = 6"="orange", "df = 7"="khaki", "df = 8"="purple", "df = 9" = "yellow")) + 
  theme(legend.position = c(0.5,0.7), legend.direction = "vertical") + 
  labs(title="(A)",x="Week", y = "Pmean (mm)") + 
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

#select optimal univariate natural spline function to describe weekly temperature
wk.temp.spline.p2 <- ggplot(wk.dlnm, aes(week,temperature)) + 
  geom_point(color="gray50") + 
  stat_smooth(method = gam, formula = y ~ ns(x, df=2), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 2")) +  
  stat_smooth(method = gam, formula = y ~ ns(x, df=3), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 3")) +  
  stat_smooth(method = gam, formula = y ~ ns(x, df=4), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 4")) +  
  stat_smooth(method = gam, formula = y ~ ns(x, df=5), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 5")) +  
  stat_smooth(method = gam, formula = y ~ ns(x, df=6), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 6")) + 
  stat_smooth(method = gam, formula = y ~ ns(x, df=7), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 7")) + 
  stat_smooth(method = gam, formula = y ~ ns(x, df=8), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 8")) + 
  stat_smooth(method = gam, formula = y ~ ns(x, df=9), se =FALSE, alpha=0.4, size=0.6, aes(color="df = 9")) + 
  scale_colour_manual(name="", values=c("df = 2"="black", "df = 3"="green", "df = 4"="blue", "df = 5"="red", "df = 6"="orange", "df = 7" = "khaki", "df = 8" = "purple", "df = 9" = "yellow")) + 
  theme(legend.position = "none", legend.direction = "vertical") + 
  labs(title="(B)",x="Week", y = "Tmean (°C)") + 
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

#plot natural model selection of weekly rainfall and temperature
plot_grid(wk.rain.spline.p1, wk.temp.spline.p2)

#============================================NTS DLNM ANALYSIS===============================================================================

#cross basis for weekly rainfall, temperature 
varknots=equalknots(wk.dlnm$rainfall, fun = "ns", df=6)
lagknots <- logknots(4, 2)
wk.cb.rain.iNTS <- crossbasis(wk.dlnm$rainfall, lag =4, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))
summary(wk.cb.rain.iNTS)

varknots=equalknots(wk.dlnm$temperature, fun = "ns", df=7)
lagknots <- logknots(4, 2)
wk.cb.temp.iNTS <- crossbasis(wk.dlnm$temperature, lag =4, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))
summary(wk.cb.temp.iNTS)

#models fitting for weekly rainfall, temperature and NTS
wk.model.iNTS0 <- glm(incidNTS ~ ns(time), family = quasipoisson(), na.action=na.exclude,  wk.dlnm)
wk.model.iNTS1 <- glm(incidNTS ~ ns(time, 7*16), family = quasipoisson(), na.action=na.exclude,  wk.dlnm)
wk.model.iNTS2 <- glm(incidNTS ~ wk.cb.rain.iNTS + ns(time, 7*16), family = quasipoisson(), na.action=na.exclude, wk.dlnm)
wk.model.iNTS3 <- glm(incidNTS ~ wk.cb.temp.iNTS + ns(time, 7*16), family = quasipoisson(), na.action=na.exclude, wk.dlnm)
wk.model.iNTS4 <- glm(incidNTS ~ wk.cb.rain.iNTS + wk.cb.temp.iNTS + ns(time, 7*16), family = quasipoisson(), na.action=na.exclude, wk.dlnm)

#prediction for weekly rainfall, temperature
wk.pred.rain.iNTS2 <- crosspred(wk.cb.rain.iNTS, wk.model.iNTS2, cen = 2.71, by=0.2)
wk.pred.rain.iNTS4 <- crosspred(wk.cb.rain.iNTS, wk.model.iNTS4, cen = 2.71, by=0.2)

wk.pred.temp.iNTS3 <- crosspred(wk.cb.temp.iNTS, wk.model.iNTS3, cen = 23.25, by=0.2)
wk.pred.temp.iNTS4 <- crosspred(wk.cb.temp.iNTS, wk.model.iNTS4, cen = 23.25, by=0.2)

#3D plots for weekly rainfall, temperature
plot(wk.pred.rain.iNTS4, xlab="Pmean (mm)", zlab="Relative risk", theta=150, phi=5, lphi=60, main="(A)")
plot(wk.pred.temp.iNTS4, xlab="Tmean (°C)", zlab="Relative risk", theta=200, phi=35, lphi=40, main="(C)")

#countour plots for weekly rainfall, temperature
plot(wk.pred.rain.iNTS4, "contour", plot.title=title("(B)", xlab ="Pmean (mm)", ylab = "Lag (weeks)"))
plot(wk.pred.temp.iNTS4, "contour", key.title=title("RR"), plot.title=title("(D)", xlab ="Tmean (°C)", ylab = "Lag (weeks)"))

#detailed analysis targeting specific rainfall, temperature and lag values guided by contour plot
summary(wk.dlnm$rainfall)
plot(wk.pred.rain.iNTS4, "slices", var=0, ci="n", col=1, ylim=c(0.5,3), lwd=3.5, main="")
for(i in 1:2) lines(wk.pred.rain.iNTS4, "slices", var=c(4,30)[i], col=i+1, lwd=3.5) > legend("topright",paste("Precipitation =",c(0,4,30)), col=1:3, lwd=3.5)
plot(wk.pred.rain.iNTS4, "slices", var=c(0, 4, 30), lag=c(0,1,4), col=3, ci.arg=list(density=90,col=grey(0.2, alpha = 0.8)), ci.level=0.95, lwd=3.5, xlab="Precipitation(mm)", ylab="Relative risk (NTS)", cex.lab=1.3, cex.axis=1.3)
plot(wk.pred.rain.iNTS4, "overall", xlab="Rainfall (mm)", ylab="Relative risk (NTS)", ylim=c(0,4), xlim=c(0,30), ci.level=0.95, ci='l', lwd=3.5)

summary(wk.dlnm$temperature)
plot(wk.pred.temp.iNTS4, "slices", var=16, ci="n", col=1, ylim=c(0.8,3), lwd=3.5, main="")
for(i in 1:2) lines(wk.pred.temp.iNTS4, "slices", var=c(25,34)[i], col=i+1, lwd=3.5) > legend("topright",paste("Temperature =",c(16,25,34)), col=1:3, lwd=3.5)
plot(wk.pred.temp.iNTS4, "slices", var=c(16, 25, 34), lag=c(0,1,4), col=4, ci.arg=list(density=90,col=grey(0.2, alpha = 0.8)), ci.level=0.95, lwd=3.5, xlab="Temperature(°C)", ylab="Relative risk (NTS)", cex.lab=1.3, cex.axis=1.3)
plot(wk.pred.temp.iNTS4, "overall", xlab="Temperature (°C)", ylab="Relative risk (NTS)", ylim=c(0,4), xlim=c(16,30), ci.level=0.95, ci='l', lwd=3.5)

#plot observed v predicted values
ggplot(wk.dlnm, aes(time)) + 
  geom_point(aes(y=wk.dlnm$incidNTS, colour="Observed"), size=1) +
  geom_point(aes(y=predict(wk.model.iNTS4, type="response"), colour="Predicted"), size=1) + 
  labs(title="", x="Week", y = "Obs or Pred NTS incidence rate") + 
  scale_colour_manual(name="", values=c("Observed"="red", "Predicted"="green")) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.8, 0.5), legend.text = element_text(size = 16))

#============================================TYPHOID DLNM ANALYSIS===============================================================================
#cross basis for weekly rainfall, temperature 
varknots=equalknots(wk.dlnm$rainfall, fun = "ns", df=6)
lagknots <- logknots(4, 2)
wk.cb.rain.typhoid <- crossbasis(wk.dlnm$rainfall, lag =4, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))
summary(wk.cb.rain.typhoid)

varknots=equalknots(wk.dlnm$temperature, fun = "ns", df=7)
lagknots <- logknots(4, 2)
wk.cb.temp.typhoid <- crossbasis(wk.dlnm$temperature, lag =4, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))
summary(wk.cb.temp.typhoid)

#models fitting for weekly rainfall, temperature and NTS
wk.model.typhoid0 <- glm(incidtyp ~ ns(time), family = quasipoisson(), na.action=na.exclude,  wk.dlnm)
wk.model.typhoid1 <- glm(incidtyp ~ ns(time, 7*16), family = quasipoisson(), na.action=na.exclude,  wk.dlnm)
wk.model.typhoid2 <- glm(incidtyp ~ wk.cb.rain.typhoid + ns(time, 7*16), family = quasipoisson(), na.action=na.exclude, wk.dlnm)
wk.model.typhoid3 <- glm(incidtyp ~ wk.cb.temp.typhoid + ns(time, 7*16), family = quasipoisson(), na.action=na.exclude, wk.dlnm)
wk.model.typhoid4 <- glm(incidtyp ~ wk.cb.rain.typhoid + wk.cb.temp.typhoid + ns(time, 7*16), family = quasipoisson(), na.action=na.exclude, wk.dlnm)

#prediction for weekly rainfall, temperature
wk.pred.rain.typhoid2 <- crosspred(wk.cb.rain.typhoid, wk.model.typhoid2, cen = 2.71, by=0.2)
wk.pred.rain.typhoid4 <- crosspred(wk.cb.rain.typhoid, wk.model.typhoid4, cen = 2.71, by=0.2)

wk.pred.temp.typhoid3 <- crosspred(wk.cb.temp.typhoid, wk.model.typhoid3, cen = 23.25, by=0.2)
wk.pred.temp.typhoid4 <- crosspred(wk.cb.temp.typhoid, wk.model.typhoid4, cen = 23.25, by=0.2)

#3D plots for weekly rainfall, temperature
plot(wk.pred.rain.typhoid4, xlab="Pmean (mm)", zlab="Relative risk", theta=50, phi=2, lphi=50, main="(A)")
plot(wk.pred.temp.typhoid4, xlab="Tmean (°C)", zlab="Relative risk", theta=150, phi=20, lphi=40, main="(C)")

#countour plots for weekly rainfall, temperature
plot(wk.pred.rain.typhoid4, "contour", key.title=title("RR"), plot.title=title("(B)", xlab ="Pmean (mm)", ylab = "Lag (weeks)"))
plot(wk.pred.temp.typhoid4, "contour", key.title=title("RR"), plot.title=title("(D)", xlab ="Tmean (°C)", ylab = "Lag (weeks)"))

#detailed analysis targeting specific rainfall, temperature and lag values guided by contour plot
summary(wk.dlnm$rainfall)
plot(wk.pred.rain.typhoid4, "slices", var=0, ci="n", col=1, ylim=c(0.5,4.5), lwd=3.5)
for(i in 1:2) lines(wk.pred.rain.typhoid4, "slices", var=c(4,30)[i], col=i+1, lwd=3.5) > legend("topright",paste("Precipitation =",c(0,4,30)), col=1:3, lwd=3.5)
plot(wk.pred.rain.typhoid4, "slices", var=c(0,4,30), lag=c(0,1,4), col=3, ci.arg=list(density=90,col=grey(0.2, alpha = 0.8)), ci.level=0.95, lwd=3.5, xlab="Precipitation(mm)", ylab="Relative risk (typhoid)", cex.lab=1.3, cex.axis=1.3)
plot(wk.pred.rain.typhoid4, "overall", xlab="Rainfall (mm)", ylab="Relative risk (typhoid)", ylim=c(0,16.5), xlim=c(0,30), ci.level=0.95, ci='l', lwd=3.5)

summary(wk.dlnm$temperature)
plot(wk.pred.temp.typhoid4, "slices", var=16, ci="n", col=1, ylim=c(0.7,1.8), lwd=3.5, main="")
for(i in 1:2) lines(wk.pred.temp.typhoid4, "slices", var=c(25,30)[i], col=i+1, lwd=3.5) > legend("topright",paste("Temperature =",c(16,25,30)), col=1:3, lwd=3.5)
plot(wk.pred.temp.typhoid4, "slices", var=c(16, 25, 30), lag=c(0,1,4), col=4, ci.arg=list(density=90,col=grey(0.2, alpha = 0.8)), ci.level=0.95, lwd=3.5, xlab="Temperature(°C)", ylab="Relative risk (typhoid)", cex.lab=1.3, cex.axis=1.3)
plot(wk.pred.temp.typhoid4, "overall", xlab="Temperature (°C)", ylab="Relative risk (typhoid)", ylim=c(0,4), xlim=c(16,30), ci.level=0.95, ci='l', lwd=3.5)

#plot observed v predicted values
ggplot(wk.dlnm, aes(time)) + 
  geom_point(aes(y=wk.dlnm$incidtyp, colour="Observed"), size=1) +
  geom_point(aes(y=predict(wk.model.typhoid4, type="response"), colour="Predicted"), size=1) + 
  labs(title="", x="Week", y = "Obs versus Pred typhoid incidence rate") + 
  scale_colour_manual(name="", values=c("Observed"="red", "Predicted"="green")) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.5, 0.5), legend.text = element_text(size = 16))

#============================================DEVIANCE ANALYSIS===============================================================================
#test for dispersion in the models above aasuming poission
wk.model.iNTS4.test <- glm(NTS ~ wk.cb.rain.iNTS + wk.cb.temp.iNTS + ns(time, 7*16), family = poisson(), wk.dlnm)
wk.model.typhoid4.test <- glm(typhi ~ wk.cb.rain.typhoid + wk.cb.temp.typhoid + ns(time, 7*16), family = poisson(), wk.dlnm)
dispersiontest(wk.model.iNTS4.test)
dispersiontest(wk.model.typhoid4.test)

#select best multivariate model fitting for NTS
summary(wk.model.iNTS0)
M0 <- data.frame(wk.dlnm,link=predict(wk.model.iNTS0, type="link"), fit=predict(wk.model.iNTS0, type="response"), pearson=residuals(wk.model.iNTS0, type="pearson"), resid=residuals(wk.model.iNTS0, type="response"), residSqr=residuals(wk.model.iNTS0, type="response")^2)
v0 <- ggplot(data=M1, aes(x=fit, y=residSqr, color="R^2 = 0.5525")) + scale_colour_manual(name="NTS, model 1",values=c("R^2 = 0.5525"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared0 <- 1-(wk.model.iNTS0$deviance/wk.model.iNTS0$null.deviance)

summary(wk.model.iNTS1)
M1 <- data.frame(wk.dlnm,link=predict(wk.model.iNTS1, type="link"), fit=predict(wk.model.iNTS1, type="response"), pearson=residuals(wk.model.iNTS1, type="pearson"), resid=residuals(wk.model.iNTS1, type="response"), residSqr=residuals(wk.model.iNTS1, type="response")^2)
v1 <- ggplot(data=M1, aes(x=fit, y=residSqr, color="R^2 = 0.8567")) + scale_colour_manual(name="NTS, model 2",values=c("R^2 = 0.8567"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared1 <- 1-(wk.model.iNTS1$deviance/wk.model.iNTS1$null.deviance)

summary(wk.model.iNTS2)
M2 <- data.frame(wk.dlnm,link=predict(wk.model.iNTS2, type="link"), fit=predict(wk.model.iNTS2, type="response"), pearson=residuals(wk.model.iNTS2, type="pearson"), resid=residuals(wk.model.iNTS2, type="response"), residSqr=residuals(wk.model.iNTS2, type="response")^2)
v2 <- ggplot(data=M2, aes(x=fit, y=residSqr, color="R^2 = 0.8625")) + scale_colour_manual(name="NTS, model 3",values=c("R^2 = 0.8625"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared2 <- 1-(wk.model.iNTS2$deviance/wk.model.iNTS2$null.deviance)

summary(wk.model.iNTS3)
M3 <- data.frame(wk.dlnm,link=predict(wk.model.iNTS3, type="link"), fit=predict(wk.model.iNTS3, type="response"), pearson=residuals(wk.model.iNTS3, type="pearson"), resid=residuals(wk.model.iNTS3, type="response"), residSqr=residuals(wk.model.iNTS3, type="response")^2)
v3 <- ggplot(data=M3, aes(x=fit, y=residSqr, color="R^2 = 0.8626")) + scale_colour_manual(name="NTS, model 4",values=c("R^2 = 0.8626"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared3 <- 1-(wk.model.iNTS3$deviance/wk.model.iNTS3$null.deviance)

summary(wk.model.iNTS4)
M4 <- data.frame(wk.dlnm,link=predict(wk.model.iNTS4, type="link"), fit=predict(wk.model.iNTS4, type="response"), pearson=residuals(wk.model.iNTS4, type="pearson"), resid=residuals(wk.model.iNTS4, type="response"), residSqr=residuals(wk.model.iNTS4, type="response")^2)
v4 <- ggplot(data=M4, aes(x=fit, y=residSqr, color="R^2 = 0.8677")) + scale_colour_manual(name="NTS, model 5",values=c("R^2 = 0.8677"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared4 <- 1-(wk.model.iNTS4$deviance/wk.model.iNTS4$null.deviance)

#select best multivariate model fitting for typhoid
summary(wk.model.typhoid0)
N0 <- data.frame(wk.dlnm,link=predict(wk.model.typhoid0, type="link"), fit=predict(wk.model.typhoid0, type="response"), pearson=residuals(wk.model.typhoid0, type="pearson"), resid=residuals(wk.model.typhoid0, type="response"), residSqr=residuals(wk.model.typhoid0, type="response")^2)
z0 <- ggplot(data=N0, aes(x=fit, y=residSqr, color="R^2 = 0.6844")) + scale_colour_manual(name="typhoid, model 1",values=c("R^2 = 0.6844"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared0 <- 1-(wk.model.typhoid0$deviance/wk.model.typhoid0$null.deviance)

summary(wk.model.typhoid1)
N1 <- data.frame(wk.dlnm,link=predict(wk.model.typhoid1, type="link"), fit=predict(wk.model.typhoid1, type="response"), pearson=residuals(wk.model.typhoid1, type="pearson"), resid=residuals(wk.model.typhoid1, type="response"), residSqr=residuals(wk.model.typhoid1, type="response")^2)
z1 <- ggplot(data=N1, aes(x=fit, y=residSqr, color="R^2 = 0.9043")) + scale_colour_manual(name="typhoid, model 2",values=c("R^2 = 0.9043"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared1 <- 1-(wk.model.typhoid1$deviance/wk.model.typhoid1$null.deviance)

summary(wk.model.typhoid2)
N2 <- data.frame(wk.dlnm,link=predict(wk.model.typhoid2, type="link"), fit=predict(wk.model.typhoid2, type="response"), pearson=residuals(wk.model.typhoid2, type="pearson"), resid=residuals(wk.model.typhoid2, type="response"), residSqr=residuals(wk.model.typhoid2, type="response")^2)
z2 <- ggplot(data=N2, aes(x=fit, y=residSqr, color="R^2 = 0.9063")) + scale_colour_manual(name="typhoid, model 3",values=c("R^2 = 0.9063"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared2 <- 1-(wk.model.typhoid2$deviance/wk.model.typhoid2$null.deviance)

summary(wk.model.typhoid3)
N3 <- data.frame(wk.dlnm,link=predict(wk.model.typhoid3, type="link"), fit=predict(wk.model.typhoid3, type="response"), pearson=residuals(wk.model.typhoid3, type="pearson"), resid=residuals(wk.model.typhoid3, type="response"), residSqr=residuals(wk.model.typhoid3, type="response")^2)
z3 <- ggplot(data=N3, aes(x=fit, y=residSqr, color="R^2 = 0.9080")) + scale_colour_manual(name="typhoid, model 4",values=c("R^2 = 0.9080"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared3 <- 1-(wk.model.typhoid3$deviance/wk.model.typhoid3$null.deviance)

summary(wk.model.typhoid4)
N4 <- data.frame(wk.dlnm,link=predict(wk.model.typhoid4, type="link"), fit=predict(wk.model.typhoid4, type="response"), pearson=residuals(wk.model.typhoid4, type="pearson"), resid=residuals(wk.model.typhoid4, type="response"), residSqr=residuals(wk.model.typhoid4, type="response")^2)
z4 <- ggplot(data=N4, aes(x=fit, y=residSqr, color="R^2 = 0.9103")) + scale_colour_manual(name="typhoid, model 5",values=c("R^2 = 0.9103"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared4 <- 1-(wk.model.typhoid4$deviance/wk.model.typhoid4$null.deviance)

#plot all squared residuals v mean fitted values
grid.arrange(grobs=list(v0,v1, v2, v3, v4, z0,z1, z2, z3, z4), ncol=5, nrow=2)

