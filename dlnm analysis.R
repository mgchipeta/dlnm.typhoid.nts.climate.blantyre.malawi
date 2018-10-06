#02/08/2018
#by Deus

#Install required packages
dlnm.analysis.packages <- c("tidyverse","lubridate","xts","ggthemes","PerformanceAnalytics","reshape2","rugarch","timetk","parallel","timeSeries","tseries","data.table","ggplot2","dlnm","broom","caret","gridExtra","splines","splines2","pspline","cowplot","mgcv","spi","chron","gridGraphics","grid","pscl","MASS", "AER", "Hmisc", "MuMIn", "VGAM","here" )

#load required packages
lapply(dlnm.analysis.packages, library, character.only=TRUE)

#========LOAD CASE DATA========

#load typhoid and NTS cases dataset.
case <- read.csv(here("Time.Series", "data", "case.csv"))
case$case_date <- dmy(case$date_s)
case$sex <-case$gender
case$age <- case$age_yrs
case$organism_type <- case$org_type
case$organism <- case$org
case$case_count <- c(1)
demog <- subset(case, select = c(case_date, age, sex, organism_type, organism, case_count))

#create separate datasets for typhi and NTS cases.
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

#convert typhi and NTS data frames to xts objects for use in time series plotting.
case_typhi.x = as.xts(case_typhi.x.w[,-1,drop = FALSE], order.by = as.Date(case_typhi.x.w[,1]))
case_iNTS.x = as.xts(case_iNTS.x.w[,-1,drop = FALSE], order.by = as.Date(case_iNTS.x.w[,1]))

#weekly aggregates for typhi and iNTS cases.
wk.typhi <- apply.weekly(case_typhi.x, FUN = sum)
wk.iNTS <- apply.weekly(case_iNTS.x, FUN = sum)

#check seasonality.
wk.typhi.seas <- as.ts(wk.typhi)
stl(wk.typhi.seas)
wk.iNTS.seas <- as.ts(wk.iNTS)
stl(wk.iNTS.seas)

#========LOAD CLIMATE DATA========

#load climate dataset.
climate <- read.csv(here("Time.Series", "data", "climate.csv"))

#calculate daily average values for temperature and rainfall.
climate$rainfall <- (climate$chil_r + climate$chic_r)/2 
climate$temperature <- ((climate$chil_mint + climate$chil_maxt)/2 + (climate$chic_mint + climate$chic_maxt)/2)/2

climate$climate_date <- dmy(climate$date)
climate <- subset(climate, select = c(climate_date, rainfall, temperature))
dlnm.climate <- climate

#create separate datasets for temperature and rainfall. these have initial daily values, no need to wrangle.
climate_rain <- subset(climate, select = c(climate_date, rainfall))
climate_temp <- subset(climate, select = c(climate_date, temperature))

#convert temperature and rainfall data frames to xts objects for use in time series plotting.
climate_rain.x = as.xts(climate_rain[,-1,drop = FALSE], order.by = as.Date(climate_rain[,1]))
climate_temp.x = as.xts(climate_temp[,-1,drop = FALSE], order.by = as.Date(climate_temp[,1]))

#weekly aggregates for rainfall and temperature
wk.rain <- apply.weekly(climate_rain.x, FUN = mean)
wk.temp <- apply.weekly(climate_temp.x, FUN = mean)

#delete unnecessary objects
rm(list = ls()[grep("^climate", ls())]) 
rm(list = ls()[grep("^case", ls())])
rm(list = ls()[grep("^alldates", ls())]) 

#=========SEASONALITY AND DEMOGRAPHICS PLOTS========

#create tibbles for easy plotting in ggplot2
wk.typhi <-tk_tbl(wk.typhi, preserve_index = TRUE, rename_index = "date") 
wk.iNTS <-tk_tbl(wk.iNTS, preserve_index = TRUE, rename_index = "date") 
wk.rain <-tk_tbl(wk.rain, preserve_index = TRUE, rename_index = "date") 
wk.temp <-tk_tbl(wk.temp, preserve_index = TRUE, rename_index = "date") 

#linear interpolation and extrapolation
census.year <-c(1998, 2008)
census.popn <-c(809397, 1022680)
census.count.1998.2008 <- approx(census.year, census.popn, n=11) 
census.count.2009.2015 <- approxExtrap(census.year, census.popn, xout=c(2009, 2010, 2011, 2012, 2013, 2014, 2015))

#incidence rates of NTS
wk.iNTS$incid <- NA
wk.iNTS2000 <- subset(wk.iNTS, year(date) == 2000)
wk.iNTS2000$incid <- wk.iNTS2000$case_count*100000/852053.6
wk.iNTS2001 <- subset(wk.iNTS, year(date) == 2001)
wk.iNTS2001$incid <- wk.iNTS2001$case_count*100000/873381.9
wk.iNTS2002 <- subset(wk.iNTS, year(date) == 2002)
wk.iNTS2002$incid <- wk.iNTS2002$case_count*100000/894710.2
wk.iNTS2003 <- subset(wk.iNTS, year(date) == 2003)
wk.iNTS2003$incid <- wk.iNTS2003$case_count*100000/916038.5
wk.iNTS2004 <- subset(wk.iNTS, year(date) == 2004)
wk.iNTS2004$incid <- wk.iNTS2004$case_count*100000/937366.8
wk.iNTS2005 <- subset(wk.iNTS, year(date) == 2005)
wk.iNTS2005$incid <- wk.iNTS2005$case_count*100000/958695.1
wk.iNTS2006 <- subset(wk.iNTS, year(date) == 2006)
wk.iNTS2006$incid <- wk.iNTS2006$case_count*100000/980023.4
wk.iNTS2007 <- subset(wk.iNTS, year(date) == 2007)
wk.iNTS2007$incid <- wk.iNTS2007$case_count*100000/1001351.7
wk.iNTS2008 <- subset(wk.iNTS, year(date) == 2008)
wk.iNTS2008$incid <- wk.iNTS2008$case_count*100000/1022680.0
wk.iNTS2009 <- subset(wk.iNTS, year(date) == 2009)
wk.iNTS2009$incid <- wk.iNTS2009$case_count*100000/1044008
wk.iNTS2010 <- subset(wk.iNTS, year(date) == 2010)
wk.iNTS2010$incid <- wk.iNTS2010$case_count*100000/1065337
wk.iNTS2011 <- subset(wk.iNTS, year(date) == 2011)
wk.iNTS2011$incid <- wk.iNTS2011$case_count*100000/1086665
wk.iNTS2012 <- subset(wk.iNTS, year(date) == 2012)
wk.iNTS2012$incid <- wk.iNTS2012$case_count*100000/1107993
wk.iNTS2013 <- subset(wk.iNTS, year(date) == 2013)
wk.iNTS2013$incid <- wk.iNTS2013$case_count*100000/1129322
wk.iNTS2014 <- subset(wk.iNTS, year(date) == 2014)
wk.iNTS2014$incid <- wk.iNTS2014$case_count*100000/1150650
wk.iNTS2015 <- subset(wk.iNTS, year(date) == 2015)
wk.iNTS2015$incid <- wk.iNTS2015$case_count*100000/1171978
wk.iNTS <- bind_rows(wk.iNTS2000,wk.iNTS2001,wk.iNTS2002,wk.iNTS2003,wk.iNTS2004,wk.iNTS2005,wk.iNTS2006,wk.iNTS2007,wk.iNTS2008,wk.iNTS2009,wk.iNTS2010,wk.iNTS2011,wk.iNTS2012,wk.iNTS2013,wk.iNTS2014,wk.iNTS2015)
rm(list = ls()[grep("^wk.iNTS20", ls())])

#incidence rates of typhoid fever
wk.typhi$incid <- NA
wk.typhi2000 <- subset(wk.typhi, year(date) == 2000)
wk.typhi2000$incid <- wk.typhi2000$case_count*100000/852053.6
wk.typhi2001 <- subset(wk.typhi, year(date) == 2001)
wk.typhi2001$incid <- wk.typhi2001$case_count*100000/873381.9
wk.typhi2002 <- subset(wk.typhi, year(date) == 2002)
wk.typhi2002$incid <- wk.typhi2002$case_count*100000/894710.2
wk.typhi2003 <- subset(wk.typhi, year(date) == 2003)
wk.typhi2003$incid <- wk.typhi2003$case_count*100000/916038.5
wk.typhi2004 <- subset(wk.typhi, year(date) == 2004)
wk.typhi2004$incid <- wk.typhi2004$case_count*100000/937366.8
wk.typhi2005 <- subset(wk.typhi, year(date) == 2005)
wk.typhi2005$incid <- wk.typhi2005$case_count*100000/958695.1
wk.typhi2006 <- subset(wk.typhi, year(date) == 2006)
wk.typhi2006$incid <- wk.typhi2006$case_count*100000/980023.4
wk.typhi2007 <- subset(wk.typhi, year(date) == 2007)
wk.typhi2007$incid <- wk.typhi2007$case_count*100000/1001351.7
wk.typhi2008 <- subset(wk.typhi, year(date) == 2008)
wk.typhi2008$incid <- wk.typhi2008$case_count*100000/1022680.0
wk.typhi2009 <- subset(wk.typhi, year(date) == 2009)
wk.typhi2009$incid <- wk.typhi2009$case_count*100000/1044008
wk.typhi2010 <- subset(wk.typhi, year(date) == 2010)
wk.typhi2010$incid <- wk.typhi2010$case_count*100000/1065337
wk.typhi2011 <- subset(wk.typhi, year(date) == 2011)
wk.typhi2011$incid <- wk.typhi2011$case_count*100000/1086665
wk.typhi2012 <- subset(wk.typhi, year(date) == 2012)
wk.typhi2012$incid <- wk.typhi2012$case_count*100000/1107993
wk.typhi2013 <- subset(wk.typhi, year(date) == 2013)
wk.typhi2013$incid <- wk.typhi2013$case_count*100000/1129322
wk.typhi2014 <- subset(wk.typhi, year(date) == 2014)
wk.typhi2014$incid <- wk.typhi2014$case_count*100000/1150650
wk.typhi2015 <- subset(wk.typhi, year(date) == 2015)
wk.typhi2015$incid <- wk.typhi2015$case_count*100000/1171978
wk.typhi <- bind_rows(wk.typhi2000,wk.typhi2001,wk.typhi2002,wk.typhi2003,wk.typhi2004,wk.typhi2005,wk.typhi2006,wk.typhi2007,wk.typhi2008,wk.typhi2009,wk.typhi2010,wk.typhi2011,wk.typhi2012,wk.typhi2013,wk.typhi2014,wk.typhi2015)
rm(list = ls()[grep("^wk.typhi20", ls())])

#weekly plots of incidence rates and climate over 16 (NTS) and 5 years (typhoid fever)
wk.p1 <- ggplot() + geom_line(data = wk.iNTS, aes(x=date, y=incid), stat = "identity", colour="tan4", size = 0.6) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + ylab('NTS incidence rate') + xlab('') + labs(title="(A)") + scale_x_date(date_breaks = "12 month", date_labels = "%b'%y") 
wk.p2 <- ggplot() + geom_line(data = wk.typhi, aes(x=date, y=incid), stat = "identity", colour="red4", size = 0.6) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + ylab('Typhoid incidence rate') + xlab('') + labs(title="(B)") + scale_x_date(date_breaks = "12 month", date_labels = "%b'%y") 
wk.p3 <- ggplot() + geom_line(data = wk.rain, aes(x=date, y=rainfall), colour = "green4", alpha = 0.5, size = 0.6) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + ylab('Precipitation (mm)') + xlab('') + labs(title="(C)") + scale_x_date(date_breaks = "12 month", date_labels ="%b'%y") 
wk.p4 <- ggplot() + geom_line(data = wk.temp, aes(x=date, y=temperature), colour = "blue4", alpha = 0.5, size = 0.6) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + ylab('Temperature (°C)') + xlab("Month'Year") + labs(title="(D)") + scale_x_date(date_breaks = "12 month", date_labels ="%b'%y") 

#combine plots.
grid.arrange(grobs=list(wk.p1, wk.p2, wk.p3, wk.p4), ncol=1, nrow=4)

#seasonal dynamics of NTS, typhoid fever and climate as boxplots
j <- seq(from=1, to=53, by=4)
k1 <-seq(from=16, to=32, by=2)
k2 <-seq(from=0, to=20, by=2)
k3 <-seq(from=0, to=5, by=1)
k4 <-seq(from=0, to=3, by=0.5)

pox1 <- ggplot(wk.temp, aes(x=week(date), y=temperature))  + 
  geom_boxplot(aes( group=week(date)), outlier.shape = NA, color="blue4", fill="blue4", alpha=0.2) + 
  labs(title="(D)",x="Week", y = "Temperature (°C)") + 
  scale_x_discrete(limits = j) + 
  scale_y_discrete(limits = k1) +
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) +
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

pox2 <- ggplot(wk.rain, aes(x=week(date), y=rainfall))  + 
  geom_boxplot(aes( group=week(date)), outlier.shape = NA, color="green4", fill="green4", alpha=0.2) + 
  labs(title="(C)",x="", y = "Precipitation (mm)") + 
  scale_x_discrete(limits = j) + 
  scale_y_discrete(limits = k2) + 
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

pox3 <- ggplot(subset(wk.iNTS, year(date)<2011), aes(x=week(date), y=incid))  + 
  geom_boxplot(aes( group=week(date)), outlier.shape = NA, color="tan4", fill="tan4", alpha=0.2) + 
  labs(title="(A)",x="", y = "NTS incidence rate") + 
  scale_x_discrete(limits = j) + 
  scale_y_discrete(limits = k3) + 
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

pox4 <- ggplot(subset(wk.typhi, year(date)>2010), aes(x=week(date), y=incid))  + 
  geom_boxplot(aes( group=week(date)), outlier.shape = NA, color="red4", fill="red", alpha=0.2) + 
  labs(title="(B)", x="", y = "Typhoid incidence rate") + 
  scale_x_discrete(limits = j) + 
  scale_y_discrete(limits = k4) + 
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

#combine boxplots.
grid.arrange(grobs=list(pox3, pox4, pox2, pox1), ncol=1, nrow=4)

#distributions of typhi and NTS cases by sex, age or both.
demog$sex[demog$sex == ""] <- NA
demog$age[demog$age == ""] <- NA
demog$date <- ymd(demog$case_date)
demog$year <- year(demog$date)

#age and sex frequency distribution
agesex.p1 <- ggplot(subset(demog, organism=="iNTS" & !is.na(age) & sex != "Unknown"), aes(x=age, color=sex)) + geom_freqpoly(position=position_dodge(width=2), binwidth=1, size=1) + theme_classic() + scale_x_continuous(limits = c(0, 91), breaks = seq(0, 91, 5)) + xlab("Age (years)") + ylab("Number of nontyphoid cases") + ggtitle("(A)") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + theme(legend.justification=c(0.5,0), legend.position = c(0.5, 0.5), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) + labs(color="Missing: 4,575 (48.1%)")
agesex.p2 <- ggplot(subset(demog, organism=="typhi" & !is.na(age) & sex != "Unknown"), aes(x=age, color=sex)) + geom_freqpoly(position=position_dodge(width=2), binwidth=1, size=1) + theme_classic() + scale_x_continuous(limits = c(0, 91), breaks = seq(0, 91, 5)) + xlab("Age (years)") + ylab("Number of typhoid cases") + ggtitle("(B)") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + theme(legend.justification=c(0.5,0), legend.position = c(0.5, 0.5), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) + labs(color="Missing: 95 (3.6%)") 

#combine plots.
grid.arrange(grobs=list(agesex.p1, agesex.p2), ncol=1, nrow=2)

#========CROSS-BASIS MODEL SELECTION FOR NTS OUTCOME========

#prepare weekly dlnm dataset
wk.dlnmN <- bind_cols(wk.iNTS, wk.rain, wk.temp, id=NULL)
wk.dlnmN$date1 <- wk.dlnmN$date2 <- wk.dlnmN$date3 <- NULL
wk.dlnmN <- rename(wk.dlnmN, NTS=case_count, incidNTS=incid)
wk.dlnmN <- subset(wk.dlnmN, year(date) < 2011)
wk.dlnmN$time <- seq.int(from = 1, to=574, by=1)
wk.dlnmN$year <- year(wk.dlnmN$date)
wk.dlnmN$month <- month(wk.dlnmN$date)
wk.dlnmN$week <- week(wk.dlnmN$date)

#test for seasonality
library(forecast)
dataN <- wk.dlnmN
dataN$incidNTS <- dataN$rainfall <- dataN$temperature <- dataN$time <- dataN$year <- dataN$month <- dataN$week <- dataN$date<- NULL
x <- ts(dataN, frequency=365.25/7)
fit <- tbats(x)
seasonal <- !is.null(fit$seasonal)

#selecting an optimal model given rainfall, temperature, lag, and seasonality and long-term trend terms
#manipulates AIC to Q-AIC
n.quasipoisson <- function(...) { 
  res <- quasipoisson(...)
  res$aic <- poisson(...)$aic 
  res
}

#defines all possible dfs for rainfall (3-7) and lag (3-5) and seasonality (3-7)
varknots=equalknots(wk.dlnmN$rainfall, fun = "ns", df=3)
lagknots <- logknots(4, fun="ns", df=4)
wk.cb.rain.nts <- crossbasis(wk.dlnmN$rainfall, lag =4, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))

sea.df <-3
nts.model <- glm(wk.dlnmN$NTS ~ wk.cb.rain.nts + ns(time, sea.df*11), family = n.quasipoisson(), na.action=na.delete, wk.dlnmN)
QAIC(nts.model, chat=summary(nts.model)$dispersion)

#defines all possible dfs for temperature (3-7) and lag (3-5) and seasonality (3-7)
varknots=equalknots(wk.dlnmN$temperature, fun = "ns", df=3)
lagknots <- logknots(4, fun="ns", df=3)
wk.cb.temp.nts <- crossbasis(wk.dlnmN$temperature, lag =4, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))

sea.df <-3
nts.model <- glm(wk.dlnmN$NTS ~ wk.cb.temp.nts + ns(time, sea.df*11), family = n.quasipoisson(), na.action=na.delete, wk.dlnmN)
QAIC(nts.model, chat=summary(nts.model)$dispersion)

#========NTS DLNM MULTIVARIATE ANALYSIS========

#final cross-basis for weekly rainfall, temperature 
sort(wk.dlnmN$rainfall, decreasing = FALSE)
varknots=equalknots(wk.dlnmN$rainfall, fun="ns", df=3)
lagknots=logknots(4, fun="ns", df=4)
wk.cb.rain.iNTS <- crossbasis(wk.dlnmN$rainfall, lag=4, argvar=list(knots=varknots), arglag=list(knots=lagknots))
summary(wk.cb.rain.iNTS)

sort(wk.dlnmN$temperature, decreasing = FALSE)
varknots=equalknots(wk.dlnmN$temperature, fun="ns", df=3)
lagknots=logknots(4, fun="ns", df=3)
wk.cb.temp.iNTS <- crossbasis(wk.dlnmN$temperature, lag=4, argvar=list(fun="ns", knots=varknots), arglag=list(knots=lagknots))
summary(wk.cb.temp.iNTS)

#models fitting for weekly rainfall, temperature and NTS (either quasipossion or negative binomial)
wk.model.iNTS4 <- glm(incidNTS ~ wk.cb.rain.iNTS + wk.cb.temp.iNTS + ns(time, (3*11)), family = quasipoisson(), na.action=na.exclude, wk.dlnmN)

#prediction for weekly rainfall, temperature
wk.pred.rain.iNTS4 <- crosspred(wk.cb.rain.iNTS, wk.model.iNTS4, cen = 0, by=0.2)
wk.pred.temp.iNTS4 <- crosspred(wk.cb.temp.iNTS, wk.model.iNTS4, cen = 23.3, by=0.2)

#countour, 3D plots and lag/predictor curves for weekly rainfall
NTS.P.3D <- function() {plot(wk.pred.rain.iNTS4, xlab="Precipitation (mm)", ylab="Lag (weeks)", zlab="Relative risk (NTS)", theta=50, phi=5, lphi=100, main="(A)", cex.lab=1.6, cex.axis=1.6)}
NTS.P.Contour <- function() {plot(wk.pred.rain.iNTS4, "contour", key.title=title("NTS"), plot.title=title("(B)", xlab ="Precipitation (mm)", ylab = "Lag (weeks)", cex.lab=1.6, cex.axis=1.6))}
plot_grid(NTS.P.3D, NTS.P.Contour)

NTS.P.slice1 <- function() {plot(wk.pred.rain.iNTS4, "slices", var=c(15,25), lag=c(1,3), col="green4", ci.arg=list(col=terrain.colors(20, alpha = 0.5)), ci.level=0.95, lwd=3.5,  col="gray0", ylab="Relative risk (NTS)", cex.lab=1.6, cex.axis=1.6)}
NTS.P.slice2 <- function() {plot(wk.pred.rain.iNTS4, "overall", xlab="Precipitation (mm)", ylab="Relative risk (NTS)", col="green4", ci.arg=list(col=terrain.colors(20, alpha = 1)), ci.level=0.95, ci='l', lwd=3.5, cex.lab=1.6, cex.axis=1.6)}
plot_grid(NTS.P.slice1, NTS.P.slice2)

#countour, 3D plots and lag/predictor curves for weekly temperature
NTS.T.3D <- function() {plot(wk.pred.temp.iNTS4, xlab="Temperature (°C)", ylab="Lag (weeks)", zlab="Relative risk (NTS)", theta=200, phi=35, lphi=40, main="(A)", cex.lab=1.6, cex.axis=1.6)}
NTS.T.Contour <- function() {plot(wk.pred.temp.iNTS4, "contour", key.title=title("NTS"), plot.title=title("(B)", xlab ="Temperature (°C)", ylab = "Lag (weeks)", cex.lab=1.6, cex.axis=1.6))}
plot_grid(NTS.T.3D, NTS.T.Contour)

NTS.T.slice1 <- function() {plot(wk.pred.temp.iNTS4, "slices", var=c(19,32), lag=c(0,1), col="red4", ci.arg=list(col=heat.colors(16, alpha = 0.5)), ci.level=0.95, lwd=3.5,  col="gray0", ylab="Relative risk (NTS)", cex.lab=1.6, cex.axis=1.6)}
NTS.T.slice2 <- function() {plot(wk.pred.temp.iNTS4, "overall", xlab="Temperature (°C)", ylab="Relative risk (NTS)", ylim=c(0,2), col="red4", ci.arg=list(col=heat.colors(20, alpha = 1)), ci.level=0.95, ci='l', lwd=3.5, cex.lab=1.6, cex.axis=1.6)}
plot_grid(NTS.T.slice1, NTS.T.slice2)

#plot observed v predicted values
ggplot(wk.dlnmN, aes(time)) + 
  geom_point(aes(y=wk.dlnmN$incidNTS, colour="Observed"), size=1) +
  geom_line(aes(y=predict(wk.model.iNTS7, type="response"), colour="Predicted7"), size=0.6) + 
  geom_line(aes(y=predict(wk.model.iNTS4, type="response"), colour="Predicted4"), size=0.6) + 
  labs(title="", x="Week", y = "Nontyphoid incidence rate") + 
  scale_colour_manual(name="", values=c("Observed"="gray", "Predicted7"="black", "Predicted4"="blue")) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.8, 0.5), legend.text = element_text(size = 16))

#model diagnostics
#plot residuals over time
residual4 <- residuals(wk.model.iNTS4,type="deviance")
plot(wk.dlnmN$date, residual4, pch=19, cex=0.4, col=grey(0.6),
     main="Residuals over time", ylab="Residuals deviance", xlab="Year")
abline(h=0,lty=2,lwd=2)

pacf(residual4,na.action=na.omit,main="From original model")

wk.model.iNTS5 <- update(wk.model.iNTS4,.~.+Lag(residual4,1))
pacf(residuals(wk.model.iNTS5,type="deviance"),na.action=na.omit,
     main="From model adjusted for residual autocorrelation")

wk.model.iNTS6 <- update(wk.model.iNTS5,.~.+Lag(residual4,2))
pacf(residuals(wk.model.iNTS6,type="deviance"),na.action=na.omit,
     main="From model adjusted for residual autocorrelation")

wk.model.iNTS7 <- update(wk.model.iNTS6,.~.+Lag(residual4,4))
pacf(residuals(wk.model.iNTS7,type="deviance"),na.action=na.omit,
     main="From model adjusted for residual autocorrelation")

#========CROSS-BASIS MODEL SELECTION WITH TYPHOID FEVER OUTCOME========

#prepare weekly typhoid dlnm dataset
wk.dlnmT <- bind_cols(wk.typhi, wk.rain, wk.temp, id=NULL)
wk.dlnmT$date1 <- wk.dlnmT$date2 <- NULL
wk.dlnmT <- rename(wk.dlnmT, typhi=case_count, incidtyp = incid)
wk.dlnmT$time <- seq.int(from = 1, to=836, by=1)
wk.dlnmT$year <- year(wk.dlnmT$date)
wk.dlnmT$month <- month(wk.dlnmT$date)
wk.dlnmT$week <- week(wk.dlnmT$date)
wk.dlnmT <- subset(wk.dlnmT, year>2010)

#test for seasonality
library(forecast)
dataT <- wk.dlnmT
dataT$typhi <- dataT$rainfall <- dataT$temperature <- dataT$time <- dataT$year <- dataT$month <- dataT$week <- dataT$date<- NULL
x <- ts(dataT, frequency=4)
fit <- tbats(x)
seasonal <- !is.null(fit$seasonal)

#manipulate the AIC so its QAIC and compare to choose an optimal model
t.quasipoisson <- function(...) { 
  res <- quasipoisson(...) 
  res$aic <- poisson(...)$aic 
  res
}

#defines all possible dfs for rainfall (3-7) and lag (3-5) and seasonality (3-7)
varknots=equalknots(wk.dlnmT$rainfall, fun = "ns", df=3)
lagknots <- logknots(4, fun="ns", df=4)
wk.cb.rain.typhi <- crossbasis(wk.dlnmT$rainfall, lag =4, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))

sea.df <-3
typhi.model <- glm(wk.dlnmT$typhi ~ wk.cb.rain.typhi + ns(time, sea.df*11), family = n.quasipoisson(), na.action=na.delete, wk.dlnmT)
QAIC(typhi.model, chat=summary(typhi.model)$dispersion)

#defines all possible dfs for temperature (3-7) and lag (3-5) and seasonality (3-7)
varknots=equalknots(wk.dlnmT$temperature, fun = "ns", df=3)
lagknots <- logknots(4, fun="ns", df=3)
wk.cb.temp.typhi <- crossbasis(wk.dlnmT$temperature, lag =4, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))

sea.df <-3
typhi.model <- glm(wk.dlnmT$typhi ~ wk.cb.temp.typhi + ns(time, sea.df*11), family = n.quasipoisson(), na.action=na.delete, wk.dlnmT)
QAIC(typhi.model, chat=summary(typhi.model)$dispersion)

#========TYPHOID DLNM MULTIVARIATE ANALYSIS========

#cross basis for weekly rainfall, temperature 
sort(wk.dlnmT$rainfall, decreasing = FALSE)
varknots=equalknots(wk.dlnmT$rainfall, fun = "ns", df=4)
lagknots <- logknots(4, df=3)
wk.cb.rain.typhoid <- crossbasis(wk.dlnmT$rainfall, lag =4, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(wk.cb.rain.typhoid)

sort(wk.dlnmT$temperature, decreasing = FALSE)
varknots=equalknots(wk.dlnmT$temperature, fun = "ns", df=3)
lagknots <- logknots(4, df=3)
wk.cb.temp.typhoid <- crossbasis(wk.dlnmT$temperature, lag =4, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(wk.cb.temp.typhoid)

#models fitting for weekly rainfall, temperature and typhoid
wk.model.typhoid4 <- glm(incidtyp ~ wk.cb.rain.typhoid + wk.cb.temp.typhoid + ns(time, 1*5), family = quasipoisson(), na.action=na.exclude, wk.dlnmT)

#prediction for weekly rainfall, temperature
wk.pred.rain.typhoid4 <- crosspred(wk.cb.rain.typhoid, wk.model.typhoid4, cen = 0, by=0.2)
wk.pred.temp.typhoid4 <- crosspred(wk.cb.temp.typhoid, wk.model.typhoid4, cen = 23.3, by=0.2)

#countour, 3D plots and lag/predictor curves for weekly rainfall
typhoid.P.3D <- function() {plot(wk.pred.rain.typhoid4, xlab="Precipitation (mm)", ylab="Lag (weeks)", zlab="Relative risk (typhoid)", theta=150, phi=5, lphi=60, main="(A)", cex.lab=1.6, cex.axis=1.6)}
typhoid.P.Contour <- function() {plot(wk.pred.rain.typhoid4, "contour", key.title=title("TYP"), plot.title=title("(B)", xlab ="Precipitation (mm)", ylab = "Lag (weeks)", cex.lab=1.6, cex.axis=1.6))}
plot_grid(typhoid.P.3D, typhoid.P.Contour)

typhoid.P.slice1 <- function() {plot(wk.pred.rain.typhoid4, "slices", var=c(15,35), lag=c(1,4), col="green4", ci.arg=list(col=terrain.colors(20, alpha = 0.5)), ci.level=0.95, lwd=3.5,  col="gray0", ylab="Relative risk (typhoid)", cex.lab=1.6, cex.axis=1.6)}
typhoid.P.slice2 <- function() {plot(wk.pred.rain.typhoid4, "overall", xlab="Precipitation (mm)", ylab="Relative risk (typhoid)", xlim=c(0,35), col="green4", ci.arg=list(col=terrain.colors(20, alpha = 1)), ci.level=0.95, ci='l', lwd=3.5, cex.lab=1.6, cex.axis=1.6)}
plot_grid(typhoid.P.slice1, typhoid.P.slice2)

#countour, 3D plots and lag/predictor curves for weekly temperature
typhoid.T.3D <- function() {plot(wk.pred.temp.typhoid4, xlab="Temperature (°C)", ylab="Lag (weeks)", zlab="Relative risk (typhoid)", theta=200, phi=35, lphi=40, main="(A)", cex.lab=1.6, cex.axis=1.6)}
typhoid.T.Contour <- function() {plot(wk.pred.temp.typhoid4, "contour", key.title=title("TYP"), plot.title=title("(B)", xlab ="Temperature (°C)", ylab = "Lag (weeks)", cex.lab=1.6, cex.axis=1.6))}
plot_grid(typhoid.T.3D, typhoid.T.Contour)

typhoid.T.slice1 <- function() {plot(wk.pred.temp.typhoid4, "slices", var=c(17,34), lag=c(0,4), col="red4", ci.arg=list(col=heat.colors(16, alpha = 0.5)), ci.level=0.95, lwd=3.5,  col="gray0", ylab="Relative risk (typhoid)", cex.lab=1.6, cex.axis=1.6)}
typhoid.T.slice2 <- function() {plot(wk.pred.temp.typhoid4, "overall", xlab="Temperature (°C)", ylab="Relative risk (typhoid)", xlim=c(17,34), col="red4", ci.arg=list(col=heat.colors(20, alpha = 1)), ci.level=0.95, ci='l', lwd=3.5, cex.lab=1.6, cex.axis=1.6)}
plot_grid(typhoid.T.slice1, typhoid.T.slice2)

#plot observed v predicted values
ggplot(wk.dlnmT, aes(time)) + 
  geom_point(aes(y=wk.dlnmT$incidtyp, colour="Observed"), size=1) +
  geom_point(aes(y=predict(wk.model.typhoid4, type="response"), colour="Predicted"), size=1) + 
  labs(title="", x="Week", y = "Typhoid incidence rate") + 
  scale_colour_manual(name="", values=c("Observed"="red", "Predicted"="green")) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.2, 0.5), legend.text = element_text(size = 16))

#model diagnostics
#plot residuals over time
residual4 <- residuals(wk.model.typhoid4,type="deviance")
plot(wk.dlnmN$date, residual4, pch=19, cex=0.4, col=grey(0.6),
     main="Residuals over time", ylab="Residuals deviance", xlab="Year")
abline(h=0,lty=2,lwd=2)

pacf(residual4,na.action=na.omit,main="From original model")

wk.model.typhoid5 <- update(wk.model.typhoid4,.~.+Lag(residual4,1))
pacf(residuals(wk.model.typhoid5,type="deviance"),na.action=na.omit,
     main="From model adjusted for residual autocorrelation")

wk.model.typhoid6 <- update(wk.model.typhoid5,.~.+Lag(residual4,2))
pacf(residuals(wk.model.iNTS6,type="deviance"),na.action=na.omit,
     main="From model adjusted for residual autocorrelation")

wk.model.typhoid7 <- update(wk.model.typhoid6,.~.+Lag(residual4,4))
pacf(residuals(wk.model.typhoid7,type="deviance"),na.action=na.omit,
     main="From model adjusted for residual autocorrelation")

#========DEVIANCE ANALYSIS FOR NTS AND TYPHOID FEVER========

#deviance and r-squared for NTS
summary(wk.model.iNTS0)
M0 <- data.frame(wk.dlnmN,link=predict(wk.model.iNTS0, type="link"), fit=predict(wk.model.iNTS0, type="response"), pearson=residuals(wk.model.iNTS0, type="pearson"), resid=residuals(wk.model.iNTS0, type="response"), residSqr=residuals(wk.model.iNTS0, type="response")^2)
v0 <- ggplot(data=M0, aes(x=fit, y=residSqr, color="R^2 = 0.6024")) + scale_colour_manual(name="NTS, model 1",values=c("R^2 = 0.6024"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared0 <- 1-(wk.model.iNTS0$deviance/wk.model.iNTS0$null.deviance)

summary(wk.model.iNTS1)
M1 <- data.frame(wk.dlnmN,link=predict(wk.model.iNTS1, type="link"), fit=predict(wk.model.iNTS1, type="response"), pearson=residuals(wk.model.iNTS1, type="pearson"), resid=residuals(wk.model.iNTS1, type="response"), residSqr=residuals(wk.model.iNTS1, type="response")^2)
v1 <- ggplot(data=M1, aes(x=fit, y=residSqr, color="R^2 = 0.8771")) + scale_colour_manual(name="NTS, model 2",values=c("R^2 = 0.8771"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared1 <- 1-(wk.model.iNTS1$deviance/wk.model.iNTS1$null.deviance)

summary(wk.model.iNTS2)
M2 <- data.frame(wk.dlnmN,link=predict(wk.model.iNTS2, type="link"), fit=predict(wk.model.iNTS2, type="response"), pearson=residuals(wk.model.iNTS2, type="pearson"), resid=residuals(wk.model.iNTS2, type="response"), residSqr=residuals(wk.model.iNTS2, type="response")^2)
v2 <- ggplot(data=M2, aes(x=fit, y=residSqr, color="R^2 = 0.8802")) + scale_colour_manual(name="NTS, model 3",values=c("R^2 = 0.8802"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared2 <- 1-(wk.model.iNTS2$deviance/wk.model.iNTS2$null.deviance)

summary(wk.model.iNTS3)
M3 <- data.frame(wk.dlnmN,link=predict(wk.model.iNTS3, type="link"), fit=predict(wk.model.iNTS3, type="response"), pearson=residuals(wk.model.iNTS3, type="pearson"), resid=residuals(wk.model.iNTS3, type="response"), residSqr=residuals(wk.model.iNTS3, type="response")^2)
v3 <- ggplot(data=M3, aes(x=fit, y=residSqr, color="R^2 = 0.8820")) + scale_colour_manual(name="NTS, model 4",values=c("R^2 = 0.8820"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared3 <- 1-(wk.model.iNTS3$deviance/wk.model.iNTS3$null.deviance)

summary(wk.model.iNTS4)
M4 <- data.frame(wk.dlnmN,link=predict(wk.model.iNTS4, type="link"), fit=predict(wk.model.iNTS4, type="response"), pearson=residuals(wk.model.iNTS4, type="pearson"), resid=residuals(wk.model.iNTS4, type="response"), residSqr=residuals(wk.model.iNTS4, type="response")^2)
v4 <- ggplot(data=M4, aes(x=fit, y=residSqr, color="R^2 = 0.8841")) + scale_colour_manual(name="NTS, model 5",values=c("R^2 = 0.8841"="blue")) + geom_point(color="tan4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
r.squared4 <- 1-(wk.model.iNTS4$deviance/wk.model.iNTS4$null.deviance)

#deviance and r-squared for typhoid
summary(wk.model.typhoid0)
N0 <- data.frame(wk.dlnmT,link=predict(wk.model.typhoid0, type="link"), fit=predict(wk.model.typhoid0, type="response"), pearson=residuals(wk.model.typhoid0, type="pearson"), resid=residuals(wk.model.typhoid0, type="response"), residSqr=residuals(wk.model.typhoid0, type="response")^2)
z0 <- ggplot(data=N0, aes(x=fit, y=residSqr, color="R^2 = 0.2237")) + scale_colour_manual(name="typhoid, model 1",values=c("R^2 = 0.2237"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared0 <- 1-(wk.model.typhoid0$deviance/wk.model.typhoid0$null.deviance)

summary(wk.model.typhoid1)
N1 <- data.frame(wk.dlnmT,link=predict(wk.model.typhoid1, type="link"), fit=predict(wk.model.typhoid1, type="response"), pearson=residuals(wk.model.typhoid1, type="pearson"), resid=residuals(wk.model.typhoid1, type="response"), residSqr=residuals(wk.model.typhoid1, type="response")^2)
z1 <- ggplot(data=N1, aes(x=fit, y=residSqr, color="R^2 = 0.8377")) + scale_colour_manual(name="typhoid, model 2",values=c("R^2 = 0.8377"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared1 <- 1-(wk.model.typhoid1$deviance/wk.model.typhoid1$null.deviance)

summary(wk.model.typhoid2)
N2 <- data.frame(wk.dlnmT,link=predict(wk.model.typhoid2, type="link"), fit=predict(wk.model.typhoid2, type="response"), pearson=residuals(wk.model.typhoid2, type="pearson"), resid=residuals(wk.model.typhoid2, type="response"), residSqr=residuals(wk.model.typhoid2, type="response")^2)
z2 <- ggplot(data=N2, aes(x=fit, y=residSqr, color="R^2 = 0.8398")) + scale_colour_manual(name="typhoid, model 3",values=c("R^2 = 0.8398"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared2 <- 1-(wk.model.typhoid2$deviance/wk.model.typhoid2$null.deviance)

summary(wk.model.typhoid3)
N3 <- data.frame(wk.dlnmT,link=predict(wk.model.typhoid3, type="link"), fit=predict(wk.model.typhoid3, type="response"), pearson=residuals(wk.model.typhoid3, type="pearson"), resid=residuals(wk.model.typhoid3, type="response"), residSqr=residuals(wk.model.typhoid3, type="response")^2)
z3 <- ggplot(data=N3, aes(x=fit, y=residSqr, color="R^2 = 0.8383")) + scale_colour_manual(name="typhoid, model 4",values=c("R^2 = 0.8383"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared3 <- 1-(wk.model.typhoid3$deviance/wk.model.typhoid3$null.deviance)

summary(wk.model.typhoid4)
N4 <- data.frame(wk.dlnmT,link=predict(wk.model.typhoid4, type="link"), fit=predict(wk.model.typhoid4, type="response"), pearson=residuals(wk.model.typhoid4, type="pearson"), resid=residuals(wk.model.typhoid4, type="response"), residSqr=residuals(wk.model.typhoid4, type="response")^2)
z4 <- ggplot(data=N4, aes(x=fit, y=residSqr, color="R^2 = 0.8445")) + scale_colour_manual(name="typhoid, model 5",values=c("R^2 = 0.8445"="blue")) + geom_point(color="red4") + stat_smooth(method="loess", se = FALSE) + theme_classic() + theme(legend.position = c(0.4,0.8))
zr.squared4 <- 1-(wk.model.typhoid4$deviance/wk.model.typhoid4$null.deviance)

#combine plots
grid.arrange(grobs=list(v0,v1, v2, v3, v4, z0,z1, z2, z3, z4), ncol=5, nrow=2)

#nontyphoid relative risk estimates and 95%CIs for potentially significant exposure and lag values (table 2)
cbind(wk.pred.rain.iNTS4$matRRfit, wk.pred.rain.iNTS4$matRRlow, wk.pred.rain.iNTS4$matRRhigh)["15",]
cbind(wk.pred.rain.iNTS4$matRRfit, wk.pred.rain.iNTS4$matRRlow, wk.pred.rain.iNTS4$matRRhigh)["35",]

cbind(wk.pred.temp.iNTS4$matRRfit, wk.pred.temp.iNTS4$matRRlow, wk.pred.temp.iNTS4$matRRhigh)["16",]
cbind(wk.pred.temp.iNTS4$matRRfit, wk.pred.temp.iNTS4$matRRlow, wk.pred.temp.iNTS4$matRRhigh)["34",]

#typhoid relative risk estimates and 95%CIs for potentially significant exposure and lag values (table2)
cbind(wk.pred.rain.typhoid4$matRRfit, wk.pred.rain.typhoid4$matRRlow, wk.pred.rain.typhoid4$matRRhigh)["15",]
cbind(wk.pred.rain.typhoid4$matRRfit, wk.pred.rain.typhoid4$matRRlow, wk.pred.rain.typhoid4$matRRhigh)["35",]

cbind(wk.pred.temp.typhoid4$matRRfit, wk.pred.temp.typhoid4$matRRlow, wk.pred.temp.typhoid4$matRRhigh)["17",]
cbind(wk.pred.temp.typhoid4$matRRfit, wk.pred.temp.typhoid4$matRRlow, wk.pred.temp.typhoid4$matRRhigh)["34",]

#plot 4-dimension between NTS, temperature, rainfall and lag
#convert predictions into a tibble
NTS.matrix.precip <- wk.pred.rain.iNTS4$matRRfit
NTS.matrix.temp <- wk.pred.temp.iNTS4$matRRfit
typhoid.matrix.precip <- wk.pred.rain.typhoid4$matRRfit
typhoid.matrix.temp <- wk.pred.temp.typhoid4$matRRfit
save(NTS.matrix.precip,file="NTS.matrix.precip.Rda")
save(NTS.matrix.temp,file="NTS.matrix.temp.Rda")
save(typhoid.matrix.precip,file="typhoid.matrix.precip.Rda")
save(typhoid.matrix.temp,file="typhoid.matrix.temp.Rda")

NTS.df.rain <-tk_tbl(wk.pred.rain.iNTS4$matRRfit, preserve_index = TRUE, rename_index = "obs")
NTS.df.temp <-tk_tbl(wk.pred.temp.iNTS4$matRRfit, preserve_index = TRUE, rename_index = "obs")
typhoid.df.rain <-tk_tbl(wk.pred.rain.typhoid4$matRRfit, preserve_index = TRUE, rename_index = "obs")
typhoid.df.temp <-tk_tbl(wk.pred.temp.typhoid4$matRRfit, preserve_index = TRUE, rename_index = "obs")

#add new descriptive variable
NTS.df.rain$exposure <- typhoid.df.rain$exposure <- "precipitation"
NTS.df.temp$exposure <- typhoid.df.temp$exposure <- "temperature" 

#combine the tibbles
NTS.df <- rbind(NTS.df.rain, NTS.df.temp)
typhoid.df <- rbind(typhoid.df.rain, typhoid.df.temp)

#reshape the dataset from wide to long
NTS.df <- melt(NTS.df, id=c("obs","exposure"))
setnames(NTS.df, "variable", "lag")
setnames(NTS.df, "value", "predicted")

typhoid.df <- melt(typhoid.df, id=c("obs","exposure"))
setnames(typhoid.df, "variable", "lag")
setnames(typhoid.df, "value", "predicted")

#4D plotting
ggplot(NTS.df, aes(x=obs, y=predicted, color=lag, lty=exposure)) + 
  geom_line(size=1) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + 
  ylab('Relative risk of NTS') + 
  xlab('Precipitation or temperature') + 
  labs(title="(A)") + 
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, 4)) +
  theme(legend.justification=c(0.5,0), legend.position = c(0.3, 0.3), legend.text = element_text(size = 13), legend.title = element_text(size = 16)) 

typhoid.4D <- ggplot(typhoid.df, aes(x=obs, y=predicted, color=lag, lty=exposure)) + 
  geom_line(size=1) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) + 
  ylab('Relative risk of typhoid') + 
  xlab('Precipitation or temperature') + 
  labs(title="(B)") + 
  scale_x_continuous(limits = c(0, 36), breaks = seq(0, 36, 4)) +
  theme(legend.position = 'none')

grid.arrange(grobs=list(NTS.4D,typhoid.4D), ncol=2, nrow=1)

#another 4D attempt
NTS.df.R <- subset(NTS.df, exposure=="precipitation")
NTS.df.T <- subset(NTS.df, exposure=="temperature")
setnames(NTS.df.R, "exposure", "precip")
setnames(NTS.df.T, "exposure", "temp")
setnames(NTS.df.R, "lag", "lag.precip")
setnames(NTS.df.T, "lag", "lag.temp")
setnames(NTS.df.R, "predicted", "predicted.precip")
setnames(NTS.df.T, "predicted", "predicted.temp")
NTS.df.RT <- merge(x = NTS.df.T, y = NTS.df.R)

lag00 <- subset(NTS.df.RT,lag.temp=="lag0" & lag.precip=="lag0")
