#02/08/2018-30/11/2018
#by Deus

#========LOAD AND ATTACH REQUIRED PACKAGES========
#Install required packages
dlnm.analysis.packages <- c("tidyverse","lubridate","xts","ggthemes","PerformanceAnalytics","reshape2","rugarch","timetk","parallel","timeSeries","tseries","data.table","ggplot2","dlnm","broom","caret","gridExtra","splines","splines2","pspline","cowplot","mgcv","spi","chron","gridGraphics","grid","pscl","MASS", "AER", "Hmisc", "MuMIn", "VGAM", "forecast", "seasonal", "plotly", "here" )
blantyre.map.packages <- c("ggmap", "rgdal", "rgeos", "tidyverse", "tmap", "maptools", "maps", "ggfortify")

#load required packages
lapply(dlnm.analysis.packages, library, character.only=TRUE)
lapply(blantyre.map.packages, library, character.only=TRUE)

#========DRAW MAP OF BLANTYRE========
#load shape file of malawi map located in directory "Time.Series/data"
malawi.map<- readOGR(dsn = "Time.Series/data", layer = "malawi_map")

#subsetting to get blantyre map only
blantyre1.map <- malawi.map@data$OBJECTID >289 & malawi.map@data$OBJECTID <297 #id from 290 to 296 
blantyre2.map <- malawi.map@data$OBJECTID >308 & malawi.map@data$OBJECTID <311 #id from 309 to 310
blantyre3.map <- malawi.map@data$OBJECTID >342  #id fom 243

#convert the shape file map into dataframe for ggplotting
blantyre.map1 <- fortify(malawi.map[blantyre1.map,])
blantyre.map2 <- fortify(malawi.map[blantyre2.map,])
blantyre.map3 <- fortify(malawi.map[blantyre3.map,])

blantyre.map <- rbind(blantyre.map1, blantyre.map2, blantyre.map3)
blantyre.map$id <- as.integer(blantyre_map$id)

#merge the blantyre map dataset with location attributes datasets
map.features <- read.csv(here("Time.Series", "data", "blantyre_features.csv"))
blantyre.demog <- read.csv(here("Time.Series", "data", "blantyre_demog.csv"))
blantyre.demog$id <- as.integer(blantyre.demog$id)
map.all <- merge(x=blantyre.map, y=blantyre.demog, by="id", x.all=TRUE)
rm(list = ls()[grep("^blantyre", ls())])

#ggplot the blantyre map with 2008 population census
ggplot() + 
  geom_polygon(data=map.all, aes(x=long, y=lat, group=group, fill=popc), colour="gray50") + 
  theme_classic() + 
  theme(axis.text.x = element_text(face="bold", size=10, color="black"), axis.text.y = element_text(face="bold", size=10, color="black")) + 
  labs(fill="(1998 - 2008) Population censuses") + xlab("Longitude") + ylab("Latitude") + 
  geom_point(data=map.features, aes(x =long, y =lat, shape=Geolocation, size=Geolocation), color="black") +
  scale_shape_manual(values=c(17, 16, 3)) +
  scale_size_manual(values=c(2,4,3)) + 
  theme(legend.key.height=unit(0.8,"line")) + 
  theme(legend.key.width=unit(0.8,"line"))

#========LOAD CASE DATA========
#load typhoid and NTS cases dataset.
case <- read.csv(here("Time.Series", "data", "case.csv"))
case$case_date <- dmy(case$case_date)
case$case_count <- c(1)

#create separate datasets for typhi and NTS cases.
case.typhi <- subset(case, organism == "typhi")
case.iNTS <- subset(case, organism == "iNTS")

#assign 0 to case_count when a date has no typhi case.
case.typhi <-aggregate(case.typhi$case_count, by=list(case.typhi$case_date), FUN=sum, na.rm=TRUE)
names(case.typhi) <- c("date", "case_count")
case.typhi <- merge(case.typhi, data.table(date=seq.Date(min(case$case_date), max(case$case_date), by="day")), by="date", all=TRUE)
case.typhi[is.na(case.typhi)] <- 0

#assign 0 to case_count when dates have no iNTS case.
case.iNTS <-aggregate(case.iNTS$case_count, by=list(case.iNTS$case_date), FUN=sum, na.rm=TRUE)
names(case.iNTS) <- c("date", "case_count")
case.iNTS <- merge(case.iNTS, data.table(date=seq.Date(min(case$case_date), max(case$case_date), by="day")), by="date", all=TRUE)
case.iNTS[is.na(case.iNTS)] <- 0

#convert data frames to xts objects for time series plotting.
case.typhi = as.xts(case.typhi[,-1,drop = FALSE], order.by = as.Date(case.typhi[,1]))
case.iNTS = as.xts(case.iNTS[,-1,drop = FALSE], order.by = as.Date(case.iNTS[,1]))

#monthly aggregates for typhi and iNTS cases.
case.typhi <- apply.monthly(case.typhi, FUN = sum)
case.iNTS <- apply.monthly(case.iNTS, FUN = sum)

#========LOAD CLIMATE DATA========
#load climate dataset.
climate <- read.csv(here("Time.Series", "data", "climate.csv"))

#calculate daily average values for temperature and rainfall.
climate$climate_date <- dmy(climate$date)
climate$rainfall <- (climate$chil_r + climate$chic_r)/2 
climate$temperature <- ((climate$chil_mint + climate$chil_maxt)/2 + (climate$chic_mint + climate$chic_maxt)/2)/2
climate <- subset(climate, select = c(climate_date, rainfall, temperature))

#create separate datasets for daily temperature and rainfall.
climate.rain <- subset(climate, select = c(climate_date, rainfall))
climate.temp <- subset(climate, select = c(climate_date, temperature))

#convert temperature and rainfall data frames to xts objects for use in time series plotting.
climate.rain = as.xts(climate.rain[,-1,drop = FALSE], order.by = as.Date(climate.rain[,1]))
climate.temp = as.xts(climate.temp[,-1,drop = FALSE], order.by = as.Date(climate.temp[,1]))

#monthly aggregates for rainfall and temperature
climate.rain <- apply.monthly(climate.rain, FUN = mean)
climate.temp <- apply.monthly(climate.temp, FUN = mean)

#=========SEASONALITY AND DEMOGRAPHICS PLOTS========
#create tibbles for adjusting seasonality
case.typhi <-tk_tbl(case.typhi, preserve_index = TRUE, rename_index = "date") 
case.iNTS <-tk_tbl(case.iNTS, preserve_index = TRUE, rename_index = "date") 
climate.rain <-tk_tbl(climate.rain, preserve_index = TRUE, rename_index = "date") 
climate.temp <-tk_tbl(climate.temp, preserve_index = TRUE, rename_index = "date") 

#seasonally-adjusted cases and climate
case.iNTS.ts <- ts(na.omit(case.iNTS$case_count), frequency = 12)
trend_n <-tk_tbl(10^(seasadj(mstl(log10(case.iNTS.ts)))), preserve_index = FALSE) #multiplicative series (log-transform)
case.iNTS$case_count_sea <- trend_n$value
case.iNTS$case_count_sea[case.iNTS$case_count_sea < 0] <- 0
setnames(case.iNTS, old="case_count", new="case_count_obs")

case.typhi.ts <- ts(na.omit(case.typhi$case_count), frequency = 12)
trend_n <-tk_tbl(10^(seasadj(mstl(log10(case.typhi.ts+1)))), preserve_index = FALSE) #multiplicative series (log-transform)
case.typhi$case_count_sea <- trend_n$value
case.typhi$case_count_sea[case.typhi$case_count_sea < 0] <- 0
setnames(case.typhi, old="case_count", new="case_count_obs")

climate.rain.ts <- ts(na.omit(climate.rain$rainfall), frequency = 12)
trend_n <-tk_tbl(abs(seasadj(mstl(climate.rain.ts))), preserve_index = FALSE) #additive series
climate.rain$rainfall_sea <- trend_n$value
climate.rain$rainfall_sea[climate.rain$rainfall_sea < 0] <- 0
setnames(climate.rain, old="rainfall", new="rainfall_obs")

climate.temp.ts <- ts(na.omit(climate.temp$temperature), frequency = 12)
trend_n <-tk_tbl(abs(seasadj(mstl(climate.temp.ts))), preserve_index = FALSE) #additive series
climate.temp$temperature_sea <- trend_n$value
climate.temp$temperature_sea[climate.temp$temperature_sea < 0] <- 0
setnames(climate.temp, old="temperature", new="temperature_obs")

#plot decomposed all series
x<-case.iNTS.ts %>% mstl() %>% ggfortify:::autoplot.ts(main="A",xlab="Years (2000-2015)",size=1,colour="orange2",is.date=FALSE) + theme_bw()
y<-case.typhi.ts %>% mstl() %>% ggfortify:::autoplot.ts(main="B",xlab="Years (2000-2015)",size=1,colour="red2",is.date=FALSE) + theme_bw()
z<-climate.rain.ts %>% mstl() %>% ggfortify:::autoplot.ts(main="C",xlab="Years (2000-2015)",size=1,colour="blue2",is.date=FALSE) + theme_bw()
v<-climate.temp.ts %>% mstl() %>% ggfortify:::autoplot.ts(main="D",xlab="Years (2000-2015)",size=1,colour="green2",is.date=FALSE) + theme_bw()
grid.arrange(grobs=list(x,y,z,v), ncol=4, nrow=1)

#plot all seasonally-adjusted series
S1<-ggplot(as.data.frame(case.iNTS)) + 
  geom_line(aes(date, case_count_obs, color="Original"), size=0.8) + 
  geom_line(aes(date, case_count_sea, color="Seasonally-adjusted"), size=0.8) + 
  scale_color_manual(values = c("Original"="black","Seasonally-adjusted"="orange2")) +
  labs(title="A", x ="", y = "iNTS cases") + 
  theme(plot.title = element_text(hjust = 0)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.7,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

S2<-ggplot(as.data.frame(case.typhi)) + 
  geom_line(aes(date, case_count_obs, color="Original"), size=0.8) + 
  geom_line(aes(date, case_count_sea, color="Seasonally-adjusted"), size=0.8) + 
  scale_color_manual(values = c("Original"="black","Seasonally-adjusted"="red2")) +
  labs(title="B", x ="", y = "Typhoid cases") + 
  theme(plot.title = element_text(hjust = 0)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.3,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

S3<-ggplot(as.data.frame(climate.rain)) + 
  geom_line(aes(date, rainfall_obs, color="Original"), size=0.8) + 
  geom_line(aes(date, rainfall_sea, color="Seasonally-adjusted"), size=0.8) + 
  scale_color_manual(values = c("Original"="black","Seasonally-adjusted"="blue2")) +
  labs(title="C", x ="Year", y = "Rainfall (mm)") + 
  theme(plot.title = element_text(hjust = 0)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.3,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

S4<-ggplot(as.data.frame(climate.temp)) + 
  geom_line(aes(date, temperature_obs, color="Original"), size=0.8) + 
  geom_line(aes(date, temperature_sea, color="Seasonally-adjusted"), size=0.8) + 
  scale_color_manual(values = c("Original"="black","Seasonally-adjusted"="green2")) +
  labs(title="D", x ="Year", y = "Temperature (°C)") + 
  theme(plot.title = element_text(hjust = 0)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.3,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

grid.arrange(grobs=list(S1, S2, S3, S4), ncol=2, nrow=2)
rm(list = ls()[grep("^trend_n", ls())])

#linear interpolation and extrapolation
census.year <-c(1998, 2008)
census.popn <-c(809397, 1022680)
census.count.1998.2008 <- approx(census.year, census.popn, n=11) 
census.count.2009.2015 <- approxExtrap(census.year, census.popn, xout=c(2009, 2010, 2011, 2012, 2013, 2014, 2015))

#incidence rates of NTS using intra+extrapolated denominators
case.iNTS$census[year(case.iNTS$date) == 2000] <-852054
case.iNTS$census[year(case.iNTS$date) == 2001] <-873382
case.iNTS$census[year(case.iNTS$date) == 2002] <-894710
case.iNTS$census[year(case.iNTS$date) == 2003] <-916039
case.iNTS$census[year(case.iNTS$date) == 2004] <-937367
case.iNTS$census[year(case.iNTS$date) == 2005] <-958695
case.iNTS$census[year(case.iNTS$date) == 2006] <-980023
case.iNTS$census[year(case.iNTS$date) == 2007] <-1001352
case.iNTS$census[year(case.iNTS$date) == 2008] <-1022680
case.iNTS$census[year(case.iNTS$date) == 2009] <-1044008
case.iNTS$census[year(case.iNTS$date) == 2010] <-1065337
case.iNTS$census[year(case.iNTS$date) == 2011] <-1086665
case.iNTS$census[year(case.iNTS$date) == 2012] <-1107993
case.iNTS$census[year(case.iNTS$date) == 2013] <-1129322
case.iNTS$census[year(case.iNTS$date) == 2014] <-1150650
case.iNTS$census[year(case.iNTS$date) == 2015] <-1171978
case.iNTS$incid_sea <-case.iNTS$case_count_sea*100000/case.iNTS$census 
case.iNTS$incid_obs <-case.iNTS$case_count_obs*100000/case.iNTS$census 

case.typhi$census[year(case.typhi$date) == 2000] <-852054
case.typhi$census[year(case.typhi$date) == 2001] <-873382
case.typhi$census[year(case.typhi$date) == 2002] <-894710
case.typhi$census[year(case.typhi$date) == 2003] <-916039
case.typhi$census[year(case.typhi$date) == 2004] <-937367
case.typhi$census[year(case.typhi$date) == 2005] <-958695
case.typhi$census[year(case.typhi$date) == 2006] <-980023
case.typhi$census[year(case.typhi$date) == 2007] <-1001352
case.typhi$census[year(case.typhi$date) == 2008] <-1022680
case.typhi$census[year(case.typhi$date) == 2009] <-1044008
case.typhi$census[year(case.typhi$date) == 2010] <-1065337
case.typhi$census[year(case.typhi$date) == 2011] <-1086665
case.typhi$census[year(case.typhi$date) == 2012] <-1107993
case.typhi$census[year(case.typhi$date) == 2013] <-1129322
case.typhi$census[year(case.typhi$date) == 2014] <-1150650
case.typhi$census[year(case.typhi$date) == 2015] <-1171978
case.typhi$incid_sea <-case.typhi$case_count_sea*100000/case.typhi$census 
case.typhi$incid_obs <-case.typhi$case_count_obs*100000/case.typhi$census 

#monthly plots of adjusted NTS, typhoid cases and climate
E1<-ggplot() + 
  geom_line(data = case.iNTS, aes(x=date, y=case_count_sea, color=" Number of iNTS cases"), stat = "identity", size = 1.0) + 
  geom_line(data = climate.rain, aes(x=date, y=rainfall_obs/0.1, color="Rainfall"), alpha=0.8, size = 0.7) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.1, name = "(mm)")) + 
  theme_bw() +
  scale_color_manual(values = c(" Number of iNTS cases"="orange2","Rainfall"="blue2")) + 
  ggtitle("A") + ylab("Cases") + xlab("Month'Year") + 
  theme(axis.title.x = element_text(size=0,face="bold"), axis.title.y = element_text(size=10, color="orange2",face="bold"), plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.text.x=element_text(face="bold", size=0), axis.title.y.right = element_text(color="blue2", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  scale_x_date(date_breaks = "24 month", date_labels ="%b'%y") +
  theme(legend.justification=c(0.5,0), legend.position = c(0.6,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) +
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

E2<-ggplot() + 
  geom_line(data = case.iNTS, aes(x=date, y=case_count_sea, color=" Number of iNTS cases"), stat = "identity", size = 1.0) + 
  geom_line(data = climate.temp, aes(x=date, y=temperature_obs/0.3, color="Temperature"), alpha=0.8, size = 0.7) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.3, name = "(°C)")) + 
  theme_bw() +
  scale_color_manual(values = c(" Number of iNTS cases"="orange2","Temperature"="green2")) + 
  ggtitle("B") + ylab("Cases") + xlab("Month'Year") + 
  theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size=10, color="orange2",face="bold"), plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.text.x=element_text(face="bold", size=8), axis.title.y.right = element_text(color="green2", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  scale_x_date(date_breaks = "24 month", date_labels ="%b'%y") +
  theme(legend.justification=c(0.5,0), legend.position = c(0.6,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) +
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

E3<-ggplot() + 
  geom_line(data = case.typhi, aes(x=date, y=case_count_sea, color=" Number of typhoid cases"), stat = "identity", size = 1.0) + 
  geom_line(data = climate.rain, aes(x=date, y=rainfall_obs/0.1, color="Rainfall"), alpha=0.8, size = 0.7) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.1, name = "(mm)")) +
  theme_bw() +
  scale_color_manual(values = c(" Number of typhoid cases"="red2","Rainfall"="blue2")) + 
  ggtitle("C") + ylab("Cases") + xlab("Month'Year") + 
  theme(axis.title.x = element_text(size=0,face="bold"), axis.title.y = element_text(size=10, color="red2",face="bold"), plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.text.x=element_text(face="bold", size=0), axis.title.y.right = element_text(color="blue2", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  scale_x_date(date_breaks = "24 month", date_labels ="%b'%y") +
  theme(legend.justification=c(0.5,0), legend.position = c(0.6,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) +
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

E4<-ggplot() + 
  geom_line(data = case.typhi, aes(x=date, y=case_count_sea, color=" Number of typhoid cases"), stat = "identity", size = 1.0) + 
  geom_line(data = climate.rain, aes(x=date, y=rainfall_obs/0.1, color="Temperature"), alpha=0.8, size = 0.7) + 
  scale_y_continuous(sec.axis = sec_axis(~.*0.1, name = "(°C)")) + 
  theme_bw() +
  scale_color_manual(values = c(" Number of typhoid cases"="red2","Temperature"="green2")) + 
  ggtitle("D") + ylab("Cases") + xlab("Month'Year") + 
  theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size=10, color="red2",face="bold"), plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.text.x=element_text(face="bold", size=8), axis.title.y.right = element_text(color="green2", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  scale_x_date(date_breaks = "24 month", date_labels ="%b'%y") +
  theme(legend.justification=c(0.5,0), legend.position = c(0.6,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) +
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

grid.arrange(grobs=list(E1, E3, E2, E4), ncol=2, nrow=2)

#monthly dynamics of NTS (11 years), typhoid (5 years) and climate as boxplots
monthNo <- (as.factor(months(subset(case.iNTS$date, year(case.iNTS$date)<2011),abbr=TRUE)))
monthTy <- (as.factor(months(subset(case.typhi$date, year(case.typhi$date)>2010),abbr=TRUE)))
monthTe <- (as.factor(months(subset(climate.temp$date,year(climate.temp$date)<2011),abbr=TRUE)))
monthTey <- (as.factor(months(subset(climate.temp$date,year(climate.temp$date)>2010),abbr=TRUE)))
monthRa <- (as.factor(months(subset(climate.rain$date,year(climate.rain$date)<2011),abbr=TRUE)))
monthRay <- (as.factor(months(subset(climate.rain$date,year(climate.rain$date)>2010),abbr=TRUE)))

pox1nts <- ggplot(subset(climate.temp, year(date)<2011), aes(x=monthTe, y=temperature_obs)) + 
  geom_boxplot(aes(group=monthTe), outlier.shape = NA, color="black", fill="green2", alpha=0.7) + 
  labs(title="D",x="Month", y = "Temperature (°C)") + 
  scale_x_discrete(limits = month.abb) + 
  theme(plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=10)) 

pox2nts <- ggplot(subset(climate.rain, year(date)<2011), aes(x=monthRa, y=rainfall_obs)) + 
  geom_boxplot(aes(group=monthRa), outlier.shape = NA, color="black", fill="blue2", alpha=0.7) + 
  labs(title="C",x="Month", y = "Rainfall (mm)") + 
  scale_x_discrete(limits = month.abb) + 
  theme(plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=10)) 

pox1typ <- ggplot(subset(climate.temp, year(date)>2010), aes(x=monthTey, y=temperature_obs)) + 
  geom_boxplot(aes(group=monthTey), outlier.shape = NA, color="black", fill="green2", alpha=0.7) + 
  labs(title="D",x="Month", y = "Temperature (°C)") + 
  scale_x_discrete(limits = month.abb) + 
  theme(plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=10)) 

pox2typ <- ggplot(subset(climate.rain, year(date)>2010), aes(x=monthRay, y=rainfall_obs)) + 
  geom_boxplot(aes(group=monthRay), outlier.shape = NA, color="black", fill="blue2", alpha=0.7) + 
  labs(title="C",x="Month", y = "Rainfall (mm)") + 
  scale_x_discrete(limits = month.abb) + 
  theme(plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=10)) 

pox3 <- ggplot(subset(case.iNTS, year(date)<2011), aes(x=monthNo, y=incid_sea))  + 
  geom_boxplot(aes(group=monthNo), outlier.shape = NA, color="black", fill="orange2", alpha=0.7) + 
  labs(title="A",x="", y = "NTS incidence rate") + 
  scale_x_discrete(limits = month.abb) + 
  theme(plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=10)) 

pox4 <- ggplot(subset(case.typhi, year(date)>2010), aes(x=monthTy, y=incid_sea))  + 
  geom_boxplot(aes( group=monthTy), outlier.shape = NA, color="black", fill="red2", alpha=0.7) + 
  labs(title="B", x="", y = "Typhoid incidence rate") + 
  scale_x_discrete(limits = month.abb) + 
  theme(plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=10)) 

grid.arrange(grobs=list(pox3,pox3, pox4,pox4, pox2nts, pox1nts,pox2typ, pox1typ), ncol=4, nrow=2)

#distributions of typhi and NTS cases by sex and age.
case$sex[case$sex == ""] <- NA
case$age[case$age == ""] <- NA
case$date <- ymd(case$case_date)
case$year <- year(case$date)

agesex.p1 <- ggplot(subset(case, organism=="iNTS" & !is.na(age) & sex != "Unknown"), aes(x=age, color=sex)) + 
  geom_freqpoly(position=position_dodge(width=1.5), binwidth=1, size=1) + 
  scale_color_manual(values=c(Male="gray40",Female="red3")) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0, 91), breaks = seq(0, 91, 5)) + 
  labs(title="A", x ="Age (years)", y = "Number of iNTS cases") + 
  theme(plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.text.x = element_text(face="bold", size=10, color="black"), axis.text.y = element_text(face="bold", size=10, color="black")) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.5, 0.6), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) + 
  labs(color="Missing sex: 2,596 (32.3%)") + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

agesex.p2<- ggplot(subset(case, organism=="typhi" & !is.na(age) & sex != "Unknown"), aes(x=age, color=sex)) + 
  geom_freqpoly(position=position_dodge(width=1.5), binwidth=1, size=1) + 
  scale_color_manual(values=c(Male="gray40",Female="red3")) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0, 91), breaks = seq(0, 91, 5)) + 
  labs(title="B", x ="Age (years)", y = "Number of typhoid cases") + 
  theme(plot.title = element_text(hjust=0,face="bold")) +
  theme(axis.text.x = element_text(face="bold", size=10, color="black"), axis.text.y = element_text(face="bold", size=10, color="black")) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.5, 0.6), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) + 
  labs(color="Missing sex: 61 (2.4%)") + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

grid.arrange(grobs=list(agesex.p1, agesex.p2), ncol=2, nrow=1)

#contour plots of seasonal dynamics for evey year
case.iNTS.spi <- subset(case.iNTS, year(case.iNTS$date)<2011, select=c(date,incid_sea)) 
case.iNTS.spi$month <- month(case.iNTS.spi$date)
case.iNTS.spi$year <- year(case.iNTS.spi$date)
case.iNTS.spi$date <- NULL
case.iNTS.spi <- spread(case.iNTS.spi, year, incid_sea)
case.iNTS.spi <- as.matrix(case.iNTS.spi)[,-1]
cnts <- plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~case.iNTS.spi, type = "contour", colorscale = 'jet', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "iNTS incidence per \n 100,000 population") %>%
layout(title="A", xaxis=list(title ="Year"), yaxis=list(title="Month"), font=list(size = 13))

case.typhi.spi <- subset(case.typhi, year(case.typhi$date)>2010, select=c(date,incid_sea)) 
case.typhi.spi$month <- month(case.typhi.spi$date)
case.typhi.spi$year <- year(case.typhi.spi$date)
case.typhi.spi$date <- NULL
case.typhi.spi <- spread(case.typhi.spi, year, incid_sea)
case.typhi.spi <- as.matrix(case.typhi.spi)[,-1]
ctyp <- plot_ly(x = c(2011,2012,2013,2014,2015), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~case.typhi.spi, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Typhoid incidence per \n 100,000 population") %>%
layout(title="D", xaxis=list(title ="Year"), yaxis=list(title="Month"), font=list(size = 13))

climate.rain.spin <- subset(climate.rain, year(climate.rain$date)<2011, select=c(date,rainfall_obs)) 
climate.rain.spin$month <- month(climate.rain.spin$date)
climate.rain.spin$year <- year(climate.rain.spin$date)
climate.rain.spin$date <- NULL
climate.rain.spin <- spread(climate.rain.spin, year, rainfall_obs)
climate.rain.spin <- as.matrix(climate.rain.spin)[,-1]
cntsrain <- plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~climate.rain.spin, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Rainfall (mm)") %>%
layout(title="B", xaxis=list(title ="Year"), yaxis=list(title="Month"), font=list(size = 13))

climate.rain.spit <- subset(climate.rain, year(climate.rain$date)>2010, select=c(date,rainfall_obs)) 
climate.rain.spit$month <- month(climate.rain.spit$date)
climate.rain.spit$year <- year(climate.rain.spit$date)
climate.rain.spit$date <- NULL
climate.rain.spit <- spread(climate.rain.spit, year, rainfall_obs)
climate.rain.spit <- as.matrix(climate.rain.spit)[,-1]
ctyprain <- plot_ly(x = c(2011,2012,2013,2014,2015), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~climate.rain.spit, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Rainfall (mm)") %>%
layout(title="E", xaxis=list(title ="Year"), yaxis=list(title="Month"), font=list(size = 13))

climate.temp.spin <- subset(climate.temp, year(climate.temp$date)<2011, select=c(date,temperature_obs)) 
climate.temp.spin$month <- month(climate.temp.spin$date)
climate.temp.spin$year <- year(climate.temp.spin$date)
climate.temp.spin$date <- NULL
climate.temp.spin <- spread(climate.temp.spin, year, temperature_obs)
climate.temp.spin <- as.matrix(climate.temp.spin)[,-1]
cntstemp <- plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~climate.temp.spin, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Temperature (°C)") %>%
layout(title="C", xaxis=list(title ="Year"), yaxis=list(title="Month"), font=list(size = 13))

climate.temp.spit <- subset(climate.temp, year(climate.temp$date)>2010, select=c(date,temperature_obs)) 
climate.temp.spit$month <- month(climate.temp.spit$date)
climate.temp.spit$year <- year(climate.temp.spit$date)
climate.temp.spit$date <- NULL
climate.temp.spit <- spread(climate.temp.spit, year, temperature_obs)
climate.temp.spit <- as.matrix(climate.temp.spit)[,-1]
ctyptemp <- plot_ly(x = c(2011,2012,2013,2014,2015), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~climate.temp.spit, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Temperature (°C)") %>%
layout(title="F", xaxis=list(title ="Year"), yaxis=list(title="Month"), font=list(size = 13))

#========CROSS-BASIS MODEL SELECTION FOR NTS OUTCOME========

#prepare weekly dlnm dataset
mo.dlnmN <- bind_cols(case.iNTS, climate.rain, climate.temp, id=NULL)
mo.dlnmN$date1 <- mo.dlnmN$date2 <- NULL
mo.dlnmN <- subset(mo.dlnmN, year(date) < 2011)
mo.dlnmN$time <- seq.int(from = 1, to=132, by=1)
mo.dlnmN$year <- year(mo.dlnmN$date)
mo.dlnmN$month <- month(mo.dlnmN$date)
mo.dlnmN$incid_seaX<-round(mo.dlnmN$incid_sea, digits = 0)

#manipulate the AIC so its QAIC and compare to choose an optimal model
n.quasipoisson <- function(...) { 
  res <- quasipoisson(...)
  res$aic <- poisson(...)$aic 
  res
}

#defines all possible dfs for rainfall and lag
lagknots <- logknots(8, fun="ns", df=3)
varknots.r=equalknots(mo.dlnmN$rainfall_obs, fun = "ns", df=3)
varknots.t=equalknots(mo.dlnmN$temperature_obs, fun = "ns", df=3)
mo.cb.rain.nts <- crossbasis(mo.dlnmN$rainfall_obs, lag =8, argvar = list(fun="ns", knots=varknots.r), arglag = list(knots=lagknots))
mo.cb.temp.nts <- crossbasis(mo.dlnmN$temperature_obs, lag =8, argvar = list(fun="ns", knots=varknots.t), arglag = list(knots=lagknots))

nts.model <- glm(mo.dlnmN$incid_seaX ~ mo.cb.rain.nts + mo.cb.temp.nts + ns(time,11), family = n.quasipoisson(), na.action=na.delete, mo.dlnmN)
QAIC(nts.model, chat=summary(nts.model)$dispersion)

#========NTS DLNM MULTIVARIATE ANALYSIS========

#cross-basis for monthly rainfall, temperature 
sort(mo.dlnmN$rainfall_obs, decreasing = FALSE)
varknots=equalknots(mo.dlnmN$rainfall_obs, fun = "ns", df=3)
lagknots <- logknots(8, fun = "ns", df=3)
mo.cb.rain.iNTS <- crossbasis(mo.dlnmN$rainfall_obs, lag =8, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.rain.iNTS)

sort(mo.dlnmN$temperature_obs, decreasing = FALSE)
varknots=equalknots(mo.dlnmN$temperature_obs, fun = "ns", df=3)
lagknots <- logknots(8, fun="ns", df=3)
mo.cb.temp.iNTS <- crossbasis(mo.dlnmN$temperature_obs, lag =8, argvar = list(fun="ns", knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.temp.iNTS)

#models fitting for monthly rainfall, temperature and NTS (either quasipossion)
mo.model.iNTS <- glm(incid_sea ~  mo.cb.rain.iNTS + mo.cb.temp.iNTS + ns(time,11), family = quasipoisson(), na.action=na.exclude, mo.dlnmN)
mo.model.iNTS1 <- glm(incid_sea ~  mo.cb.rain.iNTS + ns(time), family = quasipoisson(), na.action=na.exclude, mo.dlnmN)
mo.model.iNTS2 <- glm(incid_sea ~  mo.cb.temp.iNTS + ns(time), family = quasipoisson(), na.action=na.exclude, mo.dlnmN)

#prediction for monthly rainfall, temperature
mo.pred.rain.iNTS <- crosspred(mo.cb.rain.iNTS, mo.model.iNTS, cen = 0, by=0.2)
mo.pred.temp.iNTS <- crosspred(mo.cb.temp.iNTS, mo.model.iNTS, cen = 23, by=0.2)

#countour, 3D plots and lag/predictor curves for monthly rainfall
dev.off()
plot(mo.pred.rain.iNTS, xlab="Rainfall (mm)", ylab="Lag (month)", zlab="RR of NTS", theta=150, phi=5, lphi=100, main="A", cex.lab=1, cex.axis=1, col="gray80")
plot(mo.pred.rain.iNTS, "contour", key.title=title("NTS"), plot.title=title("(B)", xlab ="Rainfall (mm)", ylab = "Lag (month)", cex.lab=1, cex.axis=1))
plot(mo.pred.rain.iNTS, "slices", var=c(4,12), lag=c(0,3), col="blue2", ci.arg=list(col=topo.colors(20, alpha = 0.2)), ci.level=0.95, lwd=3.5,  col="gray0", ylab="RR of NTS", cex.lab=1, cex.axis=1)
plot(mo.pred.rain.iNTS, "overall", main="Overall association for first 4 months", xlab="Rainfall (mm)", ylab="RR of NTS", col="blue2", ci.arg=list(col=topo.colors(20, alpha = 1)), ci.level=0.95, ci='l', lwd=3.5, cex.lab=1, cex.axis=1)

#countour, 3D plots and lag/predictor curves for monthly temperature
dev.off()
plot(mo.pred.temp.iNTS, xlab="Temperature (°C)", ylab="Lag (month)", zlab="RR of NTS", theta=150, phi=5, lphi=100, main="(A)", cex.lab=1, cex.axis=1, col="gray80")
plot(mo.pred.temp.iNTS, "contour", key.title=title("NTS"), plot.title=title("(B)", xlab ="Temperature (°C)", ylab = "Lag (month)", cex.lab=1, cex.axis=1))
plot(mo.pred.temp.iNTS, "slices", var=c(18,30), lag=c(2,4), col="green4", ci.arg=list(col=terrain.colors(16, alpha = 0.2)), ci.level=0.95, lwd=3.5,  col="gray0", ylab="RR of NTS", cex.lab=1, cex.axis=1)
plot(mo.pred.temp.iNTS, "overall", main="Overall association for 8 months", xlab="Temperature (°C)", ylab="RR of NTS", col="green4", ci.arg=list(col=terrain.colors(20, alpha = 1)), ci.level=0.95, ci='l', lwd=3.5, cex.lab=1, cex.axis=1)

#model diagnostics diagnostics
pacf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="Autocorrelation (adjusted model)",xlim=c(0,8))
mo.model.iNTS <- update(mo.model.iNTS,.~.+Lag(residuals(mo.model.iNTS,type="deviance"),1))
mo.model.iNTS <- update(mo.model.iNTS,.~.+Lag(residuals(mo.model.iNTS,type="deviance"),5))
mo.model.iNTS <- update(mo.model.iNTS,.~.+Lag(residuals(mo.model.iNTS,type="deviance"),6))
pacf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="Autocorrelation (adjusted model)",xlim=c(0,8))
dev.off()

#model diagnostic plots
par(mfrow=c(2,2))
plot(mo.dlnmN$date, residuals(mo.model.iNTS,type="deviance"), pch=19, cex=0.4, col=grey(0.6),main="Deviance residuals", ylab="Residual deviance", xlab="Year") 
abline(h=0,lty=2,lwd=2)
pacf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="Autocorrelation (adjusted model)",xlim=c(0,8))
acf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="Correlation (adjusted model)",xlim=c(0,8))
plot(mo.dlnmN$incid_sea, pch=10, cex=0.4, col=grey(0.6), size =1 ,main="Observed v fitted", ylab="iNTS incidence rate", xlab="Month")
lines(predict(mo.model.iNTS, type="response"), col="red")
legend("topright", legend=c("observed", "fitted"), col=c("grey", "red"), lty=1:1, cex=0.8, lwd=1)
dev.off()

#model2 diagnostics diagnostics
pacf(residuals(mo.model.iNTS2,type="deviance"),na.action=na.omit,main="Autocorrelation (adjusted model)",xlim=c(0,8))
mo.model.iNTS2 <- update(mo.model.iNTS2,.~.+Lag(residuals(mo.model.iNTS2,type="deviance"),1))
mo.model.iNTS2 <- update(mo.model.iNTS2,.~.+Lag(residuals(mo.model.iNTS2,type="deviance"),5))
mo.model.iNTS2 <- update(mo.model.iNTS2,.~.+Lag(residuals(mo.model.iNTS2,type="deviance"),8))
pacf(residuals(mo.model.iNTS2,type="deviance"),na.action=na.omit,main="Autocorrelation (adjusted model)",xlim=c(0,8))
dev.off()

par(mfrow=c(2,2))
plot(mo.dlnmN$date, residuals(mo.model.iNTS,type="deviance"), pch=19, cex=0.4, col=grey(0.6),main="Deviance residuals", ylab="Residual deviance", xlab="Year") 
abline(h=0,lty=2,lwd=2)
pacf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))
acf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))
plot(mo.dlnmN$incid_sea, pch=10, cex=0.4, col=grey(0.6), size =1 ,main="Original model prediction", ylab="NTS incidence rate", xlab="Month")
lines(predict(mo.model.iNTS, type="response"), col="red")
legend("topright", legend=c("observed", "predicted"), col=c("grey", "red"), lty=1:1, cex=0.9, lwd=1)
dev.off()

#========CROSS-BASIS MODEL SELECTION WITH TYPHOID FEVER OUTCOME========

#prepare weekly typhoid dlnm dataset
mo.dlnmT <- bind_cols(case.typhi, climate.rain, climate.temp, id=NULL)
mo.dlnmT$date1 <- mo.dlnmT$date2 <- NULL
mo.dlnmT <- subset(mo.dlnmT, year(date) > 2010)
mo.dlnmT$time <- seq.int(from = 1, to=60, by=1)
mo.dlnmT$year <- year(mo.dlnmT$date)
mo.dlnmT$month <- month(mo.dlnmT$date)
mo.dlnmT$incid_seaX<-round(mo.dlnmT$incid_sea, digits = 0)

#manipulate the AIC so its QAIC and compare to choose an optimal model
t.quasipoisson <- function(...) { 
  res <- quasipoisson(...) 
  res$aic <- poisson(...)$aic 
  res
}
#defines all possible dfs for rainfall and lag
lagknots <- logknots(8, fun="ns", df=5)
varknots.r=equalknots(mo.dlnmT$rainfall_obs, fun = "ns", df=5)
varknots.t=equalknots(mo.dlnmT$temperature_obs, fun = "ns", df=5)
mo.cb.rain.typhoid <- crossbasis(mo.dlnmT$rainfall_obs, lag =8, argvar = list(fun="ns", knots=varknots.r), arglag = list(knots=lagknots))
mo.cb.temp.typhoid <- crossbasis(mo.dlnmT$temperature_obs, lag =8, argvar = list(fun="ns", knots=varknots.t), arglag = list(knots=lagknots))

typhoid.model <- glm(mo.dlnmT$incid_seaX ~ mo.cb.temp.typhoid + ns(time,5), family = t.quasipoisson(), na.action=na.delete, mo.dlnmT)
QAIC(typhoid.model, chat=summary(typhoid.model)$dispersion)

#========TYPHOID DLNM MULTIVARIATE ANALYSIS========

#cross basis for weekly rainfall, temperature 
sort(mo.dlnmT$rainfall_obs, decreasing = FALSE)
varknots=equalknots(mo.dlnmT$rainfall_obs, fun = "ns", df=3)
lagknots <- logknots(8, df=3)
mo.cb.rain.typhoid <- crossbasis(mo.dlnmT$rainfall_obs, lag =8, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.rain.typhoid)

sort(mo.dlnmT$temperature_obs, decreasing = FALSE)
varknots=equalknots(mo.dlnmT$temperature_obs, fun = "ns", df=3)
lagknots <- logknots(8, df=3)
mo.cb.temp.typhoid <- crossbasis(mo.dlnmT$temperature_obs, lag =8, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.temp.typhoid)

#models fitting for weekly rainfall, temperature and typhoid
mo.model.typhoid <- glm(incid_seaX ~ mo.cb.rain.typhoid + mo.cb.temp.typhoid + ns(time,5), family = quasipoisson(), na.action=na.exclude, mo.dlnmT)
mo.model.typhoid1 <- glm(incid_seaX ~ mo.cb.rain.typhoid + ns(time,5), family = quasipoisson(), na.action=na.exclude, mo.dlnmT)
mo.model.typhoid2 <- glm(incid_seaX ~ mo.cb.temp.typhoid + ns(time,5), family = quasipoisson(), na.action=na.exclude, mo.dlnmT)

#prediction for weekly rainfall, temperature
mo.pred.rain.typhoid <- crosspred(mo.cb.rain.typhoid, mo.model.typhoid1, cen = 0, by=0.2)
mo.pred.temp.typhoid <- crosspred(mo.cb.temp.typhoid, mo.model.typhoid2, cen = 23.3, by=0.2)

#countour, 3D plots and lag/predictor curves for weekly rainfall
dev.off()
plot(mo.pred.rain.typhoid, xlab="Rainfall (mm)", ylab="Lag (month)", zlab="RR of typhoid", theta=150, phi=5, lphi=100, main="(A)", cex.lab=1, cex.axis=1)
plot(mo.pred.rain.typhoid, "contour", key.title=title("TYP"), plot.title=title("(B)", xlab ="Rainfall (mm)", ylab = "Lag (month)", cex.lab=1, cex.axis=1))
plot(mo.pred.rain.typhoid, "slices", var=c(8,10), lag=c(1,4), col="green4", ci.arg=list(col=terrain.colors(20, alpha = 0.5)), ci.level=0.95, lwd=3.5,  col="gray0", ylab="RR of typhoid", cex.lab=1, cex.axis=1)
plot(mo.pred.rain.typhoid, "overall", main="Overall association for 8 lags", xlab="Rainfall (mm)", ylab="RR of typhoid", col="green4", ci.arg=list(col=terrain.colors(20, alpha = 1)), ci.level=0.95, ci='l', lwd=3.5, cex.lab=1, cex.axis=1)

#3D, countour, lag/predictor curves, and cumulative associations for monthly temperature
dev.off()
plot(mo.pred.temp.typhoid, xlab="Temperature (°C)", ylab="Lag (month)", zlab="RR of typhoid", theta=150, phi=5, lphi=100, main="(A)", cex.lab=1, cex.axis=1)
plot(mo.pred.temp.typhoid, "contour", key.title=title("TYP"), plot.title=title("(B)", xlab ="Temperature (°C)", ylab = "Lag (month)", cex.lab=1, cex.axis=1))
plot(mo.pred.temp.typhoid, "slices", var=c(22,30), lag=c(2,8), col="red4", ci.arg=list(col=heat.colors(16, alpha = 0.5)), ci.level=0.95, lwd=3.5,  col="gray0", ylab="RR of typhoid", cex.lab=1, cex.axis=1)
plot(mo.pred.temp.typhoid, "overall", main="Overall association for 8 lags", xlab="Temperature (°C)", ylab="RR of typhoid", col="red4", ci.arg=list(col=heat.colors(20, alpha = 1)), ci.level=0.95, ci='l', lwd=3.5, cex.lab=1, cex.axis=1)

#model diagnostics
pacf(residuals(mo.model.typhoid,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))
mo.model.typhoid <- update(mo.model.typhoid,.~.+Lag(residuals(mo.model.typhoid,type="deviance"),8))
mo.model.typhoid <- update(mo.model.typhoid,.~.+Lag(residuals(mo.model.typhoid,type="deviance"),1))
mo.model.typhoid <- update(mo.model.typhoid,.~.+Lag(residuals(mo.model.typhoid,type="deviance"),2))
mo.model.typhoid <- update(mo.model.typhoid,.~.+Lag(residuals(mo.model.typhoid,type="deviance"),6))
pacf(residuals(mo.model.typhoid,type="deviance"),na.action=na.omit,main="Autocorrelation from adjusted model",xlim=c(0,8))

#model1 diagnostics
pacf(residuals(mo.model.typhoid1,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))
mo.model.typhoid1 <- update(mo.model.typhoid1,.~.+Lag(residuals(mo.model.typhoid1,type="deviance"),1))
mo.model.typhoid1 <- update(mo.model.typhoid1,.~.+Lag(residuals(mo.model.typhoid1,type="deviance"),2))
mo.model.typhoid1 <- update(mo.model.typhoid1,.~.+Lag(residuals(mo.model.typhoid1,type="deviance"),7))
mo.model.typhoid1 <- update(mo.model.typhoid1,.~.+Lag(residuals(mo.model.typhoid1,type="deviance"),8))
pacf(residuals(mo.model.typhoid1,type="deviance"),na.action=na.omit,main="Autocorrelation from adjusted model",xlim=c(0,8))

#model2 diagnostics
pacf(residuals(mo.model.typhoid2,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))
mo.model.typhoid2 <- update(mo.model.typhoid2,.~.+Lag(residuals(mo.model.typhoid2,type="deviance"),7))
mo.model.typhoid2 <- update(mo.model.typhoid2,.~.+Lag(residuals(mo.model.typhoid2,type="deviance"),1))
mo.model.typhoid2 <- update(mo.model.typhoid2,.~.+Lag(residuals(mo.model.typhoid2,type="deviance"),8))
mo.model.typhoid2 <- update(mo.model.typhoid2,.~.+Lag(residuals(mo.model.typhoid2,type="deviance"),6))
pacf(residuals(mo.model.typhoid2,type="deviance"),na.action=na.omit,main="Autocorrelation from adjusted model",xlim=c(0,8))

#model disgnostic plots
par(mfrow=c(2,2))
plot(mo.dlnmT$date, residuals(mo.model.typhoid1,type="deviance"), pch=19, cex=0.4, col=grey(0.6),main="Deviance residuals", ylab="Residual deviance", xlab="Year") 
abline(h=0,lty=2,lwd=2)
pacf(residuals(mo.model.typhoid1,type="deviance"),na.action=na.omit,main="Autocorrelation from adjusted model",xlim=c(0,8))
acf(residuals(mo.model.typhoid1,type="deviance"),na.action=na.omit,main="Autocorrelation from adjusted model",xlim=c(0,8))
plot(mo.dlnmT$incidtyp, pch=10, cex=0.4, col=grey(0.6), size =1 ,main="Adjusted model prediction", ylab="Typhoid incidence rate", xlab="Month")
lines(predict(mo.model.typhoid1, type="response"), col="red")
legend("topleft", legend=c("observed", "predicted"), col=c("grey", "red"), lty=1:1, cex=0.8, lwd=1)
dev.off()

par(mfrow=c(2,2))
plot(mo.dlnmT$date, residuals(mo.model.typhoid2,type="deviance"), pch=19, cex=0.4, col=grey(0.6),main="Deviance residuals", ylab="Residual deviance", xlab="Year") 
abline(h=0,lty=2,lwd=2)
pacf(residuals(mo.model.typhoid2,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))
acf(residuals(mo.model.typhoid2,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))
plot(mo.dlnmT$incidtyp, pch=10, cex=0.4, col=grey(0.6), size =1 ,main="Original model prediction", ylab="Typhoid incidence rate", xlab="Month")
lines(predict(mo.model.typhoid2, type="response"), col="red")
legend("topleft", legend=c("observed", "predicted"), col=c("grey", "red"), lty=1:1, cex=0.8, lwd=1)
dev.off()

#========DLNM MULTIVARIATE ESTIMATES========

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

#========DLNM SENSITIVITY ANALYSIS========
