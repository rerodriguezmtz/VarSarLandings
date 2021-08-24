# Code used for manuscript "Spatio-temporal variability of pelagic Sargassum landings on the northern Mexican Caribbean"

Sys.setenv(TZ='GMT')
rm(list=ls(all=TRUE))
setwd("D:/GitHub/VarSarLandings")
VS= read.csv("VolSar.csv", sep = ",")

#Load libraries
library(coin)
library(doBy)
library(dplyr)
library(Hmisc)  
library(ggplot2)
library (ggpubr)
library (lattice)
library(lsr)
library (openair)
library (pgirmess) 
library(plyr)
library(reshape2)
library (rstatix)
library(scales)
library(tidyr)

# Transformations
VS$Date= as.factor(VS$Date)
VS$Year = as.factor(VS$Year)
VS$Month = as.factor(VS$Month)
VS$Site = as.factor(VS$Site)

### Order output
VS$Year=factor(VS$Year, levels = c("2018", "2019"))

# Sargasso volume per year
TableYear<- ddply(VS, c("Year"), function(d) {c(Tot = round(sum(d$Volume),2),
                                                        Min = round(min(d$m3km),2),
                                                        Max = round(max(d$m3km),2),
                                                        Mean = round(mean(d$m3km),2),
                                                        SD = round(sd(d$m3km),2))})  
TableYear

mean(VS$m3km)

# Sargasso volume per year per site
TableYear<- ddply(VS, c("Year", "Site"), function(d) {c(Tot = round(sum(d$Volume),2),
                                                Min = round(min(d$m3km),2),
                                                Max = round(max(d$m3km),2),
                                                Mean = round(mean(d$m3km),2),
                                                SD = round(sd(d$m3km),2))})  
TableYear

# Sargasso volume per year per month
TableMonth<- ddply(VS, c("Year", "Month"), function(d) {c(Tot = round(sum(d$Volume),2),
                                                        Min = round(min(d$m3km),2),
                                                        Max = round(max(d$m3km),2),
                                                        Mean = round(mean(d$m3km),2),
                                                        SD = round(sd(d$m3km),2))})  
TableMonth

# bootstrapped IC function, where x is the data set
ic=function(x) {
    bdat=numeric(10000) 
    for (i in 1:10000) bdat[i]=mean(sample(x, replace=T))
    low= round(quantile(bdat, 0.025, names=F),2)
    mean=round(quantile(bdat, 0.5, names=F), 2)
    upp= round(quantile(bdat, 0.975, names=F),2)
    cbind(mean,low,upp)
}

# Mean (CI) volume per year mean
set.seed(258)
MeanYear <- summaryBy(m3km ~ Year, FUN=ic, VS); MeanYear

# Mean (CI) volume per site per year
set.seed(258)
MeanSiteYear <- summaryBy(m3km ~ Year + Site, FUN=ic, VS); MeanSiteYear

write.table(MeanSiteYear, "TableSiteYear.csv", quote=F, sep = ",", row.names=F)

# Mean (CI) volume per month per year (all sites)
set.seed(258)
MeanYearMonth <- summaryBy(m3km ~ Year + Month, FUN=ic, VS); MeanYearMonth

write.table(MeanYearMonth, "TableYearMonth.csv", quote=F, sep = ",", row.names=F)

# Mean (CI) volume per month per year per site
set.seed(258)
MeanYearMonthSite <- summaryBy(m3km ~ Year + Site + Month, FUN=ic, VS); MeanYearMonthSite

write.table(MeanYearMonthSite, "TableSiteMonth.csv", quote=F, sep = ",", row.names=F)

# Plot volume per site per month for each year(Fig 4)
ggbarplot(VS,x = "Month", y = "m3km", 
          fill= "Year", palette = c("#39568CFF", "#FDE725FF"), facet.by = "Site",ncol=2,scales="free_y",
          position = position_dodge(0.8), legend = "bottom")+
    xlab("Month")+
    ylab(expression(paste("Volume of sargasso ", "("*~ m^3 ~Km^-1*")")))

#Statistical comparison between sites
Year1<- VS[VS$Year=="2018",]
Year2<- VS[VS$Year=="2019",]

library (pgirmess) 
kruskal.test(Year1$m3km~Year1$Site) 
Y1=compare_means(m3km ~ Site,  data = Year1)
Y1
write.table(Y1, "TableKW2018.csv", quote=F, sep = ",", row.names=F)

kruskal.test(Year2$m3km~Year2$Site) 
Y2=compare_means(m3km ~ Site,  data = Year2)
Y2
write.table(Y2, "TableKW2018.csv", quote=F, sep = ",", row.names=F)

#Statistical comparison for each site between years
SiteA<- VS[VS$Site=="A",]
SiteB<- VS[VS$Site=="B",]
SiteC<- VS[VS$Site=="C",]
SiteD<- VS[VS$Site=="D",]
SiteE<- VS[VS$Site=="E",]
SiteF<- VS[VS$Site=="F",]
SiteG<- VS[VS$Site=="G",]

coin::wilcox_test(m3km ~ Year, conf.int=T, SiteA) #do for each site

# Wind rose plot (Fig. 2)
Env=read.csv("EnvVar.csv", sep=",")

#Transformations
Env$Date= as.Date(Env$Date, format = "%d/%m/%Y", tz = "GMT")
Env$Year = as.factor(Env$Year)

windRose(Env, type = c("Year"), angle = 45)

# Daily volume (3 hotels) vs environmental factors
SED= read.csv("SarVolEnv.csv", sep = ",")

#Transformations
SED$Year = as.factor(SED$Year)
SED$Month = as.factor(SED$Month)
SED$Day = as.factor(SED$Day)
SED$date= as.Date(SED$date, format = "%d/%m/%Y", tz = "GMT")

#Sum m3 sargassum vs wind direction for 3 sites (Fig. 5)
# divide-up date by wd
wd.cut <- cut(SED[, "wd"], breaks = seq(0, 360, length = 9))

SEDtot<-aggregate(m3km ~ Year+Site+wd.cut, data = SED, FUN = sum) # to make table
write.table(SEDtot, "SummarySarWind.csv", quote=F, sep = ",", row.names=F)

# define the levels for plotting
wd <- seq(0, 360, by = 45)
levels(wd.cut) <- paste(wd[-length(wd)], "-", wd[-1], " degrees", sep = "")

# summarise by year/month and wd
summary.data <- aggregate(SED["m3km"], list(date = format(SED$date,"%Y-%m"),
                                            wd = wd.cut), sum, na.rm = TRUE)

write.table(summary.data, "SummaryWind.csv", quote=F, sep = ",", row.names=F)

# need to get into year/month/day
newdate = paste(summary.data$date,"-01", sep = "")
newdate = as.Date(newdate, format = "%Y-%m-%d")

# add to summary
summary.data <- cbind(summary.data, newdate)

# plot (Fig. 5)
xyplot(m3km ~ newdate | wd, 
       data = summary.data, 
       layout = c(4, 2),
       as.table = TRUE,
       xlab = "date",
       ylab = "Sargasso (m3 km-1)",
       panel = function(x, y) {
           panel.grid(h = -1, v = 0)
           panel.abline(v = seq(as.Date("2018/1/1"), as.Date("2019/12/31"),
                                "years"),
                        col = "grey85")
           panel.xyplot(x, y, type = "l", lwd = 2)})

# Time Plot (sargasso, SST, wind speed) (Fig. 6)
timePlot(SED, pollutant = c("m3km", "SeaTemp"), y.relation = "free") 

# Comparison of curves between Satellite vs Hotels data (Kolmogorov-Smirnov tests)-Supp table 2
# All test are done with proportion data to allow for equivalent scaling. 
SatHot= read.csv("Beach & Sea2.csv", colClasses=c(rep("factor",3), rep("numeric",2)), header=T) # OJO: LOs valores de sargazo son en m3/km para las playas y millones de toneladas para el mar 

# Are Satellite curves of 2018 = 2019?  
car.18 <- SatHot[SatHot$Site=="Sea" & SatHot$Year=="2018", 5] 
car.19 <- SatHot[SatHot$Site=="Sea" & SatHot$Year=="2019", 5]    

ks.test(car.18, car.19) # compares empirical distributions, not raw data

EEZ.18 <- SatHot[SatHot$Site=="EEZ" & SatHot$Year=="2018", 5] 
EEZ.19 <- SatHot[SatHot$Site=="EEZ" & SatHot$Year=="2019", 5]    

# Are Satellite curves = Beach curves?

site="A" # Run for each site

site.18 <- SatHot[SatHot$Site==site & SatHot$Year=="2018", 5] # proportion of m3/km playa per month per year
site.19 <- SatHot[SatHot$Site==site & SatHot$Year=="2019", 5]

ks.test(car.18, site.18); ks.test(car.19, site.19)
ks.test(EEZ.18, site.18); ks.test(EEZ.19, site.19)

