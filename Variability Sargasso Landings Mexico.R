# Code used for manuscript "Spatio-temporal variability of pelagic Sargassum landings on the northern Mexican Caribbean"

Sys.setenv(TZ='GMT')
rm(list=ls(all=TRUE))
setwd("F:/GitHub/VarSarLandings")
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

str(VS)

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

mean(VS$m3km)#overall mean

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


TableMonthSite<- ddply(VS, c("Year", "Month", "Site"), function(d) {c(Tot = round(sum(d$m3km),2))})  
TableMonthSite

write.table(TableMonthSite, "TableMonthSite.csv", quote=F, sep = ",", row.names=F)

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

# Plot volume per site per month for each year(Fig 5)
ggbarplot(VS,x = "Month", y = "m3km", 
          fill= "Year", palette = c("#39568CFF", "#FDE725FF"), facet.by = "Site",ncol=2,scales="free_y",
          position = position_dodge(0.8), legend = "bottom")+
    xlab("Month")+
    ylab(expression(paste("Volume of sargasso ", "("*~ m^3 ~Km^-1*")")))

#Statistical comparison between sites (Kruskal-Wallis and Wilcoxon)
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
Site1<- VS[VS$Site=="1",]
Site2<- VS[VS$Site=="2",]
Site3<- VS[VS$Site=="3",]
Site4<- VS[VS$Site=="4",]
Site5<- VS[VS$Site=="5",]
Site6<- VS[VS$Site=="6",]
Site7<- VS[VS$Site=="7",]

coin::wilcox_test(m3km ~ Year, conf.int=T, Site1) #do for each site

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

#Sum m3 sargassum vs wind direction for 3 sites (Fig. 6)
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

# plot (Fig. 6)
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

# Time Plot (sargasso, SST, wind speed) (Fig. 7)
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

site="1" # Run for each site

site.18 <- SatHot[SatHot$Site==site & SatHot$Year=="2018", 5] # proportion of m3/km playa per month per year
site.19 <- SatHot[SatHot$Site==site & SatHot$Year=="2019", 5]

ks.test(car.18, site.18); ks.test(car.19, site.19)
ks.test(EEZ.18, site.18); ks.test(EEZ.19, site.19)


# Correlation Sargassum landings vs SST in situ (daily for three sites)

corre <- ggscatter(SED, x = "SeaTemp", y = "m3km",
                   color = "Year", palette = "jco",
                   add = "reg.line", conf.int = TRUE, facet.by = "Year")

corre + stat_cor(aes(color = Year),method = "kendal", label.x = 22, label.y = 950)+
  scale_x_continuous(breaks = get_breaks(by = 1, from = 0),
                     limits =  c(22, 30))


# Supplementary figure 2
rm(list=ls()); gc()
f0= read.csv("Sat Sar vs SST.csv", header=T); f0

s.sea <- f0[1:24,4]; s.eez <- f0[25:48,4]
t.sea <- f0[1:24,5]; t.eez <- f0[25:48,5]

x=1:24
par(mar=c(5,4,4,5)+.1)

plot(x,s.eez, type="l", lwd=3, ylab="Sargassum biomass (10^6 tones)", xlab="Months")
lines(x,s.sea, lty=2, lwd=3)
abline(v=12, col=4, lty=3, lwd=3)
legend(13,73000, legend= c(" SW Caribbean", " EEZ"), lty= c(1, 2))

par(new=TRUE)
plot(x, t.sea, type="l", lwd=3, col=2, xaxt="n", yaxt="n", xlab="", ylab="")
lines(x, t.eez, lty=2, col=2, lwd=3)
axis(4)
mtext("SST °C", side=4,line=3, col="red")

par(mar= c(5.1, 4.1, 4.1, 2.1))

# Correlation (Kendall) Sargassum abundance (ton 10^6) satellite in SWC vs EEZ (2018= HIGHLY, 2019= SIMILAR)
#2018
cor.test(f0[1:12,4], f0[25:36,4], method="kendall") # rho= 0.24; p= 0.31

#2019
cor.test(f0[13:24,4], f0[37:48,4], method="kendall") # rho= 0.81; p< 0.00

# Correlation SST satellite vs Sargassum biomass satellite
#SWC
cor.test(f0[1:12,4], f0[1:12,5], method="kendall") # rho= 0.1212; p= 0.6384

#EEZ
cor.test(f0[13:24,4], f0[13:24,5], method="kendall") 

# Correlation Satellite data Sargassum biomass  vs SST (monthly)
car.18 <- f0[f0$Year=="2018",]
car.19 <- f0[f0$Year=="2019",]

# SST vs sargasso satellite - Change year
corre2 <- ggscatter(car.18, x = "SST", y = "Sargassum",
                    color = "Site", palette = "jco",
                    add = "reg.line", conf.int = TRUE, facet.by = "Site")

corre2 + stat_cor(aes(color = Site),method = "kendall", label.x = 26.5, label.y = 70000)+
  scale_x_continuous(breaks = get_breaks(by = 0.5, from = 0),
                     limits =  c(25, 31))

# SST vs in situ sargasso - Change year
corre2 <- ggscatter(car.19, x = "SST", y = "InSituSar",
                    color = "Site", palette = "jco",
                    add = "reg.line", conf.int = TRUE, facet.by = "Site")

corre2 + stat_cor(aes(color = Site),method = "kendall", label.x = 26.5, label.y = 40000)+
  scale_x_continuous(breaks = get_breaks(by = 0.5, from = 0),
                     limits =  c(25, 31))

