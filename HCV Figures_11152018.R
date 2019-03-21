## Figure Script: Gicquelais, Foxman, Coyle, and Eisenberg. (2019). Hepatitis C transmission in young 
## people who inject drugs: insights using a dynamic model informed by state public health surveillance. 
## Epidemics. https://doi.org/10.1016/j.epidem.2019.02.003. 

## Before running this code, run the model simulations in the code HCV_MS_11152018.m to simulate
## the hepatitis C model and output results as csv or txt files. This script creates figures 
## summarizing results and similar to those published in the article. 

################ Section 1: Packages and Functions to Load First ################### 

#note that ggplot version 2.2.1 is required, version 3.0.0 gives error
require(devtools)
install_version("ggplot2", version = "2.2.1", repos = "http://cran.us.r-project.org")
library(ggplot2);packageVersion("ggplot2")

library(ggridges)
library(grid)
library(gridExtra)
library(car)
library(plotly)
library(plyr)
library(reshape)
library(stringr)

# function for saving a legend from a figure
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# grid.arrange function for making multi-panel figures
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}



################ Section 2: Set Working Directory ################
setwd("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 11.15.2018")


############### Section 3: Importing Data Generated from Matlab ################
#### Import MDHHS Data ####

#MDHHS case counts by year
mdhhs<-read.csv("MDHHSData.csv", header=F)
colnames(mdhhs) <- c("Year",	"acute1529idu_1",	"acute1529idu_2",	"acute1529idu_3",	
                     "chronic1529idu_1",	"chronic1529idu_2",	"chronic1529idu_3")
summary(mdhhs)
mdhhs$chronic1529idu<-mdhhs$chronic1529idu_1+mdhhs$chronic1529idu_2+mdhhs$chronic1529idu_3
mdhhs$acute1529idu<-mdhhs$acute1529idu_1+mdhhs$acute1529idu_2+mdhhs$acute1529idu_3


#### Import Chronic HCV Prevalence Data (Current and Former PWID) ####

#Note that 1st row = parameter index (and becomes variable name)

#Total Chronic Cases - All Ages
LHS_ChronicPrev<-read.csv("LHS_ChronicPrev.csv", header=TRUE)
summary(LHS_ChronicPrev)
LHS_ChronicPrev$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_C<-melt(LHS_ChronicPrev, id = "Year")

#Chronic Cases - 15-19

LHS_ChronicPrev1<-read.csv("LHS_ChronicPrev1.csv", header=TRUE)
summary(LHS_ChronicPrev1)
LHS_ChronicPrev1$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS1_C<-melt(LHS_ChronicPrev1, id = "Year")

#Chronic Cases - 20-25

LHS_ChronicPrev2<-read.csv("LHS_ChronicPrev2.csv", header=TRUE)
summary(LHS_ChronicPrev2)
LHS_ChronicPrev2$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS2_C<-melt(LHS_ChronicPrev2, id = "Year")

#Chronic Cases - 26-29

LHS_ChronicPrev3<-read.csv("LHS_ChronicPrev3.csv", header=TRUE)
summary(LHS_ChronicPrev3)
LHS_ChronicPrev3$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS3_C<-melt(LHS_ChronicPrev3, id = "Year")

#Chronic Cases - 30-64

LHS_ChronicPrev4<-read.csv("LHS_ChronicPrev4.csv", header=TRUE)
summary(LHS_ChronicPrev4)
LHS_ChronicPrev4$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS4_C<-melt(LHS_ChronicPrev4, id = "Year")

#### Import Acute Cases (Current and Former PWID) ####

#Note that header row (var names) is the X + index of the parameter set (eg XParamSetNo.)

#Acute Cases - All Ages

LHS_Acute<-read.csv("LHS_Acute.csv", header=TRUE)
summary(LHS_Acute)
LHS_Acute$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)

AcuteLHS.melted<-melt(LHS_Acute, id = "Year")
summary(AcuteLHS.melted)

#Acute Cases - 15-19

LHS_Acute1<-read.csv("LHS_Acute1.csv", header=TRUE)
summary(LHS_Acute1)
LHS_Acute1$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)

keeps <- c("Year","acute1529idu_1")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute1, id = "Year")
AcuteLHS.melted_1<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_1)

#Acute Cases - 20-25

LHS_Acute2<-read.csv("LHS_Acute2.csv", header=TRUE)
summary(LHS_Acute2)
LHS_Acute2$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)

keeps <- c("Year","acute1529idu_2")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute2, id = "Year")
AcuteLHS.melted_2<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_2)

#Acute Cases - 26-29

LHS_Acute3<-read.csv("LHS_Acute3.csv", header=TRUE)
summary(LHS_Acute3)
LHS_Acute3$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)

keeps <- c("Year","acute1529idu_3")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute3, id = "Year")
AcuteLHS.melted_3<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_3)

#Acute Cases - 30-64

LHS_Acute4<-read.csv("LHS_Acute4.csv", header=TRUE)
summary(LHS_Acute4)
LHS_Acute4$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)

AcuteLHS.melted_4<-melt(LHS_Acute4, id = "Year")
summary(AcuteLHS.melted_4)

#Acute Cases - 15-29
LHS<-merge(AcuteLHS.melted_3,merge(AcuteLHS.melted_1,AcuteLHS.melted_2,by=c("variable","Year")),by=c("variable","Year"))
LHS$value1<-LHS$value+LHS$value.x+LHS$value.y
LHS$acute1529idu<-LHS$acute1529idu_1+LHS$acute1529idu_2+LHS$acute1529idu_3

keeps <- c("variable","Year","value1","acute1529idu")
AcuteLHS.melted_1529<-LHS[keeps]


#### Import LHS Parameter Data #####

params<-read.csv("LHS_ParamRSS.txt", header=F)
colnames(params) <- c("index",	"sigma1","sigma2","sigma3","sigma4","a","b",	
                      "gamma1","gamma2","gamma3","gamma4",
                      "d","eps","zeta1","zeta2","zeta3","zeta4",
                      "etan",	"etap","etaz",
                      "theta1","theta2","theta3","theta4",
                      "k1","k2","k3","k4",
                      "lambda1","lambda2","lambda3","lambda4",
                      "mu1",	"mu2", "mu3","mu4","xi","omic1","omic2","omic3","omic4",
                      "r","tau","phip","phin","psi1","psi2","psi3","psi4","psin",
                      "w","intervention_etap","Z0",
                      "c1",	"c2",	"c3",	"c4",	"c5",	"c6",	"c7",	"c8",	"c9",	"c10",
                      "c11","c12","c13","c14","c15","c16",
                      "RSS","init_b","init_theta1","init_theta2","init_theta3",
                      "init_etan","init_etap",
                      "nsduhyr","FormerInject1","FormerInject2","FormerInject3","FormerInject4",
                      "nsduhyr","PYInject1","PYInject2","PYInject3","PYInject4",
                      "Popsize1","Popsize2","Popsize3","Popsize4")
summary(params)

#Create an indicator for the best fitting 50% of parameter sets
params$best50<-ifelse(params$RSS<=median(params$RSS),1,0)
table(params$best50)

#Create an indicator for the quartile of fit
params$fitquartile<-cut(params$RSS, quantile(params$RSS, c(0,1/4,2/4,3/4,1)), include.lowest=T) 
table(params$fitquartile)

#Create a numeric quartile variable (note: change bounds as needed)
params$fitorder<-ifelse(params$fitquartile=="[290,687]",1,
                        ifelse(params$fitquartile=="(687,833]",2,
                               ifelse(params$fitquartile=="(833,1.01e+03]",3,4)))
table(params$fitorder) 

#Best 10% of paramsets
params$best10<-cut(params$RSS, quantile(params$RSS, c(0,1/10)), include.lowest=T) 
table(params$best10, useNA="ifany")

is.na(params$best10)

#Note: change as needed
params$best10<-ifelse(params$best10=="[290,565]",1,0)

params$best10<-ifelse(is.na(params$best10),0,params$best10)

table(params$best10, useNA="ifany")


#### Dataset of case counts with params ####

paramsRSS<-params
paramsRSS$variable<-paste("X", paramsRSS$index, sep="")

LHS_C <-merge(LHS_C,paramsRSS,by="variable")
AcuteLHS.melted<-merge(AcuteLHS.melted,paramsRSS,by="variable")
AcuteLHS.melted_1529<-merge(AcuteLHS.melted_1529,paramsRSS,by="variable")

LHS1_C <-merge(LHS1_C,paramsRSS,by="variable")
AcuteLHS.melted_1<-merge(AcuteLHS.melted_1,paramsRSS,by="variable")

LHS2_C <-merge(LHS2_C,paramsRSS,by="variable")
AcuteLHS.melted_2<-merge(AcuteLHS.melted_2,paramsRSS,by="variable")

LHS3_C <-merge(LHS3_C,paramsRSS,by="variable")
AcuteLHS.melted_3<-merge(AcuteLHS.melted_3,paramsRSS,by="variable")

LHS4_C <-merge(LHS4_C,paramsRSS,by="variable")
AcuteLHS.melted_4<-merge(AcuteLHS.melted_4,paramsRSS,by="variable")

#remove files no longer needed
rm(mdhhs1,keeps, paramsRSS,LHS_Acute,LHS_Acute1,LHS_Acute2,LHS_Acute3,LHS_Acute4,
   LHS_ChronicPrev,LHS_ChronicPrev1,LHS_ChronicPrev2,LHS_ChronicPrev3,LHS_ChronicPrev4,
   LHS_NewChronic,LHS_NewChronic1,LHS_NewChronic2,LHS_NewChronic3,LHS_NewChronic4)

save.image(file="C:/Users/rgicquel/Desktop/11.15.2018/Dataset Workspace_11.15.18.RData")

#### Immune - all age groups ####

#if needed, restart R here and load workspace
load("Dataset Workspace_11.15.18.RData")

#Total Immune - All Ages
LHS_Immune<-read.csv("LHS_Immune.csv", header=TRUE)
LHS_Immune$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_I<-melt(LHS_Immune, id = "Year")

#Total Immune - 15-19
LHS_Immune1<-read.csv("LHS_Immune1.csv", header=TRUE)
LHS_Immune1$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_I1<-melt(LHS_Immune1, id = "Year")

#Total Immune - 20-25
LHS_Immune2<-read.csv("LHS_Immune2.csv", header=TRUE)
LHS_Immune2$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_I2<-melt(LHS_Immune2, id = "Year")

#Total Immune - 26-29
LHS_Immune3<-read.csv("LHS_Immune3.csv", header=TRUE)
LHS_Immune3$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_I3<-melt(LHS_Immune3, id = "Year")

#Total Immune - 30+
LHS_Immune4<-read.csv("LHS_Immune4.csv", header=TRUE)
LHS_Immune4$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_I4<-melt(LHS_Immune4, id = "Year")

#### Susc - all age groups ####

#Total Susc - All Ages
LHS_Susc<-read.csv("LHS_Susc.csv", header=TRUE)
LHS_Susc$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_S<-melt(LHS_Susc, id = "Year")

#Total Susc - 15-19
LHS_Susc1<-read.csv("LHS_Susc1.csv", header=TRUE)
LHS_Susc1$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_S1<-melt(LHS_Susc1, id = "Year")

#Total Susc - 20-25
LHS_Susc2<-read.csv("LHS_Susc2.csv", header=TRUE)
LHS_Susc2$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_S2<-melt(LHS_Susc2, id = "Year")

#Total Susc - 26-29
LHS_Susc3<-read.csv("LHS_Susc3.csv", header=TRUE)
LHS_Susc3$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_S3<-melt(LHS_Susc3, id = "Year")

#Total Susc - 30+
LHS_Susc4<-read.csv("LHS_Susc4.csv", header=TRUE)
LHS_Susc4$Year[(1:151)]<-seq(from = 2000, to = 2030, by = 0.2)
LHS_S4<-melt(LHS_Susc4, id = "Year")

rm(LHS_Susc,LHS_Susc1,LHS_Susc2,LHS_Susc3,LHS_Susc4,
   LHS_Immune,LHS_Immune1,LHS_Immune2,LHS_Immune3,LHS_Immune4)

#Merge into 1 file with all PWID classes

# Acute + Chronic + Susc + Immune = Total PWID
PWID<-merge(LHS_S,AcuteLHS.melted[,1:73],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value.x"] <- "Susc"
colnames(PWID)[colnames(PWID)=="value.y"] <- "Acute"

PWID<-merge(PWID,LHS_C[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Chronic"

PWID<-merge(PWID,LHS_I[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Immune"

PWID<-merge(PWID,AcuteLHS.melted_1[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Acute_1"
PWID<-merge(PWID,AcuteLHS.melted_2[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Acute_2"
PWID<-merge(PWID,AcuteLHS.melted_3[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Acute_3"
PWID<-merge(PWID,AcuteLHS.melted_4[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Acute_4"

PWID<-merge(PWID,LHS1_C[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Chronic_1"
PWID<-merge(PWID,LHS2_C[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Chronic_2"
PWID<-merge(PWID,LHS3_C[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Chronic_3"
PWID<-merge(PWID,LHS4_C[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Chronic_4"

PWID<-merge(PWID,LHS_I1[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Immune_1"
PWID<-merge(PWID,LHS_I2[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Immune_2"
PWID<-merge(PWID,LHS_I3[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Immune_3"
PWID<-merge(PWID,LHS_I4[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Immune_4"

PWID<-merge(PWID,LHS_S1[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Susc_1"
PWID<-merge(PWID,LHS_S2[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Susc_2"
PWID<-merge(PWID,LHS_S3[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Susc_3"
PWID<-merge(PWID,LHS_S4[,1:3],by=c("variable","Year"))
colnames(PWID)[colnames(PWID)=="value"] <- "Susc_4"

rm(LHS_I,LHS_I1,LHS_I2,LHS_I3,LHS_I4,LHS_Immune,LHS_Immune1,LHS_Immune2,LHS_Immune3,LHS_Immune4,
   LHS_S,LHS_S1,LHS_S2,LHS_S3,LHS_S4,LHS_Susc,LHS_Susc1,LHS_Susc2,LHS_Susc3,LHS_Susc4,test)

save.image(file="C:/Users/rgicquel/Desktop/11.15.2018/Dataset Workspace_11.15.18.RData")


#### Save and/or Load R Workspace ####
# rm(acute1,LHS_Acute,LHS_Acute1,LHS_Acute2,LHS_Acute3,LHS_Acute4,
#     LHS_ChronicPrev,LHS_ChronicPrev1,LHS_ChronicPrev2,LHS_ChronicPrev3,LHS_ChronicPrev4,
#     LHS_NewChronic,LHS_NewChronic1,LHS_NewChronic2,LHS_NewChronic3,LHS_NewChronic4,
#     new1,prev1,paramfit,param1,df1,df2,df3,df4,df5,df6,df7,df8,int,df9,df10,df11,df12,df13,df14,df15,df16,
#     prev,acute,new)
# 
# save.image(file="/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 11.5.2018/Dataset Workspace_11.5.18.RData")

load("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 11.15.2018/Dataset Workspace_11.15.18.RData")




############### Section 4: Latin Hypercube Sampling Results ################
#### Plot of Acute Cases Fit to MDHHS Data (All, Color by RSS) #####
#Plotting Parameters
summary(AcuteLHS.melted_1529$value1*AcuteLHS.melted_1529$r)

xbreaks <- c(2000, 2010, 2020,2030)
ybreaks <- c(0,20,40)
ylabels <- c("0","20","40")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

#Dummy plots for legends:
dummy<-ggplot(AcuteLHS.melted_1529,aes(x=gamma1,y=k1,color=RSS))+
  geom_point(size=2)+#, alpha=0.5)+
  ylab(text.ylab)+
  xlab(text.xlab)+
  theme_bw(base_size=18)+
  theme(panel.border = element_rect(linetype = "solid", colour = "black"))+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_colour_gradientn(name="RSS",colours=rainbow(5))
dummy
legend1<-g_legend(dummy)


dummy1<-ggplot(AcuteLHS.melted_1529)+
  geom_point(aes(x=Year, y=acute1529idu,colour="Data"),show.legend = TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 450), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  scale_colour_manual(name=NULL, values = c("#000000"), breaks = c("Data")) +  
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Age Groups')+
  theme(plot.title = element_text(hjust = 0.5))
dummy1
Datalegend<-g_legend(dummy1)


#Plot of fit and data
p1<-ggplot(AcuteLHS.melted_1529, aes(x=Year, y=value1*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 55), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Ages: 15-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1

# title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = 0,gp = gpar(fontsize = 16))
# 
# p1 <- arrangeGrob(p1, top = title.grob)
# grid.arrange(p1)

summary(AcuteLHS.melted_1$value*AcuteLHS.melted_1$r)
summary(AcuteLHS.melted_2$value*AcuteLHS.melted_2$r)
summary(AcuteLHS.melted_3$value*AcuteLHS.melted_3$r)

ybreaks <- c(0,10,20,30)
ylabels <- c("0","10","20","30")


p2<-ggplot(AcuteLHS.melted_1, aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_1),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30.1), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  #theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  ggtitle('Ages 15-19')+
  theme(plot.title = element_text(hjust = 0.5)) 
p2
# title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = 0,gp = gpar(fontsize = 16))
# 
# p2 <- arrangeGrob(p2, top = title.grob)
# grid.arrange(p2)

p3<-ggplot(AcuteLHS.melted_2, aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_2),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30.1), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 20-25')+
  theme(plot.title = element_text(hjust = 0.5)) 
p3
# title.grob <- textGrob(label = "C)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = 0,gp = gpar(fontsize = 16))
# 
# p3 <- arrangeGrob(p3, top = title.grob)
# grid.arrange(p3)

p4<-ggplot(AcuteLHS.melted_3, aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30.1), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 26-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p4
# title.grob <- textGrob(label = "D)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = 0,gp = gpar(fontsize = 16))
# 
# p4 <- arrangeGrob(p4, top = title.grob)
# grid.arrange(p4)

summary(AcuteLHS.melted_4$value*AcuteLHS.melted_4$r)

p5<-ggplot(AcuteLHS.melted_4, aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  #geom_point(aes(x=Year, y=acute1529idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  #scale_y_continuous(limits = c(0, 120), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 30-64')+
  theme(plot.title = element_text(hjust = 0.5)) 
p5


# All Ages with title
w <- 15; h <- 4
lay<-rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7))
LHS_Acute<-grid.arrange(grid.text("All Simulations (n=10,000)",gp=gpar(fontsize=22, col="black")),
                        p2,p3,p4,p5,legend1,Datalegend,layout_matrix=lay)
ggsave(sprintf("LHS_Acute_RSSColor_1564_Title.png"), LHS_Acute, width=w, height=h)

#### Plot of Acute Cases Fit to MDHHS Data (Best 10%, Color by RSS) #####
#Plotting Parameters

summary(subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best10==1)$value1*subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best10==1)$r)

xbreaks <- c(2000, 2010, 2020, 2030)
ybreaks <- c(0,20,40,60)
ylabels <- c("0","20","40","60")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

#Dummy plots for legends:
dummy<-ggplot(data=subset(params,params$best10==1),aes(x=gamma1,y=k1,color=RSS))+
  geom_point(size=2)+#, alpha=0.5)+
  ylab(text.ylab)+
  xlab(text.xlab)+
  theme_bw(base_size=18)+
  theme(panel.border = element_rect(linetype = "solid", colour = "black"))+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_colour_gradientn(name="RSS",colours=rainbow(5))
dummy
legend1<-g_legend(dummy)


dummy1<-ggplot(AcuteLHS.melted_1529)+
  geom_point(aes(x=Year, y=acute1529idu,colour="Data"),show.legend = TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 300), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  scale_colour_manual(name=NULL, values = c("#000000"), breaks = c("Data")) +  
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Age Groups')+
  theme(plot.title = element_text(hjust = 0.5))
dummy1
Datalegend<-g_legend(dummy1)


#Plot of fit and data
p1<-ggplot(subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best10==1), 
           aes(x=Year, y=value1*r, by=variable, 
               color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 60.1), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Ages: 15-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1

# title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = 0,gp = gpar(fontsize = 16))
# 
# p1 <- arrangeGrob(p1, top = title.grob)
# grid.arrange(p1)

summary(subset(AcuteLHS.melted_1,AcuteLHS.melted_1$best10==1)$value*subset(AcuteLHS.melted_1,AcuteLHS.melted_1$best10==1)$r)
summary(subset(AcuteLHS.melted_2,AcuteLHS.melted_2$best10==1)$value*subset(AcuteLHS.melted_2,AcuteLHS.melted_2$best10==1)$r)
summary(subset(AcuteLHS.melted_3,AcuteLHS.melted_3$best10==1)$value*subset(AcuteLHS.melted_3,AcuteLHS.melted_3$best10==1)$r)

ybreaks <- c(0,10,20,30)
ylabels <- c("0","10","20","30")

p2<-ggplot(subset(AcuteLHS.melted_1,AcuteLHS.melted_1$best10==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_1),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30.1), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  #theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  ggtitle('Ages 15-19')+
  theme(plot.title = element_text(hjust = 0.5)) 
p2
# title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = 0,gp = gpar(fontsize = 16))
# 
# p2 <- arrangeGrob(p2, top = title.grob)
# grid.arrange(p2)

p3<-ggplot(subset(AcuteLHS.melted_2,AcuteLHS.melted_2$best10==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_2),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30.1), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 20-25')+
  theme(plot.title = element_text(hjust = 0.5)) 
p3
# title.grob <- textGrob(label = "C)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = 0,gp = gpar(fontsize = 16))
# 
# p3 <- arrangeGrob(p3, top = title.grob)
# grid.arrange(p3)

p4<-ggplot(subset(AcuteLHS.melted_3,AcuteLHS.melted_3$best10==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30.1), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 26-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p4
# title.grob <- textGrob(label = "D)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = 0,gp = gpar(fontsize = 16))
# 
# p4 <- arrangeGrob(p4, top = title.grob)
# grid.arrange(p4)

p5<-ggplot(subset(AcuteLHS.melted_4,AcuteLHS.melted_4$best10==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  #geom_point(aes(x=Year, y=acute1529idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  #scale_y_continuous(limits = c(0, 30), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 30-64')+
  theme(plot.title = element_text(hjust = 0.5)) 
p5

w <- 15; h <- 4
lay<-rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7))
LHS_Acute<-grid.arrange(grid.text("Best 10% of Simulations (n=1,000)",gp=gpar(fontsize=22, col="black")),
                        p2,p3,p4,p5,legend1,Datalegend,layout_matrix=lay)
ggsave(sprintf("LHS_Acute_RSSColor_Best10_1564_Title.png"), LHS_Acute, width=w, height=h)



#### Plot of Chronic HCV Prevalence (%, All) ####

# PWID$best10<-ifelse(PWID$RSS<=309.09,1,0)
# table(PWID$best10)

PWID$HCVprev<-(PWID$Acute+PWID$Chronic)/(PWID$Acute+PWID$Chronic+PWID$Susc+PWID$Immune)*100
summary(PWID$HCVprev)

PWID$HCVprev1<-(PWID$Acute_1+PWID$Chronic_1)/(PWID$Acute_1+PWID$Chronic_1+PWID$Susc_1+PWID$Immune_1)*100
summary(PWID$HCVprev1)

PWID$HCVprev2<-(PWID$Acute_2+PWID$Chronic_2)/(PWID$Acute_2+PWID$Chronic_2+PWID$Susc_2+PWID$Immune_2)*100
summary(PWID$HCVprev2)

PWID$HCVprev3<-(PWID$Acute_3+PWID$Chronic_3)/(PWID$Acute_3+PWID$Chronic_3+PWID$Susc_3+PWID$Immune_3)*100
summary(PWID$HCVprev3)

PWID$HCVprev4<-(PWID$Acute_4+PWID$Chronic_4)/(PWID$Acute_4+PWID$Chronic_4+PWID$Susc_4+PWID$Immune_4)*100
summary(PWID$HCVprev4)

PWID<-merge(params[,c(1,94)],PWID,"index",all.y=TRUE)
table(PWID$best10)

xbreaks <- c(2000, 2010, 2020,2030)
ybreaks <- c(0,25,50,75,100)
ylabels <- c("0","25","50","75","100")
text.xlab <- "Year"
text.ylab <- "HCV Prevalence (%) among PWID"
model.colors <- c('#999999', '#000000') 

#Dummy plots for legends:
dummy<-ggplot(data=subset(params,params$best10==1),aes(x=gamma1,y=k1,color=RSS))+
  geom_point(size=2)+#, alpha=0.5)+
  ylab(text.ylab)+
  xlab(text.xlab)+
  theme_bw(base_size=18)+
  theme(panel.border = element_rect(linetype = "solid", colour = "black"))+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_colour_gradientn(name="RSS",colours=rainbow(5))
dummy
legend1<-g_legend(dummy)

p1<-ggplot(subset(PWID,PWID$best10==1), aes(x=Year, y=HCVprev1, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 100), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("HCV Prevalence (%)") +
  xlab(text.xlab) +
  ggtitle('Ages 15-19')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1

p2<-ggplot(subset(PWID,PWID$best10==1), aes(x=Year, y=HCVprev2, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 100), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 20-25')+
  theme(plot.title = element_text(hjust = 0.5)) 
p2

p3<-ggplot(subset(PWID,PWID$best10==1), aes(x=Year, y=HCVprev3, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 100), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 26-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p3

p4<-ggplot(subset(PWID,PWID$best10==1), aes(x=Year, y=HCVprev4, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 100), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 30-64')+
  theme(plot.title = element_text(hjust = 0.5)) 
p4

#layout w/legend at side
w <- 15; h <- 4

#No plot title
lay<-rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5),c(1,1,1,2,2,2,3,3,3,4,4,4,5))
LHS_Chr<-grid.arrange(p1,p2,p3,p4,legend1,layout_matrix=lay)
ggsave(sprintf("LHS_ChronicPct_RSSColor_Best10_1564.pdf"), LHS_Chr, width=w, height=h)
ggsave(sprintf("LHS_ChronicPct_RSSColor_Best10_1564.png"), LHS_Chr, width=w, height=h)

#with title
w <- 15; h <- 4
lay<-rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7))

LHS_Chr<-grid.arrange(grid.text("Chronic HCV Prevalence (% of PWID)",gp=gpar(fontsize=22, col="black")),
                      p1,p2,p3,p4,legend1,layout_matrix=lay)

ggsave(sprintf("LHS_ChronicPct_RSSColor_Best10_1564_Title.png"), LHS_Chr, width=w, height=h)



#### Chronic HCV Prevalence Plots (All, Combined Plot, Case Counts) ####

#Plotting Parameters
xbreaks <- c(2000, 2005, 2010,2015)
ybreaks <- c(0, 50000,100000,150000,200000,250000,300000,350000)
ylabels <- c("0","50,000","100,000","150,000","200,000","250,000","300,000","350,000")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000')
max(LHS4_C$value)

LHS_prevalence<-ggplot(LHS_C, aes(x=Year, y=value, by=variable)) + 
  #geom_line(aes(x=Year, y=value, by=variable, colour='All Ages'), show.legend=TRUE)+
  geom_line(aes(x=LHS3_C$Year, y=LHS4_C$value, by=LHS4_C$variable, colour='30-64 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS3_C$Year, y=LHS3_C$value, by=LHS3_C$variable, colour='26-29 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS2_C$Year, y=LHS2_C$value, by=LHS2_C$variable, colour='20-25 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS1_C$Year, y=LHS1_C$value, by=LHS1_C$variable, colour='15-19 Years'), show.legend=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(LHS4_C$value+10)), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "right", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name=NULL, values = c("#CCCCCC","#999999", "#666666", "#333333"), breaks = c("30-64 Years","26-29 Years", "20-25 Years", "15-19 Years")) +
  ylab(text.ylab) +
  xlab(text.xlab) 
LHS_prevalence 

w <- 6; h <- 4
ggsave(sprintf("LHS_prevalence_AllAgesAllFits.pdf"), LHS_prevalence, width=w, height=h)
ggsave(sprintf("LHS_prevalence_AllAgesAllFits.png"), LHS_prevalence, width=w, height=h)


#Restrict to the Best-Fitting 50%
LHS_prevalence<-ggplot() + 
  #geom_line(data = subset(LHS_C,LHS_C$best50==1), aes(x=Year, y=value, by=variable, colour='All Ages'), show.legend=TRUE)+
  geom_line(data = subset(LHS4_C,LHS4_C$best50==1),aes(x=Year, y=value, by=variable, colour='30-64 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS3_C,LHS3_C$best50==1),aes(x=Year, y=value, by=variable, colour='26-29 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS2_C,LHS2_C$best50==1),aes(x=Year, y=value, by=variable, colour='20-25 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS1_C,LHS1_C$best50==1),aes(x=Year, y=value, by=variable, colour='15-19 Years'), show.legend=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(LHS4_C$value+10)), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "right", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name=NULL, values = c("#CCCCCC","#999999", "#666666", "#333333","#000000"), breaks = c("30-64 Years","26-29 Years", "20-25 Years", "15-19 Years")) +
  ylab(text.ylab) +
  xlab(text.xlab) 
LHS_prevalence

w <- 6; h <- 4
ggsave(sprintf("LHS_prevalence_Best50.pdf"), LHS_prevalence, width=w, height=h)
ggsave(sprintf("LHS_prevalence_Best50.png"), LHS_prevalence, width=w, height=h)


#### Chronic Prevalence plots among 15-29 YO (All simulations) ####
#Plot only the 15-29 YO age groups
xbreaks <- c(2000, 2005, 2010,2015)
ybreaks <- c(0,3000,6000,9000,12000)
ylabels <- c("0","3,000","6,000","9,000","12,000")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000')
max(LHS2_C$value)
max(LHS3_C$value)

LHS_prevalence<-ggplot(LHS_C, aes(x=Year, y=value, by=variable)) + 
  geom_line(aes(x=LHS3_C$Year, y=LHS3_C$value, by=LHS3_C$variable, colour='26-29 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS2_C$Year, y=LHS2_C$value, by=LHS2_C$variable, colour='20-25 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS1_C$Year, y=LHS1_C$value, by=LHS1_C$variable, colour='15-19 Years'), show.legend=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(LHS2_C$value+10)), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "right", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name=NULL, values = c("#CCCCCC", "#666666", "#000000"), breaks = c("26-29 Years", "20-25 Years", "15-19 Years")) +
  ylab(text.ylab) +
  xlab(text.xlab) 
LHS_prevalence

w <- 6; h <- 4
ggsave(sprintf("LHS_prevalence_15-29.pdf"), LHS_prevalence, width=w, height=h)
ggsave(sprintf("LHS_prevalence_15-29.png"), LHS_prevalence, width=w, height=h)




#Restrict to the Best-Fitting 50%
LHS_prevalence<-ggplot() + 
  geom_line(data = subset(LHS3_C,LHS3_C$best50==1),aes(x=Year, y=value, by=variable, colour='26-29 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS2_C,LHS2_C$best50==1),aes(x=Year, y=value, by=variable, colour='20-25 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS1_C,LHS1_C$best50==1),aes(x=Year, y=value, by=variable, colour='15-19 Years'), show.legend=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(subset(LHS2_C,LHS2_C$best50==1)$value+10)), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "right", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name=NULL, values = c("#CCCCCC","#666666", "#000000"), breaks = c("26-29 Years", "20-25 Years", "15-19 Years")) +
  ylab(text.ylab) +
  xlab(text.xlab) 
LHS_prevalence

w <- 6; h <- 4
ggsave(sprintf("LHS_prevalence_Best50_15-19.pdf"), LHS_prevalence, width=w, height=h)
ggsave(sprintf("LHS_prevalence_Best50_15-19.png"), LHS_prevalence, width=w, height=h)


#### Chronic Prevalence by Age Group (Panel Plot, Best 10% Fits) ####

#Plotting Parameters
xbreaks <- c(2000, 2010, 2020,2030)
ybreaks <- c(0,50000,100000,150000,200000)
ylabels <- c("0","50,000","100,000","150,000","200,000")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

summary(subset(LHS_C,LHS_C$best10==1)$value)


#Dummy plots for legends:
dummy<-ggplot(data=subset(params,params$best10==1),aes(x=gamma1,y=k1,color=RSS))+
  geom_point(size=2)+#, alpha=0.5)+
  ylab(text.ylab)+
  xlab(text.xlab)+
  theme_bw(base_size=18)+
  theme(panel.border = element_rect(linetype = "solid", colour = "black"))+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_colour_gradientn(name="RSS",colours=rainbow(5))
legend1<-g_legend(dummy)
dummy


#Plot of fit and data
p1<-ggplot(subset(LHS_C,LHS_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(subset(LHS_C,LHS_C$best10==1)$value+10)), 
                     breaks = c(0,200000,400000,600000), 
                     labels=c("0","200,000","400,000","600,000"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Ages: 15-64')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1

summary(subset(LHS1_C,LHS1_C$best10==1)$value)

p2<-ggplot(subset(LHS1_C,LHS1_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 2750), breaks = c(0,800,1600,2400), labels=c("0","800","1,600","2,400"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  #theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  ggtitle('Ages 15-19')+
  theme(plot.title = element_text(hjust = 0.5)) 
p2

summary(subset(LHS2_C,LHS2_C$best10==1)$value)

p3<-ggplot(subset(LHS2_C,LHS2_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 5000), breaks = c(0,2500,5000), labels=c("0","2,500","5,000"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  #theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  ggtitle('Ages 20-25')+
  theme(plot.title = element_text(hjust = 0.5)) 
p3

summary(subset(LHS3_C,LHS3_C$best10==1)$value)

p4<-ggplot(subset(LHS3_C,LHS3_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 6000), breaks = c(0,2000,4000,6000), labels=c("0","2,000","4,000","6,000"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  #theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  ggtitle('Ages 26-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p4

summary(subset(LHS4_C,LHS4_C$best10==1)$value)

p5<-ggplot(subset(LHS4_C,LHS4_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2030), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 612120), breaks = c(0,200000,400000,600000), labels=c("0","200,000","400,000","600,000"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  #theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  ggtitle('Ages 30-64')+
  theme(plot.title = element_text(hjust = 0.5)) 
p5


#layout w/legend at side 
w <- 12; h <- 6

#No plot title
lay<-rbind(c(1,1,1,1,1,1,1,2),c(1,1,1,1,1,1,1,2),c(3,3,4,4,5,5,6,6),c(3,3,4,4,5,5,6,6))
LHS_Chronic<-grid.arrange(p1,legend1,p2,p3,p4,p5,layout_matrix=lay)
ggsave(sprintf("LHS_ChronicPrev_RSSColor_Best10.pdf"), LHS_Chronic, width=w, height=h)
ggsave(sprintf("LHS_ChronicPrev_RSSColor_Best10.png"), LHS_Chronic, width=w, height=h)

w <- 15; h <- 4

#No plot title
lay<-rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5))
LHS_Chronic<-grid.arrange(p2,p3,p4,p5,legend1,layout_matrix=lay)
ggsave(sprintf("LHS_ChronicPrev_RSSColor_Best10_1564.pdf"), LHS_Chronic, width=w, height=h)
ggsave(sprintf("LHS_ChronicPrev_RSSColor_Best10_1564.png"), LHS_Chronic, width=w, height=h)

#with title
w <- 15; h <- 4
lay<-rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,6),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7),
           c(2,2,2,3,3,3,4,4,4,5,5,5,7))

LHS_Chr<-grid.arrange(grid.text("Chronic HCV Prevalence (No. of Cases)",gp=gpar(fontsize=22, col="black")),
                      p2,p3,p4,p5,legend1,layout_matrix=lay)

ggsave(sprintf("LHS_Chronic_RSSColor_Best10_1564_Title.png"), LHS_Chr, width=w, height=h)


#### 3d Plot of Injection Initiation ####

#set axis parameters
ax1 <- list(
  zeroline = TRUE,
  showline = TRUE,
  mirror = "ticks",
  gridcolor = toRGB("black"),
  gridwidth = 2,
  zerolinecolor = toRGB("black"),
  zerolinewidth = 4,
  linecolor = toRGB("black"),
  linewidth = 6,
  title='Injection Initiation: 15-19'
)
ax2 <- list(
  zeroline = TRUE,
  showline = TRUE,
  mirror = "ticks",
  gridcolor = toRGB("black"),
  gridwidth = 2,
  zerolinecolor = toRGB("black"),
  zerolinewidth = 4,
  linecolor = toRGB("black"),
  linewidth = 6,
  title='Injection Initiation: 20-25'
)
ax3 <- list(
  zeroline = TRUE,
  showline = TRUE,
  mirror = "ticks",
  gridcolor = toRGB("black"),
  gridwidth = 2,
  zerolinecolor = toRGB("black"),
  zerolinewidth = 4,
  linecolor = toRGB("black"),
  linewidth = 6,
  title='Injection Initiation: 26-29'
)

params1<-params[,c(1:81,83:94)]

#create 3d plot
p <- plot_ly(params1, x = ~theta1, y = ~theta2, z = ~theta3, color = ~RSS, opacity=1,
             marker=list(size=4),
             colors=c("red","orange","yellow","green","blue")) %>%
  add_markers() %>%
  layout(scene=list(xaxis=ax1,yaxis=ax2,zaxis=ax3))
p

p <- plot_ly(subset(params1,params$best10==1), x = ~theta1, y = ~theta2, z = ~theta3, color = ~RSS, opacity=1,
             marker=list(size=4),
             colors=c("red","orange","yellow","green","blue")) %>%
  add_markers() %>%
  layout(scene=list(xaxis=ax1,yaxis=ax2,zaxis=ax3))
p

summary(subset(params,params$best10==1)$theta3)

count<-subset(params,params$best10==1&params$theta3<=1.5)





#### LHS Parameter Histograms (Greyscale by Quartiles of Fit) #####

#Set ylabel text
text.ylab <- "# Parameter Sets"
model.colors <- c('#000000', '#666666', '#999999','#CCCCCC') #grey scheme

#Beta

#textsym <- expression(beta[1])
summary(params$b*100)
text.xlab <- bquote('% Infected'*~'Contact'^-1)
text.title <- bquote('Transmission Rate ('*beta*')')

#use to extract legend (bottom) 
p1<-ggplot(params) + 
  geom_histogram(aes(x=b*100, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,bins=30)+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits=c(0,10),breaks = c(0,5,10),labels=c("0","5","10"))+
  scale_y_continuous(limits = c(0, 1670), breaks = c(0,500,1000,1500), labels=c("0","500","1,000","1,500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p1
bottomleg<-g_legend(p1)

#use to extract legend (side)
p2<-ggplot(params) + 
  geom_histogram(aes(x=b*100, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,bins=30)+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "right", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits=c(0,10),breaks = c(0,5,10),labels=c("0","5","10"))+
  scale_y_continuous(limits = c(0, 1670), breaks = c(0,500,1000,1500), labels=c("0","500","1,000","1,500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p2
sideleg<-g_legend(p2)

#beta
p3<-ggplot(params) + 
  geom_histogram(aes(x=b*100, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,bins=30)+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits=c(0,10),breaks = c(0,5,10),labels=c("0","5","10"))+
  scale_y_continuous(limits = c(0, 1670), breaks = c(0,500,1000,1500), labels=c("0","500","1,000","1,500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p3

#sigma1
summary(params$sigma1)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 15-19 ('*sigma[1]*')')

p4<-ggplot(params) + 
  geom_histogram(aes(x=sigma1, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,breaks=seq(0,5.2,0.25))+#, 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits=c(0,5.2),breaks = c(0,2.6,5.2),labels=c("0","2.6","5.2"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p4

#sigma2
summary(params$sigma2)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 20-25 ('*sigma[2]*')')

p5<-ggplot(params) + 
  geom_histogram(aes(x=sigma2, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,breaks=seq(0,5.2,0.25))+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits=c(0,5.2),breaks = c(0,2.6,5.2),labels=c("0","2.6","5.2"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p5

#sigma3
summary(params$sigma3)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 26-29 ('*sigma[3]*')')

p6<-ggplot(params) + 
  geom_histogram(aes(x=sigma3, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,breaks=seq(0,4,0.2))+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits=c(0,4),breaks = c(0,2,4),labels=c("0.0","2.0","4.0"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p6

#sigma4
summary(params$sigma4)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 30-64 ('*sigma[4]*')')

p7<-ggplot(params) + 
  geom_histogram(aes(x=sigma4, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,breaks=seq(0,5,0.25))+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits=c(0,5),breaks = c(0,2.5,5),labels=c("0.0","2.5","5.0"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p7


#alpha
# summary(params$a)
# 
# text.xlab <- "%"
# text.title <- bquote('Sustained Virologic Response ('*alpha*')')
# 
# p8<-ggplot(params) + 
#   geom_histogram(aes(x=a*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(80,100,1), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(80,100), breaks = c(80, 85,90,95,100))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p8

#gamma1
summary(params$gamma1)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Cessation Rate: 15-19 ('*gamma[1]*')')

p9<-ggplot(params) + 
  geom_histogram(aes(x=gamma1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.05,1.2,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1.2), breaks = c(0,0.6,1.2))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p9

#gamma2
summary(params$gamma2)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Cessation Rate: 20-25 ('*gamma[2]*')')

p10<-ggplot(params) + 
  geom_histogram(aes(x=gamma2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.05,1.2,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1.2), breaks = c(0,0.6,1.2))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p10

#gamma3
summary(params$gamma3)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Cessation Rate: 26-29 ('*gamma[3]*')')

p11<-ggplot(params) + 
  geom_histogram(aes(x=gamma3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.05,1.2,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1.2), breaks = c(0,0.6,1.2))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p11


#gamma4
summary(params$gamma4)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Cessation Rate: 30-64 ('*gamma[4]*')')

p12<-ggplot(params) + 
  geom_histogram(aes(x=gamma4, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.05,1.2,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1.2), breaks = c(0,0.6,1.2))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p12

#delta
summary(params$d)

text.xlab <- bquote('%')
text.title <- bquote('Spontaneous Clearance ('*delta*')')

p13<-ggplot(params) + 
  geom_histogram(aes(x=d*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(15,50,1.5), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(15,50), breaks = c(15,32.5,50))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p13

#zeta1
summary(params$zeta1)

text.xlab <- bquote('%')
text.title <- bquote('Current PWID Prevalence: 15-19 ('*zeta[1]*')')

p14<-ggplot(params) + 
  geom_histogram(aes(x=zeta1*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.15,1.0,0.04), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.15,1.0), breaks = c(0.15,0.575,1.0))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p14

#zeta2
summary(params$zeta2)

text.xlab <- bquote('%')
text.title <- bquote('Current PWID Prevalence: 20-25 ('*zeta[2]*')')

p15<-ggplot(params) + 
  geom_histogram(aes(x=zeta2*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.26,1.0,0.035), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.26,1), breaks = c(0.3,0.65,1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p15

#zeta3
summary(params$zeta3)

text.xlab <- bquote('%')
text.title <- bquote('Current PWID Prevalence: 26-29 ('*zeta[3]*')')

p16<-ggplot(params) + 
  geom_histogram(aes(x=zeta3*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.19,1,0.04), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.19,1.0), breaks = c(0.2,0.6,1.0))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p16

#zeta4
summary(params$zeta4)

text.xlab <- bquote('%')
text.title <- bquote('Current PWID Prevalence: 30-64 ('*zeta[4]*')')

p17<-ggplot(params) + 
  geom_histogram(aes(x=zeta4*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.16,1.1,0.045), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.16,1.1), breaks = c(0.2,0.6,1.0))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p17

#etaN
summary(params$etan)

text.xlab <- "Standardized Mortality Ratio (SMR)"
text.title <- bquote('Former vs. Current PWID Mortality ('*eta[N]*')')


p18<-ggplot(params) + 
  geom_histogram(aes(x=etan, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.18,0.54,0.016), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.18,0.54), breaks = c(0.18,0.36,0.54))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p18


#etaP
summary(params$etap)

text.xlab <- "Standardized Mortality Ratio (SMR)"
text.title <- bquote('Current PWID vs. General Mortality ('*eta[P]*')')


p19<-ggplot(params) + 
  geom_histogram(aes(x=etap, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(2.5,15.3,0.6), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2.5,15.3), breaks = c(2.5,8.9,15.3))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p19

#etaZ
summary(params$etaz)

text.xlab <- "Standardized Mortality Ratio (SMR)"
text.title <- bquote('Non-PWID vs. General Mortality ('*eta[Z]*')')


p20<-ggplot(params) + 
  geom_histogram(aes(x=etaz, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(4.3,4.43,0.006), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(4.3,4.43), breaks = c(4.32,4.37,4.42))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p20


#theta1
summary(params$theta1)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 15-19 ('*theta[1]*')')


p21<-ggplot(params) + 
  geom_histogram(aes(x=theta1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.05,5,0.25), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2.5,5))+
  scale_y_continuous(limits = c(0, 5050), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p21


#theta2
summary(params$theta2)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 20-25 ('*theta[2]*')')


p22<-ggplot(params) + 
  geom_histogram(aes(x=theta2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.034,5,0.4), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2.5,5))+
  scale_y_continuous(limits = c(0, 5050), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p22

#theta3
summary(params$theta3)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 26-29 ('*theta[3]*')')


p23<-ggplot(params) + 
  geom_histogram(aes(x=theta3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.034,10,0.5), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,5,10))+
  scale_y_continuous(limits = c(0, 5050), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p23

#theta4
summary(params$theta4)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 30-64 ('*theta[4]*')')


p24<-ggplot(params) + 
  geom_histogram(aes(x=theta4, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.034,1.0,0.045), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p24


#kappa1
summary(params$k1)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Relapse Rate: 15-19 ('*kappa[1]*')')


p25<-ggplot(params) + 
  geom_histogram(aes(x=k1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,2.5,0.15), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,2.5), breaks = c(0,1.25,2.5))+
  scale_y_continuous(limits = c(0, 4200), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p25

#kappa2
summary(params$k2)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Relapse Rate: 20-25 ('*kappa[2]*')')


p26<-ggplot(params) + 
  geom_histogram(aes(x=k2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,1.66,0.1), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1.66), breaks = c(0,0.8,1.6))+
  scale_y_continuous(limits = c(0, 4200), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p26

#kappa3
summary(params$k3)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Relapse Rate: 26-29 ('*kappa[3]*')')


p27<-ggplot(params) + 
  geom_histogram(aes(x=k3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,1.43,0.1), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1.5), breaks = c(0,0.75,1.5))+
  scale_y_continuous(limits = c(0, 4200), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p27


#kappa4
summary(params$k4)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Relapse Rate: 30-64 ('*kappa[4]*')')


p28<-ggplot(params) + 
  geom_histogram(aes(x=k4, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,0.8,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,0.8), breaks = c(0,0.4,0.8))+
  scale_y_continuous(limits = c(0, 4200), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p28


#lambda1
summary(params$lambda1)

text.xlab <- bquote("%")
text.title <- bquote('HCV Prevalence: 15-19 ('*lambda[1]*')')


p29<-ggplot(params) +
  geom_histogram(aes(x=lambda1*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(10,28,0.9), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(10,28), breaks = c(10,19,28), labels=c("10.0","19.0","28.0"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p29


#lambda2
summary(params$lambda2)

text.xlab <- bquote("%")
text.title <- bquote('HCV Prevalence: 20-25 ('*lambda[2]*')')


p30<-ggplot(params) +
  geom_histogram(aes(x=lambda2*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(17,36,0.9), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(17,36), breaks = c(17,26.5,36), labels=c("17.0","26.5","36.0"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p30

#lambda3
summary(params$lambda3)

text.xlab <- bquote("%")
text.title <- bquote('HCV Prevalence: 26-29 ('*lambda[3]*')')


p31<-ggplot(params) +
  geom_histogram(aes(x=lambda3*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(28,53,1.2), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(28,53), breaks = c(28,40.5,53), labels = c("28.0","40.5","53.0"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p31


#lambda4
summary(params$lambda4)

text.xlab <- bquote("%")
text.title <- bquote('HCV Prevalence: 30-64 ('*lambda[4]*')')


p32<-ggplot(params) +
  geom_histogram(aes(x=lambda4*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(41,68,1.25), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(41,68), breaks = c(41,54.5,68), labels=c("41.0","54.5","68.0"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p32


#Crude Death Rate 1
# summary(params$mu1*100000)
# 
# text.xlab <- bquote("Deaths 100,000 Persons"^-1)
# text.title <- bquote('Mortality Rate: 15-19 ('*mu[1]*')')
# 
# p33<-ggplot(params) + 
#   geom_histogram(aes(x=mu1*100000, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(41,70,1.2), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab(text.ylab) +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(40,70), breaks = c(40,50,60,70))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p33
# 
# #Crude Death Rate 2
# summary(params$mu2*100000)
# 
# text.xlab <- bquote("Deaths 100,000 Persons"^-1)
# text.title <- bquote('Mortality Rate: 20-25 ('*mu[2]*')')
# 
# p34<-ggplot(params) + 
#   geom_histogram(aes(x=mu2*100000, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(77,99,1), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(77,100), breaks = c(80,90,100))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p34
# 
# #Crude Death Rate 3
# summary(params$mu3*100000)
# 
# text.xlab <- bquote("Deaths 100,000 Persons"^-1)
# text.title <- bquote('Mortality Rate: 26-29 ('*mu[3]*')')
# 
# p35<-ggplot(params) + 
#   geom_histogram(aes(x=mu3*100000, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(99,120,1), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(99,120), breaks = c(100,110,120))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p35
# 
# 
# #Crude Death Rate 4
# summary(params$mu4*100000)
# 
# text.xlab <- bquote("Deaths 100,000 Persons"^-1)
# text.title <- bquote('Mortality Rate: 30-64 ('*mu[4]*')')
# 
# p36<-ggplot(params) + 
#   geom_histogram(aes(x=mu4*100000, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(490,570,4), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(490,570), breaks = c(500,525,550))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p36


#xi
summary(params$xi*100)

text.xlab <- bquote("%")
text.title <- bquote('HCV Immunity to Reinfection ('*xi*')')

p37<-ggplot(params) + 
  geom_histogram(aes(x=xi*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,25,1.25), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,25), breaks = c(0,12.5,25))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p37

#omic1
summary(params$omic1*100)

text.xlab <- bquote("%")
text.title <- bquote('Former PWID Prevalence: 15-19 ('*omicron[1]*')')

p38<-ggplot(params) + 
  geom_histogram(aes(x=omic1*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.43,1.0,0.028), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.43,1), breaks = c(0.5,0.75,1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p38


#omic2
summary(params$omic2*100)

text.xlab <- bquote("%")
text.title <- bquote('Former PWID Prevalence: 20-25 ('*omicron[2]*')')

p39<-ggplot(params) + 
  geom_histogram(aes(x=omic2*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.71,1,0.014), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.7,1), breaks = c(0.7,0.85,1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p39

#omic3
summary(params$omic3*100)

text.xlab <- bquote("%")
text.title <- bquote('Former PWID Prevalence: 26-29 ('*omicron[3]*')')

p40<-ggplot(params) + 
  geom_histogram(aes(x=omic3*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.81,2.1,0.06), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.81,2.1), breaks = c(0.9,1.5,2.1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p40

#omic4
summary(params$omic4*100)

text.xlab <- bquote("%")
text.title <- bquote('Former PWID Prevalence: 30-64 ('*omicron[4]*')')

p41<-ggplot(params) + 
  geom_histogram(aes(x=omic4*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.6,2.7,0.055), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.6,2.7), breaks = c(1.6,2.15,2.7))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p41


#r
summary(params$r)
#summary(1/params$r)

text.xlab <- bquote('Surveillance Cases'*~Infections^-1)
text.title <- bquote('Reporting Rate ('*rho*')')

p42<-ggplot(params) +
  geom_histogram(aes(x=r, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.059,0.082,0.001), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.059,0.082), breaks = c(0.06,0.07,0.08))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p42


#psi1
summary(params$psi1*100)

text.xlab <- bquote('%')
text.title <- bquote('Abuse/Dependence Prevalence: 15-19 ('*psi[1]*')')

p43<-ggplot(params) + 
  geom_histogram(aes(x=psi1*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.6,2.1,0.025), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.6,2.1), breaks = c(1.6,1.85,2.1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p43

#psi2
summary(params$psi2*100)

text.xlab <- bquote('%')
text.title <- bquote('Abuse/Dependence Prevalence: 20-25 ('*psi[2]*')')

p44<-ggplot(params) + 
  geom_histogram(aes(x=psi2*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.1,1.5,0.02), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.1,1.5), breaks = c(1.1,1.3,1.5))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p44

#psi3
summary(params$psi3*100)

text.xlab <- bquote('%')
text.title <- bquote('Abuse/Dependence Prevalence: 26-29 ('*psi[3]*')')

p45<-ggplot(params) + 
  geom_histogram(aes(x=psi3*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.37,1,0.03), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.37,1), breaks = c(0.4,0.7,1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p45


#psi4
summary(params$psi4*100)

text.xlab <- bquote('%')
text.title <- bquote('Abuse/Dependence Prevalence: 30-64 ('*psi[4]*')')

p46<-ggplot(params) + 
  geom_histogram(aes(x=psi4*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.24,0.48,0.012), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.24,0.48), breaks = c(0.24,0.36,0.48))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p46


#psiN
summary(params$psin*params$Z0)

text.xlab <- bquote('Persons'*~Years^-1)
text.title <- bquote('New Abuse/Dependence ('*psi[N]*Z[0]*')')

p47<-ggplot(params) + 
  geom_histogram(aes(x=psin*Z0, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1167,4323,150), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1167,4323), breaks = c(1200,2750,4300), labels=c("1,200","2,750","4,300"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p47

#Contact Matrix: c11
summary(params$c1)
#text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 15-19 with 15-19 ('*pi["1,1"]*')')
p48<-ggplot(params) + 
  geom_histogram(aes(x=c1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,2.68,0.08), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,2.68), breaks = c(1,1.84,2.68))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p48

#Contact Matrix: c21
summary(params$c2)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 15-19 with 20-25 ('*pi["2,1"]*')')
p49<-ggplot(params) + 
  geom_histogram(aes(x=c2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,2.68,0.08), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,2.68), breaks = c(1,1.84,2.68))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p49

#Contact Matrix: c31
summary(params$c3)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 15-19 with 26-29 ('*pi["3,1"]*')')
p50<-ggplot(params) + 
  geom_histogram(aes(x=c3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.49,1,0.025), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.49,1), breaks = c(0.49,0.745,1), labels=c("0.49","0.745","1.00"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p50

#Contact Matrix: c41
summary(params$c4)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 15-19 with 30-64 ('*pi["4,1"]*')')
p51<-ggplot(params) + 
  geom_histogram(aes(x=c4, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,1,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p51

#Contact Matrix: c12
summary(params$c5)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 20-25 with 15-19 ('*pi["1,2"]*')')
p52<-ggplot(params) + 
  geom_histogram(aes(x=c5, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,2.68,0.08), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,2.68), breaks = c(1,2.84,2.68))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p52

#Contact Matrix: c22
summary(params$c6)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 20-25 with 20-25 ('*pi["2,2"]*')')
p53<-ggplot(params) + 
  geom_histogram(aes(x=c6, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,2.68,0.08), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,2.68), breaks = c(1,1.84,2.68))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p53

#Contact Matrix: c32
summary(params$c7)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 20-25 with 26-29 ('*pi["3,2"]*')')
p54<-ggplot(params) + 
  geom_histogram(aes(x=c7, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.49,1.00,0.025), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.49,1), breaks = c(0.49,0.745,1), labels=c("0.49","0.745","1.00"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p54

#Contact Matrix: c42
summary(params$c8)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 20-25 with 30-64 ('*pi["4,2"]*')')
p55<-ggplot(params) + 
  geom_histogram(aes(x=c8, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,1,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p55

#Contact Matrix: c13
summary(params$c9)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 26-29 with 15-19 ('*pi["1,3"]*')')
p56<-ggplot(params) + 
  geom_histogram(aes(x=c9, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,6.31,0.25), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,6.31), breaks = c(1,3.655,6.31), labels=c("1.0","3.655","6.31"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p56

#Contact Matrix: c23
summary(params$c10)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 26-29 with 20-25 ('*pi["2,3"]*')')
p57<-ggplot(params) + 
  geom_histogram(aes(x=c10, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,6.31,0.25), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,6.31), breaks = c(1,3.655,6.31), labels=c("1.0","3.655","6.31"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p57

#Contact Matrix: c33
summary(params$c11)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 26-29 with 26-29 ('*pi["3,3"]*')')
p58<-ggplot(params) + 
  geom_histogram(aes(x=c11, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,2.17,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,2.17), breaks = c(1,1.585,2.17), labels=c("1.0","1.585","2.17"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p58

#Contact Matrix: c43
summary(params$c12)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 26-29 with 30-64 ('*pi["4,3"]*')')
p59<-ggplot(params) + 
  geom_histogram(aes(x=c12, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,1.2,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,1.2), breaks = c(0,1,2), labels=c("0.0","1.0","2.0"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p59

#Contact Matrix: c14
summary(params$c13)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 30-64 with 15-19 ('*pi["1,4"]*')')
p60<-ggplot(params) + 
  geom_histogram(aes(x=c13, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,3.02,0.15), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,3.02), breaks = c(0,1.5,3))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p60

#Contact Matrix: c24
summary(params$c14)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 30-64 with 20-25 ('*pi["2,4"]*')')
p61<-ggplot(params) + 
  geom_histogram(aes(x=c14, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,3.02,0.15), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,3.02), breaks = c(0,1.5,3))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p61

#Contact Matrix: c34
summary(params$c15)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 30-64 with 26-29 ('*pi["3,4"]*')')
p62<-ggplot(params) + 
  geom_histogram(aes(x=c15, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.14,1.3,0.05), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.14,1.3), breaks = c(0.14,0.72,1.3))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p62

#Contact Matrix: c44
summary(params$c16)
text.xlab <- bquote("Relative Deviation from Proportionality")
text.title <- bquote('Assortativity: 30-64 with 30-64 ('*pi["4,4"]*')')
p63<-ggplot(params) + 
  geom_histogram(aes(x=c16, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.1,7.8,0.35), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.1,7.8), breaks = c(0.1,3.95,7.8))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p63


#PopSize1
# summary(params$Popsize1)
# text.xlab <- bquote("Persons")
# text.title <- bquote('Population Size: 15-19 ('*~P[1]*')')
# p80<-ggplot(params) + 
#   geom_histogram(aes(x=Popsize1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(668000,766000,4500), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab(text.ylab) +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(668000,766000), breaks = c(700000,750000),labels=c("700,000","750,000"))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p80
# 
# #PopSize2
# summary(params$Popsize2)
# text.xlab <- bquote("Persons")
# text.title <- bquote('Population Size: 20-25 ('*~P[2]*')')
# p81<-ggplot(params) + 
#   geom_histogram(aes(x=Popsize2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(771000,868000,4500), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(771000,868000), breaks = c(800000,850000),labels=c("800,000","850,000"))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p81
# 
# #PopSize3
# summary(params$Popsize3)
# text.xlab <- bquote("Persons")
# text.title <- bquote('Population Size: 26-29 ('*~P[3]*')')
# p82<-ggplot(params) + 
#   geom_histogram(aes(x=Popsize3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(467000,527000,2700), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(467000,527000), breaks = c(475000,500000,525000),labels=c("475,000","500,000","525,000"))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p82
# 
# #PopSize4
# summary(params$Popsize4)
# text.xlab <- bquote("Persons")
# text.title <- bquote('Population Size: 30-64 ('*~P[4]*')')
# p83<-ggplot(params) + 
#   geom_histogram(aes(x=Popsize4, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(4491000,4684000,9000), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(4491000,4684000), breaks = c(4500000,4575000,4650000),labels=c("4,500,000","4,575,000","4,650,000"))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p83
# 
# 
# #RSS 
summary(params$RSS)
text.xlab <- "Residual"
text.title<-"Residual Sum of Squares"
p84<-ggplot(params) +
  geom_histogram(aes(x=RSS, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', show.legend=TRUE, bins=40)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(breaks=c(300,1650,3000), labels=c("300","1,650","3,000"))+
  #scale_x_log10(breaks=c(200,2000,20000),labels=c("200","2,000","20,000"))+#limits = c(5,815), breaks = c(10,100))+
  scale_y_continuous(breaks = c(0,500,1000), labels=c("0","500","1,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p84

# p85<-ggplot(params) +
#   geom_histogram(aes(x=RSS, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', show.legend=TRUE, bins=40)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) +
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab(text.ylab) +
#   xlab(text.xlab) +
#   scale_x_log10(breaks=c(200,2000,20000),labels=c("200","2,000","20,000"))+#limits = c(5,815), breaks = c(10,100))+
#   scale_y_continuous(limits = c(0, 3000), breaks = c(0,1000,2000,3000), labels=c("0","1,000","2,000","3,000"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p85

#Saving main and supplemental figs

#Note: quartz works for saving unicode text as pdf (ggsave does not)

#All main params
w<-22; h<-12

 #with bottom legend
lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),
            c(13,13,13,13))

options(device = "quartz")
quartz()
LHS_Histograms1<-grid.arrange(p3,p13,p37,p42,
                              p21,p22,p23,p24,
                              p4,p5,p6,p7,bottomleg,layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_Params1_12_leg.pdf", "pdf",width=w,height=h)

lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),
            c(13,14,15,16),c(13,14,15,16),c(13,14,15,16),c(17,17,17,17))
options(device = "quartz")
quartz()
LHS_Histograms2<-grid.arrange(p14,p15,p16,p17,
                              p38,p39,p40,p41,
                              p43,p44,p45,p46,
                              p29,p30,p31,p32,bottomleg,layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_Params13_28_leg.pdf", "pdf",width=w,height=h)

#Assortativity of Contacts (sampled c from Polymod and sampled to 0)
lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),
            c(13,14,15,16),c(13,14,15,16),c(13,14,15,16),c(17,17,17,17))
options(device = "quartz")
quartz()
LHS_Histograms4<-grid.arrange(p48,p52,p56,p60,
                              p49,p53,p57,p61,
                              p50,p54,p58,p62,
                              p51,p55,p59,p63,bottomleg,layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_SampledContacts.pdf", "pdf",width=w,height=h)


#RSS
w<-8; h<-6
ggsave(sprintf("LHS_RSS.pdf"), p84,width=w, height=h)




############### Section 5: Acute HCV Data ################
#### Acute Data 15-29 YOs (MDHHS Data) ####
year<-seq(2000,2016,1)
acutepwid1519<-c(4,	1,	1,	1,	4,	1,	2,	7,	5,	3,	2,	2,	5,	2,	2,	1,	2*(1-0.26))
acute2025idu <- c(1,	2,	3,	7,	8,	3,	10,	9,	17,	11,	9,	7,	11,	16,	19,	18,	29*(1-0.26))
acute2629idu <- c(4,	6,	4,	5,	2,	5,	8,	4,	9,	5,	9,	3,	8,	11,	8,	13,	14*(1-0.26))

data<-data.frame(cbind(year,acutepwid1519,acute2025idu,acute2629idu))

xbreaks = c("15-19 Years", "20-25 Years","26-29 Years")

p1<-ggplot(data)+
  geom_line(aes(x=year, y=acutepwid1519),linetype="dashed",colour="#B8361B",size=1,show.legend = FALSE) +
  geom_line(aes(x=year, y=acute2025idu),linetype="dashed",colour="sea green",size=1,show.legend = FALSE)+
  geom_line(aes(x=year, y=acute2629idu),linetype="dashed",colour="#152FA5",size=1,show.legend = FALSE)+
  geom_point(aes(x=year, y=acutepwid1519,colour="15-19 Years"), size=3,show.legend = TRUE, na.rm=TRUE) + 
  geom_point(aes(x=year, y=acute2025idu,colour="20-25 Years"), size=3,show.legend = TRUE, na.rm=TRUE)+
  geom_point(aes(x=year, y=acute2629idu,colour="26-29 Years"), size=3,show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = c(2000,2003,2006,2009,2012,2015))+
  scale_y_continuous(limits = c(0, max(acute2025idu)+1))+
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_manual(name=NULL, values = c("#B8361B","sea green","#152FA5"),breaks = xbreaks,
                      #labels=str_wrap(xbreaks,15))+
                      labels=xbreaks)+
  #guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("Number of Cases") +
  xlab("Year") +
  #ggtitle('Acute and New Chronic Cases Aged 15-19')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1

p2<-ggplot(data)+
  geom_line(aes(x=year, y=acutepwid1519+acute2025idu+acute2629idu),linetype="dashed",colour="#000000",size=1,show.legend = FALSE) +
  #geom_line(aes(x=year, y=acute2025idu),colour="#000000",size=1,show.legend = FALSE)+
  #geom_line(aes(x=year, y=acute2629idu),colour="#000000",size=1,show.legend = FALSE)+
  geom_point(aes(x=year, y=acutepwid1519+acute2025idu+acute2629idu), colour="#000000",size=3,show.legend = TRUE, na.rm=TRUE) + 
  #geom_point(aes(x=year, y=acute2025idu,colour="20-25 Years"), size=1.5,show.legend = TRUE, na.rm=TRUE)+
  #geom_point(aes(x=year, y=acute2629idu,colour="26-29 Years"), size=1.5,show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = c(2000,2003,2006,2009,2012,2015))+
  scale_y_continuous(limits = c(0, max(acutepwid1519+acute2025idu+acute2629idu)+1))+
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  #scale_colour_manual(name=NULL, values = c("#B8361B","sea green","#152FA5"),breaks = xbreaks,
  #labels=str_wrap(xbreaks,15))+
  #                    labels=xbreaks)+
  #guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("Number of Reported Cases") +
  xlab("Year") +
  #ggtitle('Acute and New Chronic Cases Aged 15-19')+
  theme(plot.title = element_text(hjust = 0.5)) 
p2





