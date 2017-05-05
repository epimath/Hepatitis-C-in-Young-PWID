## Figure Script: Gicquelais et al. (2017). Hepatitis C transmission model.

## Before running this code, run the model simulations in the code HCV_MS_04202017.m to simulate
## the hepatitis C model and output results as csv or txt files. Use this script creates figures 
## summarizing results and similar to those published in the article. 

################ Section 1: Packages and Functions to Load First ################### 
library(ggplot2)
library(grid)
library(gridExtra)
library(car)
library(plotly)
library(plyr)
library(reshape)

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
setwd("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 4.20.2017")


############### Section 3: Importing Data Generated from Matlab ################
#### Import MDHHS Data ####

#MDHHS case counts by year
mdhhs<-read.csv("MDHHSData.csv", header=F)
colnames(mdhhs) <- c("Year",	"acute1530_1",	"acute1530_2",	"acute1530_3",	"acute1530idu_1",	
                          "acute1530idu_2",	"acute1530idu_3",	"chronic1530_1",	"chronic1530_2",		
                          "chronic1530_3",	"chronic1530idu_1",	"chronic1530idu_2",	"chronic1530idu_3",	
                          "unknown1530_1",	"unknown1530_2", "unknown1530_3",	
                          "chronicunk1530_1",	"chronicunk1530_2",	"chronicunk1530_3")
summary(mdhhs)
mdhhs$chronic1530idu<-mdhhs$chronic1530idu_1+mdhhs$chronic1530idu_2+mdhhs$chronic1530idu_3
mdhhs$acute1530idu<-mdhhs$acute1530idu_1+mdhhs$acute1530idu_2+mdhhs$acute1530idu_3

#### Import New Chronic Cases (Current and Former PWID) ####

#All New Chronic Cases
#note that header row (var names) is the X + index of the parameter set (eg XParamSetNo.)
LHS_NewChronic<-read.csv("LHS_NewChronic.csv", header=TRUE)
summary(LHS_NewChronic)
LHS_NewChronic$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)

keeps <- c("Year","chronic1530idu")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_NewChronic, id = "Year")
LHS.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(LHS.melted)

#New Chronic Cases - 15-19
LHS_NewChronic1<-read.csv("LHS_NewChronic1.csv", header=TRUE)
summary(LHS_NewChronic1)
LHS_NewChronic1$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)

keeps <- c("Year","chronic1530idu_1")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_NewChronic1, id = "Year")
LHS1.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(LHS1.melted)


#New Chronic Cases - 20-24
LHS_NewChronic2<-read.csv("LHS_NewChronic2.csv", header=TRUE)
summary(LHS_NewChronic2)
LHS_NewChronic2$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)

keeps <- c("Year","chronic1530idu_2")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_NewChronic2, id = "Year")
LHS2.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(LHS2.melted)

#New Chronic Cases - 25-30
LHS_NewChronic3<-read.csv("LHS_NewChronic3.csv", header=TRUE)
summary(LHS_NewChronic3)
LHS_NewChronic3$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)

keeps <- c("Year","chronic1530idu_3")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_NewChronic3, id = "Year")
LHS3.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(LHS3.melted)

#### Import Chronic HCV Prevalence Data (Current and Former PWID) ####

#Note that 1st row = parameter index (and becomes variable name)

#Total Chronic Cases - All Ages
LHS_ChronicPrev<-read.csv("LHS_ChronicPrev.csv", header=TRUE)
summary(LHS_ChronicPrev)
LHS_ChronicPrev$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)
LHS_C<-melt(LHS_ChronicPrev, id = "Year")

#Total Chronic Cases - 15-19

LHS_ChronicPrev1<-read.csv("LHS_ChronicPrev1.csv", header=TRUE)
summary(LHS_ChronicPrev1)
LHS_ChronicPrev1$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)
LHS1_C<-melt(LHS_ChronicPrev1, id = "Year")

#Total Chronic Cases - 20-24

LHS_ChronicPrev2<-read.csv("LHS_ChronicPrev2.csv", header=TRUE)
summary(LHS_ChronicPrev2)
LHS_ChronicPrev2$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)
LHS2_C<-melt(LHS_ChronicPrev2, id = "Year")

#Total Chronic Cases - 25-30

LHS_ChronicPrev3<-read.csv("LHS_ChronicPrev3.csv", header=TRUE)
summary(LHS_ChronicPrev3)
LHS_ChronicPrev3$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)
LHS3_C<-melt(LHS_ChronicPrev3, id = "Year")



#### Import Acute Cases (Current and Former PWID) ####

#Note that header row (var names) is the X + index of the parameter set (eg XParamSetNo.)

#Acute Cases - All Ages

LHS_Acute<-read.csv("LHS_Acute.csv", header=TRUE)
summary(LHS_Acute)
LHS_Acute$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)

keeps <- c("Year","acute1530idu")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute, id = "Year")
AcuteLHS.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted)

#Acute Cases - 15-19

LHS_Acute1<-read.csv("LHS_Acute1.csv", header=TRUE)
summary(LHS_Acute1)
LHS_Acute1$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)

keeps <- c("Year","acute1530idu_1")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute1, id = "Year")
AcuteLHS.melted_1<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_1)

#Acute Cases - 20-24

LHS_Acute2<-read.csv("LHS_Acute2.csv", header=TRUE)
summary(LHS_Acute2)
LHS_Acute2$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)

keeps <- c("Year","acute1530idu_2")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute2, id = "Year")
AcuteLHS.melted_2<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_2)

#Acute Cases - 25-30

LHS_Acute3<-read.csv("LHS_Acute3.csv", header=TRUE)
summary(LHS_Acute3)
LHS_Acute3$Year[(1:66)]<-seq(from = 2000, to = 2013, by = 0.2)

keeps <- c("Year","acute1530idu_3")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute3, id = "Year")
AcuteLHS.melted_3<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_3)



#### Import LHS Parameter Data #####

params<-read.csv("LHS_ParamRSS.txt", header=F)
colnames(params) <- c("index",	"b",	"phi1",	"phi2",	"phi3","d",	"k",	"eps",	"s",	"gn",	"gp",	
                          "a",	"Z0",	"r",	"lambda1",	"lambda2",	
                          "mu1",	"mu2", "mu3",	"etap",	"etan",	"intervention_etap",
                          "c1",	"c2",	"c3",	"c4",	"c5",	"c6",	"c7",	"c8",	"c9",	"RSS",
                          "init_b","init_phi1","init_phi2","init_phi3")
summary(params)

#Create an indicator for the best fitting 50% of parameter sets
params$best50<-ifelse(params$RSS<=median(params$RSS),1,0)
table(params$best50)

#Create an indicator for the quartile of fit
params$fitquartile<-cut(params$RSS, quantile(params$RSS, c(0,1/4,2/4,3/4,1)), include.lowest=T) 
table(params$fitquartile)

#Create a numeric quartile variable (note: change the bounds as needed)
params$fitorder<-ifelse(params$fitquartile=="[3.14,3.92]",1,
                 ifelse(params$fitquartile=="(3.92,4.36]",2,
                 ifelse(params$fitquartile=="(4.36,4.87]",3,4)))
table(params$fitorder) 

#### Dataset of case counts with params ####

paramsRSS<-params
paramsRSS$variable<-paste("X", paramsRSS$index, sep="")

LHS.melted<-merge(LHS.melted,paramsRSS,by="variable")
LHS_C <-merge(LHS_C,paramsRSS,by="variable")
AcuteLHS.melted<-merge(AcuteLHS.melted,paramsRSS,by="variable")

LHS1.melted<-merge(LHS1.melted,paramsRSS,by="variable")
LHS1_C <-merge(LHS1_C,paramsRSS,by="variable")
AcuteLHS.melted_1<-merge(AcuteLHS.melted_1,paramsRSS,by="variable")

LHS2.melted<-merge(LHS2.melted,paramsRSS,by="variable")
LHS2_C <-merge(LHS2_C,paramsRSS,by="variable")
AcuteLHS.melted_2<-merge(AcuteLHS.melted_2,paramsRSS,by="variable")

LHS3.melted<-merge(LHS3.melted,paramsRSS,by="variable")
LHS3_C <-merge(LHS3_C,paramsRSS,by="variable")
AcuteLHS.melted_3<-merge(AcuteLHS.melted_3,paramsRSS,by="variable")

rm(mdhhs1,keeps, paramsRSS)



#### Import intervention datasets and combine ####

#import datasets with intervention levels + all case counts at t1-t66
prev<-read.csv("ChronicPrev_SingleSeq.txt", header=F)
names(prev)[1:11] <- c("e","j",	"l",	"m",	"o",	"q",	"u",	"f",	"h",	"k",	"i")
names(prev)[12:77] <- paste("t",1:66,sep="")
prev$inttype<-ifelse(prev$e<=68,"single",ifelse(prev$e>=119,"terttoprim","primtotert"))
table(prev$inttype)

new<-read.csv("NewChronic_SingleSeq.txt", header=F)
names(new)[1:11] <- c("e","j",	"l",	"m",	"o",	"q",	"u",	"f",	"h",	"k",	"i")
names(new)[12:77] <- paste("t",1:66,sep="")
new$inttype<-ifelse(new$e<=68,"single",ifelse(new$e>=119,"terttoprim","primtotert"))
table(new$inttype)

acute<-read.csv("Acute_SingleSeq.txt", header=F)
names(acute)[1:11] <- c("e","j",	"l",	"m",	"o",	"q",	"u",	"f",	"h",	"k",	"i")
names(acute)[12:77] <- paste("t",1:66,sep="")
acute$inttype<-ifelse(acute$e<=68,"single",ifelse(acute$e>=119,"terttoprim","primtotert"))
table(acute$inttype)

#Import parameters associated with each intervention
param<-read.csv("Parameters_SingleSeq.txt", header=F)
names(param)[1:11] <- c("e","j",	"l",	"m",	"o",	"q",	"u",	"f",	"h",	"k",	"i")
names(param)[12:32] <- c("b",	"phi1",	"phi2",	"phi3","d",	"kappa",	"eps",	"s",	"gn",	"gp",	
                         "a",	"Z0",	"r",	"lambda1",	"lambda2","mu1",	"mu2",	"mu3",	"etap",	"etan","intervention_deathcurr")
names(param)[33:41] <- paste("c",1:9,sep="")
param$inttype<-ifelse(param$e<=68,"single",ifelse(param$e>=119,"terttoprim","primtotert"))
table(param$inttype)

#Create datasets with case counts at t66 summary
prev1<-prev[ ,c(1:11,77:78)]
colnames(prev1)[12]<-"ChronicPrev_End" 
prev1$index<-as.numeric(paste(prev1$e,prev1$i,sep=""))
new1<-new[ ,c(1,11,77)]
colnames(new1)[3]<-"NewChronic_End"
new1$index<-as.numeric(paste(new1$e,new1$i,sep=""))
acute1<-acute[ ,c(1,11,77)]
colnames(acute1)[3]<-"Acute_End"
acute1$index<-as.numeric(paste(acute1$e,acute1$i,sep=""))

param1<-param[ ,c(1,11:41)]
param1$index<-as.numeric(paste(param1$e,param1$i,sep=""))

merged<-merge(merge(merge(prev1, new1, by=c("e","i"),suffixes = c(".a",".b")), 
                    acute1, by=c("e","i"), suffixes = c(".c",".d")), 
              param1, by=c("e","i"), suffixes = c(".e",""))


#Create variables for intervention level
df<-merged
attach(df)
int<-j
df$gnpct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-l
df$gppct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-m
df$spct<-ifelse(int==1,"1 Year",ifelse(int==2,"24 Weeks",ifelse(int==3,"16 Weeks",ifelse(int==4,"12 Weeks","error"))))
int<-o
df$dpct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-q
df$kpct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-u
df$phipct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-f
df$conpct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-h
df$etappct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-k
df$etanpct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
detach(df)
merged<-df

table(merged$gnpct)
table(merged$gppct)
table(merged$spct)
table(merged$dpct)
table(merged$kpct)
table(merged$phipct)
table(merged$conpct)
table(merged$etappct)
table(merged$etanpct)

#Add RSS and fit information to the merged dataset
paramfit<-params[,c(1,32,37,38)]
colnames(paramfit) <- c("i","RSS","best50","fitquartile")

merged<-merge(merged,paramfit, by="i",all.x=TRUE)

#Create dataset for single interventions 
singleint<-merged[ ! merged$inttype %in% c("primtotert","terttoprim"), ]

#Create dataset for sequential interventions 
seqint<-merged[ ! merged$inttype %in% c("single"), ]

#Identify the intervention parameter in singleint
singleint$InterventionParam<-ifelse(singleint$e==1,"None",
                             ifelse(2<=singleint$e&singleint$e<=5,"gn",
                             ifelse(6<=singleint$e&singleint$e<=9,"gp",
                             ifelse(10<=singleint$e&singleint$e<=12,"s",
                             ifelse(13<=singleint$e&singleint$e<=16,"k",
                             ifelse(17<=singleint$e&singleint$e<=20,"d",
                             ifelse(21<=singleint$e&singleint$e<=24,"phi",
                             ifelse(25<=singleint$e&singleint$e<=28,"contact",
                             ifelse(29<=singleint$e&singleint$e<=32,"etap",
                             ifelse(33<=singleint$e&singleint$e<=36,"etan", 
                             ifelse(37<=singleint$e&singleint$e<=40,"gn_s",
                             ifelse(41<=singleint$e&singleint$e<=44,"gp_s",
                             ifelse(45<=singleint$e&singleint$e<=48,"k_s",
                             ifelse(49<=singleint$e&singleint$e<=52,"d_s",
                             ifelse(53<=singleint$e&singleint$e<=56,"phi_s",
                             ifelse(57<=singleint$e&singleint$e<=60,"contact_s",
                             ifelse(61<=singleint$e&singleint$e<=64,"etap_s",
                             ifelse(65<=singleint$e&singleint$e<=68,"etan_s","wrong")))))))))))))))))) 
table(singleint$InterventionParam)

#Identify the intervention level in singleint
singleint$IntLevelpct<-ifelse(singleint$e==1,"None",
                       ifelse(singleint$e %in% c(2,6,13,17,21,25,29,33,37,41,45,49,53,57,61,65),"10%",
                       ifelse(singleint$e %in% c(3,7,14,18,22,26,30,34,38,42,46,50,54,58,62,66),"20%",
                       ifelse(singleint$e %in% c(4,8,15,19,23,27,31,35,39,43,47,51,55,59,63,67),"30%",
                       ifelse(singleint$e %in% c(5,9,16,20,24,28,32,36,40,44,48,52,56,60,64,68),"40%","s")))))
table(singleint$IntLevelpct)

#Set an order variable to determine plotting order in graphs
singleint$order<-factor(ifelse(singleint$InterventionParam=="gn",1,ifelse(singleint$InterventionParam=="gn_s",2,
                  ifelse(singleint$InterventionParam=="gp",3,ifelse(singleint$InterventionParam=="gp_s",4,
                  ifelse(singleint$InterventionParam=="k",5,ifelse(singleint$InterventionParam=="d",6,
                  ifelse(singleint$InterventionParam=="contact",7, ifelse(singleint$InterventionParam=="phi",8,
                  ifelse(singleint$InterventionParam=="etap",9,ifelse(singleint$InterventionParam=="etan",10,
                  ifelse(singleint$InterventionParam=="None",11,9999))))))))))))
table(singleint$order)



#Seqint Formatting

#Create variable for intervention level
seqint$pct<-ifelse(seqint$e %in% c(69,94,119,144),"None",
            ifelse(seqint$e %in% c(70,74,78,82,86,90,95,99,103,107,111,115,
                                   120,124,128,132,136,140,145,149,153,157,161,165),"10%",
            ifelse(seqint$e %in% c(71,75,79,83,87,91,96,100,104,108,112,116,
                                   121,125,129,133,137,141,146,150,154,158,162,166),"20%",
            ifelse(seqint$e %in% c(72,76,80,84,88,92,97,101,105,109,113,117,
                                   122,126,130,134,138,142,147,151,155,159,163,167),"30%",
            ifelse(seqint$e %in% c(73,77,81,85,89,93,98,102,106,110,114,118,
                                   123,127,131,135,139,143,148,152,156,160,164,168),"40%","s")))))
table(seqint$pct)

#Create variable for order of interventions when plotting
seqint$intono<-ifelse(seqint$e %in% c(69,94,119,144),1,
               ifelse(seqint$e %in% c(70:73,95:98,120:123,145:148),2,
               ifelse(seqint$e %in% c(74:77,99:102,124:127,149:152),3,
               ifelse(seqint$e %in% c(78:81,103:106,128:131,153:156),4,
               ifelse(seqint$e %in% c(82:85,107:110,132:135,157:160),5,
               ifelse(seqint$e %in% c(86:89,111:114,136:139,161:164),6,
               ifelse(seqint$e %in% c(90:93,115:118,140:143,165:168),7,9999)))))))
table(seqint$intono)

#select the subset dataset of 'none' interventions to calculate % reduction in seqint
options(digits = 10)
noint<-subset(seqint,seqint$e %in% c(69,94))
table(noint$e)
noint <- noint[order(noint$e),] 
df<-noint[1:5000,c(1,12,15,17)]
df<-rename(df, c("ChronicPrev_End" = "ChronicPrev_e69", "Acute_End" = "Acute_e69", "NewChronic_End" = "NewChronic_e69"))
df$ChronicPrev_e94<-subset(noint,noint$e==94)$ChronicPrev_End
df$Acute_e94<-subset(noint,noint$e==94)$Acute_End
df$NewChronic_e94<-subset(noint,noint$e==94)$NewChronic_End


#merge back with intervention dataset 
seqint<-merge(seqint,df,by="i")
rm(df,noint)

# Calculate % reduction compared to appropriate none intervention
seqint$pctred_chr<-ifelse(seqint$e %in% c(69,94,119,144),0,
                   ifelse(seqint$e %in% c(70:93,120:143),100-(seqint$ChronicPrev_End/seqint$ChronicPrev_e69)*100,
                   ifelse(seqint$e %in% c(95:118,145:168),100-(seqint$ChronicPrev_End/seqint$ChronicPrev_e94)*100,99999)))

seqint$pctred_new<-ifelse(seqint$e %in% c(69,94,119,144),0,
                   ifelse(seqint$e %in% c(70:93,120:143),100-(seqint$NewChronic_End/seqint$NewChronic_e69)*100,
                   ifelse(seqint$e %in% c(95:118,145:168),100-(seqint$NewChronic_End/seqint$NewChronic_e94)*100,99999)))

seqint$pctred_acute<-ifelse(seqint$e %in% c(69,94,119,144),0,
                     ifelse(seqint$e %in% c(70:93,120:143),100-(seqint$Acute_End/seqint$Acute_e69)*100,
                     ifelse(seqint$e %in% c(95:118,145:168),100-(seqint$Acute_End/seqint$Acute_e94)*100,99999)))

summary(seqint$pctred_chr)
summary(seqint$pctred_new)
summary(seqint$pctred_acute)

#Calculate Summary Statistics of % Reduction
df1<-tapply(seqint$pctred_chr,list(seqint$e),mean)
df2<-tapply(seqint$pctred_chr,list(seqint$e),sd)
df3<-tapply(seqint$pctred_chr,list(seqint$e),min)
df4<-tapply(seqint$pctred_chr,list(seqint$e),max)
df5<-tapply(seqint$pctred_new,list(seqint$e),mean)
df6<-tapply(seqint$pctred_new,list(seqint$e),sd)
df7<-tapply(seqint$pctred_new,list(seqint$e),min)
df8<-tapply(seqint$pctred_new,list(seqint$e),max)

seqint1<-data.frame(cbind(df1,df2,df3,df4,df5,df6,df7,df8))
seqint1$e<-rownames(seqint1)
seqint1<-rename(seqint1,c("df1"="meanpctred_chr","df2"="sdpctred_chr","df3"="minpctred_chr",
                     "df4"="maxpctred_chr","df5"="meanpctred_new","df6"="sdpctred_new",
                     "df7"="minpctred_new", "df8"="maxpctred_new"))


# Subset information on each 'e' intervention
interventioninfo <- seqint[,c(2:11,13,50:58,62,63)]
interventioninfo<-unique(interventioninfo)

# Merge with seqint1
seqint1<-merge(seqint1,interventioninfo,by="e")


# % Reduction in Median
summary(seqint1$meanpctred_chr)
summary(seqint1$meanpctred_new)
summary(seqint1$minpctred_chr)
summary(seqint1$minpctred_new)
summary(seqint1$maxpctred_chr)
summary(seqint1$maxpctred_new)

#Summarize only Best 50% Fitting Parameter Sets
best50<-subset(seqint,seqint$best50==1)

#Calculate Summary Statistics of % Reduction
df1<-tapply(best50$pctred_chr,list(best50$e),mean)
df2<-tapply(best50$pctred_chr,list(best50$e),sd)
df3<-tapply(best50$pctred_chr,list(best50$e),min)
df4<-tapply(best50$pctred_chr,list(best50$e),max)
df5<-tapply(best50$pctred_new,list(best50$e),mean)
df6<-tapply(best50$pctred_new,list(best50$e),sd)
df7<-tapply(best50$pctred_new,list(best50$e),min)
df8<-tapply(best50$pctred_new,list(best50$e),max)

seqint2<-data.frame(cbind(df1,df2,df3,df4,df5,df6,df7,df8))
seqint2$e<-rownames(seqint2)
seqint2<-rename(seqint2,c("df1"="meanpctred_chr","df2"="sdpctred_chr","df3"="minpctred_chr",
                          "df4"="maxpctred_chr","df5"="meanpctred_new","df6"="sdpctred_new",
                          "df7"="minpctred_new", "df8"="maxpctred_new"))


# Subset information on each 'e' intervention
interventioninfo <- seqint[,c(2:11,13,50:58,62,63)]
interventioninfo<-unique(interventioninfo)

# Merge with seqint1
seqint2<-merge(seqint2,interventioninfo,by="e")


# % Reduction in Median
summary(seqint2$meanpctred_chr)
summary(seqint2$meanpctred_new)
summary(seqint2$minpctred_chr)
summary(seqint2$minpctred_new)
summary(seqint2$maxpctred_chr)
summary(seqint2$maxpctred_new)




############### Section 4: Latin Hypercube Sampling Results ################
#### Plot of New Chronic Cases Fit to MDHHS Data (50% Best, 1 Color) #####
#Plotting Parameters
xbreaks <- c(2000, 2005, 2010)
ybreaks <- c(0, 400, 800, 1200)
ylabels <- c("0","400","800","1,200")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

#Create a plot for all age groups combined and each age group separately

#All ages - with right legend
p1<-ggplot(subset(LHS.melted,LHS.melted$best50==1), aes(x=Year, y=value, by=variable)) + 
  geom_line(aes(x=Year, y=value, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=chronic1530idu, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 1300), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "right", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_manual(name=NULL, values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) 

#All ages - no legend 
p2<-ggplot(subset(LHS.melted,LHS.melted$best50==1), aes(x=Year, y=value, by=variable)) + 
  geom_line(aes(x=Year, y=value, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=chronic1530idu, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 1300), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name="", values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) 

title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p2 <- arrangeGrob(p2, top = title.grob)
grid.arrange(p2)

#15-19 
p3<-ggplot(subset(LHS1.melted,LHS1.melted$best50==1), aes(x=Year, y=value, by=variable)) + 
  geom_line(aes(x=Year, y=value, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=chronic1530idu_1, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 700), breaks = c(0,200,400,600))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name="", values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) 

title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p3 <- arrangeGrob(p3, top = title.grob)
grid.arrange(p3)

#20-24
p4<-ggplot(subset(LHS2.melted,LHS2.melted$best50==1), aes(x=Year, y=value, by=variable)) + 
  geom_line(aes(x=Year, y=value, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=chronic1530idu_2, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 700), breaks = c(0,200,400,600))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name="", values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab)

title.grob <- textGrob(label = "C)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p4 <- arrangeGrob(p4, top = title.grob)
grid.arrange(p4)

#25-30
p5<-ggplot(subset(LHS3.melted,LHS3.melted$best50==1), aes(x=Year, y=value, by=variable)) + 
  geom_line(aes(x=Year, y=value, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=chronic1530idu_3, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 700), breaks = c(0,200,400,600))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name="", values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) 

title.grob <- textGrob(label = "D)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p5 <- arrangeGrob(p5, top = title.grob)
grid.arrange(p5)

#compile and save as 1 figure
w <- 10; h <- 6
mylegend<-g_legend(p1)
lay<-rbind(c(1,1,1,1,1,2),c(1,1,1,1,1,2),c(1,1,1,1,1,2),c(3,3,4,4,5,5),c(3,3,4,4,5,5),c(3,3,4,4,5,5))
LHS_NewChronic<-grid.arrange(p2,mylegend,p3,p4,p5,layout_matrix=lay)
ggsave(sprintf("LHS_NewChronic.pdf"), LHS_NewChronic, width=w, height=h)


#### Plot of New Chronic Cases Fit to MDHHS Data (All, Color by RSS) #####
#Plotting Parameters
xbreaks <- c(2000, 2005, 2010)
ybreaks <- c(0, 400, 800, 1200)
ylabels <- c("0","400","800","1,200")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

#Dummy plots for legends:
dummy<-ggplot(data=params,aes(x=d,y=k,color=RSS))+
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

dummy1<-ggplot(LHS.melted)+
  geom_point(aes(x=Year, y=chronic1530idu,colour="Data"),show.legend = TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 1300), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  scale_colour_manual(name=NULL, values = c("#000000"), breaks = c("Data")) +  
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Age Groups')+
  theme(plot.title = element_text(hjust = 0.5))
Datalegend<-g_legend(dummy1)

#Plot of fit and data
p1<-ggplot(LHS.melted, aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=chronic1530idu),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 1450), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) 

title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p1 <- arrangeGrob(p1, top = title.grob)
grid.arrange(p1)

p2<-ggplot(LHS1.melted, aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=chronic1530idu_1),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 820), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) 

title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p2 <- arrangeGrob(p2, top = title.grob)
grid.arrange(p2)

p3<-ggplot(LHS2.melted, aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=chronic1530idu_2),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 820), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) 

title.grob <- textGrob(label = "C)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p3 <- arrangeGrob(p3, top = title.grob)
grid.arrange(p3)

p4<-ggplot(LHS3.melted, aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=chronic1530idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 820), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) 

title.grob <- textGrob(label = "D)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p4 <- arrangeGrob(p4, top = title.grob)
grid.arrange(p4)

#layout w/legend at side 
w <- 10; h <- 6

lay<-rbind(c(1,1,1,1,1,2),c(1,1,1,1,1,3),c(4,4,5,5,6,6),c(4,4,5,5,6,6))
LHS_NewChronic<-grid.arrange(p1,legend1,Datalegend,p2,p3,p4,layout_matrix=lay)
ggsave(sprintf("LHS_NewChronic_RSSColor.pdf"), LHS_NewChronic, width=w, height=h)

#### Plot of New Chronic Cases Fit to MDHHS Data (All, Color by Phi_i) #####
summary(params$phi1)
text.leg <- bquote('log('*phi[1]*')')
dummy1<-ggplot(data=params,aes(x=1/d,y=1/k,color=log(phi1)))+
  geom_point(size=2)+
  ylab(text.ylab)+
  xlab(text.xlab)+
  theme_bw(base_size=18)+
  theme(panel.border = element_rect(linetype = "solid", colour = "black"))+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.position = "bottom",legend.key = element_blank())+
  scale_colour_gradientn(name=text.leg,colours=rainbow(5))
legend1<-g_legend(dummy1)

summary(params$phi2)
text.leg <- bquote('log('*phi[2]*')')
dummy1<-ggplot(data=params,aes(x=1/d,y=1/k,color=log(phi2)))+
  geom_point(size=2)+
  ylab(text.ylab)+
  xlab(text.xlab)+
  theme_bw(base_size=18)+
  theme(panel.border = element_rect(linetype = "solid", colour = "black"))+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.position = "bottom",legend.key = element_blank())+
  scale_colour_gradientn(name=text.leg,colours=rainbow(5))
legend2<-g_legend(dummy1)

summary(params$phi3)
text.leg <- bquote('log('*phi[3]*')')
dummy1<-ggplot(data=params,aes(x=1/d,y=1/k,color=log(phi3)))+
  geom_point(size=2)+#, alpha=0.5)+
  ylab(text.ylab)+
  xlab(text.xlab)+
  theme_bw(base_size=18)+
  theme(panel.border = element_rect(linetype = "solid", colour = "black"))+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.position = "bottom",legend.key = element_blank())+
  scale_colour_gradientn(name=text.leg,colours=rainbow(5))
legend3<-g_legend(dummy1)

summary(LHS1.melted$value)
p2<-ggplot(LHS1.melted, aes(x=Year, y=value, by=variable, color=log(phi1))) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=chronic1530idu_1),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 150), breaks = c(0,50,100,150), labels=c("0","50","100","150"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="ks1",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('15-19 Years')+
  theme(plot.title = element_text(hjust = 0.5))

summary(LHS2.melted$value)
p3<-ggplot(LHS2.melted, aes(x=Year, y=value, by=variable, color=log(phi2))) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=chronic1530idu_2),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 600), breaks = c(0,200,400,600), labels=c("0","200","400","600"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="ks2",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('20-24 Years')+
  theme(plot.title = element_text(hjust = 0.5))

p4<-ggplot(LHS3.melted, aes(x=Year, y=value, by=variable, color=log(phi3))) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=chronic1530idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 820), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="ks3",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('25-30 Years')+
  theme(plot.title = element_text(hjust = 0.5))


w <- 15; h <- 4

lay<-rbind(c(1,2,3),c(1,2,3),c(1,2,3),c(1,2,3),c(1,2,3),c(1,2,3),c(4,5,6))
LHS_NewChronic<-grid.arrange(p2,p3,p4,legend1,legend2,legend3,layout_matrix=lay)
ggsave(sprintf("LHS_NewChronic_phiColor.pdf"), LHS_NewChronic, width=w, height=h)

#### Chronic HCV Prevalence Plots ####

#Plotting Parameters
xbreaks <- c(2000, 2005, 2010)
ybreaks <- c(0, 1000,2000,3000,4000)
ylabels <- c("0","1,000","2,000","3,000","4,000")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000')
max(LHS_C$value)

LHS_prevalence<-ggplot(LHS_C, aes(x=Year, y=value, by=variable)) + 
  geom_line(aes(x=Year, y=value, by=variable, colour='All Ages'), show.legend=TRUE)+
  geom_line(aes(x=LHS3_C$Year, y=LHS3_C$value, by=LHS3_C$variable, colour='25-30 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS2_C$Year, y=LHS2_C$value, by=LHS2_C$variable, colour='20-24 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS1_C$Year, y=LHS1_C$value, by=LHS1_C$variable, colour='15-19 Years'), show.legend=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(LHS_C$value+10)), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "right", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name=NULL, values = c("#CCCCCC","#999999", "#666666", "#000000"), breaks = c("All Ages","25-30 Years", "20-24 Years", "15-19 Years")) +
  ylab(text.ylab) +
  xlab(text.xlab) 

w <- 6; h <- 4
ggsave(sprintf("LHS_prevalence.pdf"), LHS_prevalence, width=w, height=h)

#Restrict to the Best-Fitting 50%
LHS_prevalence<-ggplot() + 
  geom_line(data = subset(LHS_C,LHS_C$best50==1), aes(x=Year, y=value, by=variable, colour='All Ages'), show.legend=TRUE)+
  geom_line(data = subset(LHS3_C,LHS3_C$best50==1),aes(x=Year, y=value, by=variable, colour='25-30 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS2_C,LHS2_C$best50==1),aes(x=Year, y=value, by=variable, colour='20-24 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS1_C,LHS1_C$best50==1),aes(x=Year, y=value, by=variable, colour='15-19 Years'), show.legend=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(LHS_C$value+10)), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "right", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name=NULL, values = c("#CCCCCC","#999999", "#666666", "#000000"), breaks = c("All Ages","25-30 Years", "20-24 Years", "15-19 Years")) +
  ylab(text.ylab) +
  xlab(text.xlab) 

w <- 6; h <- 4
ggsave(sprintf("LHS_prevalence_Best50.pdf"), LHS_prevalence, width=w, height=h)

#Create Combined Plot of Best 50% and All
p1<-ggplot(LHS_C, aes(x=Year, y=value, by=variable)) + 
  geom_line(aes(x=Year, y=value, by=variable, colour='All Ages'), show.legend=TRUE)+
  geom_line(aes(x=LHS3_C$Year, y=LHS3_C$value, by=LHS3_C$variable, colour='25-30 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS2_C$Year, y=LHS2_C$value, by=LHS2_C$variable, colour='20-24 Years'), show.legend=TRUE)+
  geom_line(aes(x=LHS1_C$Year, y=LHS1_C$value, by=LHS1_C$variable, colour='15-19 Years'), show.legend=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(LHS_C$value+10)), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name=NULL, values = c("#CCCCCC","#999999", "#666666", "#000000"), breaks = c("All Ages","25-30 Years", "20-24 Years", "15-19 Years")) +
  ylab("") +
  xlab(text.xlab) 

title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p1 <- arrangeGrob(p1, top = title.grob)
grid.arrange(p1)


#Restrict to the Best-Fitting 50%
p2<-ggplot() + 
  geom_line(data = subset(LHS_C,LHS_C$best50==1), aes(x=Year, y=value, by=variable, colour='All Ages'), show.legend=TRUE)+
  geom_line(data = subset(LHS3_C,LHS3_C$best50==1),aes(x=Year, y=value, by=variable, colour='25-30 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS2_C,LHS2_C$best50==1),aes(x=Year, y=value, by=variable, colour='20-24 Years'), show.legend=TRUE)+
  geom_line(data = subset(LHS1_C,LHS1_C$best50==1),aes(x=Year, y=value, by=variable, colour='15-19 Years'), show.legend=TRUE)+
  scale_x_continuous(limits = c(2000, 2013), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(LHS_C$value+10)), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name=NULL, values = c("#CCCCCC","#999999", "#666666", "#000000"), breaks = c("All Ages","25-30 Years", "20-24 Years", "15-19 Years")) +
  ylab(text.ylab) +
  xlab(text.xlab)

title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p2 <- arrangeGrob(p2, top = title.grob)
grid.arrange(p2)

w <- 10; h <- 4
mylegend<-g_legend(LHS_prevalence)

lay<-rbind(c(1,1,1,2,2,2,3))
options(device = "quartz")
LHS_prev<-grid.arrange(p2,p1,mylegend,layout_matrix=lay)
ggsave(sprintf("LHS_prev.pdf"), LHS_prev, width=w, height=h)

##### LHS Parameter Histograms (Black/Grey) #####

#Set ylabel text
text.ylab <- "# Parameter Sets"

#Beta

#textsym <- expression(beta[1])
summary(params$b)
text.xlab <- bquote(Days*~Years^-1*~Persons^-1)
text.title <- bquote('Transmission Rate ('*beta*')')

#use to extract legend (bottom)
p1<-ggplot(params) + 
  geom_histogram(aes(x=b, fill = "All Parameter Sets "), colour='#000000', breaks=seq(1.8e-07,2.3e-05,4e-07), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=b, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(1.8e-07,2.3e-05,4e-07), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.8e-07,2.3e-05), breaks = c(0.000001,0.00001,0.00002))+
  scale_y_continuous(limits = c(0, 300), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#use to extract legend (side)
p2<-ggplot(params) + 
  geom_histogram(aes(x=b, fill = "All Parameter Sets "), colour='#000000', breaks=seq(1.8e-07,2.3e-05,4e-07), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=b, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(1.8e-07,2.3e-05,4e-07), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "right", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.8e-07,2.3e-05), breaks = c(0.000001,0.00001,0.00002))+
  scale_y_continuous(limits = c(0, 300), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#beta: no legend
text.xlab <- bquote(Days*~Years^-1*~Persons^-1)
text.title <- bquote('Transmission Rate: ('*beta*')')

p3<-ggplot(params) + 
  geom_histogram(aes(x=b, fill = "All Parameter Sets "), colour='#000000', breaks=seq(1.8e-07,2.3e-05,4e-07), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=b, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(1.8e-07,2.3e-05,4e-07), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.8e-07,2.3e-05), breaks = c(0.000001,0.00001,0.00002))+
  scale_y_continuous(limits = c(0, 300), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#phi1
summary(params$phi1)
summary(log(params$phi1))

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: Ages 15-19 ('*phi[1]*')')
p4<-ggplot(params) + 
  geom_histogram(aes(x=phi1, fill = "All Parameter Sets "), colour='#000000', show.legend=TRUE, bins=70)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=phi1, fill = "Best 50% of Parameter Sets "), colour='#000000', na.rm=TRUE, show.legend=TRUE, bins=70)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_log10()+
  scale_y_continuous(limits=c(0,1800), breaks=c(0,500,1000,1500), labels=c("0","500","1,000","1,500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#phi2
summary(params$phi2)
text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: Ages 20-24 ('*phi[2]*')')
p5<-ggplot(params) + 
  geom_histogram(aes(x=phi2, fill = "All Parameter Sets "), colour='#000000', show.legend=TRUE,bins=70)+ 
  geom_histogram(data=subset(params,params$best50==1),aes(x=phi2, fill = "Best 50% of Parameter Sets "), colour='#000000', na.rm=TRUE, show.legend=TRUE,bins=70)+ 
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_log10(breaks=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))+
  scale_y_continuous(limits=c(0,1800), breaks=c(0,500,1000,1500), labels=c("0","500","1,000","1,500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#phi3
summary(params$phi3)
text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: Ages 25-30 ('*phi[3]*')')
p6<-ggplot(params) + 
  geom_histogram(aes(x=phi3, fill = "All Parameter Sets "), colour='#000000', show.legend=TRUE,bins=70)+ 
  geom_histogram(data=subset(params,params$best50==1),aes(x=phi3, fill = "Best 50% of Parameter Sets "), colour='#000000', na.rm=TRUE, show.legend=TRUE,bins=70)+ 
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_log10()+
  scale_y_continuous(limits=c(0,1800), breaks=c(0,500,1000,1500), labels=c("0","500","1,000","1,500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#delta
summary(params$d)
text.xlab <- bquote(Years^-1)
text.title <- bquote('Cessation Rate ('*delta*')')
p7<-ggplot(params) + 
  geom_histogram(aes(x=d, fill = "All Parameter Sets "), colour='#000000', breaks=seq(0.085,1.125,0.05), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=d, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(0.085,1.125,0.05), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.085,1.125), breaks = c(0.1, 0.5, 1.0))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#k 
summary(params$k) 
text.xlab <- bquote(Years^-1)
text.title <- bquote('Relapse Rate ('*kappa*')')
p8<-ggplot(params) + 
  geom_histogram(aes(x=k, fill = "All Parameter Sets "), colour='#000000', breaks=seq(0.1,1,0.05), show.legend=TRUE)+ 
  geom_histogram(data=subset(params,params$best50==1),aes(x=k, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(0.1,1,0.05), na.rm=TRUE, show.legend=TRUE)+ 
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.1,1), breaks = c(0.25, 0.5, 0.75,1), labels=c("0.25", "0.5", "0.75","1.0"))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#r
summary(params$r)
text.xlab <- bquote('Surveillance Cases per Infection')
text.title <- bquote('Reporting Rate ('*rho*')')
p9<-ggplot(params) + 
  geom_histogram(aes(x=r, fill = "All Parameter Sets "), colour='#000000', breaks=seq(0.059,0.082,0.0015), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=r, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(0.059,0.082,0.0015), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.059,0.082), breaks = c(0.06, 0.07, 0.08))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#eps
summary(params$eps)
text.title <- bquote('Current PWID Prevalence ('*epsilon*')')
text.xlab <- "%"
p10<-ggplot(params) + 
  geom_histogram(aes(x=eps*100, fill = "All Parameter Sets "), colour='#000000', breaks=seq(0.9,1.3,0.02), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=eps*100, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(0.9,1.3,0.02), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.9,1.3), breaks = c(1.0, 1.1, 1.2))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))



#alpha
summary(params$a)
text.xlab <- "%"
text.title <- bquote('Acute HCV Developing Chronic ('*alpha*')')
p11<-ggplot(params) + 
  geom_histogram(aes(x=a*100, fill = "All Parameter Sets "), colour='#000000', breaks=seq(75,85,0.5), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=a*100, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(75,85,0.5), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(75,85), breaks = c(75, 80, 85))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#lambda1
summary(params$lambda1)
text.xlab <- "%"
text.title <- bquote('PWID Prevalence: Ages 15-24 ('*lambda[1]*', '*lambda[2]*')')
p12<-ggplot(params) + 
  geom_histogram(aes(x=lambda1*100, fill = "All Parameter Sets "), colour='#000000', breaks=seq(0.9,1.6,0.04), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=lambda1*100, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(0.9,1.6,0.04), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.9,1.6), breaks = c(1, 1.25, 1.5))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#lambda2
summary(params$lambda2)
text.xlab <- "%"
text.title <- bquote('PWID Prevalence: Ages 25-30 ('*lambda[3]*')')
p13<-ggplot(params) + 
  geom_histogram(aes(x=lambda2*100, fill = "All Parameter Sets "), colour='#000000', breaks=seq(1.4,3.4,0.12), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=lambda2*100, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(1.4,3.4,0.12), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.4,3.4), breaks = c(1.5, 2, 2.5, 3))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Z0 (popsize of newly added 15 yo/yr)
summary(params$Z0)
text.xlab <- "Persons"
text.title <- bquote('Number of 15 Year Olds ('*Z[0]*')')
p14<-ggplot(params) + 
  geom_histogram(aes(x=Z0, fill = "All Parameter Sets "), colour='#000000', breaks=seq(133000,154000,1000), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=Z0, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(133000,154000,1000), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(133000,154000), breaks = c(133000,143000,153000), labels=c("133,000","143,000","153,000"))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c11
summary(params$c1)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 15-19 ('*theta["1,1"]*')')
p15<-ggplot(params) + 
  geom_histogram(aes(x=c1, fill = "All Parameter Sets "), colour='#000000', breaks=seq(4.1,18.1,1), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=c1, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(4.1,18.1,1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(4.1,17.6), breaks = c(5,10,15))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c12
summary(params$c2)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 20-24 ('*theta["1,2"]*', '*theta["2,1"]*')')
p16<-ggplot(params) + 
  geom_histogram(aes(x=c2, fill = "All Parameter Sets "), colour='#000000', breaks=seq(1.0,2.6,0.1), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=c2, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(1.0,2.6,0.1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.0,2.6), breaks = c(1.5,2.5))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c13
summary(params$c3)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 25-30 ('*theta["1,3"]*', '*theta["3,1"]*')')
p17<-ggplot(params) + 
  geom_histogram(aes(x=c3, fill = "All Parameter Sets "), colour='#000000', breaks=seq(0.3,1.2,0.05), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=c3, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(0.3,1.2,0.05), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.3,1.2), breaks = c(0.5,1))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c22
summary(params$c5)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-24 with 20-24 ('*theta["2,2"]*')')
p18<-ggplot(params) + 
  geom_histogram(aes(x=c5, fill = "All Parameter Sets "), colour='#000000', breaks=seq(2.4,6.4,0.2), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=c5, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(2.4,6.4,0.2), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2.4,6.4), breaks = c(3,4,5,6))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c23
summary(params$c6)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-24 with 25-30 ('*theta["2,3"]*', '*theta["3,2"]*')')
p19<-ggplot(params) + 
  geom_histogram(aes(x=c6, fill = "All Parameter Sets "), colour='#000000', breaks=seq(1.0,3.3,0.1), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=c6, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(1.0,3.3,0.1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.0,3.3), breaks = c(1.5,2.5))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c33
summary(params$c9)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 25-30 with 25-30 ('*theta["3,3"]*')')
p20<-ggplot(params) + 
  geom_histogram(aes(x=c9, fill = "All Parameter Sets "), colour='#000000', breaks=seq(1.7,3.5,0.1), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=c9, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(1.7,3.5,0.1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.7,3.5), breaks = c(2,3))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Crude Death Rate 1 (15-19)
summary(params$mu1)
text.xlab <- "Deaths per 100,000 Persons"
text.title <- bquote('Mortality Rate: Ages 15-19 ('*mu[1]*')')
p21<-ggplot(params) + 
  geom_histogram(aes(x=mu1*100000, fill = "All Parameter Sets "), colour='#000000', breaks=seq(46.3,64.3,1), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=mu1*100000, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(46.3,64.3,1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(46.3,64.3), breaks = c(50,60))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Crude Death Rate 2 (20-24)
summary(params$mu2)
text.xlab <- "Deaths per 100,000 Persons"
text.title <- bquote('Mortality Rate: Ages 20-24 ('*mu[2]*')')
p22<-ggplot(params) + 
  geom_histogram(aes(x=mu2*100000, fill = "All Parameter Sets "), colour='#000000', breaks=seq(83.5,99.2,1), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=mu2*100000, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(83.5,99.2,1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(83.5,99.2), breaks = c(85,95))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Crude Death Rate 3 (25-30)
summary(params$mu3)
text.xlab <- "Deaths per 100,000 Persons"
text.title <- bquote('Mortality Rate: Ages 25-30 ('*mu[3]*')')
p23<-ggplot(params) + 
  geom_histogram(aes(x=mu3*100000, fill = "All Parameter Sets "), colour='#000000', breaks=seq(104.9,125.9,1), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=mu3*100000, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(104.9,125.9,1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(104.9,125.9), breaks = c(85,95))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#SMR: Current PWID vs Non-PWID
summary(params$etap)
text.xlab <- "Standardized Mortality Ratio (SMR)"
text.title <- bquote('Current vs. Non-PWID Mortality ('*eta[P]*')')
p24<-ggplot(params) + 
  geom_histogram(aes(x=etap, fill = "All Parameter Sets "), colour='#000000', breaks=seq(2.5,15.5,0.5), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=etap, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(2.5,15.5,0.5), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2.5,15.5), breaks = c(5,10,15))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#SMR: Former PWID vs Current
summary(params$etan)
text.xlab <- "Standardized Mortality Ratio (SMR)"
text.title <- bquote('Former vs. Current PWID Mortality ('*eta[N]*')')
p25<-ggplot(params) + 
  geom_histogram(aes(x=etan, fill = "All Parameter Sets "), colour='#000000', breaks=seq(0.18,0.54,0.02), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=etan, fill = "Best 50% of Parameter Sets "), colour='#000000', breaks=seq(0.18,0.54,0.02), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.18,0.54), breaks = c(0.2,0.3,0.4,0.5))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#RSS 
summary(params$RSS)
text.xlab <- "Residual"
p26<-ggplot(params) + 
  geom_histogram(aes(x=RSS, fill = "All Parameter Sets "), colour='#000000', breaks=seq(3,26,0.2), show.legend=TRUE)+
  geom_histogram(data=subset(params,params$best50==1),aes(x=RSS, fill = "Best 50% of Parameter Sets "), colour='#000000',  breaks=seq(3,26,0.2), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(3,26), breaks = c(6,12,18,24))+
  ggtitle('Residual Sum of Squares')+
  theme(plot.title = element_text(hjust = 0.5))

#Saving main and supplemental figs

#Note: quartz works for saving unicode text as pdf (ggsave does not)

#Main Figure: Non-Uniformly Distributed and Contact Parameters 
w <- 20; h <- 12
mylegend<-g_legend(p1)

lay<-rbind(c(1,2,3), c(1,2,3), c(1,2,3), c(4,5,6), c(4,5,6), c(4,5,6), 
           c(7,8,9), c(7,8,9), c(7,8,9), c(10,11,12), c(10,11,12), c(10,11,12),
           c(13,13,13))
options(device = "quartz")
LHS_Histograms<-grid.arrange(p10,p12,p13,p9,p7,p8,p15,p16,p17,p18,p19,p20,mylegend,layout_matrix=lay)
quartz.save ("LHS_Histograms.pdf", "pdf",width=w,height=h)

#Non-Uniformly Distributed Params

w <- 20; h <- 9
mylegend<-g_legend(p1)

lay<-rbind(c(1,2,3), c(1,2,3), c(1,2,3), c(4,5,6), c(4,5,6), c(4,5,6), c(7,7,7))
options(device = "quartz")
LHS_Histograms_main<-grid.arrange(p10,p12,p13,p9,p7,p8,mylegend,layout_matrix=lay)
quartz.save ("LHS_Histograms_main.pdf", "pdf",width=w,height=h)

#Contact Params
w <- 20; h <- 9
lay<-rbind(c(1,2,3), c(1,2,3), c(1,2,3), c(4,5,6), c(4,5,6), c(4,5,6), c(7,7,7))
options(device = "quartz")
LHS_Histograms_supp1<-grid.arrange(p15,p16,p17,p18,p19,p20,mylegend,layout_matrix=lay)
quartz.save ("LHS_Histograms_contact.pdf", "pdf",width=w,height=h)

#Supplemental Figure: Other Parameters and RSS (RSS,other params,contact params)
w <- 20; h <- 12
mylegend<-g_legend(p1)

lay<-rbind(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4), c(5,6,7,8), c(5,6,7,8),c(5,6,7,8),c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),c(13,13,13,13))
options(device = "quartz")
LHS_Histograms_supp<-grid.arrange(p3,p4,p5,p6,p11,p14,p24,p25,p21,p22,p23,p26,mylegend,layout_matrix=lay)
quartz.save ("LHS_Histograms_supp.pdf", "pdf",width=w,height=h)

##### LHS Parameter Histograms (Greyscale by Quartiles of Fit) #####
  
#Set ylabel text
text.ylab <- "# Parameter Sets"
model.colors <- c('#000000', '#666666', '#999999','#CCCCCC') #grey scheme

#Beta

#textsym <- expression(beta[1])
summary(params$b)
text.xlab <- bquote(Days*~Years^-1*~Persons^-1)
text.title <- bquote('Transmission Rate ('*beta*')')

#use to extract legend (bottom) 
p1<-ggplot(params) + 
  geom_histogram(aes(x=b, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",breaks=seq(1.8e-07,2.3e-05,4e-07), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.8e-07,2.3e-05), breaks = c(0.000001,0.00001,0.00002))+
  scale_y_continuous(limits = c(0, 300), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#use to extract legend (side)
p2<-ggplot(params) + 
  geom_histogram(aes(x=b, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.8e-07,2.3e-05,4e-07), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "right", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 50% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.8e-07,2.3e-05), breaks = c(0.000001,0.00001,0.00002))+
  scale_y_continuous(limits = c(0, 300), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#beta: no legend
text.xlab <- bquote(Days*~Years^-1*~Persons^-1)
text.title <- bquote('Transmission Rate: ('*beta*')')

p3<-ggplot(params) + 
  geom_histogram(aes(x=b, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.8e-07,2.3e-05,4e-07), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.8e-07,2.3e-05), breaks = c(0.000001,0.00001,0.00002))+
  scale_y_continuous(limits = c(0, 300), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#phi1
summary(params$phi1)
summary(log(params$phi1))

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: Ages 15-19 ('*phi[1]*')')
p4<-ggplot(params) + 
  geom_histogram(aes(x=phi1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', show.legend=TRUE, bins=70)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_log10()+
  scale_y_continuous(limits=c(0,1800), breaks=c(0,500,1000,1500), labels=c("0","500","1,000","1,500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#phi2
summary(params$phi2)
text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: Ages 20-24 ('*phi[2]*')')
p5<-ggplot(params) + 
  geom_histogram(aes(x=phi2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', show.legend=TRUE,bins=70)+ 
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_log10(breaks=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))+
  scale_y_continuous(limits=c(0,805), breaks=c(0,200,400,600,800), labels=c("0","200","400","600","800"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#phi3
summary(params$phi3)
text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: Ages 25-30 ('*phi[3]*')')
p6<-ggplot(params) + 
  geom_histogram(aes(x=phi3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', show.legend=TRUE,bins=70)+ 
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_log10()+
  scale_y_continuous(limits=c(0,1500), breaks=c(0,500,1000,1500), labels=c("0","500","1,000","1,500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#delta
summary(params$d)
text.xlab <- bquote(Years^-1)
text.title <- bquote('Cessation Rate ('*delta*')')
p7<-ggplot(params) + 
  geom_histogram(aes(x=d, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.085,1.125,0.05), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.085,1.125), breaks = c(0.1, 0.5, 1.0))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#k 
summary(params$k) 
text.xlab <- bquote(Years^-1)
text.title <- bquote('Relapse Rate ('*kappa*')')
p8<-ggplot(params) + 
  geom_histogram(aes(x=k, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.1,1,0.05), show.legend=TRUE)+ 
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.1,1), breaks = c(0.25, 0.5, 0.75,1), labels=c("0.25", "0.5", "0.75","1.0"))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#r
summary(params$r)
text.xlab <- bquote('Surveillance Cases per Infection')
text.title <- bquote('Reporting Rate ('*rho*')')
p9<-ggplot(params) + 
  geom_histogram(aes(x=r, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.059,0.082,0.0015), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.059,0.082), breaks = c(0.06, 0.07, 0.08))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#eps
summary(params$eps)
text.title <- bquote('Current PWID Prevalence ('*epsilon*')')
text.xlab <- "%"
p10<-ggplot(params) + 
  geom_histogram(aes(x=eps*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.9,1.3,0.02), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.9,1.3), breaks = c(1.0, 1.1, 1.2))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))



#alpha
summary(params$a)
text.xlab <- "%"
text.title <- bquote('Acute HCV Developing Chronic ('*alpha*')')
p11<-ggplot(params) + 
  geom_histogram(aes(x=a*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(75,85,0.5), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(75,85), breaks = c(75, 80, 85))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#lambda1
summary(params$lambda1)
text.xlab <- "%"
text.title <- bquote('PWID Prevalence: Ages 15-24 ('*lambda[1]*', '*lambda[2]*')')
p12<-ggplot(params) + 
  geom_histogram(aes(x=lambda1*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.9,1.6,0.04), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.9,1.6), breaks = c(1, 1.25, 1.5))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#lambda2
summary(params$lambda2)
text.xlab <- "%"
text.title <- bquote('PWID Prevalence: Ages 25-30 ('*lambda[3]*')')
p13<-ggplot(params) + 
  geom_histogram(aes(x=lambda2*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.4,3.4,0.12), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.4,3.4), breaks = c(1.5, 2, 2.5, 3))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Z0 (popsize of newly added 15 yo/yr)
summary(params$Z0)
text.xlab <- "Persons"
text.title <- bquote('Number of 15 Year Olds ('*Z[0]*')')
p14<-ggplot(params) + 
  geom_histogram(aes(x=Z0, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(133000,154000,1000), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(133000,154000), breaks = c(133000,143000,153000), labels=c("133,000","143,000","153,000"))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c11
summary(params$c1)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 15-19 ('*theta["1,1"]*')')
p15<-ggplot(params) + 
  geom_histogram(aes(x=c1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(4.1,18.1,1), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(4.1,17.6), breaks = c(5,10,15))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c12
summary(params$c2)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 20-24 ('*theta["1,2"]*', '*theta["2,1"]*')')
p16<-ggplot(params) + 
  geom_histogram(aes(x=c2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.0,2.6,0.1), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.0,2.6), breaks = c(1.5,2.5))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c13
summary(params$c3)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 25-30 ('*theta["1,3"]*', '*theta["3,1"]*')')
p17<-ggplot(params) + 
  geom_histogram(aes(x=c3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.3,1.2,0.05), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.3,1.2), breaks = c(0.5,1))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c22
summary(params$c5)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-24 with 20-24 ('*theta["2,2"]*')')
p18<-ggplot(params) + 
  geom_histogram(aes(x=c5, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(2.4,6.4,0.2), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2.4,6.4), breaks = c(3,4,5,6))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c23
summary(params$c6)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-24 with 25-30 ('*theta["2,3"]*', '*theta["3,2"]*')')
p19<-ggplot(params) + 
  geom_histogram(aes(x=c6, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.0,3.3,0.1), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.0,3.3), breaks = c(1.5,2.5))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Contact Matrix: c33
summary(params$c9)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 25-30 with 25-30 ('*theta["3,3"]*')')
p20<-ggplot(params) + 
  geom_histogram(aes(x=c9, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.7,3.5,0.1), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.7,3.5), breaks = c(2,3))+
  scale_y_continuous(limits = c(0, 400), breaks = c(0,100,200,300,400))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Crude Death Rate 1 (15-19)
summary(params$mu1)
text.xlab <- "Deaths per 100,000 Persons"
text.title <- bquote('Mortality Rate: Ages 15-19 ('*mu[1]*')')
p21<-ggplot(params) + 
  geom_histogram(aes(x=mu1*100000, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(46.3,64.3,1), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(46.3,64.3), breaks = c(50,60))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Crude Death Rate 2 (20-24)
summary(params$mu2)
text.xlab <- "Deaths per 100,000 Persons"
text.title <- bquote('Mortality Rate: Ages 20-24 ('*mu[2]*')')
p22<-ggplot(params) + 
  geom_histogram(aes(x=mu2*100000, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(83.5,99.2,1), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(83.5,99.2), breaks = c(85,95))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#Crude Death Rate 3 (25-30)
summary(params$mu3)
text.xlab <- "Deaths per 100,000 Persons"
text.title <- bquote('Mortality Rate: Ages 25-30 ('*mu[3]*')')
p23<-ggplot(params) + 
  geom_histogram(aes(x=mu3*100000, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(104.9,125.9,1), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(104.9,125.9), breaks = c(85,95))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#SMR: Current PWID vs Non-PWID
summary(params$etap)
text.xlab <- "Standardized Mortality Ratio (SMR)"
text.title <- bquote('Current vs. Non-PWID Mortality ('*eta[P]*')')
p24<-ggplot(params) + 
  geom_histogram(aes(x=etap, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(2.5,15.5,0.5), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2.5,15.5), breaks = c(5,10,15))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))


#SMR: Former PWID vs Current
summary(params$etan)
text.xlab <- "Standardized Mortality Ratio (SMR)"
text.title <- bquote('Former vs. Current PWID Mortality ('*eta[N]*')')
p25<-ggplot(params) + 
  geom_histogram(aes(x=etan, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.18,0.54,0.02), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.18,0.54), breaks = c(0.2,0.3,0.4,0.5))+
  scale_y_continuous(limits = c(0, 350), breaks = c(0,100,200,300))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))

#RSS 
summary(params$RSS)
text.xlab <- "Residual"
p26<-ggplot(params) + 
  geom_histogram(aes(x=RSS, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(3,26,0.2), show.legend=TRUE)+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(3,26), breaks = c(6,12,18,24))+
  ggtitle('Residual Sum of Squares')+
  theme(plot.title = element_text(hjust = 0.5))

#Saving main and supplemental figs

#Note: quartz works for saving unicode text as pdf (ggsave does not)

#Main Figure: Non-Uniformly Distributed and Contact Parameters 
w <- 20; h <- 12
mylegend<-g_legend(p1)

lay<-rbind(c(1,2,3), c(1,2,3), c(1,2,3), c(4,5,6), c(4,5,6), c(4,5,6), 
           c(7,8,9), c(7,8,9), c(7,8,9), c(10,11,12), c(10,11,12), c(10,11,12),
           c(13,13,13))
options(device = "quartz")
LHS_Histograms<-grid.arrange(p10,p12,p13,p9,p7,p8,p15,p16,p17,p18,p19,p20,mylegend,layout_matrix=lay)
quartz.save ("LHS_Histograms_Grey.pdf", "pdf",width=w,height=h)

#Non-Uniformly Distributed Params

w <- 20; h <- 9
mylegend<-g_legend(p1)

lay<-rbind(c(1,2,3), c(1,2,3), c(1,2,3), c(4,5,6), c(4,5,6), c(4,5,6), c(7,7,7))
options(device = "quartz")
LHS_Histograms_main<-grid.arrange(p10,p12,p13,p9,p7,p8,mylegend,layout_matrix=lay)
quartz.save ("LHS_Histograms_main_Grey.pdf", "pdf",width=w,height=h)

#Contact Params
w <- 20; h <- 9
lay<-rbind(c(1,2,3), c(1,2,3), c(1,2,3), c(4,5,6), c(4,5,6), c(4,5,6), c(7,7,7))
options(device = "quartz")
LHS_Histograms_supp1<-grid.arrange(p15,p16,p17,p18,p19,p20,mylegend,layout_matrix=lay)
quartz.save ("LHS_Histograms_contact_Grey.pdf", "pdf",width=w,height=h)

#Supplemental Figure: Other Parameters and RSS (RSS,other params,contact params)
w <- 20; h <- 12
mylegend<-g_legend(p1)

lay<-rbind(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4), c(5,6,7,8), c(5,6,7,8),c(5,6,7,8),c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),c(13,13,13,13))
options(device = "quartz")
LHS_Histograms_supp<-grid.arrange(p3,p4,p5,p6,p11,p14,p24,p25,p21,p22,p23,p26,mylegend,layout_matrix=lay)
quartz.save ("LHS_Histograms_supp_Grey.pdf", "pdf",width=w,height=h)


#### Statistical testing of parameter distributions/means #### 

## test if distribution of 5,000 follows expected uniform distribution (as expected/designed)##
# Note: bounds = sampling range #
ks.test(params$d,"punif",0.085,1.123)
ks.test(params$k,"punif",0.1,1)
ks.test(params$r,"punif",(1/16.8),(1/12.3))
ks.test(params$eps,"punif",0.009,0.013)
ks.test(params$a,"punif",0.75,0.85)
ks.test(params$Z0,"punif",133231,153080)
ks.test(params$lambda1,"punif",0.009,0.016)
ks.test(params$lambda2,"punif",0.014,0.034)
ks.test(params$mu1,"punif",0.000463,0.000643)
ks.test(params$mu2,"punif",0.000835,0.000992)
ks.test(params$mu3,"punif",0.001049,0.001259)
ks.test(params$etap,"punif",2.5,15.3)
ks.test(params$etan,"punif",0.18,0.54)
ks.test(params$c1,"punif",4.19,17.53)
ks.test(params$c2,"punif",1.04,2.51)
ks.test(params$c3,"punif",0.32,1.19)
ks.test(params$c5,"punif",2.42,6.26)
ks.test(params$c6,"punif",1.08,3.23)
ks.test(params$c9,"punif",1.79,3.45)


#check equality of ranks (~means) for 50% best vs 50% worst
wilcox.test(params$b~params$best50)
wilcox.test(params$phi1~params$best50)
wilcox.test(params$phi2~params$best50)
wilcox.test(params$phi3~params$best50)
wilcox.test(params$d~params$best50)
wilcox.test(params$k~params$best50)
wilcox.test(params$r~params$best50)
wilcox.test(params$eps~params$best50)
wilcox.test(params$a~params$best50)
wilcox.test(params$Z0~params$best50)
wilcox.test(params$lambda1~params$best50)
wilcox.test(params$lambda2~params$best50)
wilcox.test(params$c1~params$best50)
wilcox.test(params$c2~params$best50)
wilcox.test(params$c3~params$best50)
wilcox.test(params$c5~params$best50)
wilcox.test(params$c6~params$best50)
wilcox.test(params$c9~params$best50)
wilcox.test(params$mu1~params$best50)
wilcox.test(params$mu2~params$best50)
wilcox.test(params$mu3~params$best50)
wilcox.test(params$etap~params$best50)
wilcox.test(params$etan~params$best50)

wilcox.test(params$RSS~params$best50)


## Do worst and best 50% come from same distribution? ##

params$var<-params$b
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$phi1
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$phi2
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$phi3
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$d
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$k
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$r
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$eps
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$a
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$Z0
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$lambda1
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$lambda1
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$mu1
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$mu2
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$mu3
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$etap
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$etan
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$c1
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$c2
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$c3
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$c5
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$c6
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$c9
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)

params$var<-params$RSS
test$worst<-ifelse(params$best50==0,params$var,NA)
test$best<-ifelse(params$best50==1,params$var,NA)
ks.test(test$worst,test$best, na.action=na.omit)


#### 3d Plot of phi1-3 ####

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
  title='phi1'
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
  title='phi2'
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
  title='phi3'
)

#create 3d plot
p <- plot_ly(params, x = ~phi1, y = ~phi2, z = ~phi3, color = ~RSS, opacity=1,
             marker=list(size=4),
             colors=c("red","orange","yellow","green","blue")) %>%
  add_markers() %>%
  layout(scene=list(xaxis=ax1,yaxis=ax2,zaxis=ax3))
p

p <- plot_ly(subset(params,params$best50==1), x = ~phi1, y = ~phi2, z = ~phi3, color = ~RSS, opacity=1,
             marker=list(size=4),
             colors=c("red","orange","yellow","green","blue")) %>%
  add_markers() %>%
  layout(scene=list(xaxis=ax1,yaxis=ax2,zaxis=ax3))
p


#### Scatter plots of RSS vs each parameter ####
likeplots<-params[,c(2:8,12:21,23:25,27:28,31:36)]
w <- 5; h <- 4

RSSscatter<-function(data,y){
  nm <- names(data)
  for (i in seq_along(nm)) {
    x<-paste("likeplots$",nm[i],sep="")
    xlab<-paste(nm[i])
    p1<-ggplot(data=likeplots,aes_string(x=x,y=y))+
    geom_point(size=2)+
    ylab("RSS")+
    xlab(xlab)+
    #ggtitle(xlab)+
    theme_bw()+
    theme(panel.border = element_rect(linetype = "solid", colour = "black"))+
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
    theme(legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
    theme(legend.key = element_blank())
  ggsave(p1, filename=paste(nm[i],"_RSS.pdf"),width=w, height=h)
  print(p1)
  }
}

RSSscatter(likeplots,likeplots$RSS)




#### Describe best-fitting 25% of parameter sets ####
summary(params$RSS)


best25<-subset(params,params$RSS<=quantile(params$RSS, 0.25) )
summary(best25[,2:36])




################ Section 5 Figures: Intervention Simulation Results ################
#### Violin Plots of Best 50% Fits of Single Interventions ####
#Set plot labels and colors
text.xlab <- "Intervention"
model.colors <- c('#CCCCCC', '#999999', '#666666', '#000000','#FFFFFF') #grey scheme


#X labels
xlabelsall <-paste0(c("Former PWID Treatment, 1 Year","Former PWID Treatment, 12 Weeks",
                      "Current PWID Treatment, 1 Year","Current PWID Treatment, 12 Weeks",
                      "Decreased Relapse","Increased Cessation", 
                      "Decreased Effective Contacts","Decreased Injection Initiation",
                      "Decreased Mortality, Current PWID","Decreased Mortality, Former PWID","None"))

#plot for legend
text.ylab <- "Chronic HCV Prevalence"
leg<-ggplot(subset(singleint,singleint$order!=9999&singleint$best50==1), 
            aes(x=order, y=ChronicPrev_End, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("","","","","","","","","","",""))+ 
  scale_y_continuous(limits=c(0,4400),breaks=c(0,1000,2000,3000,4000),labels=c("0","1,000","2,000","3,000","4,000"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") 
leg

text.ylab <- "Chronic HCV Prevalence"
p1<-ggplot(subset(singleint,singleint$order!=9999&singleint$best50==1), 
           aes(x=order, y=ChronicPrev_End, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("","","","","","","","","","",""))+ 
  scale_y_continuous(limits=c(0,4400),breaks=c(0,1000,2000,3000,4000),labels=c("0","1,000","2,000","3,000","4,000"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") 

title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p1 <- arrangeGrob(p1, top = title.grob)
grid.arrange(p1)


text.ylab <- "New Chronic HCV Cases"
p2<-ggplot(subset(singleint,singleint$order!=9999&singleint$best50==1), 
           aes(x=order, y=NewChronic_End, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=xlabelsall)+ 
  scale_y_continuous(limits=c(0,1150),breaks=c(0,300,600,900),labels=c("0","300","600","900"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab(text.xlab) 

title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p2 <- arrangeGrob(p2, top = title.grob)
grid.arrange(p2)

w <- 36; h <- 8
mylegend<-g_legend(leg)
lay<-rbind(c(1), c(1), c(1), c(2),c(2),c(2),c(3))
SingleStrategy_VP_1<-grid.arrange(p1,p2,mylegend,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_VP_1_Best50.pdf"), SingleStrategy_VP_1, width=w, height=h)

#Smaller plot excluding mortality reduction and treatment duration

#X labels
xlabelsmain <-paste0(c("Former PWID Treatment",
                       "Current PWID Treatment",
                       "Decreased Relapse","Increased Cessation", 
                       "Decreased Effective Contacts","Decreased Injection Initiation",
                       "None"))

text.ylab <- "Chronic HCV Prevalence"
p1<-ggplot(subset(singleint,singleint$best50==1&singleint$order %in% c(2,4,5,6,7,8,11)), 
           aes(x=order, y=ChronicPrev_End, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("","","","","","",""))+ 
  scale_y_continuous(limits=c(0,4400),breaks=c(0,1000,2000,3000,4000),labels=c("0","1,000","2,000","3,000","4,000"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") 

title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p1 <- arrangeGrob(p1, top = title.grob)
grid.arrange(p1)

text.ylab <- "New Chronic HCV Cases"
p2<-ggplot(subset(singleint,singleint$best50==1&singleint$order %in% c(2,4,5,6,7,8,11)), 
           aes(x=order, y=NewChronic_End, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=xlabelsmain)+ 
  scale_y_continuous(limits=c(0,1150),breaks=c(0,300,600,900),labels=c("0","300","600","900"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab(text.xlab) 

title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p2 <- arrangeGrob(p2, top = title.grob)
grid.arrange(p2)

w <- 20; h <- 8
lay<-rbind(c(1), c(1), c(1), c(2),c(2),c(2),c(3))
SingleStrategy_VP_main<-grid.arrange(p1,p2,mylegend,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_VP_main_Best50.pdf"), SingleStrategy_VP_main, width=w, height=h)

#Supplemental plot with mortality reduction, decreased treatment duration
#X labels
xlabelssupp <-paste0(c("Former PWID Treatment, 1 Year",
                       "Current PWID Treatment, 1 Year",
                       "Decreased Mortality, Current PWID","Decreased Mortality, Former PWID",
                       "None"))

text.ylab <- "Chronic HCV Prevalence"
p1<-ggplot(subset(singleint,singleint$best50==1&singleint$order %in% c(1,3,9,10,11)), 
           aes(x=order, y=ChronicPrev_End, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("","","","","","",""))+ 
  scale_y_continuous(limits=c(0,4400),breaks=c(0,1000,2000,3000,4000),labels=c("0","1,000","2,000","3,000","4,000"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") 

title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p1 <- arrangeGrob(p1, top = title.grob)
grid.arrange(p1)

text.ylab <- "New Chronic HCV Cases"
p2<-ggplot(subset(singleint,singleint$best50==1&singleint$order %in% c(1,3,9,10,11)), 
           aes(x=order, y=NewChronic_End, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=xlabelssupp)+   
  scale_y_continuous(limits=c(0,1150),breaks=c(0,300,600,900),labels=c("0","300","600","900"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab(text.xlab) 

title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p2 <- arrangeGrob(p2, top = title.grob)
grid.arrange(p2)

w <- 18; h <- 8
lay<-rbind(c(1), c(1), c(1), c(2),c(2),c(2),c(3))
SingleStrategy_VP_supp<-grid.arrange(p1,p2,mylegend,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_VP_supp_Best50.pdf"), SingleStrategy_VP_supp, width=w, height=h)


#### Line Graphs of Best 50% Primary/Tertiary Sequential Interventions ####

#legend plots: p1=points, p2=lines
p1<-ggplot(data=subset(seqint2,h==1&inttype=="terttoprim"&pct!="None"))+
  geom_point(aes(x=as.factor(intono-1),y=meanpctred_chr,colour=as.factor(pct)), size=2)+
  geom_point(aes(x=as.factor(intono-1),y=meanpctred_new,colour=as.factor(pct)), size=2)+
  scale_x_discrete(breaks = 1:6, labels=c("Treatment: Former","+Treatment: Current","+Reduced Relapse",
                                          "+Increased Cessation", "+Reduced Contacts",
                                          "+Reduced Initiation"))+ 
  scale_y_continuous(limits=c(-7,105),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  ylab("")+
  xlab("Intervention Sequence")+
  ggtitle("Tertiary to Primary")+
  scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  # scale_linetype_manual(name="Type",values=c("dotdash","solid"),labels=c("Prevalence "," New Chronic Cases"))+
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="bottom", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))

p2<-ggplot(data=subset(seqint2,h==1&inttype=="terttoprim"&pct!="None"))+
  #geom_point(aes(x=as.factor(intono-1),y=meanpctred_chr,colour=as.factor(pct)), size=2)+
  #geom_point(aes(x=as.factor(intono-1),y=meanpctred_new,colour=as.factor(pct)), size=2)+
  geom_line(aes(x=intono-1,y=meanpctred_chr, linetype="dotdash"), size=1)+
  geom_line(aes(x=intono-1,y=meanpctred_new, linetype="solid"), size=1)+
  scale_x_discrete(breaks = 1:6, labels=c("Treatment: Former","+Treatment: Current","+Reduced Relapse",
                                          "+Increased Cessation", "+Reduced Contacts",
                                          "+Reduced Initiation"))+ 
  scale_y_continuous(limits=c(-7,105),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  ylab("")+
  xlab("Intervention Sequence")+
  ggtitle("Tertiary to Primary")+
  #scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  scale_linetype_manual(name="Type",values=c("dotdash","solid"),labels=c("Prevalence "," New Chronic Cases"))+
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="bottom", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))

p3<-ggplot(data=subset(seqint2,h==1&inttype=="terttoprim"&pct!="None"))+
  geom_point(aes(x=as.factor(intono-1),y=meanpctred_chr,colour=as.factor(pct)), size=2)+
  geom_point(aes(x=as.factor(intono-1),y=meanpctred_new,colour=as.factor(pct)), size=2)+
  geom_line(aes(x=intono-1,y=meanpctred_chr,colour=as.factor(pct), linetype="dotdash"), size=1)+
  geom_line(aes(x=intono-1,y=meanpctred_new,colour=as.factor(pct), linetype="solid"), size=1)+
  scale_x_discrete(breaks = 1:6, labels=c("Treatment: Former","+Treatment: Current","+Reduced Relapse",
                                          "+Increased Cessation", "+Reduced Contacts",
                                          "+Reduced Initiation"))+ 
  scale_y_continuous(limits=c(-10,105),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  ylab("")+
  xlab("Intervention Sequence")+
  ggtitle("Tertiary to Primary")+
  scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  scale_linetype_manual(name="Type",values=c("dotdash","solid"),labels=c("Prevalence "," New Chronic Cases"))+
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))

p4<-ggplot(data=subset(seqint2,h==1&inttype=="primtotert"&pct!="None"))+
  geom_point(aes(x=as.factor(intono-1),y=meanpctred_chr,colour=as.factor(pct)), size=2)+
  geom_point(aes(x=as.factor(intono-1),y=meanpctred_new,colour=as.factor(pct)), size=2)+
  geom_line(aes(x=intono-1,y=meanpctred_chr,colour=as.factor(pct), linetype="dotdash"), size=1)+
  geom_line(aes(x=intono-1,y=meanpctred_new,colour=as.factor(pct), linetype="solid"), size=1)+
  scale_x_discrete(breaks = 1:6, labels=c("Reduced Initiation","+Reduced Contacts","+Increased Cessation",
                                          "+Reduced Relapse","+Treatment: Current","Treatment: Former"))+ 
  scale_y_continuous(limits=c(-10,105),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  ylab("Mean % Case Reduction")+
  xlab("Intervention Sequence")+
  ggtitle("Primary to Tertiary")+
  scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  scale_linetype_manual(name="Type",values=c("dotdash","solid"),labels=c("Prevalence "," New Chronic Cases"))+
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))

w <- 8; h <- 6
legendpoint<-g_legend(p1)
legendline<-g_legend(p2)
lay<-rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(3,3),c(4,4))
Multi_Sequential<-grid.arrange(p4,p3,legendpoint,legendline,layout_matrix=lay)
ggsave(sprintf("Multi_Sequential_h=k=1_Best50.pdf"), Multi_Sequential, width=w, height=h)


#### Scatterplots of Best Fitting 50% of Primary/Tertiary Sequential Interventions ####
summary(subset(seqint,seqint$best50==1)$pctred_chr)
summary(subset(seqint,seqint$best50==1)$pctred_new)

#Intervention % legend
p1<-ggplot(data=subset(seqint,h==1&inttype=="terttoprim"&pct!="None"&best50==1))+
  geom_point(aes(x=as.factor(intono-1), y=pctred_chr,
                 colour=as.factor(pct)),
             position=position_dodge(width=c(0.4,0.4)))+
  scale_x_discrete(breaks = 1:6, labels=c("Treatment: Former","+Treatment: Current","+Reduced Relapse",
                                          "+Increased Cessation", "+Reduced Contacts",
                                          "+Reduced Initiation"))+ 
  scale_y_continuous(limits=c(-70,100),breaks=c(-50,-25,0,25,50,75,100),
                     labels=c("-50","-25","0","25","50","75","100"))+
  ylab("")+
  xlab("Intervention Sequence")+
  scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="bottom", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())

p3<-ggplot(data=subset(seqint,h==1&inttype=="terttoprim"&pct!="None"&best50==1))+
  geom_point(aes(x=as.factor(intono-1), y=pctred_chr,
                 colour=as.factor(pct)),
             position=position_dodge(width=c(0.4,0.4)))+
  scale_x_discrete(breaks = 1:6, labels=c("Treatment: Former","+Treatment: Current","+Reduced Relapse",
                                          "+Increased Cessation", "+Reduced Contacts",
                                          "+Reduced Initiation"))+ 
  stat_summary(aes(x=as.factor(intono-1),y=pctred_chr,group=as.factor(pct)),
               position=position_dodge(width=c(0.4,0.4)),fun.y = median, geom="point",
               shape=23, size=2, fill="white")+
  scale_y_continuous(limits=c(-70,100),breaks=c(-50,-25,0,25,50,75,100),
                     labels=c("-50","-25","0","25","50","75","100"))+
  ylab("")+
  xlab("Intervention Sequence")+
  scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())

title.grob <- textGrob(label = "D)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p3 <- arrangeGrob(p3, top = title.grob)
grid.arrange(p3)

p4<-ggplot(data=subset(seqint,h==1&inttype=="terttoprim"&pct!="None"&best50==1))+
  geom_point(aes(x=as.factor(intono-1), y=pctred_new, 
                 colour=as.factor(pct)),
             position=position_dodge(width=c(0.4,0.4)))+
  scale_x_discrete(breaks = 1:6, labels=c("","","",
                                          "", "",
                                          ""))+ 
  stat_summary(aes(x=as.factor(intono-1),y=pctred_new,group=as.factor(pct)),
               position=position_dodge(width=c(0.4,0.4)),fun.y = median, geom="point",
               shape=23, size=2, fill="white")+
  scale_y_continuous(limits=c(-70,100),breaks=c(-50,-25,0,25,50,75,100),
                     labels=c("-50","-25","0","25","50","75","100"))+
  ylab("")+
  xlab("")+
  scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())

title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p4 <- arrangeGrob(p4, top = title.grob)
grid.arrange(p4)

p5<-ggplot(data=subset(seqint,h==1&inttype=="primtotert"&pct!="None"&best50==1))+
  geom_point(aes(x=as.factor(intono-1), y=pctred_chr,
                 colour=as.factor(pct)),
             position=position_dodge(width=c(0.4,0.4)))+
  scale_x_discrete(breaks = 1:6, labels=c("+Reduced Initiation","+Reduced Contacts",
                                          "+Increased Cessation","+Reduced Relapse",
                                          "+Treatment: Current","Treatment: Former"))+ 
  stat_summary(aes(x=as.factor(intono-1),y=pctred_chr,group=as.factor(pct)),
               position=position_dodge(width=c(0.4,0.4)),fun.y = median, geom="point",
               shape=23, size=2, fill="white")+
  scale_y_continuous(limits=c(-70,100),breaks=c(-50,-25,0,25,50,75,100),
                     labels=c("-50","-25","0","25","50","75","100"))+
  ylab("% Reduction: Prevalence")+
  xlab("Intervention Sequence")+
  scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())

title.grob <- textGrob(label = "C)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p5 <- arrangeGrob(p5, top = title.grob)
grid.arrange(p5)

p6<-ggplot(data=subset(seqint,h==1&inttype=="primtotert"&pct!="None"&best50==1))+
  geom_point(aes(x=as.factor(intono-1), y=pctred_new, 
                 colour=as.factor(pct)),
             position=position_dodge(width=c(0.4,0.4)))+
  scale_x_discrete(breaks = 1:6, labels=c("","","",
                                          "", "",
                                          ""))+ 
  stat_summary(aes(x=as.factor(intono-1),y=pctred_new,group=as.factor(pct)),
               position=position_dodge(width=c(0.4,0.4)),fun.y = median, geom="point",
               shape=23, size=2, fill="white")+
  scale_y_continuous(limits=c(-70,100),breaks=c(-50,-25,0,25,50,75,100),
                     labels=c("-50","-25","0","25","50","75","100"))+
  ylab("% Reduction: New Chronic Cases")+
  xlab("")+
  scale_colour_manual(name="Intervention %", values = c("#CCCCCC","#999999", "#666666", "#000000"), labels= c("10%","20%", "30%", "40%")) +
  theme_bw(base_size=16)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', fill = 'NA', size = 2, linetype="blank"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.key = element_blank())

title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
                       hjust = 0, vjust = 0,gp = gpar(fontsize = 16))

p6 <- arrangeGrob(p6, top = title.grob)
grid.arrange(p6)

w <- 12; h <- 12
legendpoint<-g_legend(p1)
lay<-rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(5,5))
Multi_Sequential_ranges<-grid.arrange(p6,p4,p5,p3,legendpoint,layout_matrix=lay)
ggsave(sprintf("Multi_Sequential_h=k=1_ranges_Best50.pdf"), Multi_Sequential_ranges, width=w, height=h)

#### Describe sequential results w/high impact of tert to primary on new chronic####
examine<-subset(seqint,h==1&inttype=="terttoprim"&pct!="None"&pctred_new>=10&e==c(120,121,122,123))
write.csv(examine, "examine.csv")
