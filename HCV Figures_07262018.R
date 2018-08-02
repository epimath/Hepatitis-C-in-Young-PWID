## Figure Script: Gicquelais, Foxman, Coyle, and Eisenberg. (2018). Hepatitis C transmission in young 
## people who inject drugs: insights using a dynamic model informed by state public health surveillance. 
## Submitted. 

## Before running this code, run the model simulations in the code HCV_MS_07262018.m to simulate
## the hepatitis C model and output results as csv or txt files. This script creates figures 
## summarizing results and similar to those published in the article. 

################ Section 1: Packages and Functions to Load First ################### 
library(ggplot2)
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
setwd("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 5.6.2018")


############### Section 3: Importing Data Generated from Matlab ################
#### Import MDHHS Data ####

#MDHHS case counts by year
mdhhs<-read.csv("MDHHSData.csv", header=F)
colnames(mdhhs) <- c("Year",	"acute1529idu_1",	"acute1529idu_2",	"acute1529idu_3",	
                      "chronic1529idu_1",	"chronic1529idu_2",	"chronic1529idu_3")
summary(mdhhs)
mdhhs$chronic1529idu<-mdhhs$chronic1529idu_1+mdhhs$chronic1529idu_2+mdhhs$chronic1529idu_3
mdhhs$acute1529idu<-mdhhs$acute1529idu_1+mdhhs$acute1529idu_2+mdhhs$acute1529idu_3

#### Import New Chronic Cases (Current and Former PWID) ####

#All New Chronic Cases
#note that header row (var names) is the X + index of the parameter set (eg XParamSetNo.)
LHS_NewChronic<-read.csv("LHS_NewChronic.csv", header=TRUE)
summary(LHS_NewChronic)
LHS_NewChronic$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

keeps <- c("Year","chronic1529idu")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_NewChronic, id = "Year")
LHS.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(LHS.melted)

#New Chronic Cases - 15-19
LHS_NewChronic1<-read.csv("LHS_NewChronic1.csv", header=TRUE)
summary(LHS_NewChronic1)
LHS_NewChronic1$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

keeps <- c("Year","chronic1529idu_1")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_NewChronic1, id = "Year")
LHS1.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(LHS1.melted)


#New Chronic Cases - 20-25
LHS_NewChronic2<-read.csv("LHS_NewChronic2.csv", header=TRUE)
summary(LHS_NewChronic2)
LHS_NewChronic2$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

keeps <- c("Year","chronic1529idu_2")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_NewChronic2, id = "Year")
LHS2.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(LHS2.melted)

#New Chronic Cases - 26-29
LHS_NewChronic3<-read.csv("LHS_NewChronic3.csv", header=TRUE)
summary(LHS_NewChronic3)
LHS_NewChronic3$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

keeps <- c("Year","chronic1529idu_3")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_NewChronic3, id = "Year")
LHS3.melted<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(LHS3.melted)

#New Chronic Cases - 30-64
LHS_NewChronic4<-read.csv("LHS_NewChronic4.csv", header=TRUE)
summary(LHS_NewChronic4)
LHS_NewChronic4$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

LHS4.melted<-melt(LHS_NewChronic4, id = "Year")
summary(LHS4.melted)

#New Chronic Cases - 15-29
LHS<-merge(LHS3.melted,merge(LHS1.melted,LHS2.melted,by=c("variable","Year")),by=c("variable","Year"))
LHS$value1<-LHS$value+LHS$value.x+LHS$value.y
LHS$chronic1529idu<-LHS$chronic1529idu_1+LHS$chronic1529idu_2+LHS$chronic1529idu_3

keeps <- c("variable","Year","value1","chronic1529idu")
NewChronicLHS.melted_1529<-LHS[keeps]

#### Import Chronic HCV Prevalence Data (Current and Former PWID) ####

#Note that 1st row = parameter index (and becomes variable name)

#Total Chronic Cases - All Ages
LHS_ChronicPrev<-read.csv("LHS_ChronicPrev.csv", header=TRUE)
summary(LHS_ChronicPrev)
LHS_ChronicPrev$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)
LHS_C<-melt(LHS_ChronicPrev, id = "Year")

#Chronic Cases - 15-19

LHS_ChronicPrev1<-read.csv("LHS_ChronicPrev1.csv", header=TRUE)
summary(LHS_ChronicPrev1)
LHS_ChronicPrev1$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)
LHS1_C<-melt(LHS_ChronicPrev1, id = "Year")

#Chronic Cases - 20-25

LHS_ChronicPrev2<-read.csv("LHS_ChronicPrev2.csv", header=TRUE)
summary(LHS_ChronicPrev2)
LHS_ChronicPrev2$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)
LHS2_C<-melt(LHS_ChronicPrev2, id = "Year")

#Chronic Cases - 26-29

LHS_ChronicPrev3<-read.csv("LHS_ChronicPrev3.csv", header=TRUE)
summary(LHS_ChronicPrev3)
LHS_ChronicPrev3$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)
LHS3_C<-melt(LHS_ChronicPrev3, id = "Year")

#Chronic Cases - 30-64

LHS_ChronicPrev4<-read.csv("LHS_ChronicPrev4.csv", header=TRUE)
summary(LHS_ChronicPrev4)
LHS_ChronicPrev4$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)
LHS4_C<-melt(LHS_ChronicPrev4, id = "Year")

#### Import Acute Cases (Current and Former PWID) ####

#Note that header row (var names) is the X + index of the parameter set (eg XParamSetNo.)

#Acute Cases - All Ages

LHS_Acute<-read.csv("LHS_Acute.csv", header=TRUE)
summary(LHS_Acute)
LHS_Acute$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

AcuteLHS.melted<-melt(LHS_Acute, id = "Year")
summary(AcuteLHS.melted)

#Acute Cases - 15-19

LHS_Acute1<-read.csv("LHS_Acute1.csv", header=TRUE)
summary(LHS_Acute1)
LHS_Acute1$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

keeps <- c("Year","acute1529idu_1")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute1, id = "Year")
AcuteLHS.melted_1<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_1)

#Acute Cases - 20-25

LHS_Acute2<-read.csv("LHS_Acute2.csv", header=TRUE)
summary(LHS_Acute2)
LHS_Acute2$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

keeps <- c("Year","acute1529idu_2")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute2, id = "Year")
AcuteLHS.melted_2<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_2)

#Acute Cases - 26-29

LHS_Acute3<-read.csv("LHS_Acute3.csv", header=TRUE)
summary(LHS_Acute3)
LHS_Acute3$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

keeps <- c("Year","acute1529idu_3")
mdhhs1<-mdhhs[keeps]

LHS<-melt(LHS_Acute3, id = "Year")
AcuteLHS.melted_3<-merge(LHS, mdhhs1, by ="Year", all=TRUE) 
summary(AcuteLHS.melted_3)

#Acute Cases - 30-64

LHS_Acute4<-read.csv("LHS_Acute4.csv", header=TRUE)
summary(LHS_Acute4)
LHS_Acute4$Year[(1:81)]<-seq(from = 2000, to = 2016, by = 0.2)

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

#Create a numeric quartile variable (note: change the bounds as needed)
params$fitorder<-ifelse(params$fitquartile=="[229,263]",1,
                 ifelse(params$fitquartile=="(263,297]",2,
                 ifelse(params$fitquartile=="(297,506]",3,4)))
table(params$fitorder) 

#Best 10% of paramsets
params$best10<-cut(params$RSS, quantile(params$RSS, c(0,1/10)), include.lowest=T) 
table(params$best10, useNA="ifany")

is.na(params$best10)

params$best10<-ifelse(params$best10=="[229,252]",1,0)

params$best10<-ifelse(is.na(params$best10),0,params$best10)

table(params$best10, useNA="ifany")


#### Dataset of case counts with params ####

paramsRSS<-params
paramsRSS$variable<-paste("X", paramsRSS$index, sep="")

LHS.melted<-merge(LHS.melted,paramsRSS,by="variable")
NewChronicLHS.melted_1529<-merge(NewChronicLHS.melted_1529,paramsRSS,by="variable")
LHS_C <-merge(LHS_C,paramsRSS,by="variable")
AcuteLHS.melted<-merge(AcuteLHS.melted,paramsRSS,by="variable")
AcuteLHS.melted_1529<-merge(AcuteLHS.melted_1529,paramsRSS,by="variable")

LHS1.melted<-merge(LHS1.melted,paramsRSS,by="variable")
LHS1_C <-merge(LHS1_C,paramsRSS,by="variable")
AcuteLHS.melted_1<-merge(AcuteLHS.melted_1,paramsRSS,by="variable")

LHS2.melted<-merge(LHS2.melted,paramsRSS,by="variable")
LHS2_C <-merge(LHS2_C,paramsRSS,by="variable")
AcuteLHS.melted_2<-merge(AcuteLHS.melted_2,paramsRSS,by="variable")

LHS3.melted<-merge(LHS3.melted,paramsRSS,by="variable")
LHS3_C <-merge(LHS3_C,paramsRSS,by="variable")
AcuteLHS.melted_3<-merge(AcuteLHS.melted_3,paramsRSS,by="variable")

LHS4.melted<-merge(LHS4.melted,paramsRSS,by="variable")
LHS4_C <-merge(LHS4_C,paramsRSS,by="variable")
AcuteLHS.melted_4<-merge(AcuteLHS.melted_4,paramsRSS,by="variable")

rm(mdhhs1,keeps, paramsRSS,LHS_Acute,LHS_Acute1,LHS_Acute2,LHS_Acute3,LHS_Acute4,
        LHS_ChronicPrev,LHS_ChronicPrev1,LHS_ChronicPrev2,LHS_ChronicPrev3,LHS_ChronicPrev4,
        LHS_NewChronic,LHS_NewChronic1,LHS_NewChronic2,LHS_NewChronic3,LHS_NewChronic4)

save.image(file="/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 5.6.2018/Dataset Workspace_5.6.18.RData")

#### Import intervention datasets and combine ####
load(file="/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 5.6.2018/Dataset Workspace_5.6.18.RData")

#rm(param,prev,new,acute,prev1,new1,acute1,singleint,singleint1,seqint,seqint1,merged,interventioninfo)

#import datasets with intervention levels + all case counts at t1-t66
prev<-read.csv("ChronicPrev_SingleSeq.txt", header=F)
names(prev)[1:15] <- c("e","int1",	"int2",	"int3",	"int4",	"int5",	"int6",	"int7",	
                       "int8",	"int9",	"int10","int11","i","ChronicPrev_End_All","ChronicPrev_End")
#names(prev)[14:94] <- paste("t",1:81,sep="")
prev$inttype<-ifelse(prev$e<=169,"single",ifelse(prev$e>=220,"terttoprim","primtotert"))
table(prev$inttype)

new<-read.csv("NewChronic_SingleSeq.txt", header=F)
names(new)[1:15] <- c("e","int1",	"int2",	"int3",	"int4",	"int5",	"int6",	"int7",	
                      "int8",	"int9",	"int10","int11","i","NewChronic_End_All","NewChronic_End")
#names(new)[14:94] <- paste("t",1:81,sep="")
new$inttype<-ifelse(new$e<=169,"single",ifelse(new$e>=220,"terttoprim","primtotert"))
table(new$inttype)

acute<-read.csv("Acute_SingleSeq.txt", header=F)
names(acute)[1:15] <- c("e","int1",	"int2",	"int3",	"int4",	"int5",	"int6",	"int7",	
                        "int8",	"int9",	"int10","int11","i","Acute_End_All","Acute_End")
#names(acute)[14:94] <- paste("t",1:81,sep="")
acute$inttype<-ifelse(acute$e<=169,"single",ifelse(acute$e>=220,"terttoprim","primtotert"))
table(acute$inttype)

#Import parameters associated with each intervention
param<-read.csv("Parameters_SingleSeq.txt", header=F)
names(param)[1:13] <- c("e","int1",	"int2",	"int3",	"int4",	"int5",	"int6",	"int7",	
                        "int8",	"int9",	"int10","int11","i")
names(param)[14:65] <- c("sigma1","sigma2","sigma3","sigma4","a","b",	
                       "gamma1","gamma2","gamma3","gamma4",
                       "d","eps","zeta1","zeta2","zeta3","zeta4",
                       "etan",	"etap","etaz",
                       "theta1","theta2","theta3","theta4",
                       "k1","k2","k3","k4",
                       "lambda1","lambda2","lambda3","lambda4",
                       "mu1",	"mu2", "mu3","mu4","xi","omic1","omic2","omic3","omic4",
                       "r","tau","phip","phin","psi1","psi2","psi3","psi4","psin",
                       "w","intervention_etap","Z0")
param$inttype<-ifelse(param$e<=169,"single",ifelse(param$e>=220,"terttoprim","primtotert"))
table(param$inttype)

#Create datasets with case counts at t81 (when only final data point recorded)
prev1<-prev[ ,c(1:16)]
prev1$index<-as.numeric(paste(prev1$e,prev1$i,sep=""))
new1<-new[ ,c(1,13:15)]
new1$index<-as.numeric(paste(new1$e,new1$i,sep=""))
acute1<-acute[ ,c(1,13:15)]
acute1$index<-as.numeric(paste(acute1$e,acute1$i,sep=""))

#Setup params dataset (remove int variables)
param1<-param[ ,c(1,13:65)]
param1$index<-as.numeric(paste(param1$e,param1$i,sep=""))

#merge case count summaries with params
merged<-merge(merge(merge(prev1, new1, by=c("e","i"),suffixes = c(".a",".b")), 
                    acute1, by=c("e","i"), suffixes = c(".c",".d")), 
              param1, by=c("e","i"), suffixes = c(".e",""))


#Create variables for intervention level
df<-merged
attach(df)

int<-int1
df$thetapct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-int2
df$sigmapct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-int3
df$gammapct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-int4
df$kappapct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-int5
df$phippct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-int6
df$phinpct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-int7
df$etappct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-int8
df$etazpct<-ifelse(int==1,"None",ifelse(int==2,"10%",ifelse(int==3,"20%",ifelse(int==4,"30%",ifelse(int==5,"40%","error")))))
int<-int9
df$omegapct<-ifelse(int==1,"12 Weeks",ifelse(int==2,"8 Weeks",ifelse(int==3,"16 Weeks","error")))
int<-int10
df$taupct<-ifelse(int==1,"50%",ifelse(int==2,"100%",ifelse(int==3,"75%",ifelse(int==4,"25%",ifelse(int==5,"0%","error")))))
int<-int11
df$alphapct<-ifelse(int==1,"90%",ifelse(int==2,"100%",ifelse(int==3,"80%",ifelse(int==4,"70%",ifelse(int==5,"60%","error")))))

detach(df)
merged<-df

table(merged$thetapct)
table(merged$sigmapct)
table(merged$gammapct)
table(merged$kappapct)
table(merged$phippct)
table(merged$phinpct)
table(merged$etappct)
table(merged$etazpct)
table(merged$omegapct)
table(merged$taupct)
table(merged$alphapct)


#Add RSS and fit information to the merged dataset
paramfit<-params[,c(1,70,91,92,93,94)]
colnames(paramfit) <- c("i","RSS","best50","fitquartile","fitorder","best10")

merged<-merge(merged,paramfit, by="i",all.x=TRUE)

#Create dataset for single interventions 
singleint<-merged[ ! merged$inttype %in% c("primtotert","terttoprim"), ]

#also get 'both' treatment simulation from tertoprim interventions dataset
add<-subset(merged,merged$e %in% c(250:253))

singleint<-rbind(singleint,add)

#Create dataset for sequential interventions 

seqint<-merged[ ! merged$inttype %in% c("single"), ]

#Identify the intervention parameter in singleint
singleint$InterventionParam<-ifelse(singleint$e==1,"None",
                             ifelse(2<=singleint$e&singleint$e<=5,"theta",
                             ifelse(6<=singleint$e&singleint$e<=9,"sigma",
                             ifelse(10<=singleint$e&singleint$e<=13,"gamma",
                             ifelse(14<=singleint$e&singleint$e<=17,"kappa",
                             ifelse(18<=singleint$e&singleint$e<=21,"phip",
                             ifelse(22<=singleint$e&singleint$e<=25,"phin",
                             ifelse(26<=singleint$e&singleint$e<=29,"etap",
                             ifelse(30<=singleint$e&singleint$e<=33,"etaz",
                             ifelse(singleint$e==34,"None_eta",
                             ifelse(35<=singleint$e&singleint$e<=38,"theta_eta", 
                             ifelse(39<=singleint$e&singleint$e<=42,"sigma_eta",
                             ifelse(43<=singleint$e&singleint$e<=46,"gamma_eta",
                             ifelse(47<=singleint$e&singleint$e<=50,"kappa_eta",
                             ifelse(51<=singleint$e&singleint$e<=54,"phip_eta",
                             ifelse(55<=singleint$e&singleint$e<=58,"phin_eta",
                             ifelse(59<=singleint$e&singleint$e<=61,"eta_both",
                             ifelse(62<=singleint$e&singleint$e<=65,"phi_both",
                             ifelse(66<=singleint$e&singleint$e<=69,"phip_omega",
                             ifelse(70<=singleint$e&singleint$e<=73,"phin_omega",
                             ifelse(74<=singleint$e&singleint$e<=77,"phiboth_omega",
                             ifelse(78<=singleint$e&singleint$e<=81,"phip_omega",
                             ifelse(82<=singleint$e&singleint$e<=85,"phin_omega",
                             ifelse(86<=singleint$e&singleint$e<=89,"phiboth_omega",      
                             ifelse(90<=singleint$e&singleint$e<=93,"phip_tau",    
                             ifelse(94<=singleint$e&singleint$e<=97,"phiboth_tau", 
                             ifelse(98<=singleint$e&singleint$e<=101,"phip_tau",    
                             ifelse(102<=singleint$e&singleint$e<=105,"phiboth_tau",
                             ifelse(106<=singleint$e&singleint$e<=109,"phip_tau",    
                             ifelse(110<=singleint$e&singleint$e<=113,"phiboth_tau",
                             ifelse(114<=singleint$e&singleint$e<=117,"phip_tau",    
                             ifelse(118<=singleint$e&singleint$e<=121,"phiboth_tau",
                             ifelse(122<=singleint$e&singleint$e<=125,"phip_alpha",
                             ifelse(126<=singleint$e&singleint$e<=129,"phin_alpha",
                             ifelse(130<=singleint$e&singleint$e<=133,"phiboth_alpha",
                             ifelse(134<=singleint$e&singleint$e<=137,"phip_alpha",
                             ifelse(138<=singleint$e&singleint$e<=141,"phin_alpha",
                             ifelse(142<=singleint$e&singleint$e<=145,"phiboth_alpha",
                             ifelse(146<=singleint$e&singleint$e<=149,"phip_alpha",
                             ifelse(150<=singleint$e&singleint$e<=153,"phin_alpha",
                             ifelse(154<=singleint$e&singleint$e<=157,"phiboth_alpha",
                             ifelse(158<=singleint$e&singleint$e<=161,"phip_alpha",
                             ifelse(162<=singleint$e&singleint$e<=165,"phin_alpha",
                             ifelse(166<=singleint$e&singleint$e<=169,"phiboth_alpha",
                             ifelse(250<=singleint$e&singleint$e<=253,"phiboth_eta","wrong"))))))))))))))))))))))))))))))))))))))))))))) 
table(singleint$InterventionParam)

#Identify the intervention level in singleint
singleint$IntLevelpct<-ifelse(singleint$e %in% c(1,34),"None",
                       ifelse(singleint$e %in% c(2,6,10,14,18,22,26,30,35,39,43,47,51,55,59,62,66,70,74,78,82,86,90,94,98,102,106,110,114,118,122,126,130,134,138,142,146,150,154,158,162,166,250),"10%",
                       ifelse(singleint$e %in% c(3,7,11,15,19,23,27,31,36,40,44,48,52,56,60,63,67,71,75,79,83,87,91,95,99,103,107,111,115,119,123,127,131,135,139,143,147,151,155,159,163,167,251),"20%",
                       ifelse(singleint$e %in% c(4,8,12,16,20,24,28,32,37,41,45,49,53,57,61,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,252),"30%",
                       ifelse(singleint$e %in% c(5,9,13,17,21,25,29,33,38,42,46,50,54,58,65,69,73,77,81,85,89,93,97,101,105,109,113,117,121,125,129,133,137,141,145,149,153,157,161,165,169,253),"40%","s")))))
table(singleint$IntLevelpct)

#Interventions evaluating overdose
singleint$ODint<-ifelse(singleint$e>=26&singleint$e<=61,1,
                        ifelse(singleint$e>=250,1,0))
table(singleint$ODint)
table(subset(singleint,singleint$e>=250)$ODint)

#Interventions evaluating treatment params (sensitivity analyses)
singleint$TRTsens<-ifelse(singleint$e>=62&singleint$e<=169,1,0)
table(singleint$TRTsens)

#Set an order variable to determine plotting order in graphs
singleint$order<-factor(ifelse(singleint$InterventionParam %in% c("None","None_eta"),1,
                  ifelse(singleint$InterventionParam %in% c("theta","theta_eta"),2,
                  ifelse(singleint$InterventionParam %in% c("sigma","sigma_eta"),3,
                  ifelse(singleint$InterventionParam %in% c("gamma","gamma_eta"),4,
                  ifelse(singleint$InterventionParam %in% c("kappa","kappa_eta"),5,
                  ifelse(singleint$InterventionParam %in% c("phip","phip_eta"),6,
                  ifelse(singleint$InterventionParam %in% c("phin","phin_eta"),7, 
                  ifelse(singleint$InterventionParam %in% c("phi_both","phiboth_eta"),8,
                  ifelse(singleint$InterventionParam %in% c("etap"),9,
                  ifelse(singleint$InterventionParam %in% c("etaz"),10,
                  ifelse(singleint$InterventionParam %in% c("eta_both"),11,
                  9999))))))))))))
table(singleint$order)
table(subset(singleint,singleint$order==9999)$InterventionParam)


#Calculate % reduction in acute and chronic cases

#select the subset dataset of 'none' interventions to calculate % reduction in seqint
options(digits = 10)
noint<-subset(singleint,singleint$e %in% c(1,34))
table(noint$e)
noint <- noint[order(noint$e),] 

df<-noint[1:10000,c(1,14,15,18,19,21,22)] #select e=1 (vars=i and case counts)
df<-rename(df, c("ChronicPrev_End" = "ChronicPrev_e1", "Acute_End" = "Acute_e1", "NewChronic_End" = "NewChronic_e1",
                 "ChronicPrev_End_All" = "ChronicPrev_All_e1", "Acute_End_All" = "Acute_All_e1", 
                 "NewChronic_End_All" = "NewChronic_All_e1"))

df$ChronicPrev_e34<-subset(noint,noint$e==34)$ChronicPrev_End
df$Acute_e34<-subset(noint,noint$e==34)$Acute_End
df$NewChronic_e34<-subset(noint,noint$e==34)$NewChronic_End

df$ChronicPrev_All_e34<-subset(noint,noint$e==34)$ChronicPrev_End_All
df$Acute_All_e34<-subset(noint,noint$e==34)$Acute_End_All
df$NewChronic_All_e34<-subset(noint,noint$e==34)$NewChronic_End_All


#merge back with intervention dataset 
singleint<-merge(singleint,df,by="i")
rm(df,noint,add)

# Calculate % reduction compared to appropriate none intervention
singleint$pctred_chr<-ifelse(singleint$e %in% c(1,34),0,
                          ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$ChronicPrev_End/singleint$ChronicPrev_e1)*100,
                                 ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$ChronicPrev_End/singleint$ChronicPrev_e34)*100,99999)))

singleint$pctred_new<-ifelse(singleint$e %in% c(1,34),0,
                          ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$NewChronic_End/singleint$NewChronic_e1)*100,
                                 ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$NewChronic_End/singleint$NewChronic_e34)*100,99999)))

singleint$pctred_acute<-ifelse(singleint$e %in% c(1,34),0,
                            ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$Acute_End/singleint$Acute_e1)*100,
                                   ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$Acute_End/singleint$Acute_e34)*100,99999)))

summary(singleint$pctred_chr)
summary(singleint$pctred_new)
summary(singleint$pctred_acute)

# Calculate % reduction compared to appropriate none intervention
singleint$pctred_chr_all<-ifelse(singleint$e %in% c(1,34),0,
                             ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$ChronicPrev_End_All/singleint$ChronicPrev_All_e1)*100,
                                    ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$ChronicPrev_End_All/singleint$ChronicPrev_All_e34)*100,99999)))

singleint$pctred_new_all<-ifelse(singleint$e %in% c(1,34),0,
                             ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$NewChronic_End_All/singleint$NewChronic_All_e1)*100,
                                    ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$NewChronic_End_All/singleint$NewChronic_All_e34)*100,99999)))

singleint$pctred_acute_all<-ifelse(singleint$e %in% c(1,34),0,
                               ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$Acute_End_All/singleint$Acute_All_e1)*100,
                                      ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$Acute_End_All/singleint$Acute_All_e34)*100,99999)))

summary(singleint$pctred_chr_all)
summary(singleint$pctred_new_all)
summary(singleint$pctred_acute_all)


#Calculate Summary Statistics of % Reduction
df1<-tapply(singleint$pctred_chr,list(singleint$e),mean)
df2<-tapply(singleint$pctred_chr,list(singleint$e),sd)
df3<-tapply(singleint$pctred_chr,list(singleint$e),min)
df4<-tapply(singleint$pctred_chr,list(singleint$e),max)
df5<-tapply(singleint$pctred_acute,list(singleint$e),mean)
df6<-tapply(singleint$pctred_acute,list(singleint$e),sd)
df7<-tapply(singleint$pctred_acute,list(singleint$e),min)
df8<-tapply(singleint$pctred_acute,list(singleint$e),max)
df9<-tapply(singleint$pctred_chr_all,list(singleint$e),mean)
df10<-tapply(singleint$pctred_chr_all,list(singleint$e),sd)
df11<-tapply(singleint$pctred_chr_all,list(singleint$e),min)
df12<-tapply(singleint$pctred_chr_all,list(singleint$e),max)
df13<-tapply(singleint$pctred_acute_all,list(singleint$e),mean)
df14<-tapply(singleint$pctred_acute_all,list(singleint$e),sd)
df15<-tapply(singleint$pctred_acute_all,list(singleint$e),min)
df16<-tapply(singleint$pctred_acute_all,list(singleint$e),max)

singleint1<-data.frame(cbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16))
singleint1$e<-rownames(singleint1)
singleint1<-rename(singleint1,c("df1"="meanpctred_chr","df2"="sdpctred_chr","df3"="minpctred_chr",
                          "df4"="maxpctred_chr","df5"="meanpctred_acute","df6"="sdpctred_acute",
                          "df7"="minpctred_acute", "df8"="maxpctred_acute",
                          "df9"="meanpctred_chr_all","df10"="sdpctred_chr_all","df11"="minpctred_chr_all",
                          "df12"="maxpctred_chr_all","df13"="meanpctred_acute_all","df14"="sdpctred_acute_all",
                          "df15"="minpctred_acute_all", "df16"="maxpctred_acute_all"))


# Subset information on each 'e' intervention
interventioninfo <- singleint[,c(2:13,16,77:87,93:96)] #vars=things that do not change by e
interventioninfo<-unique(interventioninfo)
table(interventioninfo$e)

# Merge with singleint1
singleint1<-merge(singleint1,interventioninfo,by="e")


# % Reduction in Median
summary(singleint1$meanpctred_chr)
summary(singleint1$meanpctred_acute)
summary(singleint1$minpctred_chr)
summary(singleint1$minpctred_acute)
summary(singleint1$maxpctred_chr)
summary(singleint1$maxpctred_acute)



#Seqint Formatting

#Create variable for intervention level
seqint$pct<-ifelse(seqint$e %in% c(170,195,220,245),"None",
            ifelse(seqint$e %in% c(171,175,179,183,187,191,196,200,204,208,212,216,
                                   221,225,229,233,237,241,246,250,254,258,262,266),"10%",
            ifelse(seqint$e %in% c(172,176,180,184,188,192,197,201,205,209,213,217,
                                   222,226,230,234,238,242,247,251,255,259,263,267),"20%",
            ifelse(seqint$e %in% c(173,177,181,185,189,193,198,202,206,210,214,218,
                                   223,227,231,235,239,243,248,252,256,260,264,268),"30%",
            ifelse(seqint$e %in% c(174,178,182,186,190,194,199,203,207,211,215,219,
                                   224,228,232,236,240,244,249,253,257,261,265,269),"40%","s")))))
table(seqint$pct)

#Create variable for order of interventions when plotting
seqint$intono<-ifelse(seqint$e %in% c(170,195,220,245),1,
               ifelse(seqint$e %in% c(171:174,196:199,221:224,246:249),2,
               ifelse(seqint$e %in% c(175:178,200:203,225:228,250:253),3,
               ifelse(seqint$e %in% c(179:182,204:207,229:232,254:257),4,
               ifelse(seqint$e %in% c(183:186,208:211,233:236,258:261),5,
               ifelse(seqint$e %in% c(187:190,212:215,237:240,262:265),6,
               ifelse(seqint$e %in% c(191:194,216:219,241:244,266:269),7,9999)))))))
table(seqint$intono)

#select the subset dataset of 'none' interventions to calculate % reduction in seqint
options(digits = 10)
noint<-subset(seqint,seqint$e %in% c(170,195))
table(noint$e)
noint <- noint[order(noint$e),] 

df<-noint[1:10000,c(1,14,15,18,19,21,22)] #select e=170 (vars=i and case counts)
df<-rename(df, c("ChronicPrev_End" = "ChronicPrev_e170", "Acute_End" = "Acute_e170", "NewChronic_End" = "NewChronic_e170",
                 "ChronicPrev_End_All" = "ChronicPrev_All_e170", "Acute_End_All" = "Acute_All_e170", "NewChronic_End_All" = "NewChronic_All_e170"))

df$ChronicPrev_e195<-subset(noint,noint$e==195)$ChronicPrev_End
df$Acute_e195<-subset(noint,noint$e==195)$Acute_End
df$NewChronic_e195<-subset(noint,noint$e==195)$NewChronic_End

df$ChronicPrev_All_e195<-subset(noint,noint$e==195)$ChronicPrev_End_All
df$Acute_All_e195<-subset(noint,noint$e==195)$Acute_End_All
df$NewChronic_All_e195<-subset(noint,noint$e==195)$NewChronic_End_All


#merge back with intervention dataset 
seqint<-merge(seqint,df,by="i")
rm(df,noint)

# Calculate % reduction compared to appropriate none intervention
seqint$pctred_chr<-ifelse(seqint$e %in% c(170,195,220,245),0,
                   ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$ChronicPrev_End/seqint$ChronicPrev_e170)*100,
                   ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$ChronicPrev_End/seqint$ChronicPrev_e195)*100,99999)))

seqint$pctred_new<-ifelse(seqint$e %in% c(170,195,220,245),0,
                   ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$NewChronic_End/seqint$NewChronic_e170)*100,
                   ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$NewChronic_End/seqint$NewChronic_e195)*100,99999)))

seqint$pctred_acute<-ifelse(seqint$e %in% c(170,195,220,245),0,
                     ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$Acute_End/seqint$Acute_e170)*100,
                     ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$Acute_End/seqint$Acute_e195)*100,99999)))

summary(seqint$pctred_chr)
summary(seqint$pctred_new)
summary(seqint$pctred_acute)

seqint$pctred_chr_all<-ifelse(seqint$e %in% c(170,195,220,245),0,
                          ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$ChronicPrev_End_All/seqint$ChronicPrev_All_e170)*100,
                                 ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$ChronicPrev_End_All/seqint$ChronicPrev_All_e195)*100,99999)))

seqint$pctred_new_all<-ifelse(seqint$e %in% c(170,195,220,245),0,
                          ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$NewChronic_End_All/seqint$NewChronic_All_e170)*100,
                                 ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$NewChronic_End_All/seqint$NewChronic_All_e195)*100,99999)))

seqint$pctred_acute_all<-ifelse(seqint$e %in% c(170,195,220,245),0,
                            ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$Acute_End_All/seqint$Acute_All_e170)*100,
                                   ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$Acute_End_All/seqint$Acute_All_e195)*100,99999)))

summary(seqint$pctred_chr_all)
summary(seqint$pctred_new_all)
summary(seqint$pctred_acute_all)

#Calculate Summary Statistics of % Reduction
df1<-tapply(seqint$pctred_chr,list(seqint$e),mean)
df2<-tapply(seqint$pctred_chr,list(seqint$e),sd)
df3<-tapply(seqint$pctred_chr,list(seqint$e),min)
df4<-tapply(seqint$pctred_chr,list(seqint$e),max)
df5<-tapply(seqint$pctred_acute,list(seqint$e),mean)
df6<-tapply(seqint$pctred_acute,list(seqint$e),sd)
df7<-tapply(seqint$pctred_acute,list(seqint$e),min)
df8<-tapply(seqint$pctred_acute,list(seqint$e),max)
df9<-tapply(seqint$pctred_chr_all,list(seqint$e),mean)
df10<-tapply(seqint$pctred_chr_all,list(seqint$e),sd)
df11<-tapply(seqint$pctred_chr_all,list(seqint$e),min)
df12<-tapply(seqint$pctred_chr_all,list(seqint$e),max)
df13<-tapply(seqint$pctred_acute_all,list(seqint$e),mean)
df14<-tapply(seqint$pctred_acute_all,list(seqint$e),sd)
df15<-tapply(seqint$pctred_acute_all,list(seqint$e),min)
df16<-tapply(seqint$pctred_acute_all,list(seqint$e),max)

seqint1<-data.frame(cbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16))
seqint1$e<-rownames(seqint1)
seqint1<-rename(seqint1,c("df1"="meanpctred_chr","df2"="sdpctred_chr","df3"="minpctred_chr",
                     "df4"="maxpctred_chr","df5"="meanpctred_acute","df6"="sdpctred_acute",
                     "df7"="minpctred_acute", "df8"="maxpctred_acute",
                     "df9"="meanpctred_chr_all","df10"="sdpctred_chr_all","df11"="minpctred_chr_all",
                     "df12"="maxpctred_chr_all","df13"="meanpctred_acute_all","df14"="sdpctred_acute_all",
                     "df15"="minpctred_acute_all", "df16"="maxpctred_acute_all"))


# Subset information on each 'e' intervention
interventioninfo <- seqint[,c(2:13,16,77:87,93:94)]
interventioninfo<-unique(interventioninfo)

# Merge with seqint1
seqint1<-merge(seqint1,interventioninfo,by="e")


# % Reduction in Median
summary(seqint1$meanpctred_chr)
summary(seqint1$meanpctred_acute)
summary(seqint1$minpctred_chr)
summary(seqint1$minpctred_acute)
summary(seqint1$maxpctred_chr)
summary(seqint1$maxpctred_acute)


#### Save and/or Load R Workspace ####
# rm(acute1,LHS_Acute,LHS_Acute1,LHS_Acute2,LHS_Acute3,LHS_Acute4,
#     LHS_ChronicPrev,LHS_ChronicPrev1,LHS_ChronicPrev2,LHS_ChronicPrev3,LHS_ChronicPrev4,
#     LHS_NewChronic,LHS_NewChronic1,LHS_NewChronic2,LHS_NewChronic3,LHS_NewChronic4,
#     new1,prev1,paramfit,param1,df1,df2,df3,df4,df5,df6,df7,df8,int,df9,df10,df11,df12,df13,df14,df15,df16,
#     prev,acute,new)
# 
# save.image(file="/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 5.6.2018/Dataset Workspace_5.6.18.RData")

load("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 5.6.2018/Dataset Workspace_5.6.18.RData")






 ############### Section 4: Latin Hypercube Sampling Results ################
#### Plot of Acute Cases Fit to MDHHS Data (50% Best, 1 Color) #####
#Plotting Parameters
xbreaks <- c(2000, 2005, 2010, 2015)
ybreaks <- c(0,25,50)
ylabels <- c("0","25","50")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

#Create a plot for all age groups combined and each age group separately

#All ages - with right legend
p1<-ggplot(subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best50==1), aes(x=Year, y=value1*r, by=variable)) + 
  geom_line(aes(x=Year, y=value1*r, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=acute1529idu, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 50), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "right", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_manual(name=NULL, values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Cases: Ages 15-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1

#All ages - no legend 
p2<-ggplot(subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best50==1), aes(x=Year, y=value1*r, by=variable)) + 
  geom_line(aes(x=Year, y=value1*r, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=acute1529idu, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 50), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name="", values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Cases: Ages 15-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p2
# title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = 0, vjust = 0,gp = gpar(fontsize = 16))
# 
# p2 <- arrangeGrob(p2, top = title.grob)
# grid.arrange(p2)

#15-19 
p3<-ggplot(subset(AcuteLHS.melted_1,AcuteLHS.melted_1$best50==1), aes(x=Year, y=value*r, by=variable)) + 
  geom_line(aes(x=Year, y=value*r, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=acute1529idu_1, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30), breaks = c(0,10,20,30))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name="", values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('Ages 15-19')+
  theme(plot.title = element_text(hjust = 0.5)) 
p3
# title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = 0, vjust = 0,gp = gpar(fontsize = 16))
# 
# p3 <- arrangeGrob(p3, top = title.grob)
# grid.arrange(p3)

#20-25
p4<-ggplot(subset(AcuteLHS.melted_2,AcuteLHS.melted_2$best50==1), aes(x=Year, y=value*r, by=variable)) + 
  geom_line(aes(x=Year, y=value*r, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=acute1529idu_2, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30), breaks = c(0,10,20,30))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name="", values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab)+
  ggtitle('Ages 20-25')+
  theme(plot.title = element_text(hjust = 0.5)) 
p4

# title.grob <- textGrob(label = "C)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = 0, vjust = 0,gp = gpar(fontsize = 16))
# 
# p4 <- arrangeGrob(p4, top = title.grob)
# grid.arrange(p4)

#26-29
p5<-ggplot(subset(AcuteLHS.melted_3,AcuteLHS.melted_3$best50==1), aes(x=Year, y=value*r, by=variable)) + 
  geom_line(aes(x=Year, y=value*r, by=variable, colour='Model'), show.legend=TRUE)+
  geom_point(aes(x=Year, y=acute1529idu_3, colour="Data"), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 30), breaks = c(0,10,20,30))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA')) +
  scale_colour_manual(name="", values = c("#000000", "#999999"), breaks = c("Model", "Data")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("") +
  xlab(text.xlab) +
  ggtitle('Ages 26-29')+
  theme(plot.title = element_text(hjust = 0.5)) 
p5

# title.grob <- textGrob(label = "D)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = 0, vjust = 0,gp = gpar(fontsize = 16))
# 
# p5 <- arrangeGrob(p5, top = title.grob)
# grid.arrange(p5)

#compile and save as 1 figure
w <- 10; h <- 6
mylegend<-g_legend(p1)
lay<-rbind(c(1,1,1,1,1,2),c(1,1,1,1,1,2),c(1,1,1,1,1,2),c(3,3,4,4,5,5),c(3,3,4,4,5,5),c(3,3,4,4,5,5))
LHS_Acute<-grid.arrange(p2,mylegend,p3,p4,p5,layout_matrix=lay)
ggsave(sprintf("LHS_Acute.pdf"), LHS_Acute, width=w, height=h)


#### Plot of Acute Cases Fit to MDHHS Data (All, Color by RSS) #####
#Plotting Parameters
summary(AcuteLHS.melted_1529$value1*AcuteLHS.melted_1529$r)

xbreaks <- c(2000, 2005, 2010, 2015)
ybreaks <- c(0,50,100,150,200,250,300)
ylabels <- c("0","50","100","150","200","250","300")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

#Dummy plots for legends:
dummy<-ggplot(data=params,aes(x=gamma1,y=k1,color=RSS))+
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



dummy1<-ggplot(AcuteLHS.melted_1529)+
  geom_point(aes(x=Year, y=acute1529idu,colour="Data"),show.legend = TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 320), breaks = ybreaks, labels=ylabels)+
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 320), breaks = ybreaks, labels=ylabels)+
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

ybreaks <- c(0,30,60,90,120)
ylabels <- c("0","30","60","90","120")

p2<-ggplot(AcuteLHS.melted_1, aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_1),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 136), breaks = ybreaks, labels=ylabels)+
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 136), breaks = ybreaks, labels=ylabels)+
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 136), breaks = ybreaks, labels=ylabels)+
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

p5<-ggplot(AcuteLHS.melted_4, aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  #geom_point(aes(x=Year, y=acute1529idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
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

#layout w/legend at side 
w <- 10; h <- 6

#No plot title
lay<-rbind(c(1,1,1,1,1,2),c(1,1,1,1,1,3),c(4,4,5,5,6,6),c(4,4,5,5,6,6))
LHS_Acute<-grid.arrange(p1,legend1,Datalegend,p2,p3,p4,layout_matrix=lay)
ggsave(sprintf("LHS_Acute_RSSColor.pdf"), LHS_Acute, width=w, height=h)
ggsave(sprintf("LHS_Acute_RSSColor.png"), LHS_Acute, width=w, height=h)

#All Ages
w <- 15; h <- 4
lay<-rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5),c(1,1,1,2,2,2,3,3,3,4,4,4,6))
LHS_Acute<-grid.arrange(p2,p3,p4,p5,legend1,Datalegend,layout_matrix=lay)
ggsave(sprintf("LHS_Acute_RSSColor_1564.pdf"), LHS_Acute, width=w, height=h)
ggsave(sprintf("LHS_Acute_RSSColor_1564.png"), LHS_Acute, width=w, height=h)

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


#Include Plot Title
# w <- 10; h <- 7
# lay<-rbind(c(1,1,1,1,1,1),c(2,2,2,2,2,3),c(2,2,2,2,2,3),c(2,2,2,2,2,4),c(2,2,2,2,2,4),
#            c(5,5,6,6,7,7),c(5,5,6,6,7,7),c(5,5,6,6,7,7),c(5,5,6,6,7,7))
# LHS_Acute<-grid.arrange(grid.text("Web Figure 1. Model Fit to Acute HCV Cases Detected by Public Health
# Surveillance in Michigan among PWID Aged 15-29 Years",gp=gpar(fontsize=18, col="black")),
#                              p1,legend1,Datalegend,p2,p3,p4,layout_matrix=lay)
# LHS_Acute
# ggsave(sprintf("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 4.16.2018/LHS_Acute_RSSColor_Title.pdf"), LHS_Acute, width=w, height=h)

#### Plot of Acute Cases Fit to MDHHS Data (Best 50%, Color by RSS) #####
#Plotting Parameters

summary(subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best50==1)$value1*subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best50==1)$r)

xbreaks <- c(2000, 2005, 2010, 2015)
ybreaks <- c(0,20,40)
ylabels <- c("0","20","40")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

#Dummy plots for legends:
dummy<-ggplot(data=subset(params,params$best50==1),aes(x=gamma1,y=k1,color=RSS))+
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 200), breaks = ybreaks, labels=ylabels)+
  theme_bw(base_size=18) +
  scale_colour_manual(name=NULL, values = c("#000000"), breaks = c("Data")) +  
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Age Groups')+
  theme(plot.title = element_text(hjust = 0.5))
dummy1
Datalegend<-g_legend(dummy1)

sub<-subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best50==1)

#Plot of fit and data
p1<-ggplot(subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best50==1), 
           aes(x=Year, y=value1*r, by=variable, 
               color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 41), breaks = ybreaks, labels=ylabels)+
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

summary(subset(AcuteLHS.melted_1,AcuteLHS.melted_1$best50==1)$value*subset(AcuteLHS.melted_1,AcuteLHS.melted_1$best50==1)$r)
summary(subset(AcuteLHS.melted_2,AcuteLHS.melted_2$best50==1)$value*subset(AcuteLHS.melted_2,AcuteLHS.melted_2$best50==1)$r)
summary(subset(AcuteLHS.melted_3,AcuteLHS.melted_3$best50==1)$value*subset(AcuteLHS.melted_3,AcuteLHS.melted_3$best50==1)$r)

ybreaks <- c(0,10,20)
ylabels <- c("0","10","20")

p2<-ggplot(subset(AcuteLHS.melted_1,AcuteLHS.melted_1$best50==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_1),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 23), breaks = ybreaks, labels=ylabels)+
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

p3<-ggplot(subset(AcuteLHS.melted_2,AcuteLHS.melted_2$best50==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_2),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 23), breaks = ybreaks, labels=ylabels)+
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

p4<-ggplot(subset(AcuteLHS.melted_3,AcuteLHS.melted_3$best50==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 23), breaks = ybreaks, labels=ylabels)+
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

p5<-ggplot(subset(AcuteLHS.melted_4,AcuteLHS.melted_4$best50==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  #geom_point(aes(x=Year, y=acute1529idu_3),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
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

#layout w/legend at side 
w <- 10; h <- 6

#No plot title
lay<-rbind(c(1,1,1,1,1,2),c(1,1,1,1,1,3),c(4,4,5,5,6,6),c(4,4,5,5,6,6))
LHS_Acute<-grid.arrange(p1,legend1,Datalegend,p2,p3,p4,layout_matrix=lay)
ggsave(sprintf("LHS_Acute_RSSColor_Best50.pdf"), LHS_Acute, width=w, height=h)
ggsave(sprintf("LHS_Acute_RSSColor_Best50.png"), LHS_Acute, width=w, height=h)

#layout w/legend at side 
w <- 15; h <- 4

#No plot title
lay<-rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5),c(1,1,1,2,2,2,3,3,3,4,4,4,6))
LHS_Acute<-grid.arrange(p2,p3,p4,p5,legend1,Datalegend,layout_matrix=lay)
ggsave(sprintf("LHS_Acute_RSSColor_Best50_1564.pdf"), LHS_Acute, width=w, height=h)
ggsave(sprintf("LHS_Acute_RSSColor_Best50_1564.png"), LHS_Acute, width=w, height=h)

# #Include Plot Title
# w <- 10; h <- 7
# lay<-rbind(c(1,1,1,1,1,1),c(2,2,2,2,2,3),c(2,2,2,2,2,3),c(2,2,2,2,2,4),c(2,2,2,2,2,4),
#            c(5,5,6,6,7,7),c(5,5,6,6,7,7),c(5,5,6,6,7,7),c(5,5,6,6,7,7))
# LHS_Acute<-grid.arrange(grid.text("Web Figure 1. Model Fit to Acute HCV Cases Detected by Public Health
#                                   Surveillance in Michigan among PWID Aged 15-29 Years",gp=gpar(fontsize=18, col="black")),
#                         p1,legend1,Datalegend,p2,p3,p4,layout_matrix=lay)
# LHS_Acute
# ggsave(sprintf("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 4.16.2018/LHS_Acute_RSSColor_Title_Best50.pdf"), LHS_Acute, width=w, height=h)



#### Plot of Acute Cases Fit to MDHHS Data (Best 10%, Color by RSS) #####
#Plotting Parameters
summary(subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best10==1)$value1*subset(AcuteLHS.melted_1529,AcuteLHS.melted_1529$best10==1)$r)

xbreaks <- c(2000, 2005, 2010, 2015)
ybreaks <- c(0,20,40)
ylabels <- c("0","20","40")
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 200), breaks = ybreaks, labels=ylabels)+
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 40), breaks = ybreaks, labels=ylabels)+
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

ybreaks <- c(0,10,20)
ylabels <- c("0","10","20")

p2<-ggplot(subset(AcuteLHS.melted_1,AcuteLHS.melted_1$best10==1), aes(x=Year, y=value*r, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  geom_point(aes(x=Year, y=acute1529idu_1),colour="#000000",show.legend = FALSE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 22), breaks = ybreaks, labels=ylabels)+
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 22), breaks = ybreaks, labels=ylabels)+
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 22), breaks = ybreaks, labels=ylabels)+
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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
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

#layout w/legend at side 
w <- 10; h <- 6

#No plot title
lay<-rbind(c(1,1,1,1,1,2),c(1,1,1,1,1,3),c(4,4,5,5,6,6),c(4,4,5,5,6,6))
LHS_Acute<-grid.arrange(p1,legend1,Datalegend,p2,p3,p4,layout_matrix=lay)
ggsave(sprintf("LHS_Acute_RSSColor_Best10.pdf"), LHS_Acute, width=w, height=h)
ggsave(sprintf("LHS_Acute_RSSColor_Best10.png"), LHS_Acute, width=w, height=h)

#layout w/legend at side 
w <- 15; h <- 4

#No plot title
lay<-rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5),c(1,1,1,2,2,2,3,3,3,4,4,4,6))
LHS_Acute<-grid.arrange(p2,p3,p4,p5,legend1,Datalegend,layout_matrix=lay)
ggsave(sprintf("LHS_Acute_RSSColor_Best10_1564.pdf"), LHS_Acute, width=w, height=h)
ggsave(sprintf("LHS_Acute_RSSColor_Best10_1564.png"), LHS_Acute, width=w, height=h)


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

# #Include Plot Title
# w <- 10; h <- 7
# lay<-rbind(c(1,1,1,1,1,1),c(2,2,2,2,2,3),c(2,2,2,2,2,3),c(2,2,2,2,2,4),c(2,2,2,2,2,4),
#            c(5,5,6,6,7,7),c(5,5,6,6,7,7),c(5,5,6,6,7,7),c(5,5,6,6,7,7))
# LHS_Acute<-grid.arrange(grid.text("Web Figure 1. Model Fit to Acute HCV Cases Detected by Public Health
#                                   Surveillance in Michigan among PWID Aged 15-29 Years",gp=gpar(fontsize=18, col="black")),
#                         p1,legend1,Datalegend,p2,p3,p4,layout_matrix=lay)
# LHS_Acute
# ggsave(sprintf("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 4.16.2018/LHS_Acute_RSSColor_Title_Best50.pdf"), LHS_Acute, width=w, height=h)



#### Chronic HCV Prevalence Plots (All, Combined Plot) ####

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


#### Chronic Prevalence by Age Group (Panel Plot, Best 50% Fits) ####
#Plotting Parameters
xbreaks <- c(2000, 2005, 2010, 2015)
ybreaks <- c(0,50000,100000,150000,200000)
ylabels <- c("0","50,000","100,000","150,000","200,000")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

summary(subset(LHS_C,LHS_C$best50==1)$value)
summary(subset(LHS1_C,LHS1_C$best50==1)$value)
summary(subset(LHS2_C,LHS2_C$best50==1)$value)
summary(subset(LHS3_C,LHS3_C$best50==1)$value)
summary(subset(LHS4_C,LHS4_C$best50==1)$value)

#Dummy plots for legends:
dummy<-ggplot(data=subset(params,params$best50==1),aes(x=gamma1,y=k1,color=RSS))+
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
p1<-ggplot(subset(LHS_C,LHS_C$best50==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(subset(LHS_C,LHS_C$best50==1)$value+10)), 
                     breaks = c(0,50000,100000,150000,200000,250000), 
                     labels=c("0","50,000","100,000","150,000","200,000","250,000"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Ages: 15-64')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1


p2<-ggplot(subset(LHS1_C,LHS1_C$best50==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 8180), breaks = c(0,2000,4000,6000,8000), labels=c("0","2,000","4,000","6,000","8,000"))+
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

p3<-ggplot(subset(LHS2_C,LHS2_C$best50==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 8180), breaks = c(0,2000,4000,6000,8000), labels=c("0","2,000","4,000","6,000","8,000"))+
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

p4<-ggplot(subset(LHS3_C,LHS3_C$best50==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 8180), breaks = c(0,2000,4000,6000,8000), labels=c("0","2,000","4,000","6,000","8,000"))+
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

p5<-ggplot(subset(LHS4_C,LHS4_C$best50==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 252000), breaks = c(0,50000,100000,150000,200000,250000), labels=c("0","50,000","100,000","150,000","200,000","250,000"))+
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
ggsave(sprintf("LHS_ChronicPrev_RSSColor_Best50.pdf"), LHS_Chronic, width=w, height=h)
ggsave(sprintf("LHS_ChronicPrev_RSSColor_Best50.png"), LHS_Chronic, width=w, height=h)

w <- 15; h <- 4

#No plot title
lay<-rbind(c(1,1,1,2,2,2,3,3,3,4,4,4,5))
LHS_Chronic<-grid.arrange(p2,p3,p4,p5,legend1,layout_matrix=lay)
ggsave(sprintf("LHS_ChronicPrev_RSSColor_Best50_1564.pdf"), LHS_Chronic, width=w, height=h)
ggsave(sprintf("LHS_ChronicPrev_RSSColor_Best50_1564.png"), LHS_Chronic, width=w, height=h)

#### Chronic Prevalence by Age Group (Panel Plot, Best 10% Fits) ####
#Plotting Parameters
xbreaks <- c(2000, 2005, 2010, 2015)
ybreaks <- c(0,50000,100000,150000,200000)
ylabels <- c("0","50,000","100,000","150,000","200,000")
text.xlab <- "Year"
text.ylab <- "Number of Cases"
model.colors <- c('#999999', '#000000') 

summary(subset(LHS_C,LHS_C$best10==1)$value)
summary(subset(LHS1_C,LHS1_C$best10==1)$value)
summary(subset(LHS2_C,LHS2_C$best10==1)$value)
summary(subset(LHS3_C,LHS3_C$best10==1)$value)
summary(subset(LHS4_C,LHS4_C$best10==1)$value)

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
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, max(subset(LHS_C,LHS_C$best10==1)$value+10)), 
                     breaks = c(0,50000,100000,150000), 
                     labels=c("0","50,000","100,000","150,000"))+
  theme_bw(base_size=18) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_gradientn(name="RSS",colours=rainbow(5))+
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  ggtitle('All Ages: 15-64')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1


p2<-ggplot(subset(LHS1_C,LHS1_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 8180), breaks = c(0,2000,4000,6000,8000), labels=c("0","2,000","4,000","6,000","8,000"))+
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

p3<-ggplot(subset(LHS2_C,LHS2_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 8180), breaks = c(0,2000,4000,6000,8000), labels=c("0","2,000","4,000","6,000","8,000"))+
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

p4<-ggplot(subset(LHS3_C,LHS3_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 8180), breaks = c(0,2000,4000,6000,8000), labels=c("0","2,000","4,000","6,000","8,000"))+
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

p5<-ggplot(subset(LHS4_C,LHS4_C$best10==1), aes(x=Year, y=value, by=variable, color=RSS)) + 
  geom_line(alpha=0.25)+
  scale_x_continuous(limits = c(2000, 2016), breaks = xbreaks)+
  scale_y_continuous(limits = c(0, 140000), breaks = c(0,30000,60000,90000,120000), labels=c("0","30,000","60,000","90,000","120,000"))+
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
  scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 2000), breaks = c(0,500,1000,1500,2000), labels=c("0","500","1,000","1,500","2,000"))+
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
  scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 2000), breaks = c(0,500,1000,1500,2000), labels=c("0","500","1,000","1,500","2,000"))+
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
  scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 2000), breaks = c(0,500,1000,1500,2000), labels=c("0","500","1,000","1,500","2,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p3

#sigma1
summary(params$sigma1)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 15-19 ('*sigma[1]*')')

p4<-ggplot(params) + 
  geom_histogram(aes(x=sigma1, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,breaks=seq(5,20,0.75))+#, 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  #scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p4

#sigma2
summary(params$sigma2)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 20-25 ('*sigma[2]*')')

p5<-ggplot(params) + 
  geom_histogram(aes(x=sigma2, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,breaks=seq(5,20,0.75))+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  #scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p5

#sigma3
summary(params$sigma3)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 26-29 ('*sigma[3]*')')

p6<-ggplot(params) + 
  geom_histogram(aes(x=sigma3, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,breaks=seq(5,20,0.75))+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  #scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p6

#sigma4
summary(params$sigma4)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 30-64 ('*sigma[4]*')')

p7<-ggplot(params) + 
  geom_histogram(aes(x=sigma4, fill = as.factor(reorder(fitorder,-fitorder))), colour="#000000",show.legend=TRUE,breaks=seq(5,20,0.75))+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  #scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
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
  geom_histogram(aes(x=theta1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.05,5,0.4), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2.5,5))+
  scale_y_continuous(limits = c(0, 5000), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p21


#theta2
summary(params$theta2)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 20-25 ('*theta[2]*')')


p22<-ggplot(params) + 
  geom_histogram(aes(x=theta2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.034,5,0.35), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2.5,5))+
  scale_y_continuous(limits = c(0, 5000), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p22

#theta3
summary(params$theta3)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 26-29 ('*theta[3]*')')


p23<-ggplot(params) + 
  geom_histogram(aes(x=theta3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.034,10,0.36), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,5,10))+
  scale_y_continuous(limits = c(0, 5000), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
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


# #lambda1
# summary(params$lambda1)
# 
# text.xlab <- bquote("%")
# text.title <- bquote('HCV Prevalence: 15-19 ('*lambda[1]*')')
# 
# 
# p29<-ggplot(params) + 
#   geom_histogram(aes(x=lambda1*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,30,1.2), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab(text.ylab) +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(0,30), breaks = c(0,10,20,30))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p29
# 
# 
# #lambda2
# summary(params$lambda2)
# 
# text.xlab <- bquote("%")
# text.title <- bquote('HCV Prevalence: 20-25 ('*lambda[2]*')')
# 
# 
# p30<-ggplot(params) + 
#   geom_histogram(aes(x=lambda2*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(3,40,1.5), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(0,40), breaks = c(0,20,40))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p30
# 
# #lambda3
# summary(params$lambda3)
# 
# text.xlab <- bquote("%")
# text.title <- bquote('HCV Prevalence: 26-29 ('*lambda[3]*')')
# 
# 
# p31<-ggplot(params) + 
#   geom_histogram(aes(x=lambda3*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(5,50,2), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(5,50), breaks = c(0,10,30,50))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p31
# 
# 
# #lambda4
# summary(params$lambda4)
# 
# text.xlab <- bquote("%")
# text.title <- bquote('HCV Prevalence: 30-64 ('*lambda[4]*')')
# 
# 
# p32<-ggplot(params) + 
#   geom_histogram(aes(x=lambda4*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(10,80,3), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(10,80), breaks = c(10,30,50,70))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p32


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
  geom_histogram(aes(x=xi*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.003,45,2), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,45), breaks = c(0,22.5,45))+
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
# summary(params$r)
# #summary(1/params$r)
# 
# text.xlab <- bquote('Surveillance Cases'*~Infections^-1)
# text.title <- bquote('Reporting Rate ('*rho*')')
# 
# p42<-ggplot(params) + 
#   geom_histogram(aes(x=r, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.02,0.082,0.003), show.legend=TRUE)+
#   theme_bw(base_size=14) +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + 
#   theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
#   scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
#   ylab("") +
#   xlab(text.xlab) +
#   scale_x_continuous(limits = c(0.02,0.082), breaks = c(0.02,0.04,0.06,0.08))+
#   scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
#   ggtitle(text.title)+
#   theme(plot.title = element_text(hjust = 0.5))
# p42


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
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 15-19 ('*pi["1,1"]*')')
p48<-ggplot(params) + 
  geom_histogram(aes(x=c1, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.5,9.5,0.35), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.5,9.5), breaks = c(1.5,5.5,9.5))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p48

#Contact Matrix: c21
summary(params$c2)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 20-25 ('*pi["2,1"]*')')
p49<-ggplot(params) + 
  geom_histogram(aes(x=c2, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.3,1.2,0.04), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.3,1.2), breaks = c(0.3,0.75,1.2))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p49

#Contact Matrix: c31
summary(params$c3)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 26-29 ('*pi["3,1"]*')')
p50<-ggplot(params) + 
  geom_histogram(aes(x=c3, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.15,0.45,0.012), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.15,0.45), breaks = c(0.15,0.3,0.45))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p50

#Contact Matrix: c41
summary(params$c4)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 30-64 ('*pi["4,1"]*')')
p51<-ggplot(params) + 
  geom_histogram(aes(x=c4, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,2.8,0.125), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,2.8), breaks = c(0,1.4,2.8))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p51

#Contact Matrix: c12
summary(params$c5)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-25 with 15-19 ('*pi["1,2"]*')')
p52<-ggplot(params) + 
  geom_histogram(aes(x=c5, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.2,1.6,0.06), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.2,1.6), breaks = c(0.2,0.9,1.6))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p52

#Contact Matrix: c22
summary(params$c6)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-25 with 20-25 ('*pi["2,2"]*')')
p53<-ggplot(params) + 
  geom_histogram(aes(x=c6, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.92,3.72,0.125), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.92,3.72), breaks = c(0.9,2.3,3.7))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p53

#Contact Matrix: c32
summary(params$c7)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-25 with 26-29 ('*pi["3,2"]*')')
p54<-ggplot(params) + 
  geom_histogram(aes(x=c7, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.45,1.75,0.06), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.45,1.75), breaks = c(0.45,1.10,1.75))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p54

#Contact Matrix: c42
summary(params$c8)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-25 with 30-64 ('*pi["4,2"]*')')
p55<-ggplot(params) + 
  geom_histogram(aes(x=c8, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,4.7,0.2), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,4.7), breaks = c(0,2.35,4.7))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p55

#Contact Matrix: c13
summary(params$c9)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 26-29 with 15-19 ('*pi["1,3"]*')')
p56<-ggplot(params) + 
  geom_histogram(aes(x=c9, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.12,0.48,0.015), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.12,0.48), breaks = c(0.12,0.3,0.48))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p56

#Contact Matrix: c23
summary(params$c10)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 26-29 with 20-25 ('*pi["2,3"]*')')
p57<-ggplot(params) + 
  geom_histogram(aes(x=c10, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.35,2.3,0.09), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.35,2.3), breaks = c(0.4,1.35,2.3))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p57

#Contact Matrix: c33
summary(params$c11)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 26-29 with 26-29 ('*pi["3,3"]*')')
p58<-ggplot(params) + 
  geom_histogram(aes(x=c11, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0.64,2.27,0.075), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.64,2.27), breaks = c(0.7,1.5,2.2))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p58

#Contact Matrix: c43
summary(params$c12)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 26-29 with 30-64 ('*pi["4,3"]*')')
p59<-ggplot(params) + 
  geom_histogram(aes(x=c12, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,5,0.2), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2.5,5))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p59

#Contact Matrix: c14
summary(params$c13)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 30-64 with 15-19 ('*pi["1,4"]*')')
p60<-ggplot(params) + 
  geom_histogram(aes(x=c13, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,3.1,0.15), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,3.1), breaks = c(0,1.5,3))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p60

#Contact Matrix: c24
summary(params$c14)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 30-64 with 20-25 ('*pi["2,4"]*')')
p61<-ggplot(params) + 
  geom_histogram(aes(x=c14, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,4.5,0.2), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,4.5), breaks = c(0,2.25,4.5))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p61

#Contact Matrix: c34
summary(params$c15)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 30-64 with 26-29 ('*pi["3,4"]*')')
p62<-ggplot(params) + 
  geom_histogram(aes(x=c15, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,6.4,0.3), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,6.4), breaks = c(0,3.2,6.4))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p62

#Contact Matrix: c44
summary(params$c16)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 30-64 with 30-64 ('*pi["4,4"]*')')
p63<-ggplot(params) + 
  geom_histogram(aes(x=c16, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(12,34,1), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(12,34), breaks = c(12,23,34))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p63

#Calculate proportions of contacts
params$pi11<-params$c1/(params$c1+params$c2+params$c3+params$c4)
params$pi21<-params$c2/(params$c1+params$c2+params$c3+params$c4)
params$pi31<-params$c3/(params$c1+params$c2+params$c3+params$c4)
params$pi41<-params$c4/(params$c1+params$c2+params$c3+params$c4)
params$pi12<-params$c5/(params$c5+params$c6+params$c7+params$c8)
params$pi22<-params$c6/(params$c5+params$c6+params$c7+params$c8)
params$pi32<-params$c7/(params$c5+params$c6+params$c7+params$c8)
params$pi42<-params$c8/(params$c5+params$c6+params$c7+params$c8)
params$pi13<-params$c9/(params$c9+params$c10+params$c11+params$c12)
params$pi23<-params$c10/(params$c9+params$c10+params$c11+params$c12)
params$pi33<-params$c11/(params$c9+params$c10+params$c11+params$c12)
params$pi43<-params$c12/(params$c9+params$c10+params$c11+params$c12)
params$pi14<-params$c13/(params$c13+params$c14+params$c15+params$c16)
params$pi24<-params$c14/(params$c13+params$c14+params$c15+params$c16)
params$pi34<-params$c15/(params$c13+params$c14+params$c15+params$c16)
params$pi44<-params$c16/(params$c13+params$c14+params$c15+params$c16)

#Contact Matrix: pi11
summary(params$pi11*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 15-19 with 15-19 ('*pi["1,1"]*'/'*sum(pi["j,1"], j==1, 4)*')')
p64<-ggplot(params) + 
  geom_histogram(aes(x=pi11*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(28,95,1.4), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(28,95), breaks = c(30,62.5,95))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p64

#Contact Matrix: pi21
summary(params$pi21*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 20-25 with 15-19 ('*pi["2,1"]*'/'*sum(pi["j,1"], j==1, 4)*')')
p65<-ggplot(params) + 
  geom_histogram(aes(x=pi21*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(2,37,0.5), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2,37), breaks = c(2,19.5,37))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p65

#Contact Matrix: pi31
summary(params$pi31*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 26-29 with 15-19 ('*pi["3,1"]*'/'*sum(pi["j,1"], j==1, 4)*')')
p66<-ggplot(params) + 
  geom_histogram(aes(x=pi31*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1.2,14.8,0.18), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.2,14.8), breaks = c(1.2,8,14.8))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p66

#Contact Matrix: pi41
summary(params$pi41*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 30-64 with 15-19 ('*pi["4,1"]*'/'*sum(pi["j,1"], j==1, 4)*')')
p67<-ggplot(params) + 
  geom_histogram(aes(x=pi41*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,56,1.2), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,56), breaks = c(0,28,56))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p67

#Contact Matrix: pi12
summary(params$pi12*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 15-19 with 20-25 ('*pi["1,2"]*'/'*sum(pi["j,2"], j==1, 4)*')')
p68<-ggplot(params) + 
  geom_histogram(aes(x=pi12*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(2,49,0.75), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2,49), breaks = c(2,25.5,49))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p68

#Contact Matrix: pi22
summary(params$pi22*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 20-25 with 20-25 ('*pi["2,2"]*'/'*sum(pi["j,2"], j==1, 4)*')')
p69<-ggplot(params) + 
  geom_histogram(aes(x=pi22*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(11,78,1.3), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(11,78), breaks = c(11,44.5,78))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p69

#Contact Matrix: pi32
summary(params$pi32*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 26-29 with 20-25 ('*pi["3,2"]*'/'*sum(pi["j,2"], j==1, 4)*')')
p70<-ggplot(params) + 
  geom_histogram(aes(x=pi32*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(5,53,0.75), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(5,53), breaks = c(5,29,53))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p70

#Contact Matrix: pi42
summary(params$pi42*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 30-64 with 20-25 ('*pi["4,2"]*'/'*sum(pi["j,2"], j==1, 4)*')')
p71<-ggplot(params) + 
  geom_histogram(aes(x=pi42*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,70,1.75), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,70), breaks = c(0,35,70))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p71

#Contact Matrix: pi13
summary(params$pi13*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 15-19 with 26-29 ('*pi["1,3"]*'/'*sum(pi["j,3"], j==1, 4)*')')
p72<-ggplot(params) + 
  geom_histogram(aes(x=pi13*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(1,27,0.27), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,27), breaks = c(1,14,27))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p72

#Contact Matrix: pi23
summary(params$pi23*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 20-25 with 26-29 ('*pi["2,3"]*'/'*sum(pi["j,3"], j==1, 4)*')')
p73<-ggplot(params) + 
  geom_histogram(aes(x=pi23*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(4.5,72,1.3), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(4.5,72), breaks = c(5,38.5,72))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p73

#Contact Matrix: pi33
summary(params$pi33*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 26-29 with 26-29 ('*pi["3,3"]*'/'*sum(pi["j,3"], j==1, 4)*')')
p74<-ggplot(params) + 
  geom_histogram(aes(x=pi33*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(8,79,1.1), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(8,79), breaks = c(8,43.5,79))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p74

#Contact Matrix: pi43
summary(params$pi43*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 30-64 with 26-29 ('*pi["4,3"]*'/'*sum(pi["j,3"], j==1, 4)*')')
p75<-ggplot(params) + 
  geom_histogram(aes(x=pi43*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,79,1.7), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,80), breaks = c(0,40,80))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p75

#Contact Matrix: pi14
summary(params$pi14*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 15-19 with 30-64 ('*pi["1,4"]*'/'*sum(pi["j,4"], j==1, 4)*')')
p76<-ggplot(params) + 
  geom_histogram(aes(x=pi14*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,18.1,0.44), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,18.1), breaks = c(0,9,18))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p76

#Contact Matrix: pi24
summary(params$pi24*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 20-25 with 30-64 ('*pi["2,4"]*'/'*sum(pi["j,4"], j==1, 4)*')')
p77<-ggplot(params) + 
  geom_histogram(aes(x=pi24*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,23.2,0.6), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,24), breaks = c(0,12,24))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p77


#Contact Matrix: pi34
summary(params$pi34*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 26-29 with 30-64 ('*pi["3,4"]*'/'*sum(pi["j,4"], j==1, 4)*')')
p78<-ggplot(params) + 
  geom_histogram(aes(x=pi34*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(0,30,0.8), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,30), breaks = c(0,15,30))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p78

#Contact Matrix: pi44
summary(params$pi44*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 30-64 with 30-64 ('*pi["4,4"]*'/'*sum(pi["j,4"], j==1, 4)*')')
p79<-ggplot(params) + 
  geom_histogram(aes(x=pi44*100, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', breaks=seq(50,100,0.9), show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(50,100), breaks = c(50,75,100))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p79

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
  scale_x_log10(breaks=c(200,2000,20000),labels=c("200","2,000","20,000"))+#limits = c(5,815), breaks = c(10,100))+
  scale_y_continuous(limits = c(0, 3000), breaks = c(0,1000,2000,3000), labels=c("0","1,000","2,000","3,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p84

p85<-ggplot(params) +
  geom_histogram(aes(x=RSS, fill = as.factor(reorder(fitorder,-fitorder))), colour='#000000', show.legend=TRUE, bins=40)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = model.colors, breaks=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))+
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_log10(breaks=c(200,2000,20000),labels=c("200","2,000","20,000"))+#limits = c(5,815), breaks = c(10,100))+
  scale_y_continuous(limits = c(0, 3000), breaks = c(0,1000,2000,3000), labels=c("0","1,000","2,000","3,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p85

#Saving main and supplemental figs

#Note: quartz works for saving unicode text as pdf (ggsave does not)

#All main params
w<-22; h<-12
lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12))
options(device = "quartz")
quartz()
LHS_Histograms1<-grid.arrange(p85,p3,p13,p37,p21,p22,p23,p24,
                              p4,p5,p6,p7,layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_Params1_12.pdf", "pdf",width=w,height=h)

options(device = "quartz")
quartz()
LHS_Histograms2<-grid.arrange(p14,p15,p16,p17,p38,p39,p40,p41,
                              p43,p44,p45,p46,layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_Params13_24.pdf", "pdf",width=w,height=h)

#with bottom legend
lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),
            c(13,13,13,13))
options(device = "quartz")
quartz()
LHS_Histograms3<-grid.arrange(p9,p10,p11,p12,p25,p26,p27,p28,
                              p19,p18,p20,p47,
                              bottomleg,
                              layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_Params25_36.pdf", "pdf",width=w,height=13.5)

options(device = "quartz")
quartz()
LHS_Histograms1<-grid.arrange(p85,p3,p13,p37,p21,p22,p23,p24,
                              p4,p5,p6,p7,bottomleg,layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_Params1_12_leg.pdf", "pdf",width=w,height=h)

options(device = "quartz")
quartz()
LHS_Histograms2<-grid.arrange(p14,p15,p16,p17,p38,p39,p40,p41,
                              p43,p44,p45,p46,bottomleg,layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_Params13_24_leg.pdf", "pdf",width=w,height=h)




#Contacts (sampled c from Polymod and sampled to 0)
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


#Proportional contacts (calculated from sampled c)
lay1<-rbind(c(1,1,2,2,3,3,4,4),c(1,1,2,2,3,3,4,4),c(1,1,2,2,3,3,4,4),
            c(5,5,6,6,7,7,8,8),c(5,5,6,6,7,7,8,8),c(5,5,6,6,7,7,8,8),
            c(9,9,10,10,11,11,12,12),c(9,9,10,10,11,11,12,12),c(9,9,10,10,11,11,12,12),
            c(13,13,14,14,15,15,16,16),c(13,13,14,14,15,15,16,16),c(13,13,14,14,15,15,16,16),
            c(17,17,17,17,17,17,17,17))
options(device = "quartz")
quartz()
LHS_Histograms5<-grid.arrange(p64,p68,p72,p76,
                              p65,p69,p73,p77,
                              p66,p70,p74,p78,
                              p67,p71,p75,p79,bottomleg,layout_matrix=lay1)
quartz.save ("LHS_Histograms_Grey_PropContacts_1.pdf", "pdf",width=w,height=h)

#RSS
w<-8; h<-6
ggsave(sprintf("LHS_RSS.pdf"), p84,width=w, height=h)



#### LHS Parameter Histograms (Black/Grey by best 10% vs rest) #####

#Set ylabel text
text.ylab <- "# Parameter Sets"


#Beta

#textsym <- expression(beta[1])
summary(params$b*100)
text.xlab <- bquote('% Infected'*~'Contact'^-1)
text.title <- bquote('Transmission Rate ('*beta*')')

#use to extract legend (bottom) 
p1<-ggplot(params) + 
  geom_histogram(aes(x=b*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 bins=25)+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  geom_histogram(data=subset(params,params$best10==1),aes(x=b*100, 
                fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 bins=25, na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 2000), breaks = c(0,500,1000,1500,2000), labels=c("0","500","1,000","1,500","2,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p1


bottomleg<-g_legend(p1)

#use to extract legend (side)
p2<-ggplot(params) + 
  geom_histogram(aes(x=b*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 bins=25)+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  geom_histogram(data=subset(params,params$best10==1),aes(x=b*100, 
                                                          fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 bins=25, na.rm=TRUE, show.legend=TRUE)+
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "right", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 2000), breaks = c(0,500,1000,1500,2000), labels=c("0","500","1,000","1,500","2,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p2
sideleg<-g_legend(p2)

#beta
p3<-ggplot(params) + 
  geom_histogram(aes(x=b*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 bins=25)+#breaks=seq(1.8e-07,2.3e-05,4e-07), 
  geom_histogram(data=subset(params,params$best10==1),aes(x=b*100, 
   fill = "Best 10% of Parameter Sets "), colour='#000000', 
   bins=25, na.rm=TRUE, show.legend=TRUE)+
                 theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 2000), breaks = c(0,500,1000,1500,2000), labels=c("0","500","1,000","1,500","2,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p3

#sigma1
summary(params$sigma1)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 15-19 ('*sigma[1]*')')

p4<-ggplot(params) + 
  geom_histogram(aes(x=sigma1, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(5,20,0.75))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=sigma1,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(5,20,0.75), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  #scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p4

#sigma2
summary(params$sigma2)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 20-25 ('*sigma[2]*')')

p5<-ggplot(params) + 
  geom_histogram(aes(x=sigma2, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(5,20,0.75))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=sigma2,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(5,20,0.75), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  #scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p5

#sigma3
summary(params$sigma3)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 26-29 ('*sigma[3]*')')

p6<-ggplot(params) + 
  geom_histogram(aes(x=sigma3, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(5,20,0.75))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=sigma3,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(5,20,0.75), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  #scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p6

#sigma4
summary(params$sigma4)

text.xlab <- bquote('Contacts')
text.title <- bquote('Total Contacts: 30-64 ('*sigma[4]*')')

p7<-ggplot(params) + 
  geom_histogram(aes(x=sigma4, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(5,20,0.75))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=sigma4,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(5,20,0.75), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  #scale_x_log10(breaks = c(0.0001,0.01,1),labels=c("0.0001","0.01","1"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p7

#gamma1
summary(params$gamma1)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Cessation Rate: 15-19 ('*gamma[1]*')')

p9<-ggplot(params) + 
  geom_histogram(aes(x=gamma1, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.05,1.2,0.05))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=gamma1,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.05,1.2,0.05), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=gamma2, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.05,1.2,0.05))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=gamma2,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.05,1.2,0.05), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=gamma3, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.05,1.2,0.05))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=gamma3,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.05,1.2,0.05), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=gamma4, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.05,1.2,0.05))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=gamma4,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.05,1.2,0.05), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=d*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(15,50,1.5))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=d*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(15,50,1.5), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=zeta1*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.15,1.0,0.04))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=zeta1*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.15,1.0,0.04), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.15,1.0), breaks = c(0.15,0.575,1))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p14

#zeta2
summary(params$zeta2)

text.xlab <- bquote('%')
text.title <- bquote('Current PWID Prevalence: 20-25 ('*zeta[2]*')')

p15<-ggplot(params) + 
  geom_histogram(aes(x=zeta2*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.26,1.0,0.035))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=zeta2*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.26,1.0,0.035), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=zeta3*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.19,1,0.04))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=zeta3*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.19,1,0.04), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=zeta4*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.16,1.1,0.045))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=zeta4*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.16,1.1,0.045), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=etan, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.18,0.54,0.016))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=etan,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.18,0.54,0.016), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=etap, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(2.5,15.3,0.6))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=etap,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(2.5,15.3,0.6), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=etaz, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(4.3,4.43,0.006))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=etaz,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(4.3,4.43,0.006), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=theta1, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.05,5,0.4))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=theta1,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.05,5,0.4), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2.5,5))+
  scale_y_continuous(limits = c(0, 5000), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p21


#theta2
summary(params$theta2)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 20-25 ('*theta[2]*')')


p22<-ggplot(params) + 
  geom_histogram(aes(x=theta2, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.034,5,0.35))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=theta2,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.034,5,0.35), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2.5,5))+
  scale_y_continuous(limits = c(0, 5000), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p22

#theta3
summary(params$theta3)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 26-29 ('*theta[3]*')')


p23<-ggplot(params) + 
  geom_histogram(aes(x=theta3, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.034,10,0.36))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=theta3,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.034,10,0.36), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,10), breaks = c(0,5,10))+
  scale_y_continuous(limits = c(0, 5000), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p23

#theta4
summary(params$theta4)

text.xlab <- bquote(Years^-1)
text.title <- bquote('Injection Initiation: 30-64 ('*theta[4]*')')


p24<-ggplot(params) + 
  geom_histogram(aes(x=theta4, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.034,1.0,0.045))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=theta4,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.034,1.0,0.045), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=k1, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,2.5,0.15))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=k1,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,2.5,0.15), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=k2, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,1.66,0.1))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=k2,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,1.66,0.1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=k3, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,1.43,0.1))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=k3,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,1.43,0.1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=k4, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,0.8,0.05))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=k4,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,0.8,0.05), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,0.8), breaks = c(0,0.4,0.8))+
  scale_y_continuous(limits = c(0, 4200), breaks = c(0,2000,4000), labels=c("0","2,000","4,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p28


#xi
summary(params$xi*100)

text.xlab <- bquote("%")
text.title <- bquote('HCV Immunity to Reinfection ('*xi*')')

p37<-ggplot(params) + 
  geom_histogram(aes(x=xi*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.003,45,2))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=xi*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.003,45,2), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,45), breaks = c(0,22.5,45))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p37

#omic1
summary(params$omic1*100)

text.xlab <- bquote("%")
text.title <- bquote('Former PWID Prevalence: 15-19 ('*omicron[1]*')')

p38<-ggplot(params) + 
  geom_histogram(aes(x=omic1*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.43,1.0,0.028))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=omic1*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.43,1.0,0.028), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=omic2*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.71,1,0.014))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=omic2*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.71,1,0.014), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=omic3*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.81,2.1,0.06))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=omic3*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.81,2.1,0.06), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=omic4*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(1.6,2.7,0.055))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=omic4*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(1.6,2.7,0.055), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.6,2.7), breaks = c(1.6,2.15,2.7))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p41

#psi1
summary(params$psi1*100)

text.xlab <- bquote('%')
text.title <- bquote('Abuse/Dependence Prevalence: 15-19 ('*psi[1]*')')

p43<-ggplot(params) + 
  geom_histogram(aes(x=psi1*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(1.6,2.1,0.025))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=psi1*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(1.6,2.1,0.025), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=psi2*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(1.1,1.5,0.02))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=psi2*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(1.1,1.5,0.02), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=psi3*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.37,1,0.03))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=psi3*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.37,1,0.03), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=psi4*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.24,0.48,0.012))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=psi4*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.24,0.48,0.012), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
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
  geom_histogram(aes(x=psin*Z0, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(1167,4323,150))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=psin*Z0,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(1167,4323,150), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1167,4323), breaks = c(1200,2750,4300), labels=c("1,200","2,750","4,300"))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p47

#Contact Matrix: c11
summary(params$c1)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 15-19 ('*pi["1,1"]*')')
p48<-ggplot(params) + 
  geom_histogram(aes(x=c1, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(1.5,9.5,0.35))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c1,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(1.5,9.5,0.35), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.5,9.5), breaks = c(1.5,5.5,9.5))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p48

#Contact Matrix: c21
summary(params$c2)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 20-25 ('*pi["2,1"]*')')
p49<-ggplot(params) + 
  geom_histogram(aes(x=c2, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.3,1.2,0.04))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c2,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.3,1.2,0.04), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.3,1.2), breaks = c(0.3,0.75,1.2))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p49

#Contact Matrix: c31
summary(params$c3)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 26-29 ('*pi["3,1"]*')')
p50<-ggplot(params) + 
  geom_histogram(aes(x=c3, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.15,0.45,0.012))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c3,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.15,0.45,0.012), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.15,0.45), breaks = c(0.15,0.3,0.45))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p50

#Contact Matrix: c41
summary(params$c4)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 15-19 with 30-64 ('*pi["4,1"]*')')
p51<-ggplot(params) + 
  geom_histogram(aes(x=c4, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,2.8,0.125))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c4,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,2.8,0.125), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,2.8), breaks = c(0,1.4,2.8))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p51

#Contact Matrix: c12
summary(params$c5)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-25 with 15-19 ('*pi["1,2"]*')')
p52<-ggplot(params) + 
  geom_histogram(aes(x=c5, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.2,1.6,0.06))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c5,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.2,1.6,0.06), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.2,1.6), breaks = c(0.2,0.9,1.6))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p52

#Contact Matrix: c22
summary(params$c6)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-25 with 20-25 ('*pi["2,2"]*')')
p53<-ggplot(params) + 
  geom_histogram(aes(x=c6, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.92,3.72,0.125))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c6,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.92,3.72,0.125), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.92,3.72), breaks = c(0.9,2.3,3.7))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p53

#Contact Matrix: c32
summary(params$c7)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-25 with 26-29 ('*pi["3,2"]*')')
p54<-ggplot(params) + 
  geom_histogram(aes(x=c7, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.45,1.75,0.06))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c7,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.45,1.75,0.06), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.45,1.75), breaks = c(0.45,1.1,1.75))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p54

#Contact Matrix: c42
summary(params$c8)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 20-25 with 30-64 ('*pi["4,2"]*')')
p55<-ggplot(params) + 
  geom_histogram(aes(x=c8, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,4.7,0.2))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c8,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,4.7,0.2), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,4.7), breaks = c(0,2.35,4.7))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p55

#Contact Matrix: c13
summary(params$c9)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 26-29 with 15-19 ('*pi["1,3"]*')')
p56<-ggplot(params) + 
  geom_histogram(aes(x=c9, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.12,0.48,0.015))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c9,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.12,0.48,0.015), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.12,0.48), breaks = c(0.12,0.3,0.48))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p56

#Contact Matrix: c23
summary(params$c10)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 26-29 with 20-25 ('*pi["2,3"]*')')
p57<-ggplot(params) + 
  geom_histogram(aes(x=c10, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.35,2.3,0.09))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c10,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.35,2.3,0.09), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.35,2.3), breaks = c(0.4,1.35,2.3))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p57

#Contact Matrix: c33
summary(params$c11)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 26-29 with 26-29 ('*pi["3,3"]*')')
p58<-ggplot(params) + 
  geom_histogram(aes(x=c11, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0.64,2.27,0.075))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c11,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0.64,2.27,0.075), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0.64,2.27), breaks = c(0.7,1.5,2.2))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p58

#Contact Matrix: c43
summary(params$c12)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 26-29 with 30-64 ('*pi["4,3"]*')')
p59<-ggplot(params) + 
  geom_histogram(aes(x=c12, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,5,0.2))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c12,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,5,0.2), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2.5,5))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p59

#Contact Matrix: c14
summary(params$c13)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 30-64 with 15-19 ('*pi["1,4"]*')')
p60<-ggplot(params) + 
  geom_histogram(aes(x=c13, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,3.1,0.15))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c13,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,3.1,0.15), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,3.1), breaks = c(0,1.5,3))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p60

#Contact Matrix: c24
summary(params$c14)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 30-64 with 20-25 ('*pi["2,4"]*')')
p61<-ggplot(params) + 
  geom_histogram(aes(x=c14, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,4.5,0.2))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c14,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,4.5,0.2), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,4.5), breaks = c(0,2.25,4.5))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p61

#Contact Matrix: c34
summary(params$c15)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 30-64 with 26-29 ('*pi["3,4"]*')')
p62<-ggplot(params) + 
  geom_histogram(aes(x=c15, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,6.4,0.3))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c15,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,6.4,0.3), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,6.4), breaks = c(0,3.2,6.4))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p62

#Contact Matrix: c44
summary(params$c16)
text.xlab <- bquote(Contacts*~Persons^-1*~Days^-1)
text.title <- bquote('Contacts: 30-64 with 30-64 ('*pi["4,4"]*')')
p63<-ggplot(params) + 
  geom_histogram(aes(x=c16, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(12,34,1))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=c16,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(12,34,1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(12,34), breaks = c(12,23,34))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p63

#Calculate proportions of contacts
params$pi11<-params$c1/(params$c1+params$c2+params$c3+params$c4)
params$pi21<-params$c2/(params$c1+params$c2+params$c3+params$c4)
params$pi31<-params$c3/(params$c1+params$c2+params$c3+params$c4)
params$pi41<-params$c4/(params$c1+params$c2+params$c3+params$c4)
params$pi12<-params$c5/(params$c5+params$c6+params$c7+params$c8)
params$pi22<-params$c6/(params$c5+params$c6+params$c7+params$c8)
params$pi32<-params$c7/(params$c5+params$c6+params$c7+params$c8)
params$pi42<-params$c8/(params$c5+params$c6+params$c7+params$c8)
params$pi13<-params$c9/(params$c9+params$c10+params$c11+params$c12)
params$pi23<-params$c10/(params$c9+params$c10+params$c11+params$c12)
params$pi33<-params$c11/(params$c9+params$c10+params$c11+params$c12)
params$pi43<-params$c12/(params$c9+params$c10+params$c11+params$c12)
params$pi14<-params$c13/(params$c13+params$c14+params$c15+params$c16)
params$pi24<-params$c14/(params$c13+params$c14+params$c15+params$c16)
params$pi34<-params$c15/(params$c13+params$c14+params$c15+params$c16)
params$pi44<-params$c16/(params$c13+params$c14+params$c15+params$c16)

#Contact Matrix: pi11
summary(params$pi11*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 15-19 with 15-19 ('*pi["1,1"]*')')
p64<-ggplot(params) + 
  geom_histogram(aes(x=pi11*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(28,95,1.4))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi11*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(28,95,1.4), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(28,95), breaks = c(30,62.5,95))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p64

#Contact Matrix: pi21
summary(params$pi21*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 20-25 with 15-19 ('*pi["2,1"]*')')
p65<-ggplot(params) + 
  geom_histogram(aes(x=pi21*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(2,37,0.5))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi21*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(2,37,0.5), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2,37), breaks = c(2,19.5,37))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p65

#Contact Matrix: pi31
summary(params$pi31*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 26-29 with 15-19 ('*pi["3,1"]*')')
p66<-ggplot(params) + 
  geom_histogram(aes(x=pi31*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(1.2,14.8,0.18))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi31*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(1.2,14.8,0.18), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1.2,14.8), breaks = c(1.2,8,14.8))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p66

#Contact Matrix: pi41
summary(params$pi41*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 30-64 with 15-19 ('*pi["4,1"]*')')
p67<-ggplot(params) + 
  geom_histogram(aes(x=pi41*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,56,1.2))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi41*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,56,1.2), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,56), breaks = c(0,28,56))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p67

#Contact Matrix: pi12
summary(params$pi12*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 15-19 with 20-25 ('*pi["1,2"]*')')
p68<-ggplot(params) + 
  geom_histogram(aes(x=pi12*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(2,49,0.75))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi12*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(2,49,0.75), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(2,49), breaks = c(2,25.5,49))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p68

#Contact Matrix: pi22
summary(params$pi22*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 20-25 with 20-25 ('*pi["2,2"]*')')
p69<-ggplot(params) + 
  geom_histogram(aes(x=pi22*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(11,78,1.3))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi22*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(11,78,1.3), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(11,78), breaks = c(11,44.5,78))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p69

#Contact Matrix: pi32
summary(params$pi32*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 26-29 with 20-25 ('*pi["3,2"]*')')
p70<-ggplot(params) + 
  geom_histogram(aes(x=pi32*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(5,53,0.75))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi32*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(5,53,0.75), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(5,53), breaks = c(5,29,53))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p70

#Contact Matrix: pi42
summary(params$pi42*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 30-64 with 20-25 ('*pi["4,2"]*')')
p71<-ggplot(params) + 
  geom_histogram(aes(x=pi42*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,70,1.75))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi42*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,70,1.75), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,70), breaks = c(0,35,70))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p71

#Contact Matrix: pi13
summary(params$pi13*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 15-19 with 26-29 ('*pi["1,3"]*')')
p72<-ggplot(params) + 
  geom_histogram(aes(x=pi13*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(1,27,0.27))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi13*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(1,27,0.27), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(1,27), breaks = c(1,14,27))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p72

#Contact Matrix: pi23
summary(params$pi23*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 20-25 with 26-29 ('*pi["2,3"]*')')
p73<-ggplot(params) + 
  geom_histogram(aes(x=pi23*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(4.5,72,1.3))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi23*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(4.5,72,1.3), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(4.5,72), breaks = c(5,38.5,72))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p73

#Contact Matrix: pi33
summary(params$pi33*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 26-29 with 26-29 ('*pi["3,3"]*')')
p74<-ggplot(params) + 
  geom_histogram(aes(x=pi33*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(8,79,1.1))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi33*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(8,79,1.1), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(8,79), breaks = c(8,43.5,79))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p74

#Contact Matrix: pi43
summary(params$pi43*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 30-64 with 26-29 ('*pi["4,3"]*')')
p75<-ggplot(params) + 
  geom_histogram(aes(x=pi43*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,79,1.7))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi43*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,79,1.7), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,80), breaks = c(0,40,80))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p75

#Contact Matrix: pi14
summary(params$pi14*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 15-19 with 30-64 ('*pi["1,4"]*')')
p76<-ggplot(params) + 
  geom_histogram(aes(x=pi14*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,18.1,0.44))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi14*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,18.1,0.44), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,18.1), breaks = c(0,9,18))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p76

#Contact Matrix: pi24
summary(params$pi24*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 20-25 with 30-64 ('*pi["2,4"]*')')
p77<-ggplot(params) + 
  geom_histogram(aes(x=pi24*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,23.2,0.6))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi24*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,23.2,0.6), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,24), breaks = c(0,12,24))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p77


#Contact Matrix: pi34
summary(params$pi34*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 26-29 with 30-64 ('*pi["3,4"]*')')
p78<-ggplot(params) + 
  geom_histogram(aes(x=pi34*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(0,30,0.8))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi34*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(0,30,0.8), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(0,30), breaks = c(0,15,30))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p78

#Contact Matrix: pi44
summary(params$pi44*100)
text.xlab <- bquote("%")
text.title <- bquote('Contacts: 30-64 with 30-64 ('*pi["4,4"]*')')
p79<-ggplot(params) + 
  geom_histogram(aes(x=pi44*100, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 breaks=seq(50,100,0.9))+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=pi44*100,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 breaks=seq(50,100,0.9), na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab("") +
  xlab(text.xlab) +
  scale_x_continuous(limits = c(50,100), breaks = c(50,75,100))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0,250,500), labels=c("0","250","500"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p79

# #RSS 
summary(params$RSS)
text.xlab <- "Residual"
text.title<-"Residual Sum of Squares"
p84<-ggplot(params) +
  geom_histogram(aes(x=RSS, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 bins=40)+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=RSS,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 bins=40, na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_log10(breaks=c(200,2000,20000),labels=c("200","2,000","20,000"))+#limits = c(5,815), breaks = c(10,100))+
  scale_y_continuous(limits = c(0, 3000), breaks = c(0,1000,2000,3000), labels=c("0","1,000","2,000","3,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p84

p85<-ggplot(params) +
  geom_histogram(aes(x=RSS, fill = "All Parameter Sets "), colour='#000000',show.legend=TRUE,
                 bins=40)+ 
  geom_histogram(data=subset(params,params$best10==1),aes(x=RSS,fill = "Best 10% of Parameter Sets "), colour='#000000', 
                 bins=40, na.rm=TRUE, show.legend=TRUE)+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  scale_fill_manual(name=NULL, values = c("#000000","#CCCCCC"), breaks = c("All Parameter Sets ","Best 10% of Parameter Sets ")) +
  ylab(text.ylab) +
  xlab(text.xlab) +
  scale_x_log10(breaks=c(200,2000,20000),labels=c("200","2,000","20,000"))+#limits = c(5,815), breaks = c(10,100))+
  scale_y_continuous(limits = c(0, 3000), breaks = c(0,1000,2000,3000), labels=c("0","1,000","2,000","3,000"))+
  ggtitle(text.title)+
  theme(plot.title = element_text(hjust = 0.5))
p85

#start here to save param figs

#Saving main and supplemental figs

#Note: quartz works for saving unicode text as pdf (ggsave does not)

#All main params
w<-22; h<-12
lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12))
options(device = "quartz")
LHS_Histograms1<-grid.arrange(p85,p3,p13,p37,p21,p22,p23,p24,
                              p4,p5,p6,p7,layout_matrix=lay1)
quartz.save ("LHS10_Histograms_Grey_Params1_12.pdf", "pdf",width=w,height=h)

LHS_Histograms2<-grid.arrange(p14,p15,p16,p17,p38,p39,p40,p41,
                              p43,p44,p45,p46,layout_matrix=lay1)
quartz.save ("LHS10_Histograms_Grey_Params13_24.pdf", "pdf",width=w,height=h)

lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),
            c(13,13,13,13))
LHS_Histograms3<-grid.arrange(p9,p10,p11,p12,p25,p26,p27,p28,
                              p19,p18,p20,p47,
                              bottomleg,
                              layout_matrix=lay1)
quartz.save ("LHS10_Histograms_Grey_Params25_36.pdf", "pdf",width=w,height=13.5)

#Contacts (sampled c from Polymod and sampled to 0)
lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),
            c(13,14,15,16),c(13,14,15,16),c(13,14,15,16),c(17,17,17,17))
LHS_Histograms4<-grid.arrange(p48,p52,p56,p60,
                              p49,p53,p57,p61,
                              p50,p54,p58,p62,
                              p51,p55,p59,p63,bottomleg,layout_matrix=lay1)
quartz.save ("LHS10_Histograms_Grey_SampledContacts.pdf", "pdf",width=w,height=h)


#Proportional contacts (calculated from sampled c)
lay1<-rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(5,6,7,8),c(5,6,7,8),c(5,6,7,8),
            c(9,10,11,12),c(9,10,11,12),c(9,10,11,12),c(13,14,15,16),c(13,14,15,16),c(13,14,15,16),c(17,17,17,17))
LHS_Histograms5<-grid.arrange(p64,p68,p72,p76,
                              p65,p69,p73,p77,
                              p66,p70,p74,p78,
                              p67,p71,p75,p79,bottomleg,layout_matrix=lay1)
quartz.save ("LHS10_Histograms_Grey_PropContacts.pdf", "pdf",width=w,height=h)

#RSS
w<-8; h<-6
ggsave(sprintf("LHS10_RSS.pdf"), p84,width=w, height=h)




################ Section 5 Figures: Intervention Simulation Results ################
#### Violin Plots of Best 10% Fits of Single Interventions ####
#Set plot labels and colors
text.xlab <- "Intervention"
#model.colors4 <- c('#CCCCCC', '#999999', '#666666', '#000000','#FFFFFF') #grey scheme
#model.colors <- c('#CCCCCC', '#666666', '#000000','#FFFFFF') #grey scheme
#model.colors <- c('#FF9999', '#CC3333','#990000','#FFFFFF') #red scheme
model.colors <- c('#99CCFF','#3366FF', '#000066','#FFFFFF') #blue scheme

table(singleint$order)
table(singleint$IntLevelpct)

summary(singleint$fitquartile)


#create new variable to force none to plot 1st in legend as white color (as actual plots do below)
singleint$plotleg<-ifelse(singleint$IntLevelpct=="None",1,ifelse(singleint$IntLevelpct=="10%",2,
                                                                 ifelse(singleint$IntLevelpct=="20%",3,4)))
singleint$type<-ifelse(singleint$order==1,"None",ifelse(singleint$order==2,"Primary",
                                                        ifelse(singleint$order %in% c(3,4,5),"Secondary",
                                                               ifelse(singleint$order %in% c(6,7,8),"Tertiary","NA"))))
  
#X labels
xlabelsmain <-paste0(c("None","Decrease Injection Initiation",
                       "Decrease Current PWID Contacts",
                       "Increase Current PWID Cessation", 
                       "Decrease Former PWID Relapse",
                       "Treat Current PWID",
                       "Treat Former PWID","Treat Current and Former PWID"))


text.ylab <- "% Reduction Chronic HCV"

leg<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                     singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%")),
            aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
 
  #scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=c("","","","","","",""))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=15) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
leg

legend1<-g_legend(leg)

allchronic<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                     singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%")),
            aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  
  #scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", 
               position=position_dodge(width=0.8)) +
  theme_bw(base_size=15) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
allchronic

allacute<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                            singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%")),
                   aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  
  #scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", 
               position=position_dodge(width=0.8)) +
  theme_bw(base_size=15) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("% Reduction Acute HCV") +
  xlab("")
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
allacute


#None, Prevalence
text.ylab <- "% Reduction Chronic HCV"
p1<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                  singleint$order %in% c(1)&singleint$IntLevelpct %in% c("None")),
                  aes(x=order, y=pctred_chr, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values="#FFFFFF", name='Intervention Level') +
  scale_x_discrete(labels="None")+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")+
  ggtitle("None")+
  theme(plot.title = element_text(hjust = 0.5))
p1


#Primary, Prevalence
text.ylab <- "% Reduction Chronic HCV"
p2<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(2)&singleint$IntLevelpct %in% c("10%","20%","40%")),
           aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap("Decrease Injection Initiation",18))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")+
  ggtitle("Primary")+
  theme(plot.title = element_text(hjust = 0.5))
p2
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

# title.grob <- textGrob(label = "  ",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = -0.5,gp = gpar(fontsize = 16))
# 
# p2 <- arrangeGrob(p2, top = title.grob)
# grid.arrange(p2)

#Secondary, Prevalence
p3<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(3:5)&singleint$IntLevelpct %in% c("10%","20%","40%")),
           aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse"),15))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("") +
  xlab("")+
  ggtitle("Secondary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p3

# p3 <- arrangeGrob(p3, top = title.grob)
# grid.arrange(p3)

#Tertiary, Prevalence
p4<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(6:8)&singleint$IntLevelpct %in% c("10%","20%","40%")),
           aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Treat Current PWID","Treat Former PWID","Treat Current & Former"),15))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("") +
  xlab("")+
  ggtitle("Tertiary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p4

# p4 <- arrangeGrob(p4, top = title.grob)
# grid.arrange(p4)

#None, Acute
text.ylab <- "% Reduction Acute HCV"
p5<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(1)&singleint$IntLevelpct %in% c("None")),
           aes(x=order, y=pctred_acute, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values="#FFFFFF", name='Intervention Level') +
  scale_x_discrete(labels="None")+ 
    scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
    stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") +
  ggtitle("None")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p5

# title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = -0.5,gp = gpar(fontsize = 16))
# 
# p5 <- arrangeGrob(p5, top = title.grob)
# grid.arrange(p5)


#Primary, Acute
text.ylab <- "% Reduction Acute HCV"
p6<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(2)&singleint$IntLevelpct %in% c("10%","20%","30%")),
           aes(x=order, y=pctred_acute, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap("Decrease Injection Initiation",18))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") +
  ggtitle("Primary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p6

# title.grob <- textGrob(label = "  ",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = -0.5,gp = gpar(fontsize = 16))
# 
# p6 <- arrangeGrob(p6, top = title.grob)
# grid.arrange(p6)

#Secondary, Acute
p7<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(3:5)&singleint$IntLevelpct %in% c("10%","20%","30%")),
           aes(x=order, y=pctred_acute, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse"),15))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("") +
  xlab("") +
  ggtitle("Secondary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p7

# p7 <- arrangeGrob(p7, top = title.grob)
# grid.arrange(p7)

#Tertiary, Acute
p8<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(6:8)&singleint$IntLevelpct %in% c("10%","20%","30%")),
           aes(x=order, y=pctred_acute, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Treat Current PWID","Treat Former PWID","Treat Current & Former"),15))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("") +
  xlab("") +
  ggtitle("Tertiary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p8

#Tertiary, Acute, Trt Both
p9<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(8)&singleint$IntLevelpct %in% c("10%","20%","30%")),
           aes(x=order, y=pctred_acute, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Treat Current & Former PWID"),15))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("") +
  xlab("") +
  ggtitle("Tertiary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p9

#Tertiary, Chronic, Trt Both
p10<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(8)&singleint$IntLevelpct %in% c("10%","20%","30%")),
           aes(x=order, y=pctred_chr, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Treat Current & Former PWID"),15))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("") +
  xlab("") +
  ggtitle("Tertiary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p10

#Tertiary, Acute, Trt Current or Former
p11<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(6,7)&singleint$IntLevelpct %in% c("10%","20%","30%")),
           aes(x=order, y=pctred_acute, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Treat Current PWID","Treat Former PWID"),15))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("") +
  xlab("") +
  ggtitle("Tertiary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p11

#Tertiary, Chronic, Trt Both
p12<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                     singleint$order %in% c(6,7)&singleint$IntLevelpct %in% c("10%","20%","30%")),
            aes(x=order, y=pctred_chr, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Treat Current PWID","Treat Former PWID"),15))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("") +
  xlab("") +
  ggtitle("Tertiary")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p12

# p8 <- arrangeGrob(p8, top = title.grob)
# grid.arrange(p8)
# 
# mylegend<-g_legend(leg)

#Main Violin Plot Figure: No Title
w <- 13.5; h <- 7

lay<-rbind(c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3), 
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3),  
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3), 
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6),
           c(7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7))
SingleStrategy_VP_main<-grid.arrange(p2,p3,p4,p6,p7,p8,legend1,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_best10_grey.pdf"), SingleStrategy_VP_main, width=w, height=h)

w <- 12; h <- 7
lay<-rbind(c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3), 
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3),  
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3), 
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6),
           c(7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7))
SingleStrategy_VP_main<-grid.arrange(p2,p3,p11,p6,p7,p12,legend1,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_best10_1_grey.pdf"), SingleStrategy_VP_main, width=w, height=h)

w <- 11; h <- 7
lay<-rbind(c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3), 
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3),  
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3), 
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6),
           c(7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7))
SingleStrategy_VP_main<-grid.arrange(p2,p3,p9,p6,p7,p10,legend1,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_best10_2_grey.pdf"), SingleStrategy_VP_main, width=w, height=h)

w <- 12; h <- 7
lay<-rbind(c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3), 
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3),  
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3), 
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6),
           c(7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7))
SingleStrategy_VP_main<-grid.arrange(p2,p3,p11,p6,p7,p12,legend1,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_best10_3_grey.pdf"), SingleStrategy_VP_main, width=w, height=h)

w <- 11; h <- 7
lay<-rbind(c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3), 
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3),  
           c(1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3), 
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6),
           c(4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6),
           c(7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7))
SingleStrategy_VP_main<-grid.arrange(p2,p3,p9,p6,p7,p10,legend1,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_best10_4_grey.pdf"), SingleStrategy_VP_main, width=w, height=h)

w <- 10.5; h <- 7
lay<-rbind(c(1), 
           c(1),  
           c(1),c(1),c(1), 
           c(2),
           c(2),
           c(2),c(2),c(2),
           c(3))
SingleStrategy_VP_main<-grid.arrange(allchronic,allacute,legend1,layout_matrix=lay)
ggsave(sprintf("SingleStrategy_best10_5_grey.pdf"), SingleStrategy_VP_main, width=w, height=h)

#### Violin Plots Primary & Tertiary (Prevalence Only, for ppt) ####

#Tertiary, Prevalence
p4<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                    singleint$order %in% c(2,6:7)&singleint$IntLevelpct %in% c("10%","20%")),
           aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#3366FF','#000066'), name='Intervention Level') +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation","Treat Current PWID for HCV","Treat Former PWID for HCV"),20))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("% Reduction in HCV Prevalence") +
  xlab("")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
p4


#### Violin Plots of Treatment Variations (Treat Current+Former PWID) % Red Best 10% ####
#X labels
# xlabelssupp <-paste0(c("None",
#                        "Sustained Virologic Response","Decrease Non-PWID Mortality",
#                        "Decrease Both"))

singleint$ordertau<-ifelse(singleint$taupct=="0%",1,
                          ifelse(singleint$taupct=="25%",2,
                                 ifelse(singleint$taupct=="50%",3,
                                        ifelse(singleint$taupct=="75%",4,
                                               ifelse(singleint$taupct=="100%",5,9999)))))


table(subset(singleint,singleint$order %in% c(8,9999)&
               singleint$InterventionParam %in% c("phi_both","phiboth_tau"))$taupct)

table(singleint$taupct)

text.ylab <- "% Reduction HCV Prevalence"
leg<-ggplot(subset(singleint,singleint$best10==1&
                     singleint$order %in% c(8,9999)&singleint$ODint==0&
                     singleint$InterventionParam %in% c("phi_both","phiboth_tau")&
                     singleint$IntLevelpct %in% c("10%","20%","40%")), 
            aes(x=as.factor(ordertau), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("0%","25%","50%","75%","100%"))+ 
  #scale_y_continuous(limits=c(9500,21000),breaks=c(10000,15000,20000),labels=c("10,000","15,000","20,000"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") 
leg

leg<-g_legend(leg)

library(stringr)

text.ylab <- "% Reduction Chronic HCV"
p1<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(8,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phi_both","phiboth_tau")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(ordertau), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("0%","25%","50%","75%","100%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(35,100),breaks=c(40,60,80,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Transmission During Treatment")+
theme(plot.title = element_text(hjust = 0.5))
p1


text.ylab <- "% Reduction Acute HCV"
p2<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(8,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phi_both","phiboth_tau")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(ordertau), y=pctred_acute, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("0%","25%","50%","75%","100%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(16,100),breaks=c(20,60,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Transmission During Treatment")+
theme(plot.title = element_text(hjust = 0.5)) 
p2


table(subset(singleint,singleint$order %in% c(8,9999)&
               singleint$InterventionParam %in% c("phi_both","phiboth_alpha"))$alphapct)

table(singleint$alphapct)

singleint$orderalpha<-ifelse(singleint$alphapct=="100%",1,
                             ifelse(singleint$alphapct=="90%",2,
                                    ifelse(singleint$alphapct=="80%",3,
                                           ifelse(singleint$alphapct=="70%",4,
                                                  ifelse(singleint$alphapct=="60%",5,9999)))))


text.ylab <- "% Reduction Chronic HCV"
p3<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(8,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phi_both","phiboth_alpha")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderalpha), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("100%","90%","80%","70%","60%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(35,100),breaks=c(40,60,80,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Sustained Virologic Response")+
theme(plot.title = element_text(hjust = 0.5))
p3

text.ylab <- "% Reduction Acute HCV"
p4<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(8,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phi_both","phiboth_alpha")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderalpha), y=pctred_acute, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("100%","90%","80%","70%","60%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(16,100),breaks=c(20,60,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Sustained Virologic Response")+
  theme(plot.title = element_text(hjust = 0.5))
p4

table(subset(singleint,singleint$order %in% c(8,9999)&
               singleint$InterventionParam %in% c("phi_both","phiboth_omega"))$omegapct)

table(singleint$omegapct)

singleint$orderomega<-ifelse(singleint$omegapct=="8 Weeks",1,
                             ifelse(singleint$omegapct=="12 Weeks",2,
                                    ifelse(singleint$omegapct=="16 Weeks",3,9999)))


text.ylab <- "% Reduction Chronic HCV"
p5<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(8,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phi_both","phiboth_omega")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderomega), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("8","12","16"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(35,100),breaks=c(40,60,80,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("Weeks")+
  ggtitle("Treatment Duration")+
  theme(plot.title = element_text(hjust = 0.5))
p5

text.ylab <- "% Reduction Acute HCV"
p6<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(8,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phi_both","phiboth_omega")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderomega), y=pctred_acute, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("8","12","16"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(16,100),breaks=c(20,60,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("Weeks")+
  ggtitle("Treatment Duration")+
  theme(plot.title = element_text(hjust = 0.5))
p6

#No title
w <- 12; h <- 12
lay<-rbind(c(1,2),c(1,2),c(1,2),c(3,4),c(3,4),c(3,4),c(5,6),c(5,6),c(5,6),c(7,7))
SingleStrategy_VP_supp<-grid.arrange(p1,p2,p3,p4,p5,p6,leg,layout_matrix=lay)
ggsave(sprintf("TrtSensVP_Both_Best10.pdf"), SingleStrategy_VP_supp, width=w, height=h)




#### Violin Plots of Treatment Variations (Treat Current PWID) % Red Best 10% ####


singleint$ordertau<-ifelse(singleint$taupct=="0%",1,
                           ifelse(singleint$taupct=="25%",2,
                                  ifelse(singleint$taupct=="50%",3,
                                         ifelse(singleint$taupct=="75%",4,
                                                ifelse(singleint$taupct=="100%",5,9999)))))


table(subset(singleint,singleint$order %in% c(6,9999)&
               singleint$InterventionParam %in% c("phip","phip_tau"))$taupct)

table(singleint$taupct)

text.ylab <- "% Reduction HCV Prevalence"
leg<-ggplot(subset(singleint,singleint$best10==1&
                     singleint$order %in% c(6,9999)&singleint$ODint==0&
                     singleint$InterventionParam %in% c("phip","phip_tau")&
                     singleint$IntLevelpct %in% c("10%","20%","40%")), 
            aes(x=as.factor(ordertau), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("0%","25%","50%","75%","100%"))+ 
  #scale_y_continuous(limits=c(9500,21000),breaks=c(10000,15000,20000),labels=c("10,000","15,000","20,000"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") 
leg

leg<-g_legend(leg)

library(stringr)

text.ylab <- "% Reduction Chronic HCV"
p1<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(6,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phip","phip_tau")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(ordertau), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("0%","25%","50%","75%","100%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(15,100),breaks=c(20,60,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Transmission During Treatment")+
  theme(plot.title = element_text(hjust = 0.5))
p1


text.ylab <- "% Reduction Acute HCV"
p2<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(6,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phip","phip_tau")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(ordertau), y=pctred_acute, by=IntLevelpct)) +
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("0%","25%","50%","75%","100%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(5,100),breaks=c(10,55,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black')) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Transmission During Treatment")+
  theme(plot.title = element_text(hjust = 0.5)) 
p2


table(subset(singleint,singleint$order %in% c(6,9999)&
               singleint$InterventionParam %in% c("phip","phip_alpha"))$alphapct)

table(singleint$alphapct)

singleint$orderalpha<-ifelse(singleint$alphapct=="100%",1,
                             ifelse(singleint$alphapct=="90%",2,
                                    ifelse(singleint$alphapct=="80%",3,
                                           ifelse(singleint$alphapct=="70%",4,
                                                  ifelse(singleint$alphapct=="60%",5,9999)))))


text.ylab <- "% Reduction Chronic HCV"
p3<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(6,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phip","phip_alpha")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderalpha), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("100%","90%","80%","70%","60%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(15,100),breaks=c(20,60,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Sustained Virologic Response")+
  theme(plot.title = element_text(hjust = 0.5))
p3

text.ylab <- "% Reduction Acute HCV"
p4<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(6,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phip","phip_alpha")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderalpha), y=pctred_acute, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("100%","90%","80%","70%","60%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(5,100),breaks=c(10,55,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Sustained Virologic Response")+
  theme(plot.title = element_text(hjust = 0.5))
p4

table(subset(singleint,singleint$order %in% c(8,9999)&
               singleint$InterventionParam %in% c("phi_both","phiboth_omega"))$omegapct)

table(singleint$omegapct)

singleint$orderomega<-ifelse(singleint$omegapct=="8 Weeks",1,
                             ifelse(singleint$omegapct=="12 Weeks",2,
                                    ifelse(singleint$omegapct=="16 Weeks",3,9999)))


text.ylab <- "% Reduction Chronic HCV"
p5<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(6,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phip","phip_omega")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderomega), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("8","12","16"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(15,100),breaks=c(20,60,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("Weeks")+
  ggtitle("Treatment Duration")+
  theme(plot.title = element_text(hjust = 0.5))
p5

text.ylab <- "% Reduction Acute HCV"
p6<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(6,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phip","phip_omega")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderomega), y=pctred_acute, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("8","12","16"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(5,100),breaks=c(10,55,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("Weeks")+
  ggtitle("Treatment Duration")+
  theme(plot.title = element_text(hjust = 0.5))
p6

#No title
w <- 12; h <- 12
lay<-rbind(c(1,2),c(1,2),c(1,2),c(3,4),c(3,4),c(3,4),c(5,6),c(5,6),c(5,6),c(7,7))
SingleStrategy_VP_supp<-grid.arrange(p1,p2,p3,p4,p5,p6,leg,layout_matrix=lay)
ggsave(sprintf("TrtSensVP_CurrPWID_Best10_grey.pdf"), SingleStrategy_VP_supp, width=w, height=h)

#### Violin Plots of % Cure on Treat Current PWID (for ppt) ####
table(subset(singleint,singleint$order %in% c(6,9999)&
               singleint$InterventionParam %in% c("phip","phip_alpha"))$alphapct)

table(singleint$alphapct)

singleint$orderalpha<-ifelse(singleint$alphapct=="100%",1,
                             ifelse(singleint$alphapct=="90%",2,
                                    ifelse(singleint$alphapct=="80%",3,
                                           ifelse(singleint$alphapct=="70%",4,
                                                  ifelse(singleint$alphapct=="60%",5,9999)))))


text.ylab <- "% Reduction in HCV Prevalence"
p3<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(6,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phip","phip_alpha")&singleint$orderalpha %in% c(2,5)&
                    singleint$IntLevelpct %in% c("10%","20%")), 
           aes(x=as.factor(orderalpha), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#3366FF','#000066'), name='Intervention Level') +
  scale_x_discrete(labels=c("90%","60%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("% Current PWID Cured by HCV Treatment")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))
p3

p3<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(6,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phip","phip_alpha")&singleint$orderalpha %in% c(2,5)&
                    singleint$IntLevelpct %in% c("10%","20%")), 
           aes(x=as.factor(orderalpha), y=pctred_chr, by=IntLevelpct)) + 
  #geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  #geom_boxplot(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  geom_line(aes(colour=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#3366FF','#000066'), name='Intervention Level') +
  scale_x_discrete(labels=c("90%","60%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("% Current PWID Cured by HCV Treatment")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))
p3

new<-subset(singleint,singleint$best10==1&
         singleint$order %in% c(6,9999)&singleint$ODint==0&
         singleint$InterventionParam %in% c("phip","phip_alpha")&singleint$orderalpha %in% c(2,5)&
         singleint$IntLevelpct %in% c("10%","20%"))

median(subset(new,new$orderalpha==2&new$IntLevelpct=="10%")$pctred_chr) #35.70323
median(subset(new,new$orderalpha==5&new$IntLevelpct=="10%")$pctred_chr) #26.27147

median(subset(new,new$orderalpha==2&new$IntLevelpct=="20%")$pctred_chr) #57.26663
median(subset(new,new$orderalpha==5&new$IntLevelpct=="20%")$pctred_chr) #44.6154

#### Violin Plots of Treatment Variations (Treat Former PWID) % Red Best 10% ####

## note: no tau simulation for former pwid (since it reflects transmission during trt)

table(subset(singleint,singleint$order %in% c(7,9999)&
               singleint$InterventionParam %in% c("phin","phin_alpha"))$alphapct)

table(singleint$alphapct)

singleint$orderalpha<-ifelse(singleint$alphapct=="100%",1,
                             ifelse(singleint$alphapct=="90%",2,
                                    ifelse(singleint$alphapct=="80%",3,
                                           ifelse(singleint$alphapct=="70%",4,
                                                  ifelse(singleint$alphapct=="60%",5,9999)))))



text.ylab <- "% Reduction HCV Prevalence"
leg<-ggplot(subset(singleint,singleint$best10==1&
                     singleint$order %in% c(7,9999)&singleint$ODint==0&
                     singleint$InterventionParam %in% c("phin","phin_alpha")&
                     singleint$IntLevelpct %in% c("10%","20%","40%")), 
            aes(x=as.factor(orderalpha), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("100%","90%","80%","70%","60%"))+#c("","","",""))+ 
  #scale_y_continuous(limits=c(9500,21000),breaks=c(10000,15000,20000),labels=c("10,000","15,000","20,000"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("") 
leg

leg<-g_legend(leg)


text.ylab <- "% Reduction Chronic HCV"
p3<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(7,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phin","phin_alpha")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderalpha), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("100%","90%","80%","70%","60%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(15,100),breaks=c(20,60,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Sustained Virologic Response")+
  theme(plot.title = element_text(hjust = 0.5))
p3

text.ylab <- "% Reduction Acute HCV"
p4<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(7,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phin","phin_alpha")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderalpha), y=pctred_acute, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("100%","90%","80%","70%","60%"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,50,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("%")+
  ggtitle("Sustained Virologic Response")+
  theme(plot.title = element_text(hjust = 0.5))
p4

table(subset(singleint,singleint$order %in% c(8,9999)&
               singleint$InterventionParam %in% c("phi_both","phiboth_omega"))$omegapct)

table(singleint$omegapct)

singleint$orderomega<-ifelse(singleint$omegapct=="8 Weeks",1,
                             ifelse(singleint$omegapct=="12 Weeks",2,
                                    ifelse(singleint$omegapct=="16 Weeks",3,9999)))


text.ylab <- "% Reduction Chronic HCV"
p5<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(7,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phin","phin_omega")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderomega), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("8","12","16"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(15,100),breaks=c(20,60,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("Weeks")+
  ggtitle("Treatment Duration")+
  theme(plot.title = element_text(hjust = 0.5))
p5

text.ylab <- "% Reduction Acute HCV"
p6<-ggplot(subset(singleint,singleint$best10==1&
                    singleint$order %in% c(7,9999)&singleint$ODint==0&
                    singleint$InterventionParam %in% c("phin","phin_omega")&
                    singleint$IntLevelpct %in% c("10%","20%","40%")), 
           aes(x=as.factor(orderomega), y=pctred_acute, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Intervention Level') +
  scale_x_discrete(labels=c("8","12","16"))+#c("","","",""))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,50,100))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("Weeks")+
  ggtitle("Treatment Duration")+
  theme(plot.title = element_text(hjust = 0.5))
p6

#No title
w <- 12; h <- 8
lay<-rbind(c(1,2),c(1,2),c(1,2),c(3,4),c(3,4),c(3,4),c(7,7))
SingleStrategy_VP_supp<-grid.arrange(p3,p4,p5,p6,leg,layout_matrix=lay)
ggsave(sprintf("TrtSensVP_FormPWID_Best10.pdf"), SingleStrategy_VP_supp, width=w, height=h)





#### Ridge Plots of Best Fitting 10% of Primary/Tertiary Sequential Interventions ####
summary(subset(seqint,seqint$best10==1)$pctred_chr)
summary(subset(seqint,seqint$best10==1)$pctred_acute)

#subset data to exclude 30% lines, select h==1, and best 50% results
int<-subset(seqint,e %in% c(170:194,220:244)&pct %in% c("10%","20%","40%")&best10==1)

#calculate median % reduction for each intervention combination
median<-cbind(melt(tapply(int$pctred_chr,list(int$intono,int$pct,int$inttype),median),id="pct"),
              melt(tapply(int$pctred_acute,list(int$intono,int$pct,int$inttype),median),id="pct"))
median<-median[,c(1,2,3,4,8)]
colnames(median) <- c("intono",	"pct",	"inttype",	"med_pctredchr", "med_pctredacute")

medianprim<-subset(median,inttype=="primtotert")
mediantert<-subset(median,inttype=="terttoprim")

#Intervention % legend
p1<-ggplot(data=subset(int,inttype=="terttoprim"),aes(y=as.factor(intono-1), x=pctred_chr))+
  geom_density_ridges(scale=1,aes(fill=pct))+
  #theme_joy()+
  geom_point(data=subset(median,inttype=="terttoprim"),
             x=subset(median,inttype=="terttoprim")$med_pctredchr,
             y=subset(median,inttype=="terttoprim")$intono-0.92,
             by=subset(median,inttype=="terttoprim")$pct,
             size=1.5,shape=23,fill="white")+
  xlab("% Reduction: Prevalence")+
  ylab("Intervention Sequence")+
  theme_bw(base_size=14)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="bottom", legend.background = element_rect(colour = 'NA', 
                                                                   fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_fill_cyclical(guide="legend",values = c("#CCCCCC", "#666666", "#000000"),
                      labels=c("10%","20%","40%"),name="Intervention %")+
  #scale_fill_cyclical(guide="legend",values = c('#99CCFF','#3366FF', '#000066'),
  #                    labels=c("10%","20%","40%"),name="Intervention %")+
  scale_y_discrete(breaks = 1:6, labels=c("Treat Former PWID","+Treat Current PWID","+Decrease Relapse  ",
                                          "+Increase Cessation ", "+Decrease Contacts ",
                                          "+Decrease Initiation  "))+ 
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around edge of figure

p1

p3<-ggplot(data=subset(int,inttype=="terttoprim"),aes(y=as.factor(intono-1), x=pctred_chr))+
  geom_density_ridges(scale=1,aes(fill=pct))+
  #theme_joy()+
  geom_point(data=subset(median,inttype=="terttoprim"),
             x=subset(median,inttype=="terttoprim")$med_pctredchr,
             y=subset(median,inttype=="terttoprim")$intono-0.92,
             by=subset(median,inttype=="terttoprim")$pct,
             size=1.5,shape=23,fill="white")+
  xlab("% Reduction: Prevalence")+
  ylab("Intervention Sequence")+
  theme_bw(base_size=14)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', 
                                                                 fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_fill_cyclical(guide="legend",values = c("#CCCCCC", "#666666", "#000000"),
                      labels=c("10%","20%","40%"),name="Intervention %")+
  #scale_fill_cyclical(guide="legend",values = c('#99CCFF','#3366FF', '#000066'),
  #                    labels=c("10%","20%","40%"),name="Intervention %")+
  scale_y_discrete(breaks = 1:6, labels=c("Treat Former PWID","+Treat Current PWID","+Decrease Relapse  ",
                                          "+Increase Cessation ", "+Decrease Contacts ",
                                          "+Decrease Initiation  "))+ 
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

p3
# title.grob <- textGrob(label = "D)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = -0.5,gp = gpar(fontsize = 16))
# 
# p3 <- arrangeGrob(p3, top = title.grob)
# grid.arrange(p3)

p4<-ggplot(data=subset(int,inttype=="terttoprim"),aes(y=as.factor(intono-1), x=pctred_acute))+
  geom_density_ridges(scale=1,aes(fill=pct))+
  #theme_joy()+
  geom_point(data=subset(median,inttype=="terttoprim"),
             x=subset(median,inttype=="terttoprim")$med_pctredacute,
             y=subset(median,inttype=="terttoprim")$intono-0.92,
             by=subset(median,inttype=="terttoprim")$pct,
             size=1.5,shape=23,fill="white")+
  ylab("Intervention Sequence")+
  xlab("% Reduction: Acute")+
  theme_bw(base_size=14)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', 
                                                                 fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_fill_cyclical(guide="legend",values = c("#CCCCCC", "#666666", "#000000"),
                      labels=c("10%","20%","40%"),name="Intervention %")+
  #scale_fill_cyclical(guide="legend",values = c('#99CCFF','#3366FF', '#000066'),
  #                    labels=c("10%","20%","40%"),name="Intervention %")+
  scale_y_discrete(breaks = 1:6, labels=c("Treat Former PWID","+Treat Current PWID","+Decrease Relapse  ",
                                          "+Increase Cessation ", "+Decrease Contacts ",
                                          "+Decrease Initiation  "))+ 
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
p4
# #+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
# 
# title.grob <- textGrob(label = "B)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = -0.5,gp = gpar(fontsize = 16))
# 
# p4 <- arrangeGrob(p4, top = title.grob)
# grid.arrange(p4)

p5<-ggplot(data=subset(int,inttype=="primtotert"),aes(y=as.factor(intono-1), x=pctred_chr))+
  geom_density_ridges(scale=1,aes(fill=pct))+
  #theme_joy()+
  geom_point(data=subset(median,inttype=="primtotert"),
             x=subset(median,inttype=="primtotert")$med_pctredchr,
             y=subset(median,inttype=="primtotert")$intono-0.92,
             by=subset(median,inttype=="primtotert")$pct,
             size=1.5,shape=23,fill="white")+
  xlab("% Reduction: Prevalence")+
  ylab("Intervention Sequence")+
  theme_bw(base_size=14)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', 
                                                                 fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_fill_cyclical(guide="legend",values = c("#CCCCCC", "#666666", "#000000"),
                      labels=c("10%","20%","40%"),name="Intervention %")+
  #scale_fill_cyclical(guide="legend",values = c('#99CCFF','#3366FF', '#000066'),
  #                    labels=c("10%","20%","40%"),name="Intervention %")+
  scale_y_discrete(breaks = 1:6, labels=c("Decrease Initiation  ","+Decrease Contacts ",
                                          "+Increase Cessation ", "+Decrease Relapse  ",
                                          "+Treat Current PWID","+Treat Former PWID"))+ 
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

p5
# title.grob <- textGrob(label = "C)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = -0.5,gp = gpar(fontsize = 16))
# 
# p5 <- arrangeGrob(p5, top = title.grob)
# grid.arrange(p5)

p6<-ggplot(data=subset(int,inttype=="primtotert"),aes(y=as.factor(intono-1), x=pctred_acute))+
  geom_density_ridges(scale=1,aes(fill=pct))+
  #theme_joy()+
  geom_point(data=subset(median,inttype=="primtotert"),
             x=subset(median,inttype=="primtotert")$med_pctredacute,
             y=subset(median,inttype=="primtotert")$intono-0.92,
             by=subset(median,inttype=="primtotert")$pct,
             size=1.5,shape=23,fill="white")+
  xlab("% Reduction: Acute")+
  ylab("Intervention Sequence")+
  theme_bw(base_size=14)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', 
                                                                 fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_fill_cyclical(guide="legend",values = c("#CCCCCC", "#666666", "#000000"),
                      labels=c("10%","20%","40%"),name="Intervention %")+
  #scale_fill_cyclical(guide="legend",values = c('#99CCFF','#3366FF', '#000066'),
  #                    labels=c("10%","20%","40%"),name="Intervention %")+
  scale_y_discrete(breaks = 1:6, labels=c("Decrease Initiation  ","+Decrease Contacts ",
                                          "+Increase Cessation ", "+Decrease Relapse  ",
                                          "+Treat Current PWID","+Treat Former PWID"))+ 
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
# 
# title.grob <- textGrob(label = "A)",x = unit(0, "lines"), y = unit(0, "lines"),
#                        hjust = -1, vjust = -0.5,gp = gpar(fontsize = 16))
# 
# p6 <- arrangeGrob(p6, top = title.grob)
# grid.arrange(p6)
p6

#no title
w <- 10; h <- 6.5
legendpoint<-g_legend(p1)
lay<-rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(5,5))
Multi_Sequential_ranges<-grid.arrange(p6,p4,p5,p3,legendpoint,layout_matrix=lay)
ggsave(sprintf("Multi_Sequential_Best10_grey.pdf"), Multi_Sequential_ranges, width=w, height=h)
#ggsave(sprintf("Multi_Sequential_Best10_Blue.pdf"), Multi_Sequential_ranges, width=w, height=h)


#### Summary Statistics for Interventions (to report in paper) ####
sub<-subset(singleint,singleint$e==62&singleint$best10==1) #trt former and current PWID at 10% per year
summary(sub$pctred_acute)
summary(sub$pctred_chr)

sub<-subset(singleint,singleint$e==22&singleint$best10==1) #trt former PWID at 10% per year
summary(sub$pctred_acute)
summary(sub$pctred_chr)

sub<-subset(singleint,singleint$e==6&singleint$best10==1) #Reduce contacts by 10%
summary(sub$pctred_acute)
summary(sub$pctred_chr)

sub<-subset(seqint,seqint$e==183&seqint$best10==1) #Reduce contacts, syr share, relapse, increase cess by 10%
summary(sub$pctred_acute)
summary(sub$pctred_chr)

sub<-subset(singleint,singleint$e==18&singleint$best10==1) #Trt SVR=90% for 10% trt level current PWID
summary(sub$pctred_acute)
summary(sub$pctred_chr)

sub<-subset(singleint,singleint$e==158&singleint$best10==1) #Trt SVR=60% for 10% trt level current PWID
summary(sub$pctred_acute)
summary(sub$pctred_chr)

sub<-subset(singleint,singleint$e==166&singleint$best10==1) #Trt SVR=60% for 10% trt level former and current PWID
summary(sub$pctred_acute)
summary(sub$pctred_chr)


#### Acute vs New Chronic 15-19 YOs (MDHHS Data) ####
year<-seq(2000,2016,1)
acutepwid1519<-c(4,	1,	1,	1,	4,	1,	2,	7,	5,	3,	2,	2,	5,	2,	2,	1,	2*(1-0.26))
newcpwid1519<-c(10,	21,	19,	40,	33,	53,	65,	57,	84,	84,	52,	46,	106,	85,	76,	76,	104*(1-0.26))

data<-data.frame(cbind(year,acutepwid1519,newcpwid1519))

xbreaks = c("Acute", "New Chronic","New Chronic, Scaled for Under-Detection")

p1<-ggplot(data)+
  geom_line(aes(x=year, y=acutepwid1519),colour="#000000",size=1,show.legend = FALSE) +
  geom_line(aes(x=year, y=newcpwid1519),colour="#000000",size=1,show.legend = FALSE)+
  geom_line(aes(x=year, y=newcpwid1519*(1/16.8)),colour="#000000",size=1,show.legend = FALSE)+
  geom_point(aes(x=year, y=acutepwid1519,colour="Acute"), size=1.5,show.legend = TRUE, na.rm=TRUE) + 
  geom_point(aes(x=year, y=newcpwid1519,colour="New Chronic"), size=1.5,show.legend = TRUE, na.rm=TRUE)+
  geom_point(aes(x=year, y=newcpwid1519*(1/16.8),colour="New Chronic, Scaled for Under-Detection"), size=1.5,show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = c(2000,2005,2010,2015))+
  scale_y_continuous(limits = c(0, max(newcpwid1519)+1))+
  theme_bw(base_size=16) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_manual(name=NULL, values = c("#B8361B","dark green","#152FA5"),breaks = xbreaks,
                      #labels=str_wrap(xbreaks,15))+
                      labels=xbreaks)+
  #guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("Number of Cases") +
  xlab("Year") +
  #ggtitle('Acute and New Chronic Cases Aged 15-19')+
  theme(plot.title = element_text(hjust = 0.5)) 
p1

leg<-g_legend(p1)

p2<-ggplot(data)+
  geom_line(aes(x=year, y=acutepwid1519),colour="#000000",size=1,show.legend = FALSE) +
  geom_line(aes(x=year, y=newcpwid1519),colour="#000000",size=1,show.legend = FALSE)+
  geom_line(aes(x=year, y=newcpwid1519*(1/16.8)),colour="#000000",size=1,show.legend = FALSE)+
  geom_point(aes(x=year, y=acutepwid1519,colour="Acute",size=1.5), show.legend = TRUE, na.rm=TRUE) + 
  geom_point(aes(x=year, y=newcpwid1519,colour="New Chronic",size=1.5), show.legend = TRUE, na.rm=TRUE)+
  geom_point(aes(x=year, y=newcpwid1519*(1/16.8),colour="New Chronic, Scaled for Under-Detection",size=1.5), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = c(2000,2005,2010,2015))+
  scale_y_continuous(limits = c(0, max(newcpwid1519)+1))+
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_manual(name=NULL, values = c("#B8361B","dark green","#152FA5"),breaks = xbreaks,
                      #labels=str_wrap(xbreaks,15))+
                      labels=xbreaks)+
  #guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("Number of Cases") +
  xlab("Year") +
  ggtitle('Acute and New Chronic Cases')+
  theme(plot.title = element_text(hjust = 0.5)) 
p2

p3<-ggplot(data)+
  geom_line(aes(x=year, y=acutepwid1519),colour="#000000",size=1,show.legend = FALSE) +
  #geom_line(aes(x=year, y=newcpwid1519),colour="#000000",size=1,show.legend = FALSE)+
  geom_line(aes(x=year, y=newcpwid1519*(1/16.8)),colour="#000000",size=1,show.legend = FALSE)+
  geom_point(aes(x=year, y=acutepwid1519,colour="Acute",size=1.5), show.legend = TRUE, na.rm=TRUE) + 
  #geom_point(aes(x=year, y=newcpwid1519,colour="New Chronic",size=1.5), show.legend = TRUE, na.rm=TRUE)+
  geom_point(aes(x=year, y=newcpwid1519*(1/16.8),colour="New Chronic, Scaled for Under-Detection",size=1.5), show.legend = TRUE, na.rm=TRUE)+
  scale_x_continuous(limits = c(2000, 2016), breaks = c(2000,2005,2010,2015))+
  scale_y_continuous(limits = c(0, max(newcpwid1519*(1/16.8))+1))+
  theme_bw(base_size=16) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'NA'), legend.background = element_rect(colour = "NA")) +
  scale_colour_manual(name=NULL, values = c("#B8361B","#152FA5"),breaks = c("Acute","New Chronic"),
                      #labels=str_wrap(xbreaks,15))+
                      labels=c("Acute","New Chronic"))+
  #guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))+
  ylab("Number of Cases") +
  xlab("Year") +
  ggtitle('Acute and New Chronic Cases Adjusted for Under-Detection')+
  theme(plot.title = element_text(hjust = 0.5)) 
p3

#no title
w <- 8; h <- 10

lay<-rbind(c(1),c(1),c(1),c(2),c(2),c(2),c(3))
chrplot<-grid.arrange(p2,p3,leg,layout_matrix=lay)
ggsave(sprintf("AcuteNewChronic1519_Data.pdf"), chrplot, width=w, height=h)


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

p1<-ggplot(data)+
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
p1





