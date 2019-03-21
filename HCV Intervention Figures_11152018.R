## Figure Script: Gicquelais, Foxman, Coyle, and Eisenberg. (2019). Hepatitis C transmission in young 
## people who inject drugs: insights using a dynamic model informed by state public health surveillance. 
## Epidemics. https://doi.org/10.1016/j.epidem.2019.02.003.

## Before running this code, run the model simulations in the code HCV_MS_11052018.m to simulate
## the hepatitis C model and output results as csv or txt files. Next, create a workspace using 
## the R Script HCV Figures_11152018.R. This script creates figures 
## summarizing results of simulated interventions similar to those published in the article. 

################ Section 1: Packages and Functions to Load First ################### 
#note that ggplot version 3.0.0 is required, version 2.2.1 gives error
require(devtools)
install_version("ggplot2", version = "3.0.0", repos = "http://cran.us.r-project.org")

library(ggplot2); packageVersion("ggplot2") #3.0.0
#require(devtools)
#install_version("ggridges", version = "0.5.0", repos = "http://cran.us.r-project.org")
library(ggridges);packageVersion("ggridges") #0.5.0
library(grid); packageVersion("grid") #3.5.1
library(gridExtra)
library(car)
library(plotly)
library(plyr); packageVersion("plyr") #1.8.4
library(reshape)
library(stringr)
library(scales); packageVersion("scales") #1.0.0
library(withr); packageVersion("withr") #2.1.2

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



################ Section 2: Set Working Directory and Create Intervention Results Workspace ################
setwd("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 11.15.2018")

################ Section 3: Import Intervention Results and Save Workspace ################

#import datasets with intervention levels + all case counts 
prev<-read.csv("ChronicPrev_SingleSeq.txt", header=F)
names(prev)[1:15] <- c("e","int1",	"int2",	"int3",	"int4",	"int5",	"int6",	"int7",	
                       "int8",	"int9",	"int10","int11","i","ChronicPrev_End_All","ChronicPrev_End")
prev$inttype<-ifelse(prev$e<=169,"single",ifelse(prev$e>=220,"terttoprim","primtotert"))
table(prev$inttype)

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

#Create datasets with case counts at t151 (final data point)
prev1<-prev[ ,c(1:16)]
prev1$index<-as.numeric(paste(prev1$e,prev1$i,sep=""))

acute1<-acute[ ,c(1,13:15)]
acute1$index<-as.numeric(paste(acute1$e,acute1$i,sep=""))

#Setup params dataset (remove int variables)
param1<-param[ ,c(1,13:65)]
param1$index<-as.numeric(paste(param1$e,param1$i,sep=""))

#merge case count summaries with params
merged<-merge(merge(prev1, acute1, by=c("e","i"),suffixes = c(".c",".d")), 
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
df$phippct<-ifelse(int==1,"None",ifelse(int==2,"3%/Year",ifelse(int==3,"6%/Year",ifelse(int==4,"12%/Year",ifelse(int==5,"24%/Year","error")))))
int<-int6
df$phinpct<-ifelse(int==1,"None",ifelse(int==2,"3%/Year",ifelse(int==3,"6%/Year",ifelse(int==4,"12%/Year",ifelse(int==5,"24%/Year","error")))))
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

rm(param,prev,acute,prev1,acute1,df,param1)


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

df<-noint[1:10000,c(1,14,15,18,19)] #select e=1 (vars=i and case counts)
df<-rename(df, c("ChronicPrev_End" = "ChronicPrev_e1", "Acute_End" = "Acute_e1",
                 "ChronicPrev_End_All" = "ChronicPrev_All_e1", "Acute_End_All" = "Acute_All_e1"))

df$ChronicPrev_e34<-subset(noint,noint$e==34)$ChronicPrev_End
df$Acute_e34<-subset(noint,noint$e==34)$Acute_End

df$ChronicPrev_All_e34<-subset(noint,noint$e==34)$ChronicPrev_End_All
df$Acute_All_e34<-subset(noint,noint$e==34)$Acute_End_All

#merge back with intervention dataset 
singleint<-merge(singleint,df,by="i")
rm(df,noint,add)

# Calculate % reduction compared to appropriate none intervention
singleint$pctred_chr<-ifelse(singleint$e %in% c(1,34),0,
                             ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$ChronicPrev_End/singleint$ChronicPrev_e1)*100,
                                    ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$ChronicPrev_End/singleint$ChronicPrev_e34)*100,99999)))

singleint$pctred_acute<-ifelse(singleint$e %in% c(1,34),0,
                               ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$Acute_End/singleint$Acute_e1)*100,
                                      ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$Acute_End/singleint$Acute_e34)*100,99999)))

summary(singleint$pctred_chr)
summary(singleint$pctred_acute)

# Calculate % reduction compared to appropriate none intervention
singleint$pctred_chr_all<-ifelse(singleint$e %in% c(1,34),0,
                                 ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$ChronicPrev_End_All/singleint$ChronicPrev_All_e1)*100,
                                        ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$ChronicPrev_End_All/singleint$ChronicPrev_All_e34)*100,99999)))

singleint$pctred_acute_all<-ifelse(singleint$e %in% c(1,34),0,
                                   ifelse(singleint$e %in% c(2:33,59:169),100-(singleint$Acute_End_All/singleint$Acute_All_e1)*100,
                                          ifelse(singleint$e %in% c(35:58,250:253),100-(singleint$Acute_End_All/singleint$Acute_All_e34)*100,99999)))

summary(singleint$pctred_chr_all)
summary(singleint$pctred_acute_all)

# Calculate # Infections Averted
singleint$averted_chr_all<-ifelse(singleint$e %in% c(1,34),0,
                                     ifelse(singleint$e %in% c(2:33,59:169),singleint$ChronicPrev_End_All-singleint$ChronicPrev_All_e1,
                                            ifelse(singleint$e %in% c(35:58,250:253),singleint$ChronicPrev_End_All-singleint$ChronicPrev_All_e34,99999)))

singleint$averted_chr<-ifelse(singleint$e %in% c(1,34),0,
                                  ifelse(singleint$e %in% c(2:33,59:169),singleint$ChronicPrev_End-singleint$ChronicPrev_e1,
                                         ifelse(singleint$e %in% c(35:58,250:253),singleint$ChronicPrev_End-singleint$ChronicPrev_e34,99999)))

singleint$averted_acute_all<-ifelse(singleint$e %in% c(1,34),0,
                                  ifelse(singleint$e %in% c(2:33,59:169),singleint$Acute_End_All-singleint$Acute_All_e1,
                                         ifelse(singleint$e %in% c(35:58,250:253),singleint$Acute_End_All-singleint$Acute_All_e34,99999)))

singleint$averted_acute<-ifelse(singleint$e %in% c(1,34),0,
                                    ifelse(singleint$e %in% c(2:33,59:169),singleint$Acute_End-singleint$Acute_e1,
                                           ifelse(singleint$e %in% c(35:58,250:253),singleint$Acute_End-singleint$Acute_e34,99999)))

summary(singleint$averted_acute,exclude=FALSE)
summary(singleint$averted_chr,exclude=FALSE)
summary(singleint$averted_acute_all,exclude=FALSE)
summary(singleint$averted_chr_all,exclude=FALSE)

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
interventioninfo <- singleint[,c(2:13,16,74:84,90:93)] #vars=things that do not change by e
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

df<-noint[1:10000,c(1,14,15,18,19)] #select e=170 (vars=i and case counts)
df<-rename(df, c("ChronicPrev_End" = "ChronicPrev_e170", "Acute_End" = "Acute_e170", 
                 "ChronicPrev_End_All" = "ChronicPrev_All_e170", "Acute_End_All" = "Acute_All_e170"))

df$ChronicPrev_e195<-subset(noint,noint$e==195)$ChronicPrev_End
df$Acute_e195<-subset(noint,noint$e==195)$Acute_End

df$ChronicPrev_All_e195<-subset(noint,noint$e==195)$ChronicPrev_End_All
df$Acute_All_e195<-subset(noint,noint$e==195)$Acute_End_All


#merge back with intervention dataset 
seqint<-merge(seqint,df,by="i")
rm(df,noint)

# Calculate % reduction compared to appropriate none intervention
seqint$pctred_chr<-ifelse(seqint$e %in% c(170,195,220,245),0,
                          ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$ChronicPrev_End/seqint$ChronicPrev_e170)*100,
                                 ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$ChronicPrev_End/seqint$ChronicPrev_e195)*100,99999)))

seqint$pctred_acute<-ifelse(seqint$e %in% c(170,195,220,245),0,
                            ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$Acute_End/seqint$Acute_e170)*100,
                                   ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$Acute_End/seqint$Acute_e195)*100,99999)))

summary(seqint$pctred_chr)
summary(seqint$pctred_acute)

seqint$pctred_chr_all<-ifelse(seqint$e %in% c(170,195,220,245),0,
                              ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$ChronicPrev_End_All/seqint$ChronicPrev_All_e170)*100,
                                     ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$ChronicPrev_End_All/seqint$ChronicPrev_All_e195)*100,99999)))

seqint$pctred_acute_all<-ifelse(seqint$e %in% c(170,195,220,245),0,
                                ifelse(seqint$e %in% c(171:194,221:244),100-(seqint$Acute_End_All/seqint$Acute_All_e170)*100,
                                       ifelse(seqint$e %in% c(196:219,246:269),100-(seqint$Acute_End_All/seqint$Acute_All_e195)*100,99999)))

summary(seqint$pctred_chr_all)
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
interventioninfo <- seqint[,c(2:13,16,74:84,90:91)]
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

#Save intervention results workspace
save.image("C:/Users/rgicquel/Desktop/11.15.2018/Dataset Workspace_11.15.18_Interv.RData")



################ Section 4: Load Intervention Simulation Results ################

load("/Users/RGicquelais/Desktop/HCV/MDCH Young HCV Modeling/Figures and Datasets/Final 11.15.2018/Dataset Workspace_11.15.18_Interv.RData")



################ Section 5: Figures of Intervention Simulation Results ################
#### Violin Plots of Best 10% Fits of Single Interventions ####
#Set plot labels and colors
text.xlab <- "Intervention"
model.colors <- c('#CCCCCC', '#666666', '#000000','#FFFFFF') #grey scheme
#model.colors <- c('#FF9999', '#CC3333','#990000','#FFFFFF') #red scheme
#model.colors <- c('#99CCFF','#3366FF', '#000066','#FFFFFF') #blue scheme

table(singleint$order)
table(singleint$IntLevelpct)

summary(singleint$fitquartile)


#create new variable to force none to plot 1st in legend as white color 
singleint$plotleg<-ifelse(singleint$IntLevelpct=="None",1,ifelse(singleint$IntLevelpct=="10%",2,
                                                                 ifelse(singleint$IntLevelpct=="20%",3,4)))
singleint$type<-ifelse(singleint$order==1,"None",ifelse(singleint$order==2,"Primary",
                                                        ifelse(singleint$order %in% c(3,4,5),"Secondary",
                                                               ifelse(singleint$order %in% c(6,7,8),"Tertiary","NA"))))

#X-axis labels
xlabelsmain <-paste0(c("None","Decrease Injection Initiation",
                       "Decrease Current PWID Contacts",
                       "Increase Current PWID Cessation", 
                       "Decrease Former PWID Relapse",
                       "Treat Current PWID",
                       "Treat Former PWID","Treat Current and Former PWID"))


text.ylab <- "% Reduction Chronic HCV"

summary(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%"))$pctred_chr)

#legend for primary/secondary
leg<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                     singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%")),
            aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Primary & Secondary Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Primary & Secondary Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=c("","","","","","",""))+ 
  scale_y_continuous(limits=c(-1,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=15) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
leg

legend1<-g_legend(leg)

#legend for treatment
leg2<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                     singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%")),
            aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Treatment Rate', labels=c("0.03/Year ","0.06/Year","0.24/Year")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Treatment Rate', labels=c("0.03/Year ","0.06/Year","0.24/Year")) +
  scale_x_discrete(labels=c("","","","","","",""))+ 
  scale_y_continuous(limits=c(-1,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=15) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
leg2

legend2<-g_legend(leg2)

allchronic<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                            singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%")),
                   aes(x=order, y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(-1,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
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
                 aes(x=order, y=pctred_acute, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
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


w <- 10.5; h <- 8
lay<-rbind(c(1), 
           c(1),  
           c(1),c(1),c(1), 
           c(2),
           c(2),
           c(2),c(2),c(2),
           c(3),c(4))
SingleStrategy_VP_main<-grid.arrange(allchronic,allacute,legend1,legend2,layout_matrix=lay)
#ggsave(sprintf("SingleStrategy_best10_5_blue.pdf"), SingleStrategy_VP_main, width=w, height=h)
ggsave(sprintf("SingleStrategy_best10_5_grey.pdf"), SingleStrategy_VP_main, width=w, height=h)

#### Violin Plots of 30-64 YO single interventions
allchronic<-ggplot(subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&
                            singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%")),
                   aes(x=order, y=pctred_chr_all, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  geom_boxplot(aes(fill=as.factor(plotleg),colour=as.factor(plotleg)), position=position_dodge(width=0.8), width=0.01) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
 scale_colour_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
   #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(-1.1,100),breaks=c(-50,-25,0,25,50,75,100),labels=c("Increase 50%","Increase 25%","Null Impact","Reduce 25%"," Reduce 50%","Reduce 75%"," Reduce 100%"))+
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
                 aes(x=order, y=pctred_acute_all, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  geom_boxplot(aes(fill=as.factor(plotleg),colour=as.factor(plotleg)), position=position_dodge(width=0.8), width=0.01) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_colour_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(-52,100),breaks=c(-50,-25,0,25,50,75,100),labels=c("Increase 50%","Increase 25%","Null Impact","Reduce 25%"," Reduce 50%","Reduce 75%"," Reduce 100%"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", 
               position=position_dodge(width=0.8)) +
  theme_bw(base_size=15) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("% Reduction Acute HCV") +
  xlab("")
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
allacute

w <- 13; h <- 8
lay<-rbind(c(1), 
           c(1),  
           c(1),c(1),c(1), 
           c(2),
           c(2),
           c(2),c(2),c(2),
           c(3),c(4))
SingleStrategy_VP_main<-grid.arrange(allchronic,allacute,legend1,legend2,layout_matrix=lay)
#ggsave(sprintf("SingleStrategy_best10_30-64YO_blue.pdf"), SingleStrategy_VP_main, width=w, height=h)
ggsave(sprintf("SingleStrategy_best10_30-64YO_grey.pdf"), SingleStrategy_VP_main, width=w, height=h)

#### Violin Plots of Treatment Variations (Treat Current PWID Sensitivity Analysis) % Red Best 10% ####


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
  scale_fill_manual(values=model.colors, name='Treatment Rate', labels=c("0.03/Year ","0.06/Year","0.24/Year")) +
  scale_x_discrete(labels=c("0%","25%","50%","75%","100%"))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
#ggsave(sprintf("TrtSensVP_CurrPWID_Best10_blue.pdf"), SingleStrategy_VP_supp, width=w, height=h)
ggsave(sprintf("TrtSensVP_CurrPWID_Best10_grey.pdf"), SingleStrategy_VP_supp, width=w, height=h)


#### Violin Plots of Treatment Variations (Treat Former PWID Sensitivity Analysis) % Red Best 10% ####

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
            aes(x=as.factor(ordertau), y=pctred_chr, by=IntLevelpct)) + 
  geom_violin(aes(fill=IntLevelpct), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=model.colors, name='Treatment Rate', labels=c("0.03/Year ","0.06/Year","0.24/Year")) +
  scale_x_discrete(labels=c("0%","25%","50%","75%","100%"))+ 
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+
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
ggsave(sprintf("TrtSensVP_FormPWID_Best10_grey.pdf"), SingleStrategy_VP_supp, width=w, height=h)



#### Ridge Plots of Best Fitting 10% of Primary/Tertiary Sequential Interventions ####
summary(subset(seqint,seqint$best10==1)$pctred_chr)
summary(subset(seqint,seqint$best10==1)$pctred_acute)

#subset data to exclude 30% lines, select h==1, and best 10% results
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
                      labels=c("10%","20%","40%"),name="Primary & Secondary Intervention Level")+
  scale_y_discrete(breaks = 1:6, labels=c("Treat Former PWID","+Treat Current PWID","+Decrease Relapse  ",
                                          "+Increase Cessation ", "+Decrease Contacts ",
                                          "+Decrease Initiation  "))+ 
  scale_x_continuous(limits=c(-1,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around edge of figure

p1

legendpoint<-g_legend(p1)

p2<-ggplot(data=subset(int,inttype=="terttoprim"),aes(y=as.factor(intono-1), x=pctred_chr))+
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
                      labels=c("0.03/Year","0.06/Year","0.24/Year"),name="Treatment Rate")+
  scale_y_discrete(breaks = 1:6, labels=c("Treat Former PWID","+Treat Current PWID","+Decrease Relapse  ",
                                          "+Increase Cessation ", "+Decrease Contacts ",
                                          "+Decrease Initiation  "))+ 
  scale_x_continuous(limits=c(-1,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around edge of figure

p2

legendtrt<-g_legend(p2)

p3<-ggplot(data=subset(int,inttype=="terttoprim"),aes(y=as.factor(intono-1), x=pctred_chr))+
  geom_density_ridges(scale=1,aes(fill=pct))+
  #theme_joy()+
  geom_point(data=subset(median,inttype=="terttoprim"),
             x=subset(median,inttype=="terttoprim")$med_pctredchr,
             y=subset(median,inttype=="terttoprim")$intono-0.92,
             by=subset(median,inttype=="terttoprim")$pct,
             size=1.5,shape=23,fill="white")+
  xlab("% Reduction: Chronic")+
  ylab("Intervention Sequence")+
  theme_bw(base_size=14)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', 
  fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_fill_cyclical(guide="legend",values = c("#CCCCCC", "#666666", "#000000"),
                      labels=c("10%","20%","40%"),name="Intervention %")+
  scale_y_discrete(breaks = 1:6, labels=c("Treat Former PWID","+Treat Current PWID","+Decrease Relapse  ",
                                          "+Increase Cessation ", "+Decrease Contacts ",
                                          "+Decrease Initiation  "))+ 
  scale_x_continuous(limits=c(-1,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

p3

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
  scale_y_discrete(breaks = 1:6, labels=c("Treat Former PWID","+Treat Current PWID","+Decrease Relapse  ",
                                          "+Increase Cessation ", "+Decrease Contacts ",
                                          "+Decrease Initiation  "))+ 
  scale_x_continuous(limits=c(-1,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
p4

p5<-ggplot(data=subset(int,inttype=="primtotert"),aes(y=as.factor(intono-1), x=pctred_chr))+
  geom_density_ridges(scale=1,aes(fill=pct))+
  #theme_joy()+
  geom_point(data=subset(median,inttype=="primtotert"),
             x=subset(median,inttype=="primtotert")$med_pctredchr,
             y=subset(median,inttype=="primtotert")$intono-0.92,
             by=subset(median,inttype=="primtotert")$pct,
             size=1.5,shape=23,fill="white")+
  xlab("% Reduction: Chronic")+
  ylab("Intervention Sequence")+
  theme_bw(base_size=14)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(legend.position="none", legend.background = element_rect(colour = 'NA', 
        fill = 'NA', size = 2, linetype="blank"))+
  theme(legend.key = element_blank())+
  scale_fill_cyclical(guide="legend",values = c("#CCCCCC", "#666666", "#000000"),
                      labels=c("10%","20%","40%"),name="Intervention %")+
  scale_y_discrete(breaks = 1:6, labels=c("Decrease Initiation  ","+Decrease Contacts ",
                                          "+Increase Cessation ", "+Decrease Relapse  ",
                                          "+Treat Current PWID","+Treat Former PWID"))+ 
  scale_x_continuous(limits=c(-1,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

p5

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
  scale_y_discrete(breaks = 1:6, labels=c("Decrease Initiation  ","+Decrease Contacts ",
                                          "+Increase Cessation ", "+Decrease Relapse  ",
                                          "+Treat Current PWID","+Treat Former PWID"))+ 
  scale_x_continuous(limits=c(-1,100),breaks=c(0,50,100),
                     labels=c("0","50","100"))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

p6

#no title
w <- 10; h <- 8

lay<-rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(5,5),c(6,6))
Multi_Sequential_ranges<-grid.arrange(p6,p4,p5,p3,legendpoint,legendtrt,layout_matrix=lay)
ggsave(sprintf("Multi_Sequential_Best10_grey.pdf"), Multi_Sequential_ranges, width=w, height=h)
#ggsave(sprintf("Multi_Sequential_Best10_Blue.pdf"), Multi_Sequential_ranges, width=w, height=h)


#### Sensitivity Analysis by Duration of Injecting (Cessation & Relapse) ####
#Set plot labels and colors
text.xlab <- "Intervention"
#model.colors4 <- c('#CCCCCC', '#999999', '#666666', '#000000','#FFFFFF') #grey scheme
model.colors <- c('#CCCCCC', '#666666', '#000000','#FFFFFF') #grey scheme
#model.colors <- c('#FF9999', '#CC3333','#990000','#FFFFFF') #red scheme
#model.colors <- c('#99CCFF','#3366FF', '#000066','#FFFFFF') #blue scheme

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

summary(1/subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%"))$gamma1)
summary(1/subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%"))$k1)

#isolate interventions of interest (single interventions)
data1<-subset(params,params$best10==1)
  
  
summary(data1$gamma1);summary(data1$gamma2);summary(data1$gamma3);summary(data1$gamma4)
summary(data1$k1);summary(data1$k2);summary(data1$k3);summary(data1$k4)

LowCess_HighRel<-subset(data1,data1$gamma1<=0.67421&data1$gamma2<=0.63125&data1$gamma3<=0.64151&data1$k1>0.1569604&data1$k2>0.371422&data1$k3>0.132667)
HighCess_LowRel<-subset(data1,data1$gamma1>0.67421&data1$gamma2>0.63125&data1$gamma3>0.64151&data1$k1<=0.1569604&data1$k2<=0.371422&data1$k3<=0.132667)

LowCess_HighRel_3064<-subset(data1,data1$gamma4<=0.55930&data1$k4>0.132667)
HighCess_LowRel_3064<-subset(data1,data1$gamma4>0.55930&data1$k4<=0.132667)

LowCess_HighRel_3064$duridu<-1/LowCess_HighRel_3064$gamma4
LowCess_HighRel_3064$durquit<-1/LowCess_HighRel_3064$k4
LowCess_HighRel_3064$LowCHighR<-1
colnames(LowCess_HighRel_3064)[colnames(LowCess_HighRel_3064)=="index"] <- "i"

HighCess_LowRel_3064$duridu<-1/HighCess_LowRel_3064$gamma4
summary(HighCess_LowRel_3064$duridu)
HighCess_LowRel_3064$durquit<-1/HighCess_LowRel_3064$k4
summary(HighCess_LowRel_3064$durquit)
HighCess_LowRel_3064$HighCLowR<-1
colnames(HighCess_LowRel_3064)[colnames(HighCess_LowRel_3064)=="index"] <- "i"

summary(params$b)

#subset only to interventions we want to plot
data1<-subset(singleint,singleint$best10==1&singleint$e<66&singleint$ODint==0&singleint$order %in% c(2:7)&singleint$IntLevelpct %in% c("10%","20%","40%"))

#merge with LowCHighRel and HighCLowRel
data2<-merge(HighCess_LowRel_3064[,c(1,11,28,95,96,97)],merge(data1,LowCess_HighRel_3064[,c(1,11,28,95,96,97)],"i",all=TRUE),"i",all=TRUE)
table(data2$HighCLowR)
table(data2$LowCHighR)

#legend for primary/secondary
leg<-ggplot(subset(data2,data2$HighCLowR==1),
            aes(x=order, y=pctred_chr_all, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Primary & Secondary Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Primary & Secondary Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=c("","","","","","",""))+ 
  scale_y_continuous(limits=c(-1,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=15) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
leg

legend1<-g_legend(leg)

#legend for treatment
leg2<-ggplot(subset(data2,data2$HighCLowR==1),
             aes(x=order, y=pctred_chr_all, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE, scale="width") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Treatment Rate', labels=c("0.03/Year ","0.06/Year","0.24/Year")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Treatment Rate', labels=c("0.03/Year ","0.06/Year","0.24/Year")) +
  scale_x_discrete(labels=c("","","","","","",""))+ 
  scale_y_continuous(limits=c(-1,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", position=position_dodge(width=0.95)) +
  theme_bw(base_size=15) +
  theme(legend.position = "bottom", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
leg2

legend2<-g_legend(leg2)

vp_highclowr<-ggplot(subset(data2,data2$HighCLowR==1),
                   aes(x=order, y=pctred_chr_all, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(-1,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", 
               position=position_dodge(width=0.8)) +
  theme_bw(base_size=15) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")+
  ggtitle("Short Duration of Injecting  (n=91 Parameter Sets)")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
vp_highclowr

vp_lowchighr<-ggplot(subset(data2,data2$LowCHighR==1),
                     aes(x=order, y=pctred_chr_all, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(-1,100),breaks=c(0,25,50,75,100),labels=c("0","25","50","75","100"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", 
               position=position_dodge(width=0.8)) +
  theme_bw(base_size=15) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab(text.ylab) +
  xlab("")+
  ggtitle("Long Duration of Injecting (n=91 Parameter Sets)")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
vp_lowchighr

vp_highclowr_a<-ggplot(subset(data2,data2$HighCLowR==1),
                     aes(x=order, y=pctred_acute_all, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(-13,100),breaks=c(-50,-25,0,25,50,75,100),labels=c("Increase 50%","Increase 25%","Null Impact","Reduce 25%"," Reduce 50%","Reduce 75%"," Reduce 100%"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", 
               position=position_dodge(width=0.8)) +
  theme_bw(base_size=15) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("% Reduction Acute HCV") +
  xlab("")+
  ggtitle("Short Duration of Injecting  (n=91 Parameter Sets)")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
vp_highclowr_a

vp_lowchighr_a<-ggplot(subset(data2,data2$LowCHighR==1),
                     aes(x=order, y=pctred_acute_all, by=IntLevelpct)) + 
  geom_violin(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.8), trim=TRUE, scale="count") +
  #geom_boxplot(aes(fill=as.factor(plotleg)), position=position_dodge(width=0.95), trim=TRUE) +
  scale_fill_manual(values=c('#CCCCCC', '#666666', '#000000'), name='Intervention Level', labels=c("10%","20%","40%")) +
  #scale_fill_manual(values=c('#99CCFF','#3366FF', '#000066'), name='Intervention Level', labels=c("10%","20%","40%")) +
  scale_x_discrete(labels=str_wrap(c("Decrease Injection Initiation",
                                     "Decrease PWID Contacts","Increase PWID Cessation","Decrease PWID Relapse",
                                     "Treat Current PWID","Treat Former PWID"),18))+ 
  scale_y_continuous(limits=c(-13,100),breaks=c(-50,-25,0,25,50,75,100),labels=c("Increase 50%","Increase 25%","Null Impact","Reduce 25%"," Reduce 50%","Reduce 75%"," Reduce 100%"))+
  stat_summary(aes(group=IntLevelpct), fun.y=median, geom="point", shape=23, size=2, fill= "white", 
               position=position_dodge(width=0.8)) +
  theme_bw(base_size=15) +
  theme(legend.position = "none", legend.key = element_rect(colour = 'black'),axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + 
  ylab("% Reduction Acute HCV") +
  xlab("")+
  ggtitle("Long Duration of Injecting (n=91 Parameter Sets)")+
  theme(plot.title = element_text(hjust = 0.5))
#+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) #adds margins around figure
vp_lowchighr_a


w <- 11; h <- 8
lay<-rbind(c(1), 
           c(1),  
           c(1),c(1),c(1), 
           c(2),
           c(2),
           c(2),c(2),c(2),
           c(3),c(4))
SingleStrategy_VP_main<-grid.arrange(vp_highclowr,vp_lowchighr,legend1,legend2,layout_matrix=lay)
#ggsave(sprintf("SingleStrategy_best10_30-64YO_blue.pdf"), SingleStrategy_VP_main, width=w, height=h)
ggsave(sprintf("SingleStrategy_best10_30-64YO_grey_chr_durinj.pdf"), SingleStrategy_VP_main, width=w, height=h)


SingleStrategy_VP_main<-grid.arrange(vp_highclowr_a,vp_lowchighr_a,legend1,legend2,layout_matrix=lay)
#ggsave(sprintf("SingleStrategy_best10_30-64YO_blue.pdf"), SingleStrategy_VP_main, width=w, height=h)
ggsave(sprintf("SingleStrategy_best10_30-64YO_grey_acute_durinj.pdf"), SingleStrategy_VP_main, width=w, height=h)


#### Summaries of Intervention % Reductions (to report in paper) ####
sub<-subset(singleint,singleint$e==62&singleint$best10==1) #trt former and current PWID at 3% per year
summary(sub$pctred_acute)
summary(sub$pctred_chr)
summary(sub$pctred_acute_all)
summary(sub$pctred_chr_all)

sub<-subset(singleint,singleint$e==22&singleint$best10==1) #trt former PWID at 3% per year
summary(sub$pctred_acute)
summary(sub$pctred_chr)
summary(sub$pctred_acute_all)
summary(sub$pctred_chr_all)

sub<-subset(singleint,singleint$e==6&singleint$best10==1) #Reduce contacts by 10%
summary(sub$pctred_acute)
summary(sub$pctred_chr)
summary(sub$pctred_acute_all)
summary(sub$pctred_chr_all)

sub<-subset(seqint,seqint$e==183&seqint$best10==1) #Reduce contacts, syr share, relapse, increase cess by 10%
summary(sub$pctred_acute)
summary(sub$pctred_chr)
summary(sub$pctred_acute_all)
summary(sub$pctred_chr_all)

sub<-subset(singleint,singleint$e==18&singleint$best10==1) #Trt SVR=90% for 3% trt level current PWID
summary(sub$pctred_acute)
summary(sub$pctred_chr)
summary(sub$pctred_acute_all)
summary(sub$pctred_chr_all)

sub<-subset(singleint,singleint$e==158&singleint$best10==1) #Trt SVR=60% for 3% trt level current PWID
summary(sub$pctred_acute)
summary(sub$pctred_chr)
summary(sub$pctred_acute_all)
summary(sub$pctred_chr_all)

sub<-subset(singleint,singleint$e==166&singleint$best10==1) #Trt SVR=60% for 3% trt level former and current PWID
summary(sub$pctred_acute)
summary(sub$pctred_chr)
summary(sub$pctred_acute_all)
summary(sub$pctred_chr_all)


