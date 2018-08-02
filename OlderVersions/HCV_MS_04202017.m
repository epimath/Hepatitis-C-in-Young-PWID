%% Main Script: Gicquelais et al. (2017). Hepatitis C tranmsission model.

% This code is the main script file to fit a hepatitis C transmission model
% to public health surveillance data from Michigan.

% Four parameters are estimated: 3 injection initiation rates (phi) and 1  
% transmission rate (beta). Latin hypercube sampling is used to
% incorporate parameter uncertainty over 5,000 simulations. A number of
% interventions are simulated. All results are exported as .csv or .txt
% files.

% Accompanying files:
    % Differential equations: HCV_DiffEq_04202017.m
    % Least squares for parameter estimation: HCV_LS_04202017.m
    % Contact matrix: Polymod.m

% For reference, the parameters are:

% Estimated parameters: 
    % param1=[b phi1 phi2 phi3]

% Sampled or set parameters: 
    % paramset=[d k eps s gn gp a Z0 r lamda1 lambda2 mu(1:3) etap etan intervention_etap];

%Param1 and paramset are combined into a single matrix called param:
    % param = [b phi1 phi2 phi3 d  k eps s  gn  gp  a  v0  r  lamda1 lambda2  mu(1:3) etap etan intervention_etap];
    % indices:[1   2    3   4   5  6  7  8  9   10  11 12  13  14       15     16-18   19   20           21];

% The inverse of age group size is in vector nu:
    % nu = [nu1,nu2,nu3];

% Compartments of the model are:
    % N[SNi; ANi; CNi; TNi; Si; Ai; Ci; Ti; Zi];
    % The first letter is the compartment type 
        % S=Susceptible [HCV Uninfected PWID], A=Acute HCV, C=Chronic HCV, 
        % T=Treated, Z=Non-PWID;
    % Former PWID have an 'N' following the compartment type while current
        % PWID and Non-PWID (Z) have no 'N';
    % i denotes the age class in years (1: 15-19, 2: 20-24, 3: 25-30;

%% HCV Surveillance Data, 2000-2013 (Michigan Department of Health and Human Services)

% Set the years data available
year = [2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013];
time = [0 1 2 3 4 5 6 7 8 9 10 11 12 13];

% Create matrices for HCV case counts by age group (15-19, 20-24, 25-30];

% All new chronic cases (regardless of history of injection drug use);
chronic1530 = [11	21	21	41	37	61	75	72	101	94	59	51	126	102;
               24	27	77	65	115	152	205	184	269	301	196	234	562	531;
               33	51	117	124	133	199	270	291	353	445	245	353	712	679]';

chronic1530total=sum(chronic1530,2);

% All cases of indeterminate state (acute or chronic);
unknown1530 = [10	8	1	0 2	6	6	1	12	1	0	2	2	5;
               15	25	8	0 4	9	6	9	22	3	4	6	2	3;
               33	56	7	0 4	8	9	10	31	4	2	9	2	3]';
           
unknown1530total=sum(unknown1530,2);

% Total new chronic and indeterminate cases;
chronicunk1530=chronic1530+unknown1530;
chronicunk1530total=sum(chronicunk1530,2);

% New chronic or indeterminate status cases (acute vs. chronic) who 
% reported yes, unknown, or missing for lifetime drug injection);
chronic1530idu = [21	29	21	41	37	60	77	61	99	89	55	48	114	97;
                  39	52	84	65	114	151	205	174	275	292	187	232	549	504;
                  66	107	124	124	134	195	265	278	356	422	233	342	675	641]';

chronic1530idutotal=sum(chronic1530idu,2);

% Total acute cases
acute1530 = [4	1	1	1	4	2	3	7	7	3	2	2	5	6;
             2	2	3	5	9	6	7	7	16	11	9	5	11	15;
             5	7	4	8	4	8	13	11	15	10	16	9	16	20]';
         
acute1530total=sum(acute1530,2);

% Total acute cases reporting yes, missing, or unknown injection of drugs
% in the 6 months before symptom onset;
acute1530idu = [4	1	1	1	4	1	2	7	6	3	2	2	5	3;
                2	2	3	5	7	3	6	7	16	10	8	5	9	12;
                5	7	4	8	3	6	13	8	11	9	10	7	10	14]';

acuteidutotal=sum(acute1530idu,2);

% Approximate the number of acute cases corrected for under-reporting using 
% Klevens RM et al. (2014). Estimating acute viral hepatitis infections from 
% nationally reported cases. American Journal of Public Health, 104(3):
% 482-487.

% Hepatitis C under-reporting factor (1/# cases per surveillance detected
% case, this estimate is between the PWID and overall rate for hepatitis C
% from Klevens et al.;
r=0.070412311;

% Multiply by acute cases observed to approximate the total number of acute
% cases;
acutecorr=(1/r)*acute1530; 
acutecorrtotal=sum(acutecorr,2);

acutecorridu=(1/r)*acute1530idu; 
acutecorrtotalidu=sum(acutecorridu,2);

% Create a dark green color
rgrn=[0 0.6 0.35];

% Plot case counts for each compartment at each time step;
figure(1)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year,chronic1530total+unknown1530total,'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,chronic1530idutotal,'Color',rgrn,'LineStyle','--','Marker','s','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,acute1530total,'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,acutecorrtotal,'rd--','LineWidth',4,'MarkerSize',10)
        ylim([0 1450])
        xlim([2000 2013])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Number of Cases')
        title('Figure 1. Michigan Hepatitis C Virus Infections Among 15-30 Year Olds -- 2000-2013')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('Total Chronic and Unconfirmed','Chronic and Unconfirmed PWID','Acute','Acute (Estimated Using CDC Correction)')

%% Age Contact Matrix
    
% Run the script entitled polymod.m to setup the contact matrices 
% used later (c, minc, and maxc).
run('Polymod.m')
%% Parameter Initial Values

% For reference, the parameters are:

% Estimated parameters: 
    % param1=[b phi1 phi2 phi3]

% Sampled or set parameters: 
    % paramset=[d k eps s gn gp a Z0 r lamda1 lambda2 mu(1:3) etap etan intervention_etap];

% Param1 and paramset are combined into a single matrix called param:
    % param = [b phi1 phi2 phi3 d  k eps s  gn  gp  a  Z0  r  lamda1 lambda2  mu(1:3) etap etan intervention_etap];
    % indices:[1   2    3   4   5  6  7  8  9   10  11 12  13  14       15     16-18   19   20         21];

% The inverse of age group size is in vector nu:
    % nu = [nu1,nu2,nu3];

%N[SNi; ANi; CNi; TNi; Si; Ai; Ci; Ti; Z]      


%param1 initial values (parameters we estimate);
    % Transmission rate
    b=0.000005;
    
    % Injection Initiation Rate Based on Steady State Values
    phi1=0.006885192;
    phi2=0.007301258;
    phi3=0.011996182;

param1=[b; phi1; phi2; phi3];
    
% paramset values (parameters we do not estimate);
    d=0.326; 
    k=0.45; 
    s=1;
    a=0.8;
    gn=0;
    gp=0;
    Z0=140166.428571429; 
    r=0.070412311;
    eps=0.011368375; 
    lambda1=0.013; 
    lambda2=0.024; 
    mu=[0.000565;0.000907;0.001139]; 
    etap=5.83; 
    etan=0.31; 
    intervention_etap=1;

paramset=[d; k; eps; s; gn; gp; a; Z0; r; lambda1; lambda2; mu(1:3); etap; etan; intervention_etap];

% combine param1 and paramset into 1 vector
param = [param1; paramset];

% nu = [nu1,nu2,nu3] as the inverse of age group size;
nu=[1/5;1/5;1/6];

%% Time Span and Age Group Population Sizes;

% Run model for time points specified;
tspan = [0 13];

% Average Michigan population size for each age group (from CDC Wonder);
popsize = [733933; 683425; 725476];

%% Initial Conditions

% Initial Condition Notes
    % Si and SNi: current PWID prevalence (Tempalski et al) and lifetime 
    % pwid prevalence (Lansky et al) less acute and chronic cases 
    
    % Ai and ANi: number of acute cases (with CDC correction) in 2000 
    % divided by 2 (acute HCV lasts 6 months);

    % Ci and CNi: (age group size * # chronic+unknown idu cases in 2000)
    % +80% of half-year of acute cases who develop chronic infection

    % Ti and TNi: No initial treated (appropriate for time period

    % Zi: Average popln size 2000-2013 for 15-30 year olds in MI (CDC Wonder)
    % less other groups 

% Param1 and paramset are combined into a single matrix called param:
    % param = [b phi1 phi2 phi3 d  k eps s  gn  gp  a  Z0  r  lamda1 lambda2  mu(1:3) etap etan intervention_deathcurr];
    % indices:[1   2    3   4   5  6  7  8  9   10  11 12  13  14       15     16-18   19   20  21];

% Determine proportion of PWID who are current PWID
propcurridu1824=param(7)/param(14);
propcurridu2534=param(7)/param(15);

% Set Initial Conditions;
N0=[0 0 0;
    0 0 0;
   (1/(r*nu(1,1)))*(1-propcurridu1824)*chronic1530idu(1,1) (1/(r*nu(2,1)))*(1-propcurridu1824)*chronic1530idu(1,2) (1/(r*nu(3,1)))*(1-propcurridu2534)*chronic1530idu(1,3);
    0 0 0;
    0 0 0;
    (1/(r*2*a))*chronic1530idu(1,1) (1/(r*2*a))*chronic1530idu(1,2) (1/(r*2*a))*chronic1530idu(1,3);
    (1/(r*nu(1,1)))*propcurridu1824*chronic1530idu(1,1) (1/(r*nu(2,1)))*propcurridu1824*chronic1530idu(1,2) (1/(r*nu(3,1)))*propcurridu2534*chronic1530idu(1,3);
    0 0 0;
    0 0 0];

% Define class Si (uninfected PWID) as popsize*eps less acutely
% and chronically infected persons;
for i=1:3
    N0(5,i) = popsize(i,1)*eps - (N0(6,i)+N0(7,i)+N0(8,i));
end

% Define class SNi (uninfected former IDU) as 
% popsize*(lifetime pwid-current pwid) less acutely
% and chronically infected persons;
N0(1,1)=popsize(1,1)*(lambda1-eps)-(N0(2,1)+N0(3,1)+N0(4,1));
N0(1,2)=popsize(2,1)*(lambda1-eps)-(N0(2,2)+N0(3,2)+N0(4,2));
N0(1,3)=popsize(3,1)*(lambda2-eps)-(N0(2,3)+N0(3,3)+N0(4,3));

% Define class Zi (uninfected non-IDU) as pop size less others defined in
% other classes;
for i=1:3
    N0(9,i) = popsize(i,1) - sum(N0(:,i));
end

% Reshape initial conditions into a vector;
N0=reshape(N0,27,1);

%% Test the ODE

[t5,x5] = ode15s(@HCV_DiffEq_04202017,tspan,N0,[],param,nu,c);

% Plot Initial Simulation Results (without parameter estimation) with Data;
for i=1:3;
    
    chronicprev=r*(x5(:,(3*i)+6*(i-1))+x5(:,(7*i)+2*(i-1)));
    acute=r*(x5(:,(2*i)+7*(i-1))+x5(:,(6*i)+3*(i-1)));
    chronicnew=r*(2*a*x5(:,(2*i)+7*(i-1))+2*a*x5(:,(6*i)+3*(i-1)));
    
figure(i+1);
hold on
plot(t5,chronicprev,'b','LineWidth',2); %chronic prevalence
hold on
plot(t5,chronicnew,'k','LineWidth',2); %new chronic PWID
hold on
plot(time,chronic1530idu(:,i),'ko','LineWidth',2) %chronic data
hold on
plot(t5,acute,'r','LineWidth',2) %acutes 
hold on
plot(time,acute1530idu(:,i),'r*','LineWidth',2)
hold on
title('Initial Conditions: HCV IDU Model')
hold on
legend('Prevalence of Chronic HCV: Model','New Chronic Cases: Model','New Chronic + Unconfirmed PWID Cases: Data','Total Acute Cases: Model','Acute (Estimated Using CDC Correction)')
end

%% Parameter Estimation (3 Injection Initiation [Phi], 1 Transmisison Rate [Beta])

y=zeros(1,3); 
[paramest RSS] = fminsearch(@(param1) HCV_LS_04202017(time,N0,param1,paramset,nu,c,chronic1530idu),param1,optimset('Display', 'iter'));
paramest=abs(paramest);

paramnew = [paramest;paramset];

RSS

% Calculate the AIC= 2(-LL+#parameters est) = 2(RSS+#parameters est)
AIC=2*(RSS+4)

%% Simulate the ODE system using the new fitted parameters

[t1,x1] = ode15s(@HCV_DiffEq_04202017,time,N0,[],paramnew,nu,c);

% Plots Results

% Plot model vs. all age groups combined;
figure(5)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*(x1(:,3)+x1(:,7)+x1(:,12)+x1(:,16)+x1(:,21)+x1(:,25)),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*(2*paramnew(11)*x1(:,2)+2*paramnew(11)*x1(:,6)+2*paramnew(11)*x1(:,11)+2*paramnew(11)*x1(:,15)+2*paramnew(11)*x1(:,20)+2*paramnew(11)*x1(:,24)),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,1)+chronic1530idu(:,2)+chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*(x1(:,2)+x1(:,6)+x1(:,11)+x1(:,15)+x1(:,20)+x1(:,24)),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,1)+acutecorridu(:,2)+acutecorridu(:,3),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,1)+acute1530idu(:,2)+acute1530idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 5. Model Fit: All Age Groups')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

%Age Group Plots

% 15-19 YO age group;
figure(6)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*(x1(:,3)+x1(:,7)),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*(2*paramnew(11)*x1(:,2)+2*paramnew(11)*x1(:,6)),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,1),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*(x1(:,2)+x1(:,6)),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,1),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,1),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 6. Model Fit: 15-19 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

% 20-24 YO age group;
figure(7)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*(x1(:,12)+x1(:,16)),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*(2*paramnew(11)*x1(:,11)+2*paramnew(11)*x1(:,15)),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,2),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*(x1(:,11)+x1(:,15)),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,2),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,2),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 7. Model Fit: 20-24 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

% 25-30 YO Age Group;
figure(8)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*(x1(:,21)+x1(:,25)),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*(2*paramnew(11)*x1(:,20)+2*paramnew(11)*x1(:,24)),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*(x1(:,20)+x1(:,24)),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,3),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 8. Model Fit: 25-30 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

%% Current PWID Compartments Only
% all age groups combined;
figure(9)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*(x1(:,7)+x1(:,16)+x1(:,25)),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*2*paramnew(11)*x1(:,6)+paramnew(13)*2*paramnew(11)*x1(:,15)+paramnew(13)*2*paramnew(11)*x1(:,24),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,1)+chronic1530idu(:,2)+chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*(x1(:,6)+x1(:,15)+x1(:,24)),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,1)+acutecorridu(:,2)+acutecorridu(:,3),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,1)+acute1530idu(:,2)+acute1530idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 9. Model Fit PWID: All Age Groups')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

% 15-19 YO age group;
figure(10)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*x1(:,7),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*2*paramnew(11)*x1(:,6),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,1),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*x1(:,6),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,1),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,1),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 10. Model Fit PWID: 15-19 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

% 20-24 YO age group;
figure(11)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*x1(:,16),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*2*paramnew(11)*x1(:,15),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,2),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*x1(:,15),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,2),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,2),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 11. Model Fit PWID: 20-24 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

% 25-30 YO Age Group;
figure(12)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*x1(:,25),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*2*paramnew(11)*x1(:,24),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*x1(:,24),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,3),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 12. Model Fit PWID: 25-30 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

%% Former PWID Compartments Only
% all age groups
figure(13)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*(x1(:,3)+x1(:,12)+x1(:,21)),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*2*paramnew(11)*x1(:,2)+paramnew(13)*2*paramnew(11)*x1(:,11)+paramnew(13)*2*paramnew(11)*x1(:,20),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,1)+chronic1530idu(:,2)+chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*(x1(:,2)+x1(:,11)+x1(:,20)),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,1)+acutecorridu(:,2)+acutecorridu(:,3),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,1)+acute1530idu(:,2)+acute1530idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 13. Model Fit Former PWID: All Age Groups')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

% 15-19 YO age group;
figure(14)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*x1(:,3),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*2*paramnew(11)*x1(:,2),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,1),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*x1(:,2),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,1),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,1),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 14. Model Fit Former PWID: 15-19 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

% 20-24 YO age group;
figure(15)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*x1(:,12),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*2*paramnew(11)*x1(:,11),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,2),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*x1(:,11),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,2),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,2),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 15. Model Fit Former PWID: 20-24 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')

% 25-30 YO Age Group;
figure(16)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
hold on
plot(t1,paramnew(13)*x1(:,21),'b','LineWidth',2);
hold on
plot(t1,paramnew(13)*2*paramnew(11)*x1(:,20),'k','LineWidth',2);
hold on
plot(time,chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
hold on
plot(t1,paramnew(13)*x1(:,20),'r','LineWidth',2)
hold on
plot(time,acutecorridu(:,3),'rd','LineWidth',2,'MarkerSize',8)
hold on
plot(time,acute1530idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
xlabel('Year')
ylabel('Number of Cases')
title('Figure 16. Model Fit Former PWID: 25-30 Age Group')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
legend('Chronic Prevalence PWID: Model','Chronic Incidence PWID: Model','New Chronic + Unconfirmed PWID Cases: Data','Acute PWID Cases: Model','Acute PWID (Estimated Using CDC Correction)','Acute PWID: Data')


%% Setup for Latin Hypercube Sample: Parameter Ranges and Scaling

% Reshape minc and maxc (from polymod.m)
minc=reshape(minc,1,9);
    
maxc=reshape(maxc,1,9);

% Set the Latin Hypercube Sample 
% lhsdesign(# simulations, # parameters sampled in paramranges)
% 5000 runs with 23 sampled parameters and 4 unsampled (intervention related)
paramscaling = lhsdesign(5000,27);

%Set the minimum and maximum sampling values for each parameter
% paramranges = [b phi1 phi2 phi3 d k eps s  gn gp a  Z0 r lambda1 lambda2 mu1:3 etap etan intervention_deathcurr c1 c2 c3 c5 c6 c9];
  paramranges = [paramest(1,1)-(0.5*paramest(1,1)) 0     0.00066 0     0.085 0.1 0.009 1 0 0 0.75 133231 1/16.8 0.009 0.014 0.000463 0.000835 0.001049 2.5  0.18 1 minc(1:3) minc(5:6) minc(9);
                 paramest(1,1)+(0.5*paramest(1,1)) 0.014 0.014   0.024 1.123 1   0.013 1 0 0 0.85 153080 1/12.3 0.016 0.034 0.000643 0.000992 0.001259 15.3 0.54 1 maxc(1:3) maxc(5:6) maxc(9)];

% time1 vector is used in recording data every 0.2 years;
time1 = [0:0.2:13];
%% Simulate the ODE Across the Latin Hypercube Sampled Space

% Setup empty matrices to record simulation results
NewChronicPWID = [];
NewChronicNPWID = [];
NewChronic = [];

NewChronicPWID1 = [];
NewChronicNPWID1 = [];
NewChronic1 = [];

NewChronicPWID2 = [];
NewChronicNPWID2 = [];
NewChronic2 = [];

NewChronicPWID3 = [];
NewChronicNPWID3 = [];
NewChronic3 = [];

ChronicPrevPWID = [];
ChronicPrevNPWID = [];
ChronicPrevRuns = [];

ChronicPrevPWID1 = [];
ChronicPrevNPWID1 = [];
ChronicPrevRuns1 = [];

ChronicPrevPWID2 = [];
ChronicPrevNPWID2 = [];
ChronicPrevRuns2 = [];

ChronicPrevPWID3 = [];
ChronicPrevNPWID3 = [];
ChronicPrevRuns3 = [];

AcutePWID = [];
AcuteNPWID = [];
AcuteRuns = [];

AcutePWID1 = [];
AcuteNPWID1 = [];
AcuteRuns1 = [];

AcutePWID2 = [];
AcuteNPWID2 = [];
AcuteRuns2 = [];

AcutePWID3 = [];
AcuteNPWID3 = [];
AcuteRuns3 = [];

SuscPWID = [];
SuscNPWID = [];
Susc = [];

SuscPWID1 = [];
SuscNPWID1 = [];
Susc1 = [];

SuscPWID2 = [];
SuscNPWID2 = [];
Susc2 = [];

SuscPWID3 = [];
SuscNPWID3 = [];
Susc3 = [];

RSSRuns = [];
ParamEstRuns=[];
CMatrixRuns = [];
InitialParamRuns=[];
    
% Loop to run model simulations;
for j = 1:length(paramscaling(:,1))

    % calculate set of params based on paramranges set above and lhs;
    fullparams = paramscaling(j,:).*(paramranges(2,:) - paramranges(1,:)) + paramranges(1,:);
    
    param1 = [fullparams(1,1:4)]';
    
    paramset = [fullparams(1,5:21)]';
    
    c = [fullparams(1,22:24); 
        0 fullparams(1,25:26);
        0 0 fullparams(1,27)];
    
    c = c+c'-diag(diag(c));
    
% param = [b phi1 phi2 phi3 d  k eps s  gn  gp  a  Z0  r  lambda1 lambda2 mu(1:3) etap etan intervention_etap];
%         [1 2    3    4    5  6  7  8  9   10  11 12  13   14       15    16-18   19   20   21];
% Number in paramset:       1  2  3  4  5   6   7  8   9    10       11    12-14   15   16   17];

% Define parameters required for initial conditions;
r=paramset(9);
a=paramset(7);
eps=paramset(3);
lambda1=paramset(10);
lambda2=paramset(11);
propcurridu1824=eps/lambda1;
propcurridu2534=eps/lambda2;

% Set Initial Conditions;
N0=[0 0 0;
    0 0 0;
   (1/(r*nu(1,1)))*(1-propcurridu1824)*chronic1530idu(1,1) (1/(r*nu(2,1)))*(1-propcurridu1824)*chronic1530idu(1,2) (1/(r*nu(3,1)))*(1-propcurridu2534)*chronic1530idu(1,3);
    0 0 0;
    0 0 0;
    (1/(r*2*a))*chronic1530idu(1,1) (1/(r*2*a))*chronic1530idu(1,2) (1/(r*2*a))*chronic1530idu(1,3);
    (1/(r*nu(1,1)))*propcurridu1824*chronic1530idu(1,1) (1/(r*nu(2,1)))*propcurridu1824*chronic1530idu(1,2) (1/(r*nu(3,1)))*propcurridu2534*chronic1530idu(1,3);
    0 0 0;
    0 0 0];

% Define class Si; 
for i=1:3
    N0(5,i) = popsize(i,1)*eps - (N0(6,i)+N0(7,i)+N0(8,i));
end

% Define class SNi;
N0(1,1)=popsize(1,1)*(lambda1-eps)-(N0(2,1)+N0(3,1)+N0(4,1));
N0(1,2)=popsize(2,1)*(lambda1-eps)-(N0(2,2)+N0(3,2)+N0(4,2));
N0(1,3)=popsize(3,1)*(lambda1-eps)-(N0(2,3)+N0(3,3)+N0(4,3));

% Define class Zi;
for i=1:3
    N0(9,i) = popsize(i,1) - sum(N0(:,i));
end

% Reshape initial conditions;
N0=reshape(N0,27,1);
    
% Parameter estimates and goodness of fit;
    [paramest RSS] = fminsearch(@(param1) HCV_LS_04202017(time,N0,param1,paramset,nu,c,chronic1530idu),param1,optimset('Display', 'iter'));
    
    paramest = abs(paramest);
    paramnew = [paramest;paramset];
    
% Run model with sampled and estimated parameters
       
    [t,x] = ode15s(@HCV_DiffEq_04202017,time1,N0,[],paramnew,nu,c);
            
% Save Output 

% Total Cases;
 Y_newchronic_pwid = paramnew(13)*(2*paramnew(11)*x(:,6)+2*paramnew(11)*x(:,15)+2*paramnew(11)*x(:,24));
 Y_newchronic_npwid = paramnew(13)*(2*paramnew(11)*x(:,2)+2*paramnew(11)*x(:,11)+2*paramnew(11)*x(:,20));
 Y_newchronic = paramnew(13)*(2*paramnew(11)*x(:,2)+2*paramnew(11)*x(:,6)+2*paramnew(11)*x(:,11)+2*paramnew(11)*x(:,15)+2*paramnew(11)*x(:,20)+2*paramnew(11)*x(:,24));

 Y_chronicprev_pwid = paramnew(13)*(x(:,7)+x(:,16)+x(:,25));
 Y_chronicprev_npwid = paramnew(13)*(x(:,3)+x(:,12)+x(:,21));
 Y_chronicprev = paramnew(13)*(x(:,3)+x(:,7)+x(:,12)+x(:,16)+x(:,21)+x(:,25));

 Y_acute_pwid = paramnew(13)*(x(:,6)+x(:,15)+x(:,24));
 Y_acute_npwid = paramnew(13)*(x(:,2)+x(:,11)+x(:,20));
 Y_acute = paramnew(13)*(x(:,2)+x(:,6)+x(:,11)+x(:,15)+x(:,20)+x(:,24));

 Y_S_pwid = x(:,5)+x(:,14)+x(:,23);
 Y_S_npwid = x(:,1)+x(:,10)+x(:,19);
 Y_S = x(:,5)+x(:,14)+x(:,23)+x(:,1)+x(:,10)+x(:,19);

 NewChronicPWID = [NewChronicPWID, Y_newchronic_pwid];
 NewChronicNPWID = [NewChronicNPWID, Y_newchronic_npwid];
 NewChronic = [NewChronic, Y_newchronic];

 ChronicPrevPWID = [ChronicPrevPWID, Y_chronicprev_pwid];
 ChronicPrevNPWID = [ChronicPrevNPWID, Y_chronicprev_npwid];
 ChronicPrevRuns = [ChronicPrevRuns, Y_chronicprev];

 AcutePWID = [AcutePWID, Y_acute_pwid];
 AcuteNPWID = [AcuteNPWID, Y_acute_npwid];
 AcuteRuns = [AcuteRuns, Y_acute];

 SuscPWID = [SuscPWID, Y_S_pwid];
 SuscNPWID = [SuscNPWID, Y_S_npwid];
 Susc = [Susc, Y_S];

% 15-19 Year Olds;
 Y_newchronic_pwid1 = paramnew(13)*(2*paramnew(11)*x(:,6));
 Y_newchronic_npwid1 = paramnew(13)*(2*paramnew(11)*x(:,2));
 Y_newchronic1 = paramnew(13)*(2*paramnew(11)*x(:,2)+2*paramnew(11)*x(:,6));

 Y_chronicprev_pwid1 = paramnew(13)*x(:,7);
 Y_chronicprev_npwid1 = paramnew(13)*x(:,3);
 Y_chronicprev1 = paramnew(13)*(x(:,3)+x(:,7));

 Y_acute_pwid1 = paramnew(13)*x(:,6);
 Y_acute_npwid1 = paramnew(13)*x(:,2);
 Y_acute1 = paramnew(13)*(x(:,2)+x(:,6));

 Y_S_pwid1 = x(:,5);
 Y_S_npwid1 = x(:,1);
 Y_S1 = x(:,5)+x(:,1);

 NewChronicPWID1 = [NewChronicPWID1, Y_newchronic_pwid1];
 NewChronicNPWID1 = [NewChronicNPWID1, Y_newchronic_npwid1];
 NewChronic1 = [NewChronic1, Y_newchronic1];

 ChronicPrevPWID1 = [ChronicPrevPWID1, Y_chronicprev_pwid1];
 ChronicPrevNPWID1 = [ChronicPrevNPWID1, Y_chronicprev_npwid1];
 ChronicPrevRuns1 = [ChronicPrevRuns1, Y_chronicprev1];

 AcutePWID1 = [AcutePWID1, Y_acute_pwid1];
 AcuteNPWID1 = [AcuteNPWID1, Y_acute_npwid1];
 AcuteRuns1 = [AcuteRuns1, Y_acute1];

 SuscPWID1 = [SuscPWID1, Y_S_pwid1];
 SuscNPWID1 = [SuscNPWID1, Y_S_npwid1];
 Susc1 = [Susc1, Y_S1];

% 20-24 Year Olds;
 Y_newchronic_pwid2 = paramnew(13)*(2*paramnew(11)*x(:,15));
 Y_newchronic_npwid2 = paramnew(13)*(2*paramnew(11)*x(:,11));
 Y_newchronic2 = paramnew(13)*(2*paramnew(11)*x(:,11)+2*paramnew(11)*x(:,15));

 Y_chronicprev_pwid2 = paramnew(13)*x(:,16);
 Y_chronicprev_npwid2 = paramnew(13)*x(:,12);
 Y_chronicprev2 = paramnew(13)*(x(:,12)+x(:,16));

 Y_acute_pwid2 = paramnew(13)*x(:,15);
 Y_acute_npwid2 = paramnew(13)*x(:,11);
 Y_acute2 = paramnew(13)*(x(:,11)+x(:,15));

 Y_S_pwid2 = x(:,14);
 Y_S_npwid2 = x(:,10);
 Y_S2 = x(:,14)+x(:,10);

 NewChronicPWID2 = [NewChronicPWID2, Y_newchronic_pwid2];
 NewChronicNPWID2 = [NewChronicNPWID2, Y_newchronic_npwid2];
 NewChronic2 = [NewChronic2, Y_newchronic2];

 ChronicPrevPWID2 = [ChronicPrevPWID2, Y_chronicprev_pwid2];
 ChronicPrevNPWID2 = [ChronicPrevNPWID2, Y_chronicprev_npwid2];
 ChronicPrevRuns2 = [ChronicPrevRuns2, Y_chronicprev2];

 AcutePWID2 = [AcutePWID2, Y_acute_pwid2];
 AcuteNPWID2 = [AcuteNPWID2, Y_acute_npwid2];
 AcuteRuns2 = [AcuteRuns2, Y_acute2];

 SuscPWID2 = [SuscPWID2, Y_S_pwid2];
 SuscNPWID2 = [SuscNPWID2, Y_S_npwid2];
 Susc2 = [Susc2, Y_S2];

% 25-30 Year Olds;
 Y_newchronic_pwid3 = paramnew(13)*(2*paramnew(11)*x(:,24));
 Y_newchronic_npwid3 = paramnew(13)*(2*paramnew(11)*x(:,20));
 Y_newchronic3 = paramnew(13)*(2*paramnew(11)*x(:,20)+2*paramnew(11)*x(:,24));

 Y_chronicprev_pwid3 = paramnew(13)*x(:,25);
 Y_chronicprev_npwid3 = paramnew(13)*x(:,21);
 Y_chronicprev3 = paramnew(13)*(x(:,21)+x(:,25));

 Y_acute_pwid3 = paramnew(13)*x(:,24);
 Y_acute_npwid3 = paramnew(13)*x(:,20);
 Y_acute3 = paramnew(13)*(x(:,20)+x(:,24));

 Y_S_pwid3 = x(:,23);
 Y_S_npwid3 = x(:,19);
 Y_S3 = x(:,23)+x(:,19);

 NewChronicPWID3 = [NewChronicPWID3, Y_newchronic_pwid3];
 NewChronicNPWID3 = [NewChronicNPWID3, Y_newchronic_npwid3];
 NewChronic3 = [NewChronic3, Y_newchronic3];

 ChronicPrevNPWID3 = [ChronicPrevNPWID3, Y_chronicprev_npwid3];
 ChronicPrevRuns3 = [ChronicPrevRuns3, Y_chronicprev3];

 AcutePWID3 = [AcutePWID3, Y_acute_pwid3];
 AcuteNPWID3 = [AcuteNPWID3, Y_acute_npwid3];
 AcuteRuns3 = [AcuteRuns3, Y_acute3];

 SuscPWID3 = [SuscPWID3, Y_S_pwid3];
 SuscNPWID3 = [SuscNPWID3, Y_S_npwid3];
 Susc3 = [Susc3, Y_S3];

 % RSS;
 RSSRuns = [RSSRuns, RSS];

 % Parameters;
 ParamEstRuns = [ParamEstRuns, paramnew];

 % Contact Matrix;
 CMatrixRuns = [CMatrixRuns; reshape(c,1,9)];

 %Initial Values of Estimated Params;
 InitialParamRuns=[InitialParamRuns, param1];

end

%Save the Workspace
save('HCVLHS04202017.mat')

%% Load the Latin Hypercube Sampling Results

load('HCVLHS04202017.mat')

%% Identify the Best Fitting Parameter Set

% output run index with minimum residual sum of squares
[minRSS,I]=min(RSSRuns)
Best=[minRSS,I];
Best=Best(1,2);

paramscaling(Best,:);
bestp=ParamEstRuns(:,Best);

[maxRSS,J]=max(RSSRuns)
medianRSS=median(RSSRuns)
meanRSS=mean(RSSRuns)
sdRSS=std(RSSRuns)


%% Parameter and Fit Figures

% RSS
figure(17)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(RSSRuns)
hold on
ylabel('Number of Runs')
hold on
xlabel('Residual Sum of Squares')
hold on
title('Residual Sum of Squares')

% Beta Estimates
figure(18)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(ParamEstRuns(1,:))
hold on
plot(ParamEstRuns(1,Best),1,'r*','LineWidth',2,'MarkerSize',8)
hold on
ylabel('Number of Runs')
hold on
xlabel('Estimate')
hold on
title('Transmission Rate')
legend('Beta','Best Fit')

% Phi1 Estimates
figure(19)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(ParamEstRuns(2,:))
hold on
plot(ParamEstRuns(2,Best),1,'r*','LineWidth',2,'MarkerSize',8)
hold on
ylabel('Number of Runs')
hold on
xlabel('Estimate')
hold on
title('Injection Initiation: 15-19 Year Old Age Class')
legend('phi1','Best Fit')

% Phi2 Estimates
figure(20)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(ParamEstRuns(3,:))
hold on
plot(ParamEstRuns(3,Best),1,'r*','LineWidth',2,'MarkerSize',8)
hold on
ylabel('Number of Runs')
hold on
xlabel('Estimate')
hold on
title('Injection Initiation: 20-24 Year Old Age Class')
legend('Phi2','Best Fit')

% Phi3 Estimates
figure(21)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(ParamEstRuns(4,:))
hold on
plot(ParamEstRuns(4,Best),1,'r*','LineWidth',2,'MarkerSize',8)
hold on
ylabel('Number of Runs')
hold on
xlabel('Estimate')
hold on
title('Injection Initiation: 25-30 Year Old Age Class')
legend('phi3','Best Fit')

%% LHS Case Count Figures: All Age Groups (Combined Results)

% Chronic and Acute Current and Former PWID;
figure(22)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic,'Color',[0.7 0.7 0.7],'LineWidth',2);
        hold on
    p2=plot(time1,NewChronic(:,Best),'k','LineWidth',2);
        hold on
    p3=plot(time1,AcuteRuns,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,AcuteRuns(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
    p5=plot(time,chronic1530idu(:,1)+chronic1530idu(:,2)+chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
        hold on
    p6=plot(time,acutecorr(:,1)+acutecorr(:,2)+acutecorr(:,3),'rd','LineWidth',2,'MarkerSize',8)
        hold on
    p7=plot(time,acute1530(:,1)+acute1530(:,2)+acute1530(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS: Acute and Chronic')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
%saveas(gcf,'AcuteNewChronic_All','fig')


% Chronic Current and Former PWID
figure(23)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic,'Color',[0.7 0.7 0.7],'LineWidth',2);
    hold on
    p2=plot(time1,NewChronic(:,Best),'k','LineWidth',2);
    hold on
    p3=plot(time,chronic1530idu(:,1)+chronic1530idu(:,2)+chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Fit to Chronic PWID')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','New Chronic + Unconfirmed PWID Cases: Data')
%saveas(gcf,'NewChronic_All','fig')

% Chronic Prevalence (Current and Former PWID)
figure(24)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,ChronicPrevRuns,'Color',[0.7 0.7 0.85],'LineWidth',2);
    hold on
    p2=plot(time1,ChronicPrevRuns(:,Best),'Color',[0 0 1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Chronic Prevalence')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
       legend([p1(1) p2(1)], 'Chronic Prevalence: Model','Chronic Prevalence: Best Fit')
%saveas(gcf,'ChronicPrev_All','fig')

% Susceptible Current and Former PWID
figure(25)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,Susc,'Color',[0.8 0.9 0.8],'LineWidth',2);
    hold on
    p2=plot(time1,Susc(:,Best),'Color',[0 0.4 0.1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Susceptible PWID')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1)], 'Susceptible PWID: Model','Susceptible PWID: Best Fit')
%saveas(gcf,'Susc_All','fig')

%% LHS Case Count Figures: 15-19 Age Group

% Chronic and Acute PWID;
figure(26)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic1,'Color',[0.7 0.7 0.7],'LineWidth',2);
        hold on
    p2=plot(time1,NewChronic1(:,Best),'k','LineWidth',2);
        hold on
    p3=plot(time1,AcuteRuns1,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,AcuteRuns1(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
    p5=plot(time,chronic1530idu(:,1),'ko','LineWidth',2,'MarkerSize',8)
        hold on
    p6=plot(time,acutecorr(:,1),'rd','LineWidth',2,'MarkerSize',8)
        hold on
    p7=plot(time,acute1530(:,1),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS: Acute and Chronic 15-19 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
%saveas(gcf,'AcuteNewChronic_1519','fig')

% Chronic PWID
figure(27)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic1,'Color',[0.7 0.7 0.7],'LineWidth',2);
    hold on
    p2=plot(time1,NewChronic1(:,Best),'k','LineWidth',2);
    hold on
    p3=plot(time,chronic1530idu(:,1),'ko','LineWidth',2,'MarkerSize',8)
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Fit to Chronic PWID 15-19 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','New Chronic + Unconfirmed PWID Cases: Data')
%saveas(gcf,'NewChronic_1519','fig')

% Chronic Prevalence
figure(28)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,ChronicPrevRuns1,'Color',[0.7 0.7 0.85],'LineWidth',2);
    hold on
    p2=plot(time1,ChronicPrevRuns1(:,424),'Color',[0 0 1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Chronic Prevalence 15-19 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
       legend([p1(1) p2(1)], 'Chronic Prevalence: Model','Chronic Prevalence: Best Fit')
%saveas(gcf,'ChronicPrev_1519','fig')

% Susceptible PWID
figure(29)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,Susc1,'Color',[0.8 0.9 0.8],'LineWidth',2);
    hold on
    p2=plot(time1,Susc1(:,Best),'Color',[0 0.4 0.1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Susceptible PWID 15-19 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1)], 'Susceptible PWID: Model','Susceptible PWID: Best Fit')
%saveas(gcf,'Susc_1519','fig')

%% LHS Case Count Figures: 20-24 Age Group

% Chronic and Acute PWID;
figure(30)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic2,'Color',[0.7 0.7 0.7],'LineWidth',2);
        hold on
    p2=plot(time1,NewChronic2(:,Best),'k','LineWidth',2);
        hold on
    p3=plot(time1,AcuteRuns2,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,AcuteRuns2(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
    p5=plot(time,chronic1530idu(:,2),'ko','LineWidth',2,'MarkerSize',8)
        hold on
    p6=plot(time,acutecorr(:,2),'rd','LineWidth',2,'MarkerSize',8)
        hold on
    p7=plot(time,acute1530(:,2),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS: Acute and Chronic 20-24 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
%saveas(gcf,'AcuteNewChronic_2024','fig')


% Chronic PWID
figure(31)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic2,'Color',[0.7 0.7 0.7],'LineWidth',2);
    hold on
    p2=plot(time1,NewChronic2(:,Best),'k','LineWidth',2);
    hold on
    p3=plot(time,chronic1530idu(:,2),'ko','LineWidth',2,'MarkerSize',8)
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Fit to Chronic PWID 20-24 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','New Chronic + Unconfirmed PWID Cases: Data')
%saveas(gcf,'NewChronic_2024','fig')

% Chronic Prevalence
figure(32)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,ChronicPrevRuns2,'Color',[0.7 0.7 0.85],'LineWidth',2);
    hold on
    p2=plot(time1,ChronicPrevRuns2(:,Best),'Color',[0 0 1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Chronic Prevalence 20-24 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
       legend([p1(1) p2(1)], 'Chronic Prevalence: Model','Chronic Prevalence: Best Fit')
%saveas(gcf,'ChronicPrev_2024','fig')

% Susceptible PWID
figure(33)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,Susc2,'Color',[0.8 0.9 0.8],'LineWidth',2);
    hold on
    p2=plot(time1,Susc2(:,Best),'Color',[0 0.4 0.1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Susceptible PWID 20-24 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1)], 'Susceptible PWID: Model','Susceptible PWID: Best Fit')
%saveas(gcf,'Susc_2024','fig')

%% LHS Case Count Figures: 25-30 Age Group

% Chronic and Acute PWID;
figure(34)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic3,'Color',[0.7 0.7 0.7],'LineWidth',2);
        hold on
    p2=plot(time1,NewChronic3(:,Best),'k','LineWidth',2);
        hold on
    p3=plot(time1,AcuteRuns3,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,AcuteRuns3(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
    p5=plot(time,chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
        hold on
    p6=plot(time,acutecorr(:,3),'rd','LineWidth',2,'MarkerSize',8)
        hold on
    p7=plot(time,acute1530(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS: Acute and Chronic 25-30 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
%saveas(gcf,'AcuteNewChronic_2530','fig')

% Chronic PWID
figure(35)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic3,'Color',[0.7 0.7 0.7],'LineWidth',2);
    hold on
    p2=plot(time1,NewChronic3(:,Best),'k','LineWidth',2);
    hold on
    p3=plot(time,chronic1530idu(:,3),'ko','LineWidth',2,'MarkerSize',8)
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Fit to Chronic PWID 25-30 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','New Chronic + Unconfirmed PWID Cases: Data')
%saveas(gcf,'NewChronic_2530','fig')

% Chronic Prevalence
figure(36)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,ChronicPrevRuns3,'Color',[0.7 0.7 0.85],'LineWidth',2);
    hold on
    p2=plot(time1,ChronicPrevRuns3(:,Best),'Color',[0 0 1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Chronic Prevalence 25-30 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
       legend([p1(1) p2(1)], 'Chronic Prevalence: Model','Chronic Prevalence: Best Fit')
%saveas(gcf,'ChronicPrev_2530','fig')

% Susceptible PWID
figure(37)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,Susc3,'Color',[0.8 0.9 0.8],'LineWidth',2);
    hold on
    p2=plot(time1,Susc3(:,Best),'Color',[0 0.4 0.1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 13])
    ylabel('Number of Cases')
    title('LHS Susceptible PWID 25-30 Age Class')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1)], 'Susceptible PWID: Model','Susceptible PWID: Best Fit')
%saveas(gcf,'Susc_2530','fig')

%% Compile parameters, contact matrix, RSS, and initial parameters 
LHS_ParamRSS = [];
LHS_ParamRSS(:,1) = 1:5000;
LHS_ParamRSS(:,2:22) = ParamEstRuns';
LHS_ParamRSS(:,23:31) = CMatrixRuns;
LHS_ParamRSS(:,32) = RSSRuns';
LHS_ParamRSS(:,33:36) = InitialParamRuns';

%% Save Results as TXT or CSV files

dlmwrite('LHS_ParamRSS.txt', [LHS_ParamRSS], 'precision', 10); 

dlmwrite('LHS_AcuteNPWID.csv',[LHS_ParamRSS(:,1)'; AcuteNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcuteNPWID1.csv',[LHS_ParamRSS(:,1)';AcuteNPWID1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcuteNPWID2.csv',[LHS_ParamRSS(:,1)';AcuteNPWID2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcuteNPWID3.csv',[LHS_ParamRSS(:,1)';AcuteNPWID3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_AcutePWID.csv',[LHS_ParamRSS(:,1)';AcutePWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcutePWID1.csv',[LHS_ParamRSS(:,1)';AcutePWID1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcutePWID2.csv',[LHS_ParamRSS(:,1)';AcutePWID2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcutePWID3.csv',[LHS_ParamRSS(:,1)';AcutePWID3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_Acute.csv',[LHS_ParamRSS(:,1)';AcuteRuns], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Acute1.csv',[LHS_ParamRSS(:,1)';AcuteRuns1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Acute2.csv',[LHS_ParamRSS(:,1)';AcuteRuns2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Acute3.csv',[LHS_ParamRSS(:,1)';AcuteRuns3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_ChronicPrevNPWID.csv',[LHS_ParamRSS(:,1)';ChronicPrevNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevNPWID1.csv',[LHS_ParamRSS(:,1)';ChronicPrevNPWID1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevNPWID2.csv',[LHS_ParamRSS(:,1)';ChronicPrevNPWID2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevNPWID3.csv',[LHS_ParamRSS(:,1)';ChronicPrevNPWID3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_ChronicPrevPWID.csv',[LHS_ParamRSS(:,1)';ChronicPrevPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevPWID1.csv',[LHS_ParamRSS(:,1)';ChronicPrevPWID1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevPWID2.csv',[LHS_ParamRSS(:,1)';ChronicPrevPWID2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevPWID3.csv',[LHS_ParamRSS(:,1)';ChronicPrevPWID3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_ChronicPrev.csv',[LHS_ParamRSS(:,1)';ChronicPrevRuns], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrev1.csv',[LHS_ParamRSS(:,1)';ChronicPrevRuns1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrev2.csv',[LHS_ParamRSS(:,1)';ChronicPrevRuns2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrev3.csv',[LHS_ParamRSS(:,1)';ChronicPrevRuns3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_NewChronic.csv',[LHS_ParamRSS(:,1)';NewChronic], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronic1.csv',[LHS_ParamRSS(:,1)';NewChronic1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronic2.csv',[LHS_ParamRSS(:,1)';NewChronic2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronic3.csv',[LHS_ParamRSS(:,1)';NewChronic3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_NewChronicNPWID.csv',[LHS_ParamRSS(:,1)';NewChronicNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicNPWID1.csv',[LHS_ParamRSS(:,1)';NewChronicNPWID1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicNPWID2.csv',[LHS_ParamRSS(:,1)';NewChronicNPWID2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicNPWID3.csv',[LHS_ParamRSS(:,1)';NewChronicNPWID3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_NewChronicPWID.csv',[LHS_ParamRSS(:,1)';NewChronicPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicPWID1.csv',[LHS_ParamRSS(:,1)';NewChronicPWID1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicPWID2.csv',[LHS_ParamRSS(:,1)';NewChronicPWID2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicPWID3.csv',[LHS_ParamRSS(:,1)';NewChronicPWID3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_Susc.csv',[LHS_ParamRSS(:,1)';Susc], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Susc1.csv',[LHS_ParamRSS(:,1)';Susc1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Susc2.csv',[LHS_ParamRSS(:,1)';Susc2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Susc3.csv',[LHS_ParamRSS(:,1)';Susc3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_SuscNPWID.csv',[LHS_ParamRSS(:,1)';SuscNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscNPWID1.csv',[LHS_ParamRSS(:,1)';SuscNPWID1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscNPWID2.csv',[LHS_ParamRSS(:,1)';SuscNPWID2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscNPWID3.csv',[LHS_ParamRSS(:,1)';SuscNPWID3], 'delimiter', ',', 'precision', 10);
dlmwrite('LHS_SuscPWID.csv',[LHS_ParamRSS(:,1)';SuscPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscPWID1.csv',[LHS_ParamRSS(:,1)';SuscPWID1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscPWID2.csv',[LHS_ParamRSS(:,1)';SuscPWID2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscPWID3.csv',[LHS_ParamRSS(:,1)';SuscPWID3], 'delimiter', ',', 'precision', 10);
                                     

MDHHSData = [];
MDHHSData(:,1) = year'; 
MDHHSData(:,2:4) = acute1530;
MDHHSData(:,5:7) = acute1530idu; 
MDHHSData(:,8:10) = chronic1530; 
MDHHSData(:,11:13) = chronic1530idu;
MDHHSData(:,14:16) = unknown1530;
MDHHSData(:,17:19) = chronicunk1530;
dlmwrite('MDHHSData.csv',[MDHHSData], 'delimiter', ',', 'precision', 10);

%% Setup Subsets of Interventions

load('HCVLHS04202017.mat')

% Set intervention levels

% Treatment in former PWID [index j];
intervention_gn = [0 0.1 0.2 0.3 0.4];
 
% Treatment in current PWID [index l];
intervention_gp = [0 0.1 0.2 0.3 0.4];
 
% 1/duration of treatment (1 year, 24 weeks, 16 weeks, 12 weeks) [index m];
intervention_s = [1 1/(168/365) 1/(112/365) 1/(84/365)]; 
 
% Relapse rate (k) [index o];
intervention_k = [1 0.9 0.8 0.7 0.6];
 
% Cessation rate (d) [index q];
intervention_d = [1 1.1 1.2 1.3 1.4];
 
% Injection initiation (phi) [index u];
intervention_phi = [1 0.9 0.8 0.7 0.6];
 
% Contact matrix [index f];
intervention_contact = [1 0.9 0.8 0.7 0.6];

% Mortality rate: current PWID [index h];
intervention_etap = [1 0.9 0.8 0.7 0.6];

% Mortality rate: former PWID [index k];
intervention_etan = [1 0.9 0.8 0.7 0.6];

% Setup intervention combinations we want results from using the
    % intervention indices from above. These are also listed in
    % listed in file: InterventionCombos.xlsx, which assigns an index e
    % for single and sequential runs;
    
               %j   l   m   o   q   u   f   h   k    
interventionrun=[1	1	1	1	1	1	1	1	1	;
                2	1	1	1	1	1	1	1	1	;
                3	1	1	1	1	1	1	1	1	;
                4	1	1	1	1	1	1	1	1	;
                5	1	1	1	1	1	1	1	1	;
                1	2	1	1	1	1	1	1	1	;
                1	3	1	1	1	1	1	1	1	;
                1	4	1	1	1	1	1	1	1	;
                1	5	1	1	1	1	1	1	1	;
                1	1	2	1	1	1	1	1	1	;
                1	1	3	1	1	1	1	1	1	;
                1	1	4	1	1	1	1	1	1	;
                1	1	1	2	1	1	1	1	1	;
                1	1	1	3	1	1	1	1	1	;
                1	1	1	4	1	1	1	1	1	;
                1	1	1	5	1	1	1	1	1	;
                1	1	1	1	2	1	1	1	1	;
                1	1	1	1	3	1	1	1	1	;
                1	1	1	1	4	1	1	1	1	;
                1	1	1	1	5	1	1	1	1	;
                1	1	1	1	1	2	1	1	1	;
                1	1	1	1	1	3	1	1	1	;
                1	1	1	1	1	4	1	1	1	;
                1	1	1	1	1	5	1	1	1	;
                1	1	1	1	1	1	2	1	1	;
                1	1	1	1	1	1	3	1	1	;
                1	1	1	1	1	1	4	1	1	;
                1	1	1	1	1	1	5	1	1	;
                1	1	1	1	1	1	1	2	1	;
                1	1	1	1	1	1	1	3	1	;
                1	1	1	1	1	1	1	4	1	;
                1	1	1	1	1	1	1	5	1	;
                1	1	1	1	1	1	1	1	2	;
                1	1	1	1	1	1	1	1	3	;
                1	1	1	1	1	1	1	1	4	;
                1	1	1	1	1	1	1	1	5	;
                2	1	4	1	1	1	1	1	1	;
                3	1	4	1	1	1	1	1	1	;
                4	1	4	1	1	1	1	1	1	;
                5	1	4	1	1	1	1	1	1	;
                1	2	4	1	1	1	1	1	1	;
                1	3	4	1	1	1	1	1	1	;
                1	4	4	1	1	1	1	1	1	;
                1	5	4	1	1	1	1	1	1	;
                1	1	4	2	1	1	1	1	1	;
                1	1	4	3	1	1	1	1	1	;
                1	1	4	4	1	1	1	1	1	;
                1	1	4	5	1	1	1	1	1	;
                1	1	4	1	2	1	1	1	1	;
                1	1	4	1	3	1	1	1	1	;
                1	1	4	1	4	1	1	1	1	;
                1	1	4	1	5	1	1	1	1	;
                1	1	4	1	1	2	1	1	1	;
                1	1	4	1	1	3	1	1	1	;
                1	1	4	1	1	4	1	1	1	;
                1	1	4	1	1	5	1	1	1	;
                1	1	4	1	1	1	2	1	1	;
                1	1	4	1	1	1	3	1	1	;
                1	1	4	1	1	1	4	1	1	;
                1	1	4	1	1	1	5	1	1	;
                1	1	4	1	1	1	1	2	1	;
                1	1	4	1	1	1	1	3	1	;
                1	1	4	1	1	1	1	4	1	;
                1	1	4	1	1	1	1	5	1	;
                1	1	4	1	1	1	1	1	2	;
                1	1	4	1	1	1	1	1	3	;
                1	1	4	1	1	1	1	1	4	;
                1	1	4	1	1	1	1	1	5	;
                1	1	4	1	1	1	1	1	1	;
                1	1	4	1	1	2	1	1	1	;
                1	1	4	1	1	3	1	1	1	;
                1	1	4	1	1	4	1	1	1	;
                1	1	4	1	1	5	1	1	1	;
                1	1	4	1	1	2	2	1	1	;
                1	1	4	1	1	3	3	1	1	;
                1	1	4	1	1	4	4	1	1	;
                1	1	4	1	1	5	5	1	1	;
                1	1	4	1	2	2	2	1	1	;
                1	1	4	1	3	3	3	1	1	;
                1	1	4	1	4	4	4	1	1	;
                1	1	4	1	5	5	5	1	1	;
                1	1	4	2	2	2	2	1	1	;
                1	1	4	3	3	3	3	1	1	;
                1	1	4	4	4	4	4	1	1	;
                1	1	4	5	5	5	5	1	1	;
                1	2	4	2	2	2	2	1	1	;
                1	3	4	3	3	3	3	1	1	;
                1	4	4	4	4	4	4	1	1	;
                1	5	4	5	5	5	5	1	1	;
                2	2	4	2	2	2	2	1	1	;
                3	3	4	3	3	3	3	1	1	;
                4	4	4	4	4	4	4	1	1	;
                5	5	4	5	5	5	5	1	1	;
                1	1	4	1	1	1	1	5	5	;
                1	1	4	1	1	2	1	5	5	;
                1	1	4	1	1	3	1	5	5	;
                1	1	4	1	1	4	1	5	5	;
                1	1	4	1	1	5	1	5	5	;
                1	1	4	1	1	2	2	5	5	;
                1	1	4	1	1	3	3	5	5	;
                1	1	4	1	1	4	4	5	5	;
                1	1	4	1	1	5	5	5	5	;
                1	1	4	1	2	2	2	5	5	;
                1	1	4	1	3	3	3	5	5	;
                1	1	4	1	4	4	4	5	5	;
                1	1	4	1	5	5	5	5	5	;
                1	1	4	2	2	2	2	5	5	;
                1	1	4	3	3	3	3	5	5	;
                1	1	4	4	4	4	4	5	5	;
                1	1	4	5	5	5	5	5	5	;
                1	2	4	2	2	2	2	5	5	;
                1	3	4	3	3	3	3	5	5	;
                1	4	4	4	4	4	4	5	5	;
                1	5	4	5	5	5	5	5	5	;
                2	2	4	2	2	2	2	5	5	;
                3	3	4	3	3	3	3	5	5	;
                4	4	4	4	4	4	4	5	5	;
                5	5	4	5	5	5	5	5	5	;
                1	1	4	1	1	1	1	1	1	;
                2	1	4	1	1	1	1	1	1	;
                3	1	4	1	1	1	1	1	1	;
                4	1	4	1	1	1	1	1	1	;
                5	1	4	1	1	1	1	1	1	;
                2	2	4	1	1	1	1	1	1	;
                3	3	4	1	1	1	1	1	1	;
                4	4	4	1	1	1	1	1	1	;
                5	5	4	1	1	1	1	1	1	;
                2	2	4	2	1	1	1	1	1	;
                3	3	4	3	1	1	1	1	1	;
                4	4	4	4	1	1	1	1	1	;
                5	5	4	5	1	1	1	1	1	;
                2	2	4	2	2	1	1	1	1	;
                3	3	4	3	3	1	1	1	1	;
                4	4	4	4	4	1	1	1	1	;
                5	5	4	5	5	1	1	1	1	;
                2	2	4	2	2	1	2	1	1	;
                3	3	4	3	3	1	3	1	1	;
                4	4	4	4	4	1	4	1	1	;
                5	5	4	5	5	1	5	1	1	;
                2	2	4	2	2	2	2	1	1	;
                3	3	4	3	3	3	3	1	1	;
                4	4	4	4	4	4	4	1	1	;
                5	5	4	5	5	5	5	1	1	;
                1	1	4	1	1	1	1	5	5	;
                2	1	4	1	1	1	1	5	5	;
                3	1	4	1	1	1	1	5	5	;
                4	1	4	1	1	1	1	5	5	;
                5	1	4	1	1	1	1	5	5	;
                2	2	4	1	1	1	1	5	5	;
                3	3	4	1	1	1	1	5	5	;
                4	4	4	1	1	1	1	5	5	;
                5	5	4	1	1	1	1	5	5	;
                2	2	4	2	1	1	1	5	5	;
                3	3	4	3	1	1	1	5	5	;
                4	4	4	4	1	1	1	5	5	;
                5	5	4	5	1	1	1	5	5	;
                2	2	4	2	2	1	1	5	5	;
                3	3	4	3	3	1	1	5	5	;
                4	4	4	4	4	1	1	5	5	;
                5	5	4	5	5	1	1	5	5	;
                2	2	4	2	2	1	2	5	5	;
                3	3	4	3	3	1	3	5	5	;
                4	4	4	4	4	1	4	5	5	;
                5	5	4	5	5	1	5	5	5	;
                2	2	4	2	2	2	2	5	5	;
                3	3	4	3	3	3	3	5	5	;
                4	4	4	4	4	4	4	5	5	;
                5	5	4	5	5	5	5	5	5	];

%% Clear Unnecessary parts of workspace and Save
clearvars -except ParamEstRuns CMatrixRuns RSSRuns chronic1530idu nu popsize time1 interventionrun intervention_gn intervention_gp intervention_s intervention_k intervention_d intervention_phi intervention_contact intervention_deathcurr intervention_deathform

save('HCV_IntSetup.mat')

%% Simulate Interventions

% Interventions are simulated in groups of 500 parameter sets at a time
for indexnum=1:10;
    
    indexmin=((indexnum-1)*500)+1; 
    indexmax=indexnum*500;

% Load intervention workspace    
load('HCV_IntSetup.mat') 

% Setup matrices to save results
chronicprev_results=[];
newchronic_results=[];
acute_results=[];
parameters_results=[];

% Loop over each intervention
for e=1:168;

intlev=interventionrun(e,1:9); 

% Assign the appropriate level to each intervention parameter  
j=intlev(1); 
l=intlev(2);
m=intlev(3);
o=intlev(4);
q=intlev(5);
u=intlev(6);
f=intlev(7);
h=intlev(8);
k=intlev(9);

% Select which block of 500 parameter sets will be simulated
for i=indexmin:indexmax;

% Select a parameter set i 
param=ParamEstRuns(:,i);

% Reshape contact matrix as 3x3 and scale by intervention
c=CMatrixRuns(i,:);
c=reshape(c,3,3)*intervention_contact(f); 
 
% Param1 and paramset are combined into a single matrix called param:
    % param = [b phi1 phi2 phi3 d  k eps s  gn  gp  a  v0  r  lamda1 lambda2  mu(1:3) etap etan intervention_etap];
    % indices:[1   2    3   4   5  6  7  8  9   10  11 12  13  14       15     16-18   19   20           21];

% Assign intervention parameters
% Note: etan and etap are further transformed in the ODE file
param(2)=intervention_phi(u)*ParamEstRuns(2,i);
param(3)=intervention_phi(u)*ParamEstRuns(3,i);
param(4)=intervention_phi(u)*ParamEstRuns(4,i);
param(5)=intervention_d(q)*ParamEstRuns(5,i);
param(6)=intervention_k(o)*ParamEstRuns(6,i);
param(8)=intervention_s(m);
param(9)=intervention_gn(j);
param(10)=intervention_gp(l);
param(19)=ParamEstRuns(19,i);
param(20)=intervention_etan(k)*ParamEstRuns(20,i);
param(21)=intervention_etap(h);

% Set initial conditions
eps=param(7);
a=param(11);
r=param(13);
lambda1=param(14);
lambda2=param(15);
propcurridu1824=eps/lambda1;
propcurridu2534=eps/lambda2;

%Initial Conditions;
N0=[0 0 0;
    0 0 0;
   (1/(r*nu(1,1)))*(1-propcurridu1824)*chronic1530idu(1,1) (1/(r*nu(2,1)))*(1-propcurridu1824)*chronic1530idu(1,2) (1/(r*nu(3,1)))*(1-propcurridu2534)*chronic1530idu(1,3);
    0 0 0;
    0 0 0;
    (1/(r*2*a))*chronic1530idu(1,1) (1/(r*2*a))*chronic1530idu(1,2) (1/(r*2*a))*chronic1530idu(1,3);
    (1/(r*nu(1,1)))*propcurridu1824*chronic1530idu(1,1) (1/(r*nu(2,1)))*propcurridu1824*chronic1530idu(1,2) (1/(r*nu(3,1)))*propcurridu2534*chronic1530idu(1,3);
    0 0 0;
    0 0 0];

%Define class Si;
for s=1:3
    N0(5,s) = popsize(s,1)*eps - (N0(6,s)+N0(7,s)+N0(8,s));
end

%Define class SNi;
N0(1,1)=popsize(1,1)*(lambda1-eps)-(N0(2,1)+N0(3,1)+N0(4,1));
N0(1,2)=popsize(2,1)*(lambda1-eps)-(N0(2,2)+N0(3,2)+N0(4,2));
N0(1,3)=popsize(3,1)*(lambda2-eps)-(N0(2,3)+N0(3,3)+N0(4,3));

%Defined class Zi;
for s=1:3
    N0(9,s) = popsize(s,1) - sum(N0(:,s));
end

% Reshape initial conditions;
N0=reshape(N0,27,1);

% Run ODE for ith parameter set
[t,x] = ode15s(@HCV_DiffEq_04202017,time1,N0,[],param,nu,c);
 
% Collect total (all age) data from each run
    Y_newchronic = param(13)*(2*param(11)*x(:,2)+2*param(11)*x(:,6)+2*param(11)*x(:,11)+2*param(11)*x(:,15)+2*param(11)*x(:,20)+2*param(11)*x(:,24));
    
    Y_chronicprev = param(13)*(x(:,3)+x(:,7)+x(:,12)+x(:,16)+x(:,21)+x(:,25));
    
    Y_acute = param(13)*(x(:,2)+x(:,6)+x(:,11)+x(:,15)+x(:,20)+x(:,24));

% Add results to results matrices   
chronicprev_results = [chronicprev_results; e,j,l,m,o,q,u,f,h,k,i,Y_chronicprev'];
newchronic_results = [newchronic_results; e,j,l,m,o,q,u,f,h,k,i,Y_newchronic'];
acute_results = [acute_results; e,j,l,m,o,q,u,f,h,k,i,Y_acute'];
parameters_results = [parameters_results; e,j,l,m,o,q,u,f,h,k,i,param', reshape(c,1,9)];
      
    end;
    
% Display the interation (e) number in command window;
chronicprev_results(end,1)

    end;

% Save the intervention workspace;
save(['HCVInterventions04202017' num2str(indexnum) '.mat'])

% Display the index number in command window;
indexnum

% Clear workspace and results to run the next set of 500 parameter sets
clear

end;

% Combine the 10 Intervention Workspaces

chronicprev_results1 = [];
newchronic_results1 = [];
acute_results1 = [];
parameters_results1 = [];
  
for indexnum=1:10;
    
    load(['HCVInterventions04202017' num2str(indexnum) '.mat']);
                            
    chronicprev_results1 = [chronicprev_results1;chronicprev_results];
    newchronic_results1 = [newchronic_results1;newchronic_results];
    acute_results1 = [acute_results1;acute_results];
    parameters_results1 = [parameters_results1;parameters_results];

end;

% Save as a new workspace
save('InterventionResults_04202017.mat','-v7.3')

% Save intervention simulation results as txt files
dlmwrite('ChronicPrev_SingleSeq.txt',[chronicprev_results1], 'precision', 10);
dlmwrite('NewChronic_SingleSeq.txt',[newchronic_results1], 'precision', 10);
dlmwrite('Acute_SingleSeq.txt',[acute_results1], 'precision', 10);
dlmwrite('Parameters_SingleSeq.txt',[parameters_results1], 'precision', 10);

%% View the intervention results workspace
load('InterventionResults_04202017.mat')
