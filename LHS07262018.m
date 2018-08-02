%% Latin Hypercube Sample: Parameter Ranges and Scaling

load('HCVLHS_Setup_04162018.mat');

% Rescale and reshape contact matrix
    
    %Set 30-64 yo contacts with other age groups to 0 sampling minimum;
    minc(4,1:3)=0;
    minc(1:3,4)=0;

    %Reshape contact matrix 
    minc=reshape(minc,1,16);
    maxc=reshape(maxc,1,16);

% Setup the Latin Hypercube Sample 

% lhsdesign(# simulations, # parameters sampled in paramranges)
% 10,000 runs with 73 sampled parameters (some not actually sampled but are intervention related)
paramscaling = lhsdesign(10000,73);

% Set parameter upper and lower sampling bounds
%paramranges#s sigma1 sigma2 sigma3	 sigma4	alpha  beta        gamma1	gamma2	gamma3	gamma4	delta  epsilon	zeta1	zeta2  zeta3  zeta4	etaN	etaP   etaZ	  theta1  theta2  theta3  theta4   kappa  kappa	kappa  kappa  lambda1  labmda2 lambda3 lambda4 mu1	   mu2	    mu3	     mu4	   xi	omic1  omic2  omic3	 omic4	rho 	tau	 phiP phiN	psi1  psi2	psi3   psi4	  psiN	 omega	intervention_etaP	Z0     contacts popsize1:4                    nsduhyr];                   
%paramranges#s 1      2       3      4      5      6           7       8       9       10      11     12       13      14     15     16       17     18     19     20      21      22      23       24     25    26     27      28      29      30      31      32       33       34       35       36   37     38     39     40     41      42   43   44    45    46    47     48     49     50         51              52     53:68    69     70     71     72       73  ];                   
  paramranges=[5      5       5      5      0.8   0.0000035   0.07    0.07    0.07    0.07    0.15   1/0.5    0.0015 0.0026 0.0019 0.0016  0.18    2.5    4.3    0.05    0.034   0.034   0.034    0.37   0.37  0.16   0.16    0.1    0.2      0.3     0.552     0.000553 0.000915 0.00119  0.00493  0    0.0043 0.0071 0.0081 0.016  1/50    0    0    0     0.016 0.011 0.0037 0.0024 0.0081 365.25/56  1               144097 minc     723319 815922 483070 4605115  1;
               20     20      20     20     1     0.0000035   1.2     1.2     1.2     1.2     0.5    1/0.5    0.010  0.010  0.010  0.0110  0.54    15.3   4.43   1.5     1.5     1.5     1        2.6    2.6   0.92   0.92    0.1    0.2      0.3     0.552     0.000553 0.000915 0.00119  0.00493  0.45 0.010  0.010  0.021  0.027  1/50    0    0    0     0.021 0.015 0.010  0.0048 0.03   365.25/112 1               144097 maxc     723319 815922 483070 4605115  15];

% time1 vector is used to record data every 0.2 years from time 0 to 16;
time1 = [0:0.2:16];

% set upper and lower bounds for param1

x1=[0   0.05   0.034   0.034]' %lower bounds for param1 (beta, theta1, theta2, theta3);
x2=[0.1 5      5      10]' %upper bounds for param1 (beta, theta1, theta2, theta3);


%Simulate the ODE Across the Latin Hypercube Sampled Space

% Setup empty matrices to record simulation results
NewChronicPWID = [];
NewChronicNPWID = [];
NewChronic = [];

        NewChronic_Age1PWID=[];
        NewChronic_Age1NPWID=[];
        NewChronic_Age1=[];

        NewChronic_Age2PWID=[];
        NewChronic_Age2NPWID=[];
        NewChronic_Age2=[];

        NewChronic_Age3PWID=[];
        NewChronic_Age3NPWID=[];
        NewChronic_Age3=[];

        NewChronic_Age4PWID=[];
        NewChronic_Age4NPWID=[];
        NewChronic_Age4=[];

ChronicPrevPWID = [];
ChronicPrevNPWID = [];
ChronicPrevRuns = [];

        ChronicPrev_Age1PWID = [];
        ChronicPrev_Age1NPWID = [];
        ChronicPrev_Age1 = [];

        ChronicPrev_Age2PWID = [];
        ChronicPrev_Age2NPWID = [];
        ChronicPrev_Age2 = [];

        ChronicPrev_Age3PWID = [];
        ChronicPrev_Age3NPWID = [];
        ChronicPrev_Age3 = [];

        ChronicPrev_Age4PWID = [];
        ChronicPrev_Age4NPWID = [];
        ChronicPrev_Age4 = [];

AcutePWID = [];
AcuteNPWID = [];
AcuteRuns = [];

        Acute_Age1PWID = [];
        Acute_Age1NPWID = [];
        Acute_Age1 = [];

        Acute_Age2PWID = [];
        Acute_Age2NPWID = [];
        Acute_Age2 = [];

        Acute_Age3PWID = [];
        Acute_Age3NPWID = [];
        Acute_Age3 = [];

        Acute_Age4PWID = [];
        Acute_Age4NPWID = [];
        Acute_Age4 = [];

SuscPWID = [];
SuscNPWID = [];
Susc = [];

        Susc_Age1PWID = [];
        Susc_Age1NPWID = [];
        Susc_Age1 = [];

        Susc_Age2PWID = [];
        Susc_Age2NPWID = [];
        Susc_Age2 = [];

        Susc_Age3PWID = [];
        Susc_Age3NPWID = [];
        Susc_Age3 = [];

        Susc_Age4PWID = [];
        Susc_Age4NPWID = [];
        Susc_Age4 = [];

ImmunePWID = [];
ImmuneNPWID = [];
Immune = [];

        Immune_Age1PWID = [];
        Immune_Age1NPWID = [];
        Immune_Age1 = [];

        Immune_Age2PWID = [];
        Immune_Age2NPWID = [];
        Immune_Age2 = [];

        Immune_Age3PWID = [];
        Immune_Age3NPWID = [];
        Immune_Age3 = [];

        Immune_Age4PWID = [];
        Immune_Age4NPWID = [];
        Immune_Age4 = [];

AbDep=[];
        AbDep_Age1=[];
        AbDep_Age2=[];
        AbDep_Age3=[];
        AbDep_Age4=[];

RSSRuns = [];
ParamEstRuns=[];
CMatrixRuns = [];
InitialParamRuns=[];
PopSizeRuns=[];
FormerInjectRuns=[];
PYInjectRuns=[];

% Parfor Loop to run model simulations over all sampled parameter sets;
% second number (8) specifies number of workers;
parfor (j = 1:length(paramscaling(:,1)),8)
    
        %alternatively can use a for loop with 1 worker
        %for j = 1:length(paramscaling(:,1)) 
        
    % select set of params based on paramranges set above and lhs;
    fullparams = paramscaling(j,:).*(paramranges(2,:) - paramranges(1,:)) + paramranges(1,:);
    
%paramranges#s sigma1 sigma2 sigma3	 sigma4	alpha  beta	gamma1	gamma2	gamma3	gamma4	delta  epsilon	zeta1	zeta2  zeta3   zeta4	etaN	etaP   etaZ	  theta1  theta2  theta3  theta4   kappa  kappa	kappa  kappa  lambda1  labmda2 lambda3 lambda4 mu1	   mu2	   mu3	   mu4	   xi	omic1  omic2  omic3	 omic4	rho   	tau	 phiP phiN	psi1  psi2	psi3   psi4	  psiN	 omega	intervention_etaP	Z0     contacts popsize1:4                    nsduhyr];                   
%paramranges#s 1      2       3      4      5      6    7       8       9       10      11     12       13      14     15      16       17      18     19     20      21      22      23       24     25    26     27      28      29      30      31      32      33      34      35      36   37     38     39     40     41      42   43   44    45    46    47     48     49     50         51              52     53:68    69     70     71     72       73  ];                   
    
    param1 = [fullparams(1,6) fullparams(20:22)]';
    
    paramset = [fullparams(1,1:5) fullparams(1,7:19) fullparams(23:52)]';
    
    c = reshape(fullparams(1,53:68),4,4)

    %Calculate proportion of contacts in each age group by column
    pi=zeros(4,4);
    for i=1:4 %rows
        for l=1:4 %columns
        pi(i,l)=c(i,l)./sum(c(1:4,l));
        end
    end

    pi=reshape(pi,16,1);
    
    popsize=fullparams(69:72)';
    
    %select pyinject and formerinject prevalences by rounding nsduhyr to
    %nearest integer
    nsduhyr=round(fullparams(73),0);
    
%set values of Ci and Fi for the calculation of kappa;
F=formerinject(nsduhyr,:);
C=pyinject(nsduhyr,:);

%define paramset values required to calculate kappa
g1=paramset(6);
g2=paramset(7);
g3=paramset(8);
g4=paramset(9);
mu1=paramset(28);
mu2=paramset(29);
mu3=paramset(30);
mu4=paramset(31);
etap=paramset(17);
etan=paramset(16);

%Rescale etan to 1/etap in the case that etan*etap<1
if etan*etap<1 
   paramset(16) = 1/etap;
end;
    
   
% Calculate k based on selected params using steady state eqn for Fi with
% values sampled for this paramset
omic1=paramset(33);
omic2=paramset(34);
omic3=paramset(35);
omic4=paramset(36);
zeta1=paramset(12);
zeta2=paramset(13);
zeta3=paramset(14);
zeta4=paramset(15);

paramset(20)=g1*zeta1*popsize(1,1)/(omic1*popsize(1,1))-mu1*etap*etan-nu(1,1)
paramset(21)=(g2*zeta2*popsize(2,1)+nu(1,1)*omic1*popsize(1,1))/(omic2*popsize(2,1))-mu2*etap*etan-nu(2,1)
paramset(22)=(g3*zeta3*popsize(3,1)+nu(2,1)*omic2*popsize(2,1))/(omic3*popsize(3,1))-mu3*etap*etan-nu(3,1)
paramset(23)=(g4*zeta4*popsize(4,1)+nu(3,1)*omic3*popsize(3,1))/(omic4*popsize(4,1))-mu4*etap*etan-nu(4,1)
 
for i=1:4
    if paramset(i+19)<0
        paramset(i+19)=0.1
    end
end
    
% Define parameters required for initial conditions;
lambda1=paramset(24);
lambda2=paramset(25);
lambda3=paramset(26);
lambda4=paramset(27);
d=paramset(10);
xi=paramset(32);
r=paramset(37);
eps=paramset(11);
psi1=paramset(41);
psi2=paramset(42);
psi3=paramset(43);
psi4=paramset(44);
omic1=paramset(33);
omic2=paramset(34);
omic3=paramset(35);
omic4=paramset(36);
zeta1=paramset(12);
zeta2=paramset(13);
zeta3=paramset(14);
zeta4=paramset(15);

   
% Set Initial Conditions;
N0=[0 0 0 0;
    0 0 0 0;
    omic1*lambda1*(1-d)*popsize(1,1) omic2*lambda2*(1-d)*popsize(2,1) omic3*lambda3*(1-d)*popsize(3,1) omic4*lambda4*(1-d)*popsize(4,1);
    omic1*d*xi*lambda1*popsize(1,1)  omic2*d*xi*lambda2*popsize(2,1)  omic3*d*xi*lambda3*popsize(3,1)  omic4*d*xi*lambda4*popsize(4,1);
    0 0 0 0;
    0 0 0 0; 
    (1/r)*acute1529idu(1,1) (1/r)*acute1529idu(1,2) (1/r)*acute1529idu(1,3) 4*(1/r)*acute1529idu(1,3);
    zeta1*lambda1*(1-d)*popsize(1,1) zeta2*lambda2*(1-d)*popsize(2,1)  zeta3*lambda3*(1-d)*popsize(3,1) zeta4*lambda4*(1-d)*popsize(4,1);
    zeta1*d*xi*lambda1*popsize(1,1)  zeta2*d*xi*lambda2*popsize(2,1)   zeta3*d*xi*lambda3*popsize(3,1)  zeta4*d*xi*lambda4*popsize(4,1);
    0 0 0 0;
    psi1*popsize(1,1) psi2*popsize(2,1) psi3*popsize(3,1) psi4*popsize(4,1)];

% Define class Si (uninfected PWID) as popsize*pyinject less sum of current PWID classes;

zeta=[zeta1 zeta2 zeta3 zeta4];
for i=1:4
    N0(6,i) = popsize(i,1)*zeta(1,i) - (N0(7,i)+N0(8,i)+N0(9,i)+N0(10,i));
end

% Define class SNi (uninfected former PWID) as popsize*formerinject less sum of former PWID classes;

omic=[omic1 omic2 omic3 omic4];
for i=1:4
    N0(1,i) = popsize(i,1)*omic(1,i) - (N0(2,i)+N0(3,i)+N0(4,i)+N0(5,i));
end

% Reshape initial conditions into a vector;
N0=reshape(N0,44,1);

%Estimate Params
y=zeros(1,3); 

[paramest RSS] = fminsearchbnd(@(param1) HCV_LS_04072018(time,N0,param1,paramset,nu,pi,popsize,acute1529idu),param1,x1,x2,optimset('MaxFunEvals',5000));

% save results where model converged (an RSS value was provided)
if RSS>0
    
paramest=abs(paramest);
  
paramnew = [paramset(1:5);paramest(1); paramset(6:18); paramest(2:4); paramset(19:48)];
    
% Run model with sampled and estimated parameters

[t,x] = ode15s(@HCV_DiffEq_03272018,time1,N0,[],paramnew,nu,pi,pop_pyabdepnoinject);
            
% Save Output 

% Total Cases;
Y_newchronic_pwid = eps*(1-d)*(x(:,11*1-4)+x(:,11*2-4)+x(:,11*3-4)+x(:,11*4-4)); 
Y_newchronic_npwid = eps*(1-d)*(x(:,11*1-9)+x(:,11*2-9)+x(:,11*3-9)+x(:,11*4-9));
 
Y_chronicprev_pwid = (x(:,11*1-3)+x(:,11*2-3)+x(:,11*3-3)+x(:,11*4-3));
Y_chronicprev_npwid = (x(:,11*1-8)+x(:,11*2-8)+x(:,11*3-8)+x(:,11*4-8));

Y_acute_pwid = (x(:,11*1-4)+x(:,11*2-4)+x(:,11*3-4)+x(:,11*4-4));
Y_acute_npwid = (x(:,11*1-9)+x(:,11*2-9)+x(:,11*3-9)+x(:,11*4-9));

Y_S_pwid = x(:,11*1-5)+x(:,11*2-5)+x(:,11*3-5)+x(:,11*4-5);
Y_S_npwid = x(:,11*1-10)+x(:,11*2-10)+x(:,11*3-10)+x(:,11*4-10);
 
Y_I_pwid = x(:,11*1-2)+x(:,11*2-2)+x(:,11*3-2)+x(:,11*4-2);
Y_I_npwid = x(:,11*1-7)+x(:,11*2-7)+x(:,11*3-7)+x(:,11*4-7);
 
Y_Z=x(:,11*1)+x(:,11*2)+x(:,11*3)+x(:,11*4);

 NewChronicPWID(:,j) = Y_newchronic_pwid;
 NewChronicNPWID(:,j) = Y_newchronic_npwid;
 NewChronic(:,j) = Y_newchronic_pwid+Y_newchronic_npwid;

 ChronicPrevPWID(:,j) = Y_chronicprev_pwid;
 ChronicPrevNPWID(:,j) = Y_chronicprev_npwid;
 ChronicPrevRuns(:,j) = Y_chronicprev_pwid+Y_chronicprev_npwid;

 AcutePWID(:,j) = Y_acute_pwid;
 AcuteNPWID(:,j) = Y_acute_npwid;
 AcuteRuns(:,j) = Y_acute_pwid+Y_acute_npwid;

 SuscPWID(:,j) = Y_S_pwid;
 SuscNPWID(:,j) = Y_S_npwid;
 Susc(:,j) = Y_S_pwid+Y_S_npwid;
 
 ImmunePWID(:,j) = Y_I_pwid;
 ImmuneNPWID(:,j) = Y_I_npwid;
 Immune(:,j) = Y_I_pwid+Y_I_npwid;
 
 AbDep(:,j)= Y_Z;
 

%Save 15-19 results
i=1;

Y_newchronic_pwid = eps*(1-d)*x(:,11*i-4);
Y_newchronic_npwid = eps*(1-d)*x(:,11*i-9);

        NewChronic_Age1PWID(:,j)=Y_newchronic_pwid;
        NewChronic_Age1NPWID(:,j)=Y_newchronic_npwid;
        NewChronic_Age1(:,j)=Y_newchronic_pwid+Y_newchronic_npwid;

Y_chronicprev_pwid = x(:,11*i-3);
Y_chronicprev_npwid = x(:,11*i-8);

        ChronicPrev_Age1PWID(:,j)=Y_chronicprev_pwid;
        ChronicPrev_Age1NPWID(:,j)=Y_chronicprev_npwid;
        ChronicPrev_Age1(:,j)=Y_chronicprev_pwid+Y_chronicprev_npwid;

Y_acute_pwid = x(:,11*i-4);
Y_acute_npwid = x(:,11*i-9);

        Acute_Age1PWID(:,j)=Y_acute_pwid;
        Acute_Age1NPWID(:,j)=Y_acute_npwid;
        Acute_Age1(:,j)=Y_acute_pwid+Y_acute_npwid;

Y_S_pwid = x(:,11*i-5);
Y_S_npwid = x(:,11*i-10);

        Susc_Age1PWID(:,j)=Y_S_pwid;
        Susc_Age1NPWID(:,j)=Y_S_npwid;
        Susc_Age1(:,j)=Y_S_pwid+Y_S_npwid;

Y_I_pwid = x(:,11*i-2);
Y_I_npwid = x(:,11*i-7);

        Immune_Age1PWID(:,j)=Y_I_pwid;
        Immune_Age1NPWID(:,j)=Y_I_npwid;
        Immune_Age1(:,j)=Y_I_pwid+Y_I_npwid;

Y_Z=x(:,11*i);

        AbDep_Age1(:,j)=Y_Z;
        
%Save 20-25 results
i=2;

Y_newchronic_pwid = eps*(1-d)*x(:,11*i-4);
Y_newchronic_npwid = eps*(1-d)*x(:,11*i-9);

        NewChronic_Age2PWID(:,j)=Y_newchronic_pwid;
        NewChronic_Age2NPWID(:,j)=Y_newchronic_npwid;
        NewChronic_Age2(:,j)=Y_newchronic_pwid+Y_newchronic_npwid;

Y_chronicprev_pwid = x(:,11*i-3);
Y_chronicprev_npwid = x(:,11*i-8);

        ChronicPrev_Age2PWID(:,j)=Y_chronicprev_pwid;
        ChronicPrev_Age2NPWID(:,j)=Y_chronicprev_npwid;
        ChronicPrev_Age2(:,j)=Y_chronicprev_pwid+Y_chronicprev_npwid;

Y_acute_pwid = x(:,11*i-4);
Y_acute_npwid = x(:,11*i-9);

        Acute_Age2PWID(:,j)=Y_acute_pwid;
        Acute_Age2NPWID(:,j)=Y_acute_npwid;
        Acute_Age2(:,j)=Y_acute_pwid+Y_acute_npwid;

Y_S_pwid = x(:,11*i-5);
Y_S_npwid = x(:,11*i-10);

        Susc_Age2PWID(:,j)=Y_S_pwid;
        Susc_Age2NPWID(:,j)=Y_S_npwid;
        Susc_Age2(:,j)=Y_S_pwid+Y_S_npwid;

Y_I_pwid = x(:,11*i-2);
Y_I_npwid = x(:,11*i-7);

        Immune_Age2PWID(:,j)=Y_I_pwid;
        Immune_Age2NPWID(:,j)=Y_I_npwid;
        Immune_Age2(:,j)=Y_I_pwid+Y_I_npwid;

Y_Z=x(:,11*i);

        AbDep_Age2(:,j)=Y_Z;

%Save 26-29 results
i=3;

Y_newchronic_pwid = eps*(1-d)*x(:,11*i-4);
Y_newchronic_npwid = eps*(1-d)*x(:,11*i-9);

        NewChronic_Age3PWID(:,j)=Y_newchronic_pwid;
        NewChronic_Age3NPWID(:,j)=Y_newchronic_npwid;
        NewChronic_Age3(:,j)=Y_newchronic_pwid+Y_newchronic_npwid;

Y_chronicprev_pwid = x(:,11*i-3);
Y_chronicprev_npwid = x(:,11*i-8);

        ChronicPrev_Age3PWID(:,j)=Y_chronicprev_pwid;
        ChronicPrev_Age3NPWID(:,j)=Y_chronicprev_npwid;
        ChronicPrev_Age3(:,j)=Y_chronicprev_pwid+Y_chronicprev_npwid;

Y_acute_pwid = x(:,11*i-4);
Y_acute_npwid = x(:,11*i-9);

        Acute_Age3PWID(:,j)=Y_acute_pwid;
        Acute_Age3NPWID(:,j)=Y_acute_npwid;
        Acute_Age3(:,j)=Y_acute_pwid+Y_acute_npwid;

Y_S_pwid = x(:,11*i-5);
Y_S_npwid = x(:,11*i-10);

        Susc_Age3PWID(:,j)=Y_S_pwid;
        Susc_Age3NPWID(:,j)=Y_S_npwid;
        Susc_Age3(:,j)=Y_S_pwid+Y_S_npwid;

Y_I_pwid = x(:,11*i-2);
Y_I_npwid = x(:,11*i-7);

        Immune_Age3PWID(:,j)=Y_I_pwid;
        Immune_Age3NPWID(:,j)=Y_I_npwid;
        Immune_Age3(:,j)=Y_I_pwid+Y_I_npwid;

Y_Z=x(:,11*i);

        AbDep_Age3(:,j)=Y_Z;

%Save 30-64 results
i=4;

Y_newchronic_pwid = eps*(1-d)*x(:,11*i-4);
Y_newchronic_npwid = eps*(1-d)*x(:,11*i-9);

        NewChronic_Age4PWID(:,j)=Y_newchronic_pwid;
        NewChronic_Age4NPWID(:,j)=Y_newchronic_npwid;
        NewChronic_Age4(:,j)=Y_newchronic_pwid+Y_newchronic_npwid;

Y_chronicprev_pwid = x(:,11*i-3);
Y_chronicprev_npwid = x(:,11*i-8);

        ChronicPrev_Age4PWID(:,j)=Y_chronicprev_pwid;
        ChronicPrev_Age4NPWID(:,j)=Y_chronicprev_npwid;
        ChronicPrev_Age4(:,j)=Y_chronicprev_pwid+Y_chronicprev_npwid;

Y_acute_pwid = x(:,11*i-4);
Y_acute_npwid = x(:,11*i-9);

        Acute_Age4PWID(:,j)=Y_acute_pwid;
        Acute_Age4NPWID(:,j)=Y_acute_npwid;
        Acute_Age4(:,j)=Y_acute_pwid+Y_acute_npwid;

Y_S_pwid = x(:,11*i-5);
Y_S_npwid = x(:,11*i-10);

        Susc_Age4PWID(:,j)=Y_S_pwid;
        Susc_Age4NPWID(:,j)=Y_S_npwid;
        Susc_Age4(:,j)=Y_S_pwid+Y_S_npwid;

Y_I_pwid = x(:,11*i-2);
Y_I_npwid = x(:,11*i-7);

        Immune_Age4PWID(:,j)=Y_I_pwid;
        Immune_Age4NPWID(:,j)=Y_I_npwid;
        Immune_Age4(:,j)=Y_I_pwid+Y_I_npwid;

Y_Z=x(:,11*i);

        AbDep_Age4(:,j)=Y_Z;
        
 % RSS;
 RSSRuns(1,j) = RSS;

 % Parameters;
 ParamEstRuns(:,j) = paramnew;

 % Contact Matrix;
 CMatrixRuns(j,:) = reshape(c,1,16);
 
 % Population Size;
 PopSizeRuns(j,:) = popsize';


 % formerinject;
 FormerInjectRuns(j,:) = [nsduhyr F(1,:)];
  
 % pyinject;
 PYInjectRuns(j,:) = [nsduhyr C(1,:)];
 
  %Initial Values of Estimated Params and initial values of etan and etap; 
 % Note that etan gets scaled to 1/etap if etan*etap<1;
 InitialParamRuns(j,:) = [param1' fullparams(17) fullparams(18)];

 j

%save results when model does not converge (no RSS value provided)
else
    
 % RSS;
 RSSRuns(1,j) = nan(1,1);

 % Parameters;
 ParamEstRuns(:,j) = [paramset(1:5);param1(1); paramset(6:18); param1(2:4); paramset(19:48)];

 % Contact Matrix;
 CMatrixRuns(j,:) = reshape(c,1,16);
 
 % Population Size;
 PopSizeRuns(j,:) = popsize';

 % formerinject;
 FormerInjectRuns(j,:) = [nsduhyr F(1,:)];
 
 % pyinject;
 PYInjectRuns(j,:) = [nsduhyr C(1,:)];
 
 %Initial Values of Estimated Params and initial values of etan and etap; 
 % Note that etan gets scaled to 1/etap if etan*etap<1;
 InitialParamRuns(j,:) = [param1' fullparams(17) fullparams(18)];
 
 NewChronicPWID(:,j) = nan(81,1);
 NewChronicNPWID(:,j) = nan(81,1);
 NewChronic(:,j) = nan(81,1);

 ChronicPrevPWID(:,j) = nan(81,1);
 ChronicPrevNPWID(:,j) = nan(81,1);
 ChronicPrevRuns(:,j) = nan(81,1);

 AcutePWID(:,j) = nan(81,1);
 AcuteNPWID(:,j) = nan(81,1);
 AcuteRuns(:,j) = nan(81,1);

 SuscPWID(:,j) = nan(81,1);
 SuscNPWID(:,j) = nan(81,1);
 Susc(:,j) = nan(81,1);
 
 ImmunePWID(:,j) = nan(81,1);
 ImmuneNPWID(:,j) = nan(81,1);
 Immune(:,j) = nan(81,1);
 
 AbDep(:,j)= nan(81,1);
 

NewChronic_Age1PWID(:,j)=nan(81,1);
NewChronic_Age1NPWID(:,j)=nan(81,1);
NewChronic_Age1(:,j)=nan(81,1);

NewChronic_Age2PWID(:,j)=nan(81,1);
NewChronic_Age2NPWID(:,j)=nan(81,1);
NewChronic_Age2(:,j)=nan(81,1);

NewChronic_Age3PWID(:,j)=nan(81,1);
NewChronic_Age3NPWID(:,j)=nan(81,1);
NewChronic_Age3(:,j)=nan(81,1);

NewChronic_Age4PWID(:,j)=nan(81,1);
NewChronic_Age4NPWID(:,j)=nan(81,1);
NewChronic_Age4(:,j)=nan(81,1);

ChronicPrev_Age1PWID(:,j)=nan(81,1);
ChronicPrev_Age1NPWID(:,j)=nan(81,1);
ChronicPrev_Age1(:,j)=nan(81,1);

ChronicPrev_Age2PWID(:,j)=nan(81,1);
ChronicPrev_Age2NPWID(:,j)=nan(81,1);
ChronicPrev_Age2(:,j)=nan(81,1);

ChronicPrev_Age3PWID(:,j)=nan(81,1);
ChronicPrev_Age3NPWID(:,j)=nan(81,1);
ChronicPrev_Age3(:,j)=nan(81,1);

ChronicPrev_Age4PWID(:,j)=nan(81,1);
ChronicPrev_Age4NPWID(:,j)=nan(81,1);
ChronicPrev_Age4(:,j)=nan(81,1);

Acute_Age1PWID(:,j)=nan(81,1);
Acute_Age1NPWID(:,j)=nan(81,1);
Acute_Age1(:,j)=nan(81,1);

Acute_Age2PWID(:,j)=nan(81,1);
Acute_Age2NPWID(:,j)=nan(81,1);
Acute_Age2(:,j)=nan(81,1);

Acute_Age3PWID(:,j)=nan(81,1);
Acute_Age3NPWID(:,j)=nan(81,1);
Acute_Age3(:,j)=nan(81,1);

Acute_Age4PWID(:,j)=nan(81,1);
Acute_Age4NPWID(:,j)=nan(81,1);
Acute_Age4(:,j)=nan(81,1);

Susc_Age1PWID(:,j)=nan(81,1);
Susc_Age1NPWID(:,j)=nan(81,1);
Susc_Age1(:,j)=nan(81,1);

Susc_Age2PWID(:,j)=nan(81,1);
Susc_Age2NPWID(:,j)=nan(81,1);
Susc_Age2(:,j)=nan(81,1);

Susc_Age3PWID(:,j)=nan(81,1);
Susc_Age3NPWID(:,j)=nan(81,1);
Susc_Age3(:,j)=nan(81,1);

Susc_Age4PWID(:,j)=nan(81,1);
Susc_Age4NPWID(:,j)=nan(81,1);
Susc_Age4(:,j)=nan(81,1);

Immune_Age1PWID(:,j)=nan(81,1);
Immune_Age1NPWID(:,j)=nan(81,1);
Immune_Age1(:,j)=nan(81,1);

Immune_Age2PWID(:,j)=nan(81,1);
Immune_Age2NPWID(:,j)=nan(81,1);
Immune_Age2(:,j)=nan(81,1);

Immune_Age3PWID(:,j)=nan(81,1);
Immune_Age3NPWID(:,j)=nan(81,1);
Immune_Age3(:,j)=nan(81,1);

Immune_Age4PWID(:,j)=nan(81,1);
Immune_Age4NPWID(:,j)=nan(81,1);
Immune_Age4(:,j)=nan(81,1);


AbDep_Age1(:,j)=nan(81,1);
AbDep_Age2(:,j)=nan(81,1);
AbDep_Age3(:,j)=nan(81,1);
AbDep_Age4(:,j)=nan(81,1);
    
end

fprintf(['LHS number %i.\n'],j)

end

disp('LHScomplete') 

%Save the Workspace
save('HCVLHS05062018.mat')