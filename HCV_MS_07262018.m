%% Main Script: Gicquelais et al. (2018). Hepatitis C tranmsission model.

% This code is the main script file to fit a hepatitis C transmission model
% to public health surveillance data from Michigan.

% Four parameters are estimated: 3 injection initiation rates (theta_i) and
% 1  transmission rate (beta). Latin hypercube sampling is used to
% incorporate parameter uncertainty over 10,000 simulations. Several
% interventions are simulated. All results are exported as .csv or .txt
% files.

% Accompanying files:
    % Differential equations: HCV_DiffEq_07262018.m
    % Least squares for parameter estimation: HCV_LS_07262018.m
    % Contact matrix: PolymodPhysical_02182018.m
    % Latin hypercube sampling: LHS07262018.m

% For reference, the parameters are:

% Estimated parameters: 
    % param1= [b theta1 theta2 theta3];

% Sampled or set parameters: 
    % paramset = [sigma1 sigma2 sigma3 sigma4 a g1 g2 g3 g4 d eps 
    %             zeta1 zeta2 zeta3 zeta4 etan etap etaz theta4 
    %             k1 k2 k3 k4 lambda1 lambda2 lambda3 lambda4 
    %             mu1 mu2 mu3 mu4 xi omic1 omic2 omic3 omic4 r tau 
    %             phip phin psi1 psi2 psi3 psi4 psin w intervention_etap Z0]

%Param1 and paramset are combined into a single matrix called param:
%     param  =[sigma1 sigma2 sigma3 sigma4 
%              1      2     3       4          
%
%          a b  g(1:4) d  eps zeta(1:4) etan etap etaz theta(1:4) k(1:4) lambda(1:4)
%          5 6  7:10   11 12  13:16     17   18   19   20:23      24:27   28:31
%
%          mu(1:4) xi omic(1:4) r  tau phip phin psi(1:4) psin w  intervention_etap Z0];
%          32:35   36 37:40     41 42  43   44   45:48    49   50          51       52  


% The inverse of age group size is in vector nu:
    % nu = [nu1,nu2,nu3,nu4];

% Compartments of the model are:
    % N[SNi; ANi; CNi; INi; TNi; Si; Ai; Ci; Ii; Ti; Zi];
    % The first letter is the compartment type 
        % S=Susceptible [HCV Uninfected PWID], A=Acute HCV, C=Chronic HCV, 
        % I=Immune, T=Treated, Z=Non-PWID;
    % Former PWID have an 'N' following the compartment type while current
        % PWID and Non-PWID (Z) have no 'N';
    % i denotes the age class in years (1: 15-19, 2: 20-25, 3: 26-29, 
        % 4: 30-64);

%% HCV Surveillance Data, 2000-2016 (Michigan Department of Health and Human Services)

% Set the years data available
year = [2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016];
time = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];

% Create matrices for HCV case counts by age group (15-19, 20-25, 26-29];

% Total acute cases reporting yes, missing, or unknown injection of drugs
% in the 6 months before symptom onset;
acute1529idu = [4	1	1	1	4	1	2	7	5	3	2	2	5	2	2	1	2*(1-0.26)
                1	2	3	7	8	3	10	9	17	11	9	7	11	16	19	18	29*(1-0.26)
                4	6	4	5	2	5	8	4	9	5	9	3	8	11	8	13	14*(1-0.26)]';

acute1529idutotal=sum(acute1529idu,2);

% Approximate the number of acute cases corrected for under-reporting using 
% Klevens RM et al. (2014). Estimating acute viral hepatitis infections from 
% nationally reported cases. American Journal of Public Health, 104(3):
% 482-487.

% Hepatitis C under-reporting factor (1/# cases per surveillance detected
% case, this estimate is between the PWID and overall rate for hepatitis C
% from Klevens et al. and Onofrey et al.;
r=1/16.8;
rlow=1/12.3;
rhigh=1/168;

% Multiply by acute cases observed to approximate the total number of acute
% cases;
acutecorridu=(1/r)*acute1529idu; 
acutecorrtotalidu=sum(acutecorridu,2);

% New chronic or indeterminate status cases (acute vs. chronic) who 
% reported yes, unknown, or missing for lifetime drug injection);
chronic1529idu = [10	21	19	40	33	53	65	57	84	84	52	46	106	85	76	76	104*(1-0.26)
                  30	29	95	79	132	165	247	204	308	356	213	288	667	612	686	648	902
                  21	35	73	82	81	122	157	188	221	271	162	208	435	419	592	639	915]';

chronic1529idutotal=sum(chronic1529idu,2);

% Cases of indeterminate chronic vs acute status who reported yes, 
% unknown, or missing lifetime drug injections;
unknown1529idu=[9	8	1	0	1	0	2	0	4	0	0	0	0	0	0	0	1*(1-0.26)
                19	30	9	0	1	5	2	2	8	0	1	0	0	0	0	1	1
                25	34	5	0	1	3	2	4	13	0	1	3	1	0	0	0	1]';

unknown1529idutotal=sum(unknown1529idu,2);

% Sum unknown and new chronic cases;
chronic1529idu_model=chronic1529idu+unknown1529idu;
chronic1529idu_modeltotal=sum(chronic1529idu_model,2);


% Create a dark green color
rgrn=[0 0.6 0.35];

% Plot case counts for each compartment at each time step;
figure(1)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year,chronic1529idu_modeltotal,'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    %plot(year,chronic1529idutotal+unknown1529idutotal,'Color',rgrn,'LineStyle','--','Marker','s','LineWidth',4,'MarkerSize',10)
    %    hold on
    plot(year,acute1529idutotal,'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,acutecorrtotalidu,'rd--','LineWidth',4,'MarkerSize',10)
        ylim([0 2000])
        xlim([2000 2016])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Number of Cases')
        title('Figure 1. Michigan HCV Infections Among 15-29 Year Olds -- 2000-2016')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('Chronic and Unconfirmed PWID','Acute PWID','Acute PWID (Estimated Using CDC Correction)')
%saveas(gcf,'Chronic Data','fig')

figure(2)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year,chronic1529idu_model(:,1),'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,chronic1529idu_model(:,2),'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,chronic1529idu_model(:,3),'rd--','LineWidth',4,'MarkerSize',10)
        ylim([0 1000])
        xlim([2000 2016])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Number of Cases')
        title('Figure 2. Michigan New Chronic HCV Infections Among 15-29 Year Olds -- 2000-2016')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('15-19','20-25','26-29')
%saveas(gcf,'Chronic Data by Age','fig')

figure(3)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year,(1/r)*acute1529idu(:,1),'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,(1/r)*acute1529idu(:,2),'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,(1/r)*acute1529idu(:,3),'rd--','LineWidth',4,'MarkerSize',10)
        ylim([0 600])
        xlim([2000 2016])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Number of Cases')
        title('Figure 3. Michigan Acute HCV Infections Among 15-29 Year Olds -- 2000-2016')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('15-19','20-25','26-29')
%saveas(gcf,'Acute Data by Age','fig')

figure(4)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
    plot(year,(1/rlow)*acute1529idu(:,1),'r--','LineWidth',2,'MarkerSize',10)
        hold on
    plot(year,(1/r)*acute1529idu(:,1),'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,(1/rhigh)*acute1529idu(:,1),'r--','LineWidth',2,'MarkerSize',10)
        hold on
    plot(year,chronic1529idu(:,1),'b*--','LineWidth',4,'MarkerSize',10)
        hold on
        ylim([0 200])
        xlim([2000 2016])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Number of Cases')
        title('Figure 4. Michigan Acute & NC 15-19 Year Olds -- 2000-2016')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('acute: low under-report','acute: avg underreport','acute: high under-report','newchronic')
        
figure(5)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
    plot(year,(1/rlow)*acute1529idu(:,2),'r--','LineWidth',2,'MarkerSize',10)
        hold on
    plot(year,(1/r)*acute1529idu(:,2),'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year,(1/rhigh)*acute1529idu(:,2),'r--','LineWidth',2,'MarkerSize',10)
        hold on
    plot(year,chronic1529idu(:,2),'b*--','LineWidth',4,'MarkerSize',10)
        hold on
        ylim([0 1000])
        xlim([2000 2016])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Number of Cases')
        title('Figure 5. Michigan Acute & NC 20-25 Year Olds -- 2000-2016')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('acute: low under-report','acute: avg underreport','acute: high under-report','newchronic')
        
  figure(6)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
    plot(year,(1/rlow)*acute1529idu(:,3),'r--','LineWidth',2,'MarkerSize',10)
        hold on
    plot(year,(1/r)*acute1529idu(:,3),'ko--','LineWidth',4,'MarkerSize',10)
        hold on
   plot(year,(1/rhigh)*acute1529idu(:,3),'r--','LineWidth',2,'MarkerSize',10)
        hold on
    plot(year,chronic1529idu(:,3),'b*--','LineWidth',4,'MarkerSize',10)
        hold on
        ylim([0 1000])
        xlim([2000 2016])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Number of Cases')
        title('Figure 6. Michigan Acute & NC 26-29 Year Olds -- 2000-2016')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('acute: low under-report','acute: avg underreport','acute: high under-report','newchronic')

%% National Survey on Drug Use and Health (NSDUH) Data, 2000-2014

% Precision of all estimates to two significant digits;

% Prevalence of past year abuse and dependence without injection drug use by age (15-19, 20-25, 26-29, 30-64);
pyabdepnoinject = [ 0.018	0.016	0.021	0.022	0.026	0.025	0.024	0.022	0.024	0.015	0.020	0.019	0.016	0.011	0.011
                    0.013	0.014	0.023	0.020	0.025	0.027	0.022	0.029	0.027	0.025	0.020	0.020	0.024	0.019	0.016
                    0.0063	0.011	0.019	0.022	0.020	0.018	0.025	0.018	0.019	0.012	0.013	0.012	0.019	0.015	0.018
                    0.0035	0.0066	0.009	0.0084	0.0085	0.0075	0.0082	0.008	0.0077	0.0076	0.0071	0.0072	0.0081	0.0075	0.0094]';
               
% Prevalence of past year injection drug use by age (15-19, 20-25, 26-29, 30-64);
pyinject = [0.0015	0.0037	0.0025	0.0021	0.0016	0.0020	0.0013	0.0018	0.0013	0.0015	0.0018	0.0017	0.0031	0.0011	0.0014
            0.0026	0.0035	0.002	0.0036	0.0037	0.0052	0.0028	0.0039	0.0028	0.0027	0.0029	0.0040	0.0064	0.0063	0.0061
            0.0019	0.0023	0.0006	0.0023	0.0062	0.0044	0.0067	0.0017	0.0035	0.0032	0.0051	0.0026	0.0045	0.0051	0.0087
            0.0016	0.0018	0.0023	0.0014	0.0020	0.0017	0.0018	0.0014	0.0016	0.0025	0.0018	0.0018	0.0018	0.0024	0.0016]';
        
% Prevalence of former injection drug use (lifetime - past year idu) by age (15-19, 20-25, 26-29, 30-64);
formerinject = [0.0043	0.0051	0.0042	0.0032	0.0037	0.0045	0.0046	0.0030	0.0027	0.0026	0.0034	0.0019	0.0021	0.0020	0.0026
                0.0071	0.013	0.011	0.011	0.012	0.010	0.011	0.011	0.0079	0.013	0.010	0.0097	0.0096	0.011	0.011
                0.0081	0.0099	0.0084	0.018	0.011	0.0096	0.016	0.017	0.019	0.016	0.015	0.015	0.016	0.012	0.015
                0.016	0.018	0.023	0.021	0.018	0.019	0.021	0.022	0.023	0.019	0.018	0.019	0.021	0.022	0.020]';

% Prevalence of past year use without injection drug use by age (15-19, 20-25, 26-29, 30-64); 
pyusenoinject = [0.11	0.13	0.15	0.15	0.15	0.14	0.14	0.14	0.12	0.12	0.12	0.11	0.1	0.088	0.089
                0.092	0.12	0.15	0.15	0.15	0.16	0.16	0.15	0.15	0.15	0.14	0.13	0.13	0.12	0.12
                0.053	0.077	0.10	0.11	0.10	0.089	0.11	0.10	0.10	0.10	0.10	0.082	0.11	0.11	0.088
                0.027	0.031	0.046	0.048	0.043	0.043	0.049	0.047	0.045	0.044	0.048	0.041	0.049	0.043	0.045]';

%95% CI bounds to sample former pwid;
pyinjectmin=[0.00086	0.0026	0.0014	0.0011	0.0010	0.0011	0.00076	0.00093	0.00074	0.00086	0.00091	0.00089	0.0018	0.00059	0.00064
            0.0016	0.0020	0.0011	0.0024	0.0025	0.0037	0.0019	0.0026	0.0017	0.0018	0.0019	0.0025	0.0046	0.0045	0.0044
            0.00057	0.00056	0.00001	0.00058	0.0029	0.0024	0.0033	0.00057	0.0015	0.00092	0.00175	0.00096	0.0019	0.0023	0.0053
            0.00069	0.0011	0.0014	0.00076	0.0012	0.00093	0.00088	0.00076	0.00086	0.0015	0.0011	0.00096	0.001	0.0015	0.0010]';

pyinjectmax=[0.0025	0.0050	0.0041	0.0035	0.0025	0.0033	0.0020	0.0030	0.0020	0.0024	0.0033	0.0028	0.0050	0.002	0.0028
            0.0039	0.0056	0.0033	0.0051	0.0052	0.0070	0.0039	0.0057	0.0043	0.0039	0.0041	0.0061	0.0087	0.0085	0.0083
            0.0045	0.0063	0.0036	0.0059	0.012	0.0076	0.012	0.0039	0.0069	0.0079	0.012	0.0056	0.0090	0.0098	0.013
            0.0032	0.0028	0.0034	0.0022	0.0032	0.0029	0.0032	0.0023	0.0027	0.0040	0.0029	0.0032	0.0029	0.0035	0.0023]';

%95% CI bounds to sample former pwid;
formerinjectmin=[0.0032	0.0036	0.0029	0.0020	0.0024	0.0029	0.0030	0.0021	0.0017	0.0016	0.0019	0.0011	0.0013	0.00086	0.0015
                0.0054	0.0102	0.0084	0.0092	0.0095	0.0083	0.0088	0.0088	0.006	0.011	0.0079	0.0077	0.0075	0.0078	0.0083
                0.0052	0.0049	0.0045	0.013	0.0075	0.0064	0.011	0.011	0.013	0.011	0.011	0.011	0.0098	0.0081	0.011
                0.014	0.016	0.020	0.017	0.015	0.016	0.018	0.018	0.020	0.016	0.015	0.015	0.017	0.018	0.017]';

formerinjectmax=[0.0055	0.0069	0.0058	0.005	0.0054	0.0067	0.0067	0.0043	0.0041	0.0040	0.0055	0.0030	0.0033	0.0038	0.0043
                0.0093	0.016	0.013	0.014	0.015	0.012	0.015	0.014	0.010	0.016	0.013	0.012	0.012	0.014	0.014
                0.012	0.018	0.014	0.024	0.016	0.014	0.022	0.024	0.026	0.023	0.021	0.021	0.023	0.018	0.020
                0.019	0.021	0.026	0.025	0.021	0.022	0.024	0.025	0.026	0.022	0.022	0.023	0.024	0.026	0.023]';

           
figure(7)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year(1,1:15),pyabdepnoinject(:,1)*100,'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),pyabdepnoinject(:,2)*100,'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),pyabdepnoinject(:,3)*100,'rd--','LineWidth',4,'MarkerSize',10)
        hold on 
    plot(year(1,1:15),pyabdepnoinject(:,4)*100,'Color',rgrn,'LineStyle','--','Marker','s','LineWidth',4,'MarkerSize',10)
        ylim([0 3])
        xlim([2000 2014])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Prevalence (%)')
        title('Figure 7. Past Year Abuse or Dependence Prevalence (Excl IDU), US -- 2000-2014')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('15-19','20-25','26-29','30-64')
        
figure(8)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year(1,1:15),pyinject(:,1)*100,'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),pyinject(:,2)*100,'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),pyinject(:,3)*100,'rd--','LineWidth',4,'MarkerSize',10)
        hold on 
    plot(year(1,1:15),pyinject(:,4)*100,'Color',rgrn,'LineStyle','--','Marker','s','LineWidth',4,'MarkerSize',10)
        ylim([0 1])
        xlim([2000 2014])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Prevalence (%)')
        title('Figure 8. Past Year Injection Drug Use, US -- 2000-2014')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('15-19','20-25','26-29','30-64')
        
figure(9)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year(1,1:15),formerinject(:,1)*100,'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),formerinject(:,2)*100,'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),formerinject(:,3)*100,'rd--','LineWidth',4,'MarkerSize',10)
        hold on 
    plot(year(1,1:15),formerinject(:,4)*100,'Color',rgrn,'LineStyle','--','Marker','s','LineWidth',4,'MarkerSize',10)
        ylim([0 3])
        xlim([2000 2014])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Prevalence (%)')
        title('Figure 9. Former Injection Drug Use, US -- 2000-2014')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('15-19','20-25','26-29','30-64')
        
figure(10)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year(1,1:15),pyusenoinject(:,1)*100,'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),pyusenoinject(:,2)*100,'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),pyusenoinject(:,3)*100,'rd--','LineWidth',4,'MarkerSize',10)
        hold on 
    plot(year(1,1:15),pyusenoinject(:,4)*100,'Color',rgrn,'LineStyle','--','Marker','s','LineWidth',4,'MarkerSize',10)
        ylim([0 20])
        xlim([2000 2014])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Prevalence (%)')
        title('Figure 10. Past Year Use (No IDU), US -- 2000-2014')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('15-19','20-25','26-29','30-64')

%% Examine Ratio of Former vs Current PWID

ratio=formerinject./pyinject

figure(11)
    set(gca,'LineWidth',2,'FontSize',24,'FontName','Calibri')
        hold on
    plot(year(1,1:15),ratio(:,1),'ko--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),ratio(:,2),'b*--','LineWidth',4,'MarkerSize',10)
        hold on
    plot(year(1,1:15),ratio(:,3),'rd--','LineWidth',4,'MarkerSize',10)
        hold on 
    plot(year(1,1:15),ratio(:,4),'Color',rgrn,'LineStyle','--','Marker','s','LineWidth',4,'MarkerSize',10)
        ylim([0 20])
        xlim([2000 2014])
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        xlabel('Year')
        ylabel('Ratio')
        title('Figure 11. Ratio of Former:Past Year IDU by Age, US -- 2000-2014')
        set(gca,'LineWidth',4,'FontSize',24,'FontName','Calibri')
        legend('15-19','20-25','26-29','30-64')

%% Age Contact Matrix from Polymod (Mossong et al.)
    
% Run the script entitled PolymodPhysical_02182018.m to setup the contact matrices 
% used later (c, minc, and maxc).
run('PolymodPhysical_02182018.m')

%% Time Span and Age Group Population Sizes;

% Run model for time points specified;
tspan = [0 16];

% Average Michigan population size for each age group (from CDC Wonder);
popsize = [723319 815922 483070 4605115]';


%% Parameter Initial Values

% For reference, the parameters are:

% Estimated parameters: 
    % param1=[sigma1 sigma2 sigma3 sigma4]

% Sampled or set parameters: 
% param   =[sigma1 sigma2 sigma3 sigma4 
%              1      2     3       4          

%          a b  g(1:4) d  eps zeta(1:4) etan etap etaz theta(1:4) k(1:4) lambda(1:4)
%          5 6  7:10   11 12  13:16     17   18   19   20:23      24:27   28:31
%
%          mu(1:4) xi omic(1:4) r  tau phip phin psi(1:4) psin w  intervention_etap Z0];
%          32:35   36 37:40     41 42  43   44   45:48    49   50          51       52          
    
% The inverse of age group size is in vector nu:
    % nu = [nu1,nu2,nu3,nu4];

% Polymod-based contact matrix input:
    % pi=[pi11 pi12 pi13 pi14
    %     pi21 pi22 pi23 pi24
    %     pi31 pi32 pi33 pi34
    %     pi41 pi42 pi43 pi44];
    
%N[SNi; ANi; CNi; INi; TNi; Si; Ai; Ci; Ii; Ti; Z]      


%param1 initial values (parameters we estimate);
    b=0.000003;%3; %probability of infection given contact;
    theta1=0.05; %injection initiation rate (by age);
    theta2=0.034;
    theta3=0.034;

param1=[b theta1 theta2 theta3]';

%paramset values (set parameters)
    sigma1=5; %# syringe sharing partners
    sigma2=5;
    sigma3=5;
    sigma4=5;
    a=0.894; %proportion treated who achieve SVR;
    g1=1.2; %cessation rate;
    g2=1.2;
    g3=1.2;
    g4=1.2;
    d=0.3; %proportion of acutes who spontaneously resolve infection;
    eps=1/0.5; %1/duration of acute infection;
    zeta1=0.01; %prevalence of past year drug injection, year 2000;
    zeta2=0.01;
    zeta3=0.01;
    zeta4=0.01;
    etan=0.37; %SMR: former vs current PWID mortality;
    etap=9.8; %SMR: current PWID mortality vs general popln;
    etaz=4.36; %SMR: drug user mortality vs general popln;
    theta4=0.034; %injection initiation: 30-64;
    lambda1=0.1; %Initial HCV prevalence;
    lambda2=0.2;
    lambda3=0.3;
    lambda4=0.552; 
    mu1=0.000553; %Mortality rate; 
    mu2=0.000915; 
    mu3=0.00119; 
    mu4=0.00493; 
    xi=0.2; %proprtion of acutes who resolve infection that develop immunity;
    omic1=0.01%prevalence of past year drug injection, year 2000;
    omic2=0.01;
    omic3=0.021;
    omic4=0.027;
    r=1/50; %reporting rate (surveillance-detected case per true cases;
    tau=0; %proportion of treated individuals who contribute to transmission;
    phip=0; %treatment rate, current PWID;
    phin=0; %treatment rate, former PWID;
    psi1=0.018; %prevalence of past year dependence abuse, year 2000;
    psi2=0.013;
    psi3=0.0063;
    psi4=0.0035;
    psin=0.019; %average prevalence of past year dependence abuse among 15-19 yos;
    w=365/84; %1/duration of treatment;
    intervention_etap=1;
    Z0=144097; %Average size of 15 yo population for 2000-2016;
    

% nu = [nu1,nu2,nu3,nu4] as the inverse of age group size;
nu=[1/5;1/6;1/4;1/35];

% mid-period values for current and former IDU
F=formerinject(8,:);
C=pyinject(8,:);
    
%Calculate k based on other params using staedy state eqn for Fi
k1=g1*zeta1*popsize(1,1)/(omic1*popsize(1,1))-mu1*etap*etan-nu(1,1)
k2=(g2*zeta2*popsize(2,1)+nu(1,1)*omic1*popsize(1,1))/(omic2*popsize(2,1))-mu2*etap*etan-nu(2,1)
k3=(g3*zeta3*popsize(3,1)+nu(2,1)*omic2*popsize(2,1))/(omic3*popsize(3,1))-mu3*etap*etan-nu(3,1)
k4=(g4*zeta4*popsize(4,1)+nu(3,1)*omic3*popsize(3,1))/(omic4*popsize(4,1))-mu4*etap*etan-nu(4,1)
 
%Reset ki to 0.1 if negative from steady state equations above;
if k1<0
   k1=0.1
end

if k2<0
   k2=0.1
end

if k3<0
   k3=0.1
end

if k4<0
   k4=0.1
end

% combine non-estimated parameters
paramset=[sigma1 sigma2 sigma3 sigma4 a g1 g2 g3 g4 d eps zeta1 zeta2 zeta3 zeta4 etan etap etaz theta4 k1 k2 k3 k4 lambda1 lambda2 lambda3 lambda4 mu1 mu2 mu3 mu4 xi omic1 omic2 omic3 omic4 r tau phip phin psi1 psi2 psi3 psi4 psin w intervention_etap Z0]';
 

% combine param1 and paramset into 1 vector
param = [paramset(1:5);param1(1);paramset(6:18);param1(2:4);paramset(19:end)];



% check that etan*etap>1 (if <1 then former PWID mortality will be lower
% than mortality for general population (if etan*etap<1 then rescale to
% etan=1/etap so that etan*etap=1 and moratlity of former PWID same as
% general population

etan*etap

%% replace acute 15-19 with new chronics

% replacing acute 15-19 with scaled version of newchronic15-19 

acute1529idu(:,1)
acute1529idu(:,1)=(1/16.8)*chronic1529idu_model(:,1);
acute1529idu(:,1)

%% Initial Conditions

% Initial Condition Notes
    % SNi: former PWID prevalence - infecteds in N classes
    %      omici*popsizei-sum(ANi+CNi+INi+TNi);   
    % ANi: 0 (assume no acutes are former PWID at start;
    % CNi: former PWID prevalence * HCV prevalence * proportion of cases
    %      that become chronic
    %      omici*(1-d)*lambdai*popsizei
    % INi: former PWID prevalence * HCV prevalence * proportion of cases
    %      that resolve * proportion of cases that have sterilizing
    %      immunity
    %      omici*d*xi*lambdai*popsizei
    % TNi: 0 (no treatment initially)
    % Si: current PWID prevalence - infecteds in N classes
    %      zetai*popsizei-sum(Ai+Ci+Ii+Ti);
    % Ai: number of acute cases in 2000 corrected for under-reporting
    %     (1/r)*acute1529idu 
    % Ci: current PWID prevalence * HCV prevalence * proportion of cases
    %      that become chronic
    %      zetai*(1-d)*lambdai*popsizei
    % Ii: current PWID prevalence * HCV prevalence * proportion of cases
    %      that resolve * proportion of cases that have sterilizing
    %      immunity
    %      zetai*d*xi*lambdai*popsizei
    % Ti: 0 (no treatment at initially)
    % Zi: Past Year Abuse/Dependence Reporting No IDU year 2000
   
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

%% Transform contact matrix to proportions

% Calculate proportion of total contacts made with each age group
% (proportion of each column total for each row)

pi=zeros(4,4);
for i=1:4 %rows
    for j=1:4 %columns
    pi(i,j)=c(i,j)./sum(c(1:4,j));
    end
end

pi=reshape(pi,16,1);

%% Test the ODE

[t5,x5] = ode15s(@HCV_DiffEq_07262018,tspan,N0,[],param,nu,pi);

% Plot Initial Simulation Results (without parameter estimation) with Data;
for i=1:3;
    
    chronicnew=(eps*(1-d))*(x5(:,(11*i-4))+x5(:,11*i-9));
    acute=r*(x5(:,(11*i-4))+x5(:,11*i-9));
    chronicprev=(x5(:,11*i-3)+x5(:,11*i-8));
    formerpwid=x5(:,11*i-10)+x5(:,11*i-9)+x5(:,11*i-8)+x5(:,11*i-7)+x5(:,11*i-6);
    currentpwid=x5(:,11*i-5)+x5(:,11*i-4)+x5(:,11*i-3)+x5(:,11*i-2)+x5(:,11*i-1);
    nsduh=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14]; %nsduh time series
        
figure(i+11);
hold on
plot(t5,chronicprev,'b','LineWidth',2); %chronic prevalence
hold on
plot(t5,chronicnew,'k','LineWidth',2); %new chronic PWID
hold on
plot(time,chronic1529idu_model(:,i),'ko','LineWidth',2) %chronic data
hold on
title('Initial Conditions: HCV IDU Model')
hold on
legend('Prevalence of Chronic HCV: Model','New Chronic Cases: Model','New Chronic + Unconfirmed PWID Cases: Data')

figure(i+14)
hold on
plot(t5,chronicprev,'b','LineWidth',2); %chronic hcv prevalence
hold on
plot(t5,currentpwid,'k','LineWidth',2); %current pwid prevalence
hold on
plot(nsduh,pyinject(:,i)*popsize(i,1),'ko','LineWidth',2)
hold on
plot(t5,formerpwid,'r','LineWidth',2) %former pwid prevalence
hold on
plot(nsduh,formerinject(:,i)*popsize(i,1),'r*','LineWidth',2)
hold on
title('IC: PWID and HCV prevalence')
hold on
legend('HCV Prevalence: Model','Current PWID: Model','Current PWID: NSDUH','Former PWID: Model','Former PWID: NSDUH')

figure(i+17);
hold on
plot(t5,acute,'k','LineWidth',2) %acutes 
hold on
plot(time,acute1529idu(:,i),'k*','LineWidth',2)
hold on
title('Initial Conditions: HCV IDU Model')
hold on
legend('Total Acute Cases: Model','Acute: Data')

end

%final compartment (30-64)
    chronicnew=(eps*(1-d))*(x5(:,(11*4-4))+x5(:,11*4-9));
    acute=r*(x5(:,(11*4-4))+x5(:,11*4-9));
    chronicprev=(x5(:,11*4-3)+x5(:,11*4-8));
    formerpwid=x5(:,11*4-10)+x5(:,11*4-9)+x5(:,11*4-8)+x5(:,11*4-7)+x5(:,11*4-6);
    currentpwid=x5(:,11*4-5)+x5(:,11*4-4)+x5(:,11*4-3)+x5(:,11*4-2)+x5(:,11*4-1);
    
figure(21);
hold on
plot(t5,chronicprev,'b','LineWidth',2); %chronic prevalence
hold on
plot(t5,chronicnew,'k','LineWidth',2); %new chronic PWID
hold on
title('Initial Conditions: HCV IDU Model')
hold on
legend('Prevalence of Chronic HCV: Model','New Chronic Cases: Model')

figure(22)
hold on
plot(t5,chronicprev,'b','LineWidth',2); %chronic hcv prevalence
hold on
plot(t5,currentpwid,'k','LineWidth',2); %current pwid prevalence
hold on
plot(nsduh,pyinject(:,4)*popsize(4,1),'ko','LineWidth',2)
hold on
plot(t5,formerpwid,'r','LineWidth',2) %former pwid prevalence
hold on
plot(nsduh,formerinject(:,4)*popsize(4,1),'r*','LineWidth',2)
hold on
title('IC: PWID and HCV prevalence')
hold on
legend('HCV Prevalence: Model','Current PWID: Model','Current PWID: NSDUH','Former PWID: Model','Former PWID: NSDUH')

figure(23);
hold on
plot(t5,acute,'k','LineWidth',2) %acutes 
hold on
title('Initial Conditions: HCV IDU Model')
hold on
legend('Total Acute Cases: Model')

%% One Parameter Estimation for Base Params Set Above
% Estimate 3 Injection Initiation (Theta), 1 Transmisison Rate (Beta)

y=zeros(1,3);


% set upper and lower bounds for param1 estimated parameters
x1=[0   0.05   0.034   0.034]' %lower bounds for param1 (beta, theta1, theta2, theta3);
x2=[0.1 5      5      10]' %upper bounds for param1 (beta, theta1, theta2, theta3);

% Fit to data using bounds set above as lower/upper limits for est params
[paramest RSS] = fminsearchbnd(@(param1) HCV_LS_07262018(time,N0,param1,paramset,nu,pi,popsize,acute1529idu),param1,x1,x2,optimset('Display', 'iter','MaxFunEvals',5000));

% Fit to data without restricting estimated parameter upper/lower bounds
% [paramest RSS] = fminsearch(@(param1) HCV_LS_07262018(time,N0,param1,paramset,nu,pi,popsize,acute1529idu),param1,optimset('Display', 'iter','MaxFunEvals',3000));

paramest=abs(paramest);
  
%create new parameter matrix with estimated paramest
paramnew = [paramset(1:5);paramest(1); paramset(6:18); paramest(2:4); paramset(19:end)];

RSS

% Calculate the AIC= 2(-LL+#parameters est) = 2(RSS+#parameters est)
AIC=2*(RSS+4)

%% Simulate the ODE using the new fitted parameters

[t1,x1] = ode15s(@HCV_DiffEq_07262018,tspan,N0,[],paramnew,nu,pi);

% Plot Initial Simulation Results (without parameter estimation) with Data;
for i=1:3;
    
    chronicnew=(eps*(1-d))*(x1(:,(11*i-4))+x1(:,11*i-9));
    acute=paramnew(41)*(x1(:,(11*i-4))+x1(:,11*i-9));
    chronicprev=(x1(:,11*i-3)+x1(:,11*i-8));
    formerpwid=x1(:,11*i-10)+x1(:,11*i-9)+x1(:,11*i-8)+x1(:,11*i-7)+x1(:,11*i-6);
    currentpwid=x1(:,11*i-5)+x1(:,11*i-4)+x1(:,11*i-3)+x1(:,11*i-2)+x1(:,11*i-1);
    nsduh=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14]; %nsduh time series
        
figure(i+23);
hold on
plot(t1,chronicprev,'b','LineWidth',2); %chronic prevalence
hold on
plot(t1,chronicnew,'k','LineWidth',2); %new chronic PWID
hold on
plot(time,chronic1529idu_model(:,i),'ko','LineWidth',2) %chronic data
hold on
title('Chronics')
hold on
legend('Prevalence of Chronic HCV: Model','New Chronic Cases: Model','New Chronic + Unconfirmed PWID Cases: Data')

figure(i+25)
hold on
plot(t1,chronicprev,'b','LineWidth',2); %chronic hcv prevalence 
hold on
plot(t1,currentpwid,'k','LineWidth',2); %current pwid prevalence
hold on
plot(nsduh,pyinject(:,i)*popsize(i,1),'ko','LineWidth',2)
hold on
plot(t1,formerpwid,'r','LineWidth',2) %former pwid prevalence
hold on
plot(nsduh,formerinject(:,i)*popsize(i,1),'r*','LineWidth',2)
hold on
title('PWID and HCV prevalence')
hold on
legend('HCV Prevalence: Model','Current PWID: Model','Current PWID: NSDUH','Former PWID: Model','Former PWID: NSDUH')

figure(i+28);
hold on
plot(t1,acute,'k','LineWidth',2) %acutes 
hold on
plot(time,acute1529idu(:,i),'k*','LineWidth',2)
hold on
title('Acutes')
hold on
legend('Total Acute Cases: Model','Acute: Data')
end

%final compartment
    chronicnew=(eps*(1-d))*(x1(:,(11*4-4))+x1(:,11*4-9));
    acute=paramnew(41)*(x1(:,(11*4-4))+x1(:,11*4-9));
    chronicprev=(x1(:,11*4-3)+x1(:,11*4-8));
    formerpwid=x1(:,11*4-10)+x1(:,11*4-9)+x1(:,11*4-8)+x1(:,11*4-7)+x1(:,11*4-6);
    currentpwid=x1(:,11*4-5)+x1(:,11*4-4)+x1(:,11*4-3)+x1(:,11*4-2)+x1(:,11*4-1);
    
figure(32);
hold on
plot(t1,chronicprev,'b','LineWidth',2); %chronic prevalence
hold on
plot(t1,chronicnew,'k','LineWidth',2); %new chronic PWID
hold on
title('Initial Conditions: HCV IDU Model')
hold on
legend('Prevalence of Chronic HCV: Model','New Chronic Cases: Model')

figure(33)
hold on
plot(t1,chronicprev,'b','LineWidth',2); %chronic hcv prevalence 
hold on
plot(t1,currentpwid,'k','LineWidth',2); %current pwid prevalence
hold on
plot(nsduh,pyinject(:,4)*popsize(4,1),'ko','LineWidth',2)
hold on
plot(t1,formerpwid,'r','LineWidth',2) %former pwid prevalence
hold on
plot(nsduh,formerinject(:,4)*popsize(4,1),'r*','LineWidth',2)
hold on
title('IC: PWID and HCV prevalence')
hold on
legend('HCV Prevalence: Model','Current PWID: Model','Current PWID: NSDUH','Former PWID: Model','Former PWID: NSDUH')

figure(34);
hold on
plot(t1,acute,'k','LineWidth',2) %acutes 
hold on
title('Acutes')
hold on
legend('Total Acute Cases: Model','Acute: Data')

%% save workspace for Latin Hypercube Sampling (LHS)

save('HCVLHS_Setup_04162018.mat','chronic1529idu_model','minc','maxc','time','nu','pyinjectmin','pyinjectmax','formerinjectmin','formerinjectmax','nsduh','acute1529idu','pyinject','formerinject')

%% Run Latin Hypercube Sample: Parameter Ranges and Scaling

% note that this script implements the latin hypercube sampling, currently
% set to re-run the ODE and fit to data for 10,000 simulations using 8
% nodes/workers (run on a high performance computing cluster for best results)

run('LHS07262018.m')

%% Load the Latin Hypercube Sampling Results

load('HCVLHS05062018.mat')

%% Identify the Best Fitting Parameter Set

ParamNA=[];

ParamEstRuns1=[];
RSSRuns1=[];

%Identify Na runs that did not converge
for j=1:10000

if isnan(RSSRuns(1,j))
    ParamNA=[ParamNA, ParamEstRuns(:,j)];

else 
    ParamEstRuns1=[ParamEstRuns1, ParamEstRuns(:,j)];
    RSSRuns1=[RSSRuns1, RSSRuns(1,j)];
    
end
end

% output run index with minimum residual sum of squares
[minRSS,I]=min(RSSRuns1)
Best=[minRSS,I];
Best=Best(1,2);

paramscaling(Best,:);
bestp=ParamEstRuns(:,Best);
bestc=reshape(CMatrixRuns(Best,:),4,4);

[maxRSS,J]=max(RSSRuns1)
medianRSS=median(RSSRuns1)
meanRSS=mean(RSSRuns1)
sdRSS=std(RSSRuns1)

%% Examine estimated parameter values
summaryPE=[];

for i=1:52
    i
    [maxPE,J]=max(ParamEstRuns1(i,:))
    [minPE,R]=min(ParamEstRuns1(i,:))

    summaryPE(i,1)=max(ParamEstRuns1(i,:));
    summaryPE(i,2)=min(ParamEstRuns1(i,:));
    summaryPE(i,3)=median(ParamEstRuns1(i,:));
    summaryPE(i,4)=mean(ParamEstRuns1(i,:));
    summaryPE(i,5)=std(ParamEstRuns1(i,:));

end

% Contact matrix

for i=1:16
    summaryPE(i+52,1)=max(CMatrixRuns(:,i));
    summaryPE(i+52,2)=min(CMatrixRuns(:,i));
    summaryPE(i+52,3)=median(CMatrixRuns(:,i));
    summaryPE(i+52,4)=mean(CMatrixRuns(:,i));
    summaryPE(i+56,5)=std(CMatrixRuns(:,i));

end

%% Examine parameter sets where model did not converge
%Note: not applicable for current simulation, all sets converged
summaryNA=[];

for i=1:52
    i
    [maxPE,J]=max(ParamNA(i,:))
    [minPE,R]=min(ParamNA(i,:))

    summaryNA(i,1)=max(ParamNA(i,:));
    summaryNA(i,2)=min(ParamNA(i,:));
    summaryNA(i,3)=median(ParamNA(i,:));
    summaryNA(i,4)=mean(ParamNA(i,:));
    summaryNA(i,5)=std(ParamNA(i,:));

end


%% Parameter and Fit Figures

% RSS
figure(35)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(RSSRuns)
hold on
ylabel('Number of Runs')
hold on
xlabel('Residual Sum of Squares')
hold on
title('Residual Sum of Squares')

% theta1 Estimates
figure(36)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(ParamEstRuns(20,:))
hold on
plot(ParamEstRuns(20,Best),1,'r*','LineWidth',2,'MarkerSize',8)
hold on
ylabel('Number of Runs')
hold on
xlabel('Estimate')
hold on
title('Injection Initiation: 15-19')
legend('Theta1','Best Fit')

% Theta2 Estimates
figure(37)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(ParamEstRuns(21,:))
hold on
plot(ParamEstRuns(21,Best),1,'r*','LineWidth',2,'MarkerSize',8)
hold on
ylabel('Number of Runs')
hold on
xlabel('Estimate')
hold on
title('Injection Initiation: 20-25')
legend('Theta2','Best Fit')

% Theta3 Estimates
figure(38)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(ParamEstRuns(22,:))
hold on
plot(ParamEstRuns(22,Best),1,'r*','LineWidth',2,'MarkerSize',8)
hold on
ylabel('Number of Runs')
hold on
xlabel('Estimate')
hold on
title('Injection Initiation: 26-29')
legend('Theta3','Best Fit')

% Beta

figure(39)
set(gca,'LineWidth',2,'FontSize',20,'FontName','Calibri')
    hold on
histogram(ParamEstRuns1(6,:),'Normalization','probability','BinWidth',1e-04)
hold on
ylabel('Number of Runs')
hold on
xlabel('Beta')
hold on
title('Beta')
legend('Est')

%% LHS Case Count Figures: All Age Groups (Combined Results)
time1 = [0:0.2:16];
r=ParamEstRuns(41,:);


acuteprev=[]
acuteprevformer=[]
acuteprevcurrent=[]

for j=1:length(AcuteRuns)
    
    if isnan(RSSRuns(j))
    acuteprev(:,j)=nan;
    acuteprevformer(:,j)=nan;
    acuteprevcurrent(:,j)=nan;  
    else
    acuteprev(:,j)=r(1,j)*(Acute_Age1(:,j)+Acute_Age2(:,j)+Acute_Age3(:,j));
    acuteprevformer(:,j)=r(1,j)*(Acute_Age1NPWID(:,j)+Acute_Age2NPWID(:,j)+Acute_Age3NPWID(:,j));
    acuteprevcurrent(:,j)=r(1,j)*(Acute_Age1PWID(:,j)+Acute_Age2PWID(:,j)+Acute_Age3PWID(:,j));
    end
end

% Chronic and Acute Current and Former PWID;
figure(40)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
%      p1=plot(time1,NewChronic_Age1+NewChronic_Age2+NewChronic_Age3,'Color',[0.7 0.7 0.7],'LineWidth',2);
%          hold on
%      p2=plot(time1,NewChronic_Age1(:,Best)+NewChronic_Age2(:,Best)+NewChronic_Age3(:,Best),'k','LineWidth',2);
%         hold on
    p3=plot(time1,acuteprev,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,acuteprev(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
%     p5=plot(time,chronic1529idu_model(:,1)+chronic1529idu_model(:,2)+chronic1529idu_model(:,3),'ko','LineWidth',2,'MarkerSize',8)
%         hold on
    p6=plot(time,acute1529idu(:,1)+acute1529idu(:,2)+acute1529idu(:,3),'rd','LineWidth',2,'MarkerSize',8)
%         hold on
%     p7=plot(time,acute1529idu(:,1)+acute1529idu(:,2)+acute1529idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 16])
    ylabel('Number of Cases')
    title('LHS: Acute ')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    %legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
    legend([p3(1) p4(1) p6(1)], 'Acute: Model','Acute: Best Fit','Acute: Data')
    %saveas(gcf,'AcuteNewChronic_All','fig')

% Chronic Current and Former PWID
figure(41)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,NewChronic_Age1+NewChronic_Age2+NewChronic_Age3,'Color',[0.7 0.7 0.7],'LineWidth',2);
    hold on
    p2=plot(time1,NewChronic_Age1(:,Best)+NewChronic_Age2(:,Best)+NewChronic_Age3(:,Best),'k','LineWidth',2);
    hold on
    p3=plot(time,chronic1529idu_model(:,1)+chronic1529idu_model(:,2)+chronic1529idu_model(:,3),'ko','LineWidth',2,'MarkerSize',8)
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 16])
    ylabel('Number of Cases')
    title('LHS Fit to Chronic PWID')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1) p3(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','New Chronic + Unconfirmed PWID Cases: Data')
%saveas(gcf,'NewChronic_All','fig')

% Chronic Prevalence (Current and Former PWID)
figure(42)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,ChronicPrevRuns,'Color',[0.7 0.7 0.85],'LineWidth',2);
    hold on
    p2=plot(time1,ChronicPrevRuns(:,Best),'Color',[0 0 1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 16])
    ylabel('Number of Cases')
    title('LHS Chronic Prevalence')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
       legend([p1(1) p2(1)], 'Chronic Prevalence: Model','Chronic Prevalence: Best Fit')
%saveas(gcf,'ChronicPrev_All','fig')

% Susceptible Current and Former PWID
figure(43)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
    p1=plot(time1,Susc,'Color',[0.8 0.9 0.8],'LineWidth',2);
    hold on
    p2=plot(time1,Susc(:,Best),'Color',[0 0.4 0.1],'LineWidth',2);
    hold on
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 16])
    ylabel('Number of Cases')
    title('LHS Susceptible PWID')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    legend([p1(1) p2(1)], 'Susceptible PWID: Model','Susceptible PWID: Best Fit')
%saveas(gcf,'Susc_All','fig')

%% Age 15-19 Results
r=ParamEstRuns(41,:)
nsduh=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14]; %nsduh time series

    chronicprev=[];
    chronicprevformer=[];
    chronicprevcurrent=[];
    acuteprev=[];
    acuteprevformer=[];
    acuteprevcurrent=[];
    formerpwid=[];
    currentpwid=[];
    FormerInjectPopRuns=formerinject;
    PYInjectPopRuns=pyinject;   
%15-19;
i=1;
for j=1:length(ChronicPrev_Age1);
    
    chronicprev(:,j)=ChronicPrev_Age1(:,j);
    chronicprevformer(:,j)=ChronicPrev_Age1NPWID(:,j);
    chronicprevcurrent(:,j)=ChronicPrev_Age1PWID(:,j);
    acuteprev(:,j)=r(1,j)*Acute_Age1(:,j);
    acuteprevformer(:,j)=r(1,j)*Acute_Age1NPWID(:,j);
    acuteprevcurrent(:,j)=r(1,j)*Acute_Age1PWID(:,j);
    formerpwid=Susc_Age1NPWID(:,j)+Acute_Age1NPWID(:,j)+chronicprevformer+Immune_Age1NPWID(:,j);
    currentpwid=Susc_Age1PWID(:,j)+Acute_Age1PWID(:,j)+chronicprevcurrent+Immune_Age1PWID(:,j);

    
end
    
figure(i+43)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
%      p1=plot(time1,NewChronic_Age1+NewChronic_Age2+NewChronic_Age3,'Color',[0.7 0.7 0.7],'LineWidth',2);
%          hold on
%      p2=plot(time1,NewChronic_Age1(:,Best)+NewChronic_Age2(:,Best)+NewChronic_Age3(:,Best),'k','LineWidth',2);
%         hold on
    p3=plot(time1,acuteprev,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,acuteprev(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
%     p5=plot(time,chronic1529idu_model(:,1)+chronic1529idu_model(:,2)+chronic1529idu_model(:,3),'ko','LineWidth',2,'MarkerSize',8)
%         hold on
    p6=plot(time,acute1529idu(:,1),'rd','LineWidth',2,'MarkerSize',8)
%         hold on
%     p7=plot(time,acute1529idu(:,1)+acute1529idu(:,2)+acute1529idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 16])
    ylabel('Number of Cases')
    title('LHS: Acute 15-19 ')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    %legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
    legend([p3(1) p4(1) p6(1)], 'Acute: Model','Acute: Best Fit','Acute: Data')
    %saveas(gcf,'AcuteNewChronic_All','fig')

figure(i+44)
hold on
p1=plot(time1,chronicprev,'b','LineWidth',2); %chronic hcv prevalence 
hold on
p2=plot(time1,currentpwid,'Color',[0.7 0.7 0.7],'LineWidth',2); %current pwid prevalence
hold on
p3=plot(time1,formerpwid,'Color',[1 0.8 0.8],'LineWidth',2); %former pwid prevalence
hold on
 p4=plot(nsduh,formerinject(:,i)*PopSizeRuns(Best,i),'r*','LineWidth',2)
 hold on
p5=plot(nsduh,pyinject(:,i)*PopSizeRuns(Best,i),'ko','LineWidth',2)
hold on
title('PWID and HCV prevalence: 15-19');
hold on
legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'HCV Prevalence: Model','Current PWID: Model','Former PWID: Model','Former PWID: NSDUH','Current PWID: NSDUH');

figure(i+45)
hold on
p1=plot(time1,chronicprev(:,Best),'b','LineWidth',2); %chronic hcv prevalence 
hold on
p2=plot(time1,currentpwid(:,Best),'k','LineWidth',2); %current pwid prevalence
hold on
p3=plot(time1,formerpwid(:,Best),'r','LineWidth',2) %former pwid prevalence
hold on
 p4=plot(nsduh,formerinject(:,i)*PopSizeRuns(Best,i),'r*','LineWidth',2)
 hold on
p5=plot(nsduh,pyinject(:,i)*PopSizeRuns(Best,i),'ko','LineWidth',2)
hold on
title('PWID and HCV prevalence: 15-19')
hold on
legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'HCV Prevalence: Model','Current PWID: Model','Former PWID: Model','Former PWID: NSDUH','Current PWID: NSDUH')

%% Age 20-25 Results
%20-25;
i=2;
for j=1:length(ChronicPrev_Age1)
    
    chronicprev(:,j)=ChronicPrev_Age2(:,j);
    chronicprevformer(:,j)=ChronicPrev_Age2NPWID(:,j);
    chronicprevcurrent(:,j)=ChronicPrev_Age2PWID(:,j);
    acuteprev(:,j)=r(1,j)*Acute_Age2(:,j);
    acuteprevformer(:,j)=r(1,j)*Acute_Age2NPWID(:,j);
    acuteprevcurrent(:,j)=r(1,j)*Acute_Age2PWID(:,j);
    formerpwid=Susc_Age2NPWID(:,j)+Acute_Age2NPWID(:,j)+chronicprevformer+Immune_Age2NPWID(:,j);
    currentpwid=Susc_Age2PWID(:,j)+Acute_Age2PWID(:,j)+chronicprevcurrent+Immune_Age2PWID(:,j);

    
end
    
figure(i+43)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
%      p1=plot(time1,NewChronic_Age1+NewChronic_Age2+NewChronic_Age3,'Color',[0.7 0.7 0.7],'LineWidth',2);
%          hold on
%      p2=plot(time1,NewChronic_Age1(:,Best)+NewChronic_Age2(:,Best)+NewChronic_Age3(:,Best),'k','LineWidth',2);
%         hold on
    p3=plot(time1,acuteprev,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,acuteprev(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
%     p5=plot(time,chronic1529idu_model(:,1)+chronic1529idu_model(:,2)+chronic1529idu_model(:,3),'ko','LineWidth',2,'MarkerSize',8)
%         hold on
    p6=plot(time,acute1529idu(:,2),'rd','LineWidth',2,'MarkerSize',8)
%         hold on
%     p7=plot(time,acute1529idu(:,1)+acute1529idu(:,2)+acute1529idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 16])
    ylabel('Number of Cases')
    title('LHS: Acute 20-25 ')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    %legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
    legend([p3(1) p4(1) p6(1)], 'Acute: Model','Acute: Best Fit','Acute: Data')
    %saveas(gcf,'AcuteNewChronic_All','fig')


figure(i+44)
hold on
p1=plot(time1,chronicprev,'b','LineWidth',2); %chronic hcv prevalence 
hold on
p2=plot(time1,currentpwid,'Color',[0.7 0.7 0.7],'LineWidth',2); %current pwid prevalence
hold on
p3=plot(time1,formerpwid,'Color',[1 0.8 0.8],'LineWidth',2); %former pwid prevalence
hold on
 p4=plot(nsduh,formerinject(:,i)*PopSizeRuns(Best,i),'r*','LineWidth',2)
 hold on
p5=plot(nsduh,pyinject(:,i)*PopSizeRuns(Best,i),'ko','LineWidth',2)
hold on
title('PWID and HCV prevalence: 20-25');
hold on
legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'HCV Prevalence: Model','Current PWID: Model','Former PWID: Model','Former PWID: NSDUH','Current PWID: NSDUH')

figure(i+45)
hold on
p1=plot(time1,chronicprev(:,Best),'b','LineWidth',2); %chronic hcv prevalence 
hold on
p2=plot(time1,currentpwid(:,Best),'k','LineWidth',2); %current pwid prevalence
hold on
p3=plot(time1,formerpwid(:,Best),'r','LineWidth',2) %former pwid prevalence
hold on
 p4=plot(nsduh,formerinject(:,i)*PopSizeRuns(Best,i),'r*','LineWidth',2)
 hold on
p5=plot(nsduh,pyinject(:,i)*PopSizeRuns(Best,i),'ko','LineWidth',2)
hold on
title('PWID and HCV prevalence: 20-25')
hold on
legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'HCV Prevalence: Model','Current PWID: Model','Former PWID: Model','Former PWID: NSDUH','Current PWID: NSDUH')
%% Age 26-29 Results
%26-29;
i=3;
for j=1:length(ChronicPrev_Age1)
    
    chronicprev(:,j)=ChronicPrev_Age3(:,j);
    chronicprevformer(:,j)=ChronicPrev_Age3NPWID(:,j);
    chronicprevcurrent(:,j)=ChronicPrev_Age3PWID(:,j);
    acuteprev(:,j)=r(1,j)*Acute_Age3(:,j);
    acuteprevformer(:,j)=r(1,j)*Acute_Age3NPWID(:,j);
    acuteprevcurrent(:,j)=r(1,j)*Acute_Age3PWID(:,j);
    formerpwid=Susc_Age3NPWID(:,j)+Acute_Age3NPWID(:,j)+chronicprevformer+Immune_Age3NPWID(:,j);
    currentpwid=Susc_Age3PWID(:,j)+Acute_Age3PWID(:,j)+chronicprevcurrent+Immune_Age3PWID(:,j);

end
    
figure(i+43)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
%      p1=plot(time1,NewChronic_Age1+NewChronic_Age2+NewChronic_Age3,'Color',[0.7 0.7 0.7],'LineWidth',2);
%          hold on
%      p2=plot(time1,NewChronic_Age1(:,Best)+NewChronic_Age2(:,Best)+NewChronic_Age3(:,Best),'k','LineWidth',2);
%         hold on
    p3=plot(time1,acuteprev,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,acuteprev(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
%     p5=plot(time,chronic1529idu_model(:,1)+chronic1529idu_model(:,2)+chronic1529idu_model(:,3),'ko','LineWidth',2,'MarkerSize',8)
%         hold on
    p6=plot(time,acute1529idu(:,3),'rd','LineWidth',2,'MarkerSize',8)
%         hold on
%     p7=plot(time,acute1529idu(:,1)+acute1529idu(:,2)+acute1529idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 16])
    ylabel('Number of Cases')
    title('LHS: Acute 26-29 ')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    %legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
    legend([p3(1) p4(1) p6(1)], 'Acute: Model','Acute: Best Fit','Acute: Data')
    %saveas(gcf,'AcuteNewChronic_All','fig')

figure(i+44)
hold on
p1=plot(time1,chronicprev,'b','LineWidth',2); %chronic hcv prevalence 
hold on
p2=plot(time1,currentpwid,'Color',[0.7 0.7 0.7],'LineWidth',2); %current pwid prevalence
hold on
p3=plot(time1,formerpwid,'Color',[1 0.8 0.8],'LineWidth',2); %former pwid prevalence
hold on
 p4=plot(nsduh,formerinject(:,i)*PopSizeRuns(Best,i),'r*','LineWidth',2)
 hold on
p5=plot(nsduh,pyinject(:,i)*PopSizeRuns(Best,i),'ko','LineWidth',2)
hold on
title('PWID and HCV prevalence: 26-29');
hold on
legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'HCV Prevalence: Model','Current PWID: Model','Former PWID: Model','Former PWID: NSDUH','Current PWID: NSDUH')

figure(i+45)
hold on
p1=plot(time1,chronicprev(:,Best),'b','LineWidth',2); %chronic hcv prevalence 
hold on
p2=plot(time1,currentpwid(:,Best),'k','LineWidth',2); %current pwid prevalence
hold on
p3=plot(time1,formerpwid(:,Best),'r','LineWidth',2) %former pwid prevalence
hold on
 p4=plot(nsduh,formerinject(:,i)*PopSizeRuns(Best,i),'r*','LineWidth',2)
 hold on
p5=plot(nsduh,pyinject(:,i)*PopSizeRuns(Best,i),'ko','LineWidth',2)
hold on
title('PWID and HCV prevalence: 26-29')
hold on
legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'HCV Prevalence: Model','Current PWID: Model','Former PWID: Model','Former PWID: NSDUH','Current PWID: NSDUH')
%% Age 30-64 Results
%30-64;
i=4;
for j=1:length(ChronicPrev_Age1)

    chronicprev(:,j)=ChronicPrev_Age4(:,j);
    chronicprevformer(:,j)=ChronicPrev_Age4NPWID(:,j);
    chronicprevcurrent(:,j)=ChronicPrev_Age4PWID(:,j);
    acuteprev(:,j)=Acute_Age4(:,j);
    acuteprevformer(:,j)=Acute_Age3NPWID(:,j);
    acuteprevcurrent(:,j)=Acute_Age3PWID(:,j);
    formerpwid=Susc_Age4NPWID(:,j)+Acute_Age4NPWID(:,j)+chronicprevformer+Immune_Age4NPWID(:,j);
    currentpwid=Susc_Age4PWID(:,j)+Acute_Age4PWID(:,j)+chronicprevcurrent+Immune_Age4PWID(:,j);
    
end
    
figure(i+43)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    hold on
%      p1=plot(time1,NewChronic_Age1+NewChronic_Age2+NewChronic_Age3,'Color',[0.7 0.7 0.7],'LineWidth',2);
%          hold on
%      p2=plot(time1,NewChronic_Age1(:,Best)+NewChronic_Age2(:,Best)+NewChronic_Age3(:,Best),'k','LineWidth',2);
%         hold on
    p3=plot(time1,acuteprev,'Color',[1 0.8 0.8],'LineWidth',2);
        hold on
    p4=plot(time1,acuteprev(:,Best),'Color',[1 0 0],'LineWidth',2);
        hold on
%     p5=plot(time,chronic1529idu_model(:,1)+chronic1529idu_model(:,2)+chronic1529idu_model(:,3),'ko','LineWidth',2,'MarkerSize',8)
%         hold on
    %p6=plot(time,acute1529idu(:,3),'rd','LineWidth',2,'MarkerSize',8)
%         hold on
%     p7=plot(time,acute1529idu(:,1)+acute1529idu(:,2)+acute1529idu(:,3),'b*','LineWidth',2,'MarkerSize',8)
set(gca,'LineWidth',2,'FontSize',18,'FontName','Calibri')
    xlabel('Year')
    xlim([0 16])
    ylabel('Number of Cases')
    title('LHS: Acute 30-64 ')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Calibri')
    %legend([p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1)], 'Chronic Incidence: Model','Chronic Incidence: Best Fit','Acute: Model','Acute: Best Fit','New Chronic + Unconfirmed PWID Cases: Data','Acute: CDC Correction','Acute: Data')
    legend([p3(1) p4(1) ], 'Acute: Model','Acute: Best Fit')
    %saveas(gcf,'AcuteNewChronic_All','fig')

figure(i+44)
hold on
p1=plot(time1,chronicprev,'b','LineWidth',2); %chronic hcv prevalence 
hold on
p2=plot(time1,currentpwid,'Color',[0.7 0.7 0.7],'LineWidth',2); %current pwid prevalence
hold on
p3=plot(time1,formerpwid,'Color',[1 0.8 0.8],'LineWidth',2); %former pwid prevalence
hold on
 p4=plot(nsduh,formerinject(:,i)*PopSizeRuns(Best,i),'r*','LineWidth',2)
 hold on
p5=plot(nsduh,pyinject(:,i)*PopSizeRuns(Best,i),'ko','LineWidth',2)
hold on
title('PWID and HCV prevalence: 30-64');
hold on
legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'HCV Prevalence: Model','Current PWID: Model','Former PWID: Model','Former PWID: NSDUH','Current PWID: NSDUH')

figure(i+45)
hold on
p1=plot(time1,chronicprev(:,Best),'b','LineWidth',2); %chronic hcv prevalence 
hold on
p2=plot(time1,currentpwid(:,Best),'k','LineWidth',2); %current pwid prevalence
hold on
p3=plot(time1,formerpwid(:,Best),'r','LineWidth',2) %former pwid prevalence
hold on
 p4=plot(nsduh,formerinject(:,i)*PopSizeRuns(Best,i),'r*','LineWidth',2)
 hold on
p5=plot(nsduh,pyinject(:,i)*PopSizeRuns(Best,i),'ko','LineWidth',2)
hold on
title('PWID and HCV prevalence: 30-64')
hold on
legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'HCV Prevalence: Model','Current PWID: Model','Former PWID: Model','Former PWID: NSDUH','Current PWID: NSDUH')


%% Compile parameters, contact matrix, RSS, and initial parameters 
LHS_ParamRSS = [];
LHS_ParamRSS(:,1) = 1:10000; %create a param indexing variable
LHS_ParamRSS(:,2:53) = ParamEstRuns1';
LHS_ParamRSS(:,54:69) = CMatrixRuns;
LHS_ParamRSS(:,70) = RSSRuns1';
LHS_ParamRSS(:,71:76) = InitialParamRuns;
LHS_ParamRSS(:,77:81) = FormerInjectRuns;
LHS_ParamRSS(:,82:86) = PYInjectRuns;
LHS_ParamRSS(:,87:90) = PopSizeRuns;

%% Save Results as TXT or CSV files

dlmwrite('LHS_ParamRSS.txt', [LHS_ParamRSS], 'precision', 10); 

    dlmwrite('LHS_AcuteNPWID.csv',[LHS_ParamRSS(:,1)'; AcuteNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcuteNPWID1.csv',[LHS_ParamRSS(:,1)';Acute_Age1NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcuteNPWID2.csv',[LHS_ParamRSS(:,1)';Acute_Age2NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcuteNPWID3.csv',[LHS_ParamRSS(:,1)';Acute_Age3NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcuteNPWID4.csv',[LHS_ParamRSS(:,1)';Acute_Age4NPWID], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_AcutePWID.csv',[LHS_ParamRSS(:,1)';AcutePWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcutePWID1.csv',[LHS_ParamRSS(:,1)';Acute_Age1PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcutePWID2.csv',[LHS_ParamRSS(:,1)';Acute_Age2PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcutePWID3.csv',[LHS_ParamRSS(:,1)';Acute_Age3PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_AcutePWID4.csv',[LHS_ParamRSS(:,1)';Acute_Age4PWID], 'delimiter', ',', 'precision', 10);
    
dlmwrite('LHS_Acute.csv',[LHS_ParamRSS(:,1)';AcuteRuns], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Acute1.csv',[LHS_ParamRSS(:,1)';Acute_Age1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Acute2.csv',[LHS_ParamRSS(:,1)';Acute_Age2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Acute3.csv',[LHS_ParamRSS(:,1)';Acute_Age3], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Acute4.csv',[LHS_ParamRSS(:,1)';Acute_Age4], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_ChronicPrevNPWID.csv',[LHS_ParamRSS(:,1)';ChronicPrevNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevNPWID1.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age1NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevNPWID2.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age2NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevNPWID3.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age3NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevNPWID4.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age4NPWID], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_ChronicPrevPWID.csv',[LHS_ParamRSS(:,1)';ChronicPrevPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevPWID1.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age1PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevPWID2.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age2PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevPWID3.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age3PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrevPWID4.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age4PWID], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_ChronicPrev.csv',[LHS_ParamRSS(:,1)';ChronicPrevRuns], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrev1.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrev2.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrev3.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age3], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ChronicPrev4.csv',[LHS_ParamRSS(:,1)';ChronicPrev_Age4], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_NewChronic.csv',[LHS_ParamRSS(:,1)';NewChronic], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronic1.csv',[LHS_ParamRSS(:,1)';NewChronic_Age1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronic2.csv',[LHS_ParamRSS(:,1)';NewChronic_Age2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronic3.csv',[LHS_ParamRSS(:,1)';NewChronic_Age3], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronic4.csv',[LHS_ParamRSS(:,1)';NewChronic_Age4], 'delimiter', ',', 'precision', 10);
    
    dlmwrite('LHS_NewChronicNPWID.csv',[LHS_ParamRSS(:,1)';NewChronicNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicNPWID1.csv',[LHS_ParamRSS(:,1)';NewChronic_Age1NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicNPWID2.csv',[LHS_ParamRSS(:,1)';NewChronic_Age2NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicNPWID3.csv',[LHS_ParamRSS(:,1)';NewChronic_Age3NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicNPWID4.csv',[LHS_ParamRSS(:,1)';NewChronic_Age4NPWID], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_NewChronicPWID.csv',[LHS_ParamRSS(:,1)';NewChronicPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicPWID1.csv',[LHS_ParamRSS(:,1)';NewChronic_Age1PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicPWID2.csv',[LHS_ParamRSS(:,1)';NewChronic_Age2PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicPWID3.csv',[LHS_ParamRSS(:,1)';NewChronic_Age3PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_NewChronicPWID4.csv',[LHS_ParamRSS(:,1)';NewChronic_Age4PWID], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_Susc.csv',[LHS_ParamRSS(:,1)';Susc], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Susc1.csv',[LHS_ParamRSS(:,1)';Susc_Age1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Susc2.csv',[LHS_ParamRSS(:,1)';Susc_Age2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Susc3.csv',[LHS_ParamRSS(:,1)';Susc_Age3], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Susc4.csv',[LHS_ParamRSS(:,1)';Susc_Age4], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_SuscNPWID.csv',[LHS_ParamRSS(:,1)';SuscNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscNPWID1.csv',[LHS_ParamRSS(:,1)';Susc_Age1NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscNPWID2.csv',[LHS_ParamRSS(:,1)';Susc_Age2NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscNPWID3.csv',[LHS_ParamRSS(:,1)';Susc_Age3NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscNPWID4.csv',[LHS_ParamRSS(:,1)';Susc_Age4NPWID], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_SuscPWID.csv',[LHS_ParamRSS(:,1)';SuscPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscPWID1.csv',[LHS_ParamRSS(:,1)';Susc_Age1PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscPWID2.csv',[LHS_ParamRSS(:,1)';Susc_Age2PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscPWID3.csv',[LHS_ParamRSS(:,1)';Susc_Age3PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_SuscPWID4.csv',[LHS_ParamRSS(:,1)';Susc_Age4PWID], 'delimiter', ',', 'precision', 10);
   
    dlmwrite('LHS_Immune.csv',[LHS_ParamRSS(:,1)';Immune], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Immune1.csv',[LHS_ParamRSS(:,1)';Immune_Age1], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Immune2.csv',[LHS_ParamRSS(:,1)';Immune_Age2], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Immune3.csv',[LHS_ParamRSS(:,1)';Immune_Age3], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_Immune4.csv',[LHS_ParamRSS(:,1)';Immune_Age4], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_ImmuneNPWID.csv',[LHS_ParamRSS(:,1)';ImmuneNPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ImmuneNPWID1.csv',[LHS_ParamRSS(:,1)';Immune_Age1NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ImmuneNPWID2.csv',[LHS_ParamRSS(:,1)';Immune_Age2NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ImmuneNPWID3.csv',[LHS_ParamRSS(:,1)';Immune_Age3NPWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ImmuneNPWID4.csv',[LHS_ParamRSS(:,1)';Immune_Age4NPWID], 'delimiter', ',', 'precision', 10);

    dlmwrite('LHS_ImmunePWID.csv',[LHS_ParamRSS(:,1)';ImmunePWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ImmunePWID1.csv',[LHS_ParamRSS(:,1)';Immune_Age1PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ImmunePWID2.csv',[LHS_ParamRSS(:,1)';Immune_Age2PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ImmunePWID3.csv',[LHS_ParamRSS(:,1)';Immune_Age3PWID], 'delimiter', ',', 'precision', 10);
    dlmwrite('LHS_ImmunePWID4.csv',[LHS_ParamRSS(:,1)';Immune_Age4PWID], 'delimiter', ',', 'precision', 10);

year = [2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016];

MDHHSData = [];
MDHHSData(:,1) = year'; 
MDHHSData(:,2:4) = acute1529idu; 
MDHHSData(:,5:7) = chronic1529idu_model;
dlmwrite('MDHHSData.csv',[MDHHSData], 'delimiter', ',', 'precision', 10);

%% Setup Subsets of Interventions

load('HCVLHS05062018.mat')

% Set intervention levels

% Injection Initiation [index int1];
intervention_theta = [1 0.9 0.8 0.7 0.6];
 
% PWID contacts [index int2];
intervention_sigma = [1 0.9 0.8 0.7 0.6];
 
% Cessation [index int3];
intervention_gamma = [1 1.1 1.2 1.3 1.4]; 
 
% Relapse rate (k) [index int4];
intervention_kappa = [1 0.9 0.8 0.7 0.6];
 
% Treatment Current PWID [index int5];
intervention_phip = [0 0.1 0.2 0.3 0.4];
 
% Treatment Former PWID [index int6];
intervention_phin = [0 0.1 0.2 0.3 0.4];
 
% Death Current PWID [index int7];
intervention_etap = [1 0.9 0.8 0.7 0.6];

% Death Non-PWID [index int8];
intervention_etaz = [1 0.9 0.8 0.7 0.6];

% Treatment Duration [index int9];
intervention_omega = [365.25/84 365.25/56 365.25/112];
                    %12 weeks, 8 weeks, 16 weeks

% Proportion of Treated who Transmit [index k];
intervention_tau = [0.5 1 0.75 0.25 0];

% SVR [index k];
intervention_alpha = [0.9 1 0.8 0.7 0.6];

% Setup intervention combinations we want results from using the
    % intervention values from above and indices below. 
    % These are also listed in file: InterventionCombos_04162018.xlsx, 
    % which assigns an index e for single and combination interventions;
    
               %int1 int2 int3 int4 int5 int6 int7 int8 int9 int10 int11 
interventionrun=[1	1	1	1	1	1	1	1	1	1	1		;
                2	1	1	1	1	1	1	1	1	1	1		;
                3	1	1	1	1	1	1	1	1	1	1		;
                4	1	1	1	1	1	1	1	1	1	1		;
                5	1	1	1	1	1	1	1	1	1	1		;
                1	2	1	1	1	1	1	1	1	1	1		;
                1	3	1	1	1	1	1	1	1	1	1		;
                1	4	1	1	1	1	1	1	1	1	1		;
                1	5	1	1	1	1	1	1	1	1	1		;
                1	1	2	1	1	1	1	1	1	1	1		;
                1	1	3	1	1	1	1	1	1	1	1		;
                1	1	4	1	1	1	1	1	1	1	1		;
                1	1	5	1	1	1	1	1	1	1	1		;
                1	1	1	2	1	1	1	1	1	1	1		;
                1	1	1	3	1	1	1	1	1	1	1		;
                1	1	1	4	1	1	1	1	1	1	1		;
                1	1	1	5	1	1	1	1	1	1	1		;
                1	1	1	1	2	1	1	1	1	1	1		;
                1	1	1	1	3	1	1	1	1	1	1		;
                1	1	1	1	4	1	1	1	1	1	1		;
                1	1	1	1	5	1	1	1	1	1	1		;
                1	1	1	1	1	2	1	1	1	1	1		;
                1	1	1	1	1	3	1	1	1	1	1		;
                1	1	1	1	1	4	1	1	1	1	1		;
                1	1	1	1	1	5	1	1	1	1	1		;
                1	1	1	1	1	1	2	1	1	1	1		;
                1	1	1	1	1	1	3	1	1	1	1		;
                1	1	1	1	1	1	4	1	1	1	1		;
                1	1	1	1	1	1	5	1	1	1	1		;
                1	1	1	1	1	1	1	2	1	1	1		;
                1	1	1	1	1	1	1	3	1	1	1		;
                1	1	1	1	1	1	1	4	1	1	1		;
                1	1	1	1	1	1	1	5	1	1	1		;
                1	1	1	1	1	1	5	5	1	1	1		;
                2	1	1	1	1	1	5	5	1	1	1		;
                3	1	1	1	1	1	5	5	1	1	1		;
                4	1	1	1	1	1	5	5	1	1	1		;
                5	1	1	1	1	1	5	5	1	1	1		;
                1	2	1	1	1	1	5	5	1	1	1		;
                1	3	1	1	1	1	5	5	1	1	1		;
                1	4	1	1	1	1	5	5	1	1	1		;
                1	5	1	1	1	1	5	5	1	1	1		;
                1	1	2	1	1	1	5	5	1	1	1		;
                1	1	3	1	1	1	5	5	1	1	1		;
                1	1	4	1	1	1	5	5	1	1	1		;
                1	1	5	1	1	1	5	5	1	1	1		;
                1	1	1	2	1	1	5	5	1	1	1		;
                1	1	1	3	1	1	5	5	1	1	1		;
                1	1	1	4	1	1	5	5	1	1	1		;
                1	1	1	5	1	1	5	5	1	1	1		;
                1	1	1	1	2	1	5	5	1	1	1		;
                1	1	1	1	3	1	5	5	1	1	1		;
                1	1	1	1	4	1	5	5	1	1	1		;
                1	1	1	1	5	1	5	5	1	1	1		;
                1	1	1	1	1	2	5	5	1	1	1		;
                1	1	1	1	1	3	5	5	1	1	1		;
                1	1	1	1	1	4	5	5	1	1	1		;
                1	1	1	1	1	5	5	5	1	1	1		;
                1	1	1	1	1	1	2	2	1	1	1		;
                1	1	1	1	1	1	3	3	1	1	1		;
                1	1	1	1	1	1	4	4	1	1	1		;
                1	1	1	1	2	2	1	1	1	1	1		;
                1	1	1	1	3	3	1	1	1	1	1		;
                1	1	1	1	4	4	1	1	1	1	1		;
                1	1	1	1	5	5	1	1	1	1	1		;
                1	1	1	1	2	1	1	1	2	1	1		;
                1	1	1	1	3	1	1	1	2	1	1		;
                1	1	1	1	4	1	1	1	2	1	1		;
                1	1	1	1	5	1	1	1	2	1	1		;
                1	1	1	1	1	2	1	1	2	1	1		;
                1	1	1	1	1	3	1	1	2	1	1		;
                1	1	1	1	1	4	1	1	2	1	1		;
                1	1	1	1	1	5	1	1	2	1	1		;
                1	1	1	1	2	2	1	1	2	1	1		;
                1	1	1	1	3	3	1	1	2	1	1		;
                1	1	1	1	4	4	1	1	2	1	1		;
                1	1	1	1	5	5	1	1	2	1	1		;
                1	1	1	1	2	1	1	1	3	1	1		;
                1	1	1	1	3	1	1	1	3	1	1		;
                1	1	1	1	4	1	1	1	3	1	1		;
                1	1	1	1	5	1	1	1	3	1	1		;
                1	1	1	1	1	2	1	1	3	1	1		;
                1	1	1	1	1	3	1	1	3	1	1		;
                1	1	1	1	1	4	1	1	3	1	1		;
                1	1	1	1	1	5	1	1	3	1	1		;
                1	1	1	1	2	2	1	1	3	1	1		;
                1	1	1	1	3	3	1	1	3	1	1		;
                1	1	1	1	4	4	1	1	3	1	1		;
                1	1	1	1	5	5	1	1	3	1	1		;
                1	1	1	1	2	1	1	1	1	2	1		;
                1	1	1	1	3	1	1	1	1	2	1		;
                1	1	1	1	4	1	1	1	1	2	1		;
                1	1	1	1	5	1	1	1	1	2	1		;
                1	1	1	1	2	2	1	1	1	2	1		;
                1	1	1	1	3	3	1	1	1	2	1		;
                1	1	1	1	4	4	1	1	1	2	1		;
                1	1	1	1	5	5	1	1	1	2	1		;
                1	1	1	1	2	1	1	1	1	3	1		;
                1	1	1	1	3	1	1	1	1	3	1		;
                1	1	1	1	4	1	1	1	1	3	1		;
                1	1	1	1	5	1	1	1	1	3	1		;
                1	1	1	1	2	2	1	1	1	3	1		;
                1	1	1	1	3	3	1	1	1	3	1		;
                1	1	1	1	4	4	1	1	1	3	1		;
                1	1	1	1	5	5	1	1	1	3	1		;
                1	1	1	1	2	1	1	1	1	4	1		;
                1	1	1	1	3	1	1	1	1	4	1		;
                1	1	1	1	4	1	1	1	1	4	1		;
                1	1	1	1	5	1	1	1	1	4	1		;
                1	1	1	1	2	2	1	1	1	4	1		;
                1	1	1	1	3	3	1	1	1	4	1		;
                1	1	1	1	4	4	1	1	1	4	1		;
                1	1	1	1	5	5	1	1	1	4	1		;
                1	1	1	1	2	1	1	1	1	5	1		;
                1	1	1	1	3	1	1	1	1	5	1		;
                1	1	1	1	4	1	1	1	1	5	1		;
                1	1	1	1	5	1	1	1	1	5	1		;
                1	1	1	1	2	2	1	1	1	5	1		;
                1	1	1	1	3	3	1	1	1	5	1		;
                1	1	1	1	4	4	1	1	1	5	1		;
                1	1	1	1	5	5	1	1	1	5	1		;
                1	1	1	1	2	1	1	1	1	1	2		;
                1	1	1	1	3	1	1	1	1	1	2		;
                1	1	1	1	4	1	1	1	1	1	2		;
                1	1	1	1	5	1	1	1	1	1	2		;
                1	1	1	1	1	2	1	1	1	1	2		;
                1	1	1	1	1	3	1	1	1	1	2		;
                1	1	1	1	1	4	1	1	1	1	2		;
                1	1	1	1	1	5	1	1	1	1	2		;
                1	1	1	1	2	2	1	1	1	1	2		;
                1	1	1	1	3	3	1	1	1	1	2		;
                1	1	1	1	4	4	1	1	1	1	2		;
                1	1	1	1	5	5	1	1	1	1	2		;
                1	1	1	1	2	1	1	1	1	1	3		;
                1	1	1	1	3	1	1	1	1	1	3		;
                1	1	1	1	4	1	1	1	1	1	3		;
                1	1	1	1	5	1	1	1	1	1	3		;
                1	1	1	1	1	2	1	1	1	1	3		;
                1	1	1	1	1	3	1	1	1	1	3		;
                1	1	1	1	1	4	1	1	1	1	3		;
                1	1	1	1	1	5	1	1	1	1	3		;
                1	1	1	1	2	2	1	1	1	1	3		;
                1	1	1	1	3	3	1	1	1	1	3		;
                1	1	1	1	4	4	1	1	1	1	3		;
                1	1	1	1	5	5	1	1	1	1	3		;
                1	1	1	1	2	1	1	1	1	1	4		;
                1	1	1	1	3	1	1	1	1	1	4		;
                1	1	1	1	4	1	1	1	1	1	4		;
                1	1	1	1	5	1	1	1	1	1	4		;
                1	1	1	1	1	2	1	1	1	1	4		;
                1	1	1	1	1	3	1	1	1	1	4		;
                1	1	1	1	1	4	1	1	1	1	4		;
                1	1	1	1	1	5	1	1	1	1	4		;
                1	1	1	1	2	2	1	1	1	1	4		;
                1	1	1	1	3	3	1	1	1	1	4		;
                1	1	1	1	4	4	1	1	1	1	4		;
                1	1	1	1	5	5	1	1	1	1	4		;
                1	1	1	1	2	1	1	1	1	1	5		;
                1	1	1	1	3	1	1	1	1	1	5		;
                1	1	1	1	4	1	1	1	1	1	5		;
                1	1	1	1	5	1	1	1	1	1	5		;
                1	1	1	1	1	2	1	1	1	1	5		;
                1	1	1	1	1	3	1	1	1	1	5		;
                1	1	1	1	1	4	1	1	1	1	5		;
                1	1	1	1	1	5	1	1	1	1	5		;
                1	1	1	1	2	2	1	1	1	1	5		;
                1	1	1	1	3	3	1	1	1	1	5		;
                1	1	1	1	4	4	1	1	1	1	5		;
                1	1	1	1	5	5	1	1	1	1	5		;
                1	1	1	1	1	1	1	1	1	1	1	;
                2	1	1	1	1	1	1	1	1	1	1	;
                3	1	1	1	1	1	1	1	1	1	1	;
                4	1	1	1	1	1	1	1	1	1	1	;
                5	1	1	1	1	1	1	1	1	1	1	;
                2	2	1	1	1	1	1	1	1	1	1	;
                3	3	1	1	1	1	1	1	1	1	1	;
                4	4	1	1	1	1	1	1	1	1	1	;
                5	5	1	1	1	1	1	1	1	1	1	;
                2	2	2	1	1	1	1	1	1	1	1	;
                3	3	3	1	1	1	1	1	1	1	1	;
                4	4	4	1	1	1	1	1	1	1	1	;
                5	5	5	1	1	1	1	1	1	1	1	;
                2	2	2	2	1	1	1	1	1	1	1	;
                3	3	3	3	1	1	1	1	1	1	1	;
                4	4	4	4	1	1	1	1	1	1	1	;
                5	5	5	5	1	1	1	1	1	1	1	;
                2	2	2	2	2	1	1	1	1	1	1	;
                3	3	3	3	3	1	1	1	1	1	1	;
                4	4	4	4	4	1	1	1	1	1	1	;
                5	5	5	5	5	1	1	1	1	1	1	;
                2	2	2	2	2	2	1	1	1	1	1	;
                3	3	3	3	3	3	1	1	1	1	1	;
                4	4	4	4	4	4	1	1	1	1	1	;
                5	5	5	5	5	5	1	1	1	1	1	;
                1	1	1	1	1	1	5	5	1	1	1	;
                2	1	1	1	1	1	5	5	1	1	1	;
                3	1	1	1	1	1	5	5	1	1	1	;
                4	1	1	1	1	1	5	5	1	1	1	;
                5	1	1	1	1	1	5	5	1	1	1	;
                2	2	1	1	1	1	5	5	1	1	1	;
                3	3	1	1	1	1	5	5	1	1	1	;
                4	4	1	1	1	1	5	5	1	1	1	;
                5	5	1	1	1	1	5	5	1	1	1	;
                2	2	2	1	1	1	5	5	1	1	1	;
                3	3	3	1	1	1	5	5	1	1	1	;
                4	4	4	1	1	1	5	5	1	1	1	;
                5	5	5	1	1	1	5	5	1	1	1	;
                2	2	2	2	1	1	5	5	1	1	1	;
                3	3	3	3	1	1	5	5	1	1	1	;
                4	4	4	4	1	1	5	5	1	1	1	;
                5	5	5	5	1	1	5	5	1	1	1	;
                2	2	2	2	2	1	5	5	1	1	1	;
                3	3	3	3	3	1	5	5	1	1	1	;
                4	4	4	4	4	1	5	5	1	1	1	;
                5	5	5	5	5	1	5	5	1	1	1	;
                2	2	2	2	2	2	5	5	1	1	1	;
                3	3	3	3	3	3	5	5	1	1	1	;
                4	4	4	4	4	4	5	5	1	1	1	;
                5	5	5	5	5	5	5	5	1	1	1	;
                1	1	1	1	1	1	1	1	1	1	1	;
                1	1	1	1	1	2	1	1	1	1	1	;
                1	1	1	1	1	3	1	1	1	1	1	;
                1	1	1	1	1	4	1	1	1	1	1	;
                1	1	1	1	1	5	1	1	1	1	1	;
                1	1	1	1	2	2	1	1	1	1	1	;
                1	1	1	1	3	3	1	1	1	1	1	;
                1	1	1	1	4	4	1	1	1	1	1	;
                1	1	1	1	5	5	1	1	1	1	1	;
                1	1	1	2	2	2	1	1	1	1	1	;
                1	1	1	3	3	3	1	1	1	1	1	;
                1	1	1	4	4	4	1	1	1	1	1	;
                1	1	1	5	5	5	1	1	1	1	1	;
                1	1	2	2	2	2	1	1	1	1	1	;
                1	1	3	3	3	3	1	1	1	1	1	;
                1	1	4	4	4	4	1	1	1	1	1	;
                1	1	5	5	5	5	1	1	1	1	1	;
                1	2	2	2	2	2	1	1	1	1	1	;
                1	3	3	3	3	3	1	1	1	1	1	;
                1	4	4	4	4	4	1	1	1	1	1	;
                1	5	5	5	5	5	1	1	1	1	1	;
                2	2	2	2	2	2	1	1	1	1	1	;
                3	3	3	3	3	3	1	1	1	1	1	;
                4	4	4	4	4	4	1	1	1	1	1	;
                5	5	5	5	5	5	1	1	1	1	1	;
                1	1	1	1	1	1	5	5	1	1	1	;
                1	1	1	1	1	2	5	5	1	1	1	;
                1	1	1	1	1	3	5	5	1	1	1	;
                1	1	1	1	1	4	5	5	1	1	1	;
                1	1	1	1	1	5	5	5	1	1	1	;
                1	1	1	1	2	2	5	5	1	1	1	;
                1	1	1	1	3	3	5	5	1	1	1	;
                1	1	1	1	4	4	5	5	1	1	1	;
                1	1	1	1	5	5	5	5	1	1	1	;
                1	1	1	2	2	2	5	5	1	1	1	;
                1	1	1	3	3	3	5	5	1	1	1	;
                1	1	1	4	4	4	5	5	1	1	1	;
                1	1	1	5	5	5	5	5	1	1	1	;
                1	1	2	2	2	2	5	5	1	1	1	;
                1	1	3	3	3	3	5	5	1	1	1	;
                1	1	4	4	4	4	5	5	1	1	1	;
                1	1	5	5	5	5	5	5	1	1	1	;
                1	2	2	2	2	2	5	5	1	1	1	;
                1	3	3	3	3	3	5	5	1	1	1	;
                1	4	4	4	4	4	5	5	1	1	1	;
                1	5	5	5	5	5	5	5	1	1	1	;
                2	2	2	2	2	2	5	5	1	1	1	;
                3	3	3	3	3	3	5	5	1	1	1	;
                4	4	4	4	4	4	5	5	1	1	1	;
                5	5	5	5	5	5	5	5	1	1	1	];


%% Clear Unnecessary parts of workspace and Save

clearvars -except time1 nu ParamEstRuns CMatrixRuns PopSizeRuns acute1529idu chronic1529idu_model interventionrun intervention_alpha intervention_etap intervention_etaz intervention_gamma intervention_kappa intervention_omega intervention_phin intervention_phip intervention_sigma intervention_tau intervention_theta

save('HCV_IntSetup_05026018.mat')

%% Simulate Interventions

%Note: Recommend to simulate code below using high performance computing
%resources.

% Interventions are simulated in groups of 1000 parameter sets at a time
for indexnum=1:10; %# groups of paramsets
    
    indexmin=((indexnum-1)*1000)+1; 
    indexmax=indexnum*1000;

% Load intervention workspace    
load('HCV_IntSetup_05062018.mat') 

% Setup matrices to save results
chronicprev_results=[];
newchronic_results=[];
acute_results=[];
acute_results_report=[];
parameters_results=[];

% Loop over each intervention
for e=1:269;

intlev=interventionrun(e,1:11); 

% Assign the appropriate level to each intervention parameter  
ParamEstRuns1=ParamEstRuns;

int1=intlev(1,1); %theta (inj init);
    ParamEstRuns1(20,:)=ParamEstRuns(20,:)*intervention_theta(int1);
    ParamEstRuns1(21,:)=ParamEstRuns(21,:)*intervention_theta(int1);
    ParamEstRuns1(22,:)=ParamEstRuns(22,:)*intervention_theta(int1);
    ParamEstRuns1(23,:)=ParamEstRuns(23,:)*intervention_theta(int1);
int2=intlev(1,2); %sigma (contacts);
    ParamEstRuns1(1,:)=ParamEstRuns(1,:)*intervention_sigma(int2);
    ParamEstRuns1(2,:)=ParamEstRuns(2,:)*intervention_sigma(int2);
    ParamEstRuns1(3,:)=ParamEstRuns(3,:)*intervention_sigma(int2);
    ParamEstRuns1(4,:)=ParamEstRuns(4,:)*intervention_sigma(int2);
int3=intlev(1,3); %gamma (cessation);
    ParamEstRuns1(7,:)= ParamEstRuns(7,:)*intervention_gamma(int3);
    ParamEstRuns1(8,:)= ParamEstRuns(8,:)*intervention_gamma(int3);
    ParamEstRuns1(9,:)= ParamEstRuns(9,:)*intervention_gamma(int3);
    ParamEstRuns1(10,:)= ParamEstRuns(10,:)*intervention_gamma(int3);
int4=intlev(1,4); %kappa (relapse);
    ParamEstRuns1(24,:)=ParamEstRuns(24,:)*intervention_kappa(int4);
    ParamEstRuns1(25,:)=ParamEstRuns(25,:)*intervention_kappa(int4);
    ParamEstRuns1(26,:)=ParamEstRuns(26,:)*intervention_kappa(int4);
    ParamEstRuns1(27,:)=ParamEstRuns(27,:)*intervention_kappa(int4);
int5=intlev(1,5); %phip (current PWID trt);
    ParamEstRuns1(43,:)=intervention_phip(int5);
int6=intlev(1,6); %phin (former PWID trt);
    ParamEstRuns1(44,:)=intervention_phin(int6);
int7=intlev(1,7); %etap (current PWID mortality);
    ParamEstRuns1(51,:)=ParamEstRuns(51,:)*intervention_etap(int7);
int8=intlev(1,8); %etaz (non-PWID mortality);
    ParamEstRuns1(19,:)=ParamEstRuns(19,:)*intervention_etaz(int8);
int9=intlev(1,9); %omega (treatment duration);
    ParamEstRuns1(50,:)=intervention_omega(int9);
int10=intlev(1,10); %tau (proportion treated who transmit);
    ParamEstRuns1(42,:)=intervention_tau(int10);
int11=intlev(1,11); %alpha (treatment SVR);
    ParamEstRuns1(5,:)=intervention_alpha(int11);

% Select which block of 1250 parameter sets will be simulated
for i=indexmin:indexmax;

% Select a parameter set i 
param=ParamEstRuns1(:,i);

%contact matrix
c = reshape(CMatrixRuns(i,1:16),4,4);

    %Calculate proportion of contacts in each age group by column
    pi=zeros(4,4);
    for y=1:4 %rows
        for l=1:4 %columns
        pi(y,l)=c(y,l)./sum(c(1:4,l));
        end
    end

    pi=reshape(pi,16,1);
    
    popsize=PopSizeRuns(i,1:4)';
    
% Define parameters required for initial conditions;
lambda1=param(28);
lambda2=param(29);
lambda3=param(30);
lambda4=param(31);
d=param(11);
xi=param(36);
r=param(41);
eps=param(12);
psi1=param(45);
psi2=param(46);
psi3=param(47);
psi4=param(48);
omic1=param(37);
omic2=param(38);
omic3=param(39);
omic4=param(40);
zeta1=param(13);
zeta2=param(14);
zeta3=param(15);
zeta4=param(16);

   
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
for y=1:4
    N0(6,y) = popsize(y,1)*zeta(1,y) - (N0(7,y)+N0(8,y)+N0(9,y)+N0(10,y));
end

% Define class SNi (uninfected former PWID) as popsize*formerinject less sum of former PWID classes;

omic=[omic1 omic2 omic3 omic4];
for y=1:4
    N0(1,y) = popsize(y,1)*omic(1,y) - (N0(2,y)+N0(3,y)+N0(4,y)+N0(5,y));
end

% Reshape initial conditions into a vector;
N0=reshape(N0,44,1);

%Run ODE
[t,x] = ode15s(@HCV_DiffEq_07262018,time1,N0,[],param,nu,pi);

% Collect total (all age) data from each run
Y_newchronic = eps*(1-d)*(x(81,11*1-4)+x(81,11*2-4)+x(81,11*3-4)+x(81,11*4-4)+x(81,11*1-9)+x(81,11*2-9)+x(81,11*3-9)+x(81,11*4-9)); 
Y_newchronic1 = eps*(1-d)*(x(81,11*1-4)+x(81,11*2-4)+x(81,11*3-4)+x(81,11*1-9)+x(81,11*2-9)+x(81,11*3-9)); 
 
Y_chronicprev = x(81,11*1-3)+x(81,11*2-3)+x(81,11*3-3)+x(81,11*4-3)+x(81,11*1-8)+x(81,11*2-8)+x(81,11*3-8)+x(81,11*4-8);
Y_chronicprev1 = x(81,11*1-3)+x(81,11*2-3)+x(81,11*3-3)+x(81,11*1-8)+x(81,11*2-8)+x(81,11*3-8);

Y_acute = x(81,11*1-4)+x(81,11*2-4)+x(81,11*3-4)+x(81,11*4-4)+x(81,11*1-9)+x(81,11*2-9)+x(81,11*3-9)+x(81,11*4-9);
Y_acute1 = x(81,11*1-4)+x(81,11*2-4)+x(81,11*3-4)+x(81,11*1-9)+x(81,11*2-9)+x(81,11*3-9);

Y_acute_report = r*(x(81,11*1-4)+x(81,11*2-4)+x(81,11*3-4)+x(81,11*4-4)+x(81,11*1-9)+x(81,11*2-9)+x(81,11*3-9)+x(81,11*4-9));
Y_acute_report1 = r*(x(81,11*1-4)+x(81,11*2-4)+x(81,11*3-4)+x(81,11*1-9)+x(81,11*2-9)+x(81,11*3-9));

% Add results to results matrices   
chronicprev_results = [chronicprev_results; e,int1,int2,int3,int4,int5,int6,int7,int8,int9,int10,int11,i,Y_chronicprev,Y_chronicprev1];
newchronic_results = [newchronic_results; e,int1,int2,int3,int4,int5,int6,int7,int8,int9,int10,int11,i,Y_newchronic,Y_newchronic1];
acute_results = [acute_results; e,int1,int2,int3,int4,int5,int6,int7,int8,int9,int10,int11,i,Y_acute,Y_acute1];
acute_results_report = [acute_results_report; e,int1,int2,int3,int4,int5,int6,int7,int8,int9,int10,int11,i,Y_acute_report,Y_acute_report1];
parameters_results = [parameters_results; e,int1,int2,int3,int4,int5,int6,int7,int8,int9,int10,int11,i,param'];

    end;
    
% Display the interation (e) number in command window;
chronicprev_results(end,1)

    end;

% Save the intervention workspace;
save(['HCVInterventions05062018' num2str(indexnum) '.mat'])

% Display the index number in command window;
indexnum

% Clear workspace and results to run the next set of 1000 parameter sets
clear

end;

%% Combine the 10 Intervention Workspaces

chronicprev_results1 = [];
newchronic_results1 = [];
acute_results1 = [];
acute_results_report1=[];
parameters_results1 = [];
  
for indexnum=1:10;
    
    load(['HCVInterventions05062018' num2str(indexnum) '.mat']);
                            
    chronicprev_results1 = [chronicprev_results1;chronicprev_results];
    newchronic_results1 = [newchronic_results1;newchronic_results];
    acute_results1 = [acute_results1;acute_results];
    acute_results_report1 = [acute_results_report1;acute_results_report];
    parameters_results1 = [parameters_results1;parameters_results];

end;

% Save as a new workspace
save('InterventionResults_05062018.mat','-v7.3')

% Save intervention simulation results as txt files
dlmwrite('ChronicPrev_SingleSeq.txt',[chronicprev_results1], 'precision', 10);
dlmwrite('NewChronic_SingleSeq.txt',[newchronic_results1], 'precision', 10);
dlmwrite('Acute_SingleSeq.txt',[acute_results1], 'precision', 10);
dlmwrite('AcuteReport_SingleSeq.txt',[acute_results_report1], 'precision', 10);
dlmwrite('Parameters_SingleSeq.txt',[parameters_results1], 'precision', 10);

%% View the intervention results workspace
load('InterventionResults_05062018.mat')