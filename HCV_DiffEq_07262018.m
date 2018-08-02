%% Ordinary Differential Equations: Gicquelais et al. (2018). Hepatitis C tranmsission model.

function hcv = HCV_DiffEq_07262018(t,N,param,nu,pi);

hcv = zeros(44,1);

N=reshape(N,11,4);

pi=reshape(pi,4,4);

% For reference, the parameters are:

% param   =[sigma1 sigma2 sigma3 sigma4 
%              1      2     3       4          

%          a b  g(1:4) d  eps zeta(1:4) etan etap etaz theta(1:4) k(1:4) lambda(1:4)
%          5 6  7:10   11 12  13:16     17   18   19   20:23      24:27   28:31
%
%          mu(1:4) xi omic(1:4) r  tau phip phin psi(1:4) psin w  intervention_etap Z0];
%          32:35   36 37:40     41 42  43   44   45:48    49   50          51       52  
%         
% The inverse of age group size is in vector nu:
    % nu = [nu1,nu2,nu3,nu4];

% The compartments are:
% N[SNi; ANi; CNi; INi; TNi; Si; Ai; Ci; INi; Ti; Zi]

% Define parameters
sigma=[param(1) param(2) param(3) param(4)];

a=param(5);
b=param(6);
g=param(7:10);
d=param(11);
eps=param(12);
zeta=param(13:16);
etan=param(17);
etap=param(18); 
etaz=param(19);
theta=param(20:23);
k=param(24:27);
lambda=param(28:31);
mu=param(32:35);
xi=param(36);
omic=param(37:40);
r=param(41);
tau=param(42);
phip=param(43);
phin=param(44);
psi=(45:48);
psin=param(49);
w=param(50);
intervention_etap=param(51);
Z0=param(52);

death_former=etan*etap;

etap=etap*intervention_etap; %etap*any intervention that scales it

% Contacts with infectious people by age group:
% Number of infectious people contacted in each age group is the 
% product of infectious people in each age compartment and pi, the 
% proportion of contacts in each age group

% 1x4 matrix of infecteds who contribute to transmission
inf=zeros(1,4);
for j=1:4;
    inf(1,j) = N(7,j)+N(8,j)+tau*N(10,j);
end;

% 1x4 matrix of infecteds who contribute to transmission and make contact
% with a susceptible of each age group
con=zeros(1,4);
con=inf*pi;

% Differential equations for the 15-19 age group
v=1;

dSNdt(v) = w*a*N(5,v)+d*(1-xi)*eps*N(2,v)+g(v)*N(6,v)-k(v)*N(1,v)-mu(v)*death_former*N(1,v)-nu(v)*N(1,v);
dANdt(v) = -eps*N(2,v)+g(v)*N(7,v)-k(v)*N(2,v)-mu(v)*death_former*N(2,v)-nu(v)*N(2,v);
dCNdt(v) = (1-d)*eps*N(2,v)+w*(1-a)*N(5,v)-phin*N(3,v)+g(v)*N(8,v)-k(v)*N(3,v)-mu(v)*death_former*N(3,v)-nu(v)*N(3,v);
dINdt(v) = d*xi*eps*N(2,v)+g(v)*N(9,v)-k(v)*N(4,v)-mu(v)*death_former*N(4,v)-nu(v)*N(4,v);
dTNdt(v) = phin*N(3,v)-w*N(5,v)+g(v)*N(10,v)-k(v)*N(5,v)-mu(v)*death_former*N(5,v)-nu(v)*N(5,v);

dSdt(v) = theta(v)*N(11,v)-b*N(6,v)*sigma(v)*con(1,v)+w*a*N(10,v)+d*(1-xi)*eps*N(7,v)+k(v)*N(1,v)-g(v)*N(6,v)-mu(v)*etap*N(6,v)-nu(v)*N(6,v);
dAdt(v) = b*N(6,v)*sigma(v)*con(1,v)-eps*N(7,v)+k(v)*N(2,v)-g(v)*N(7,v)-mu(v)*etap*N(7,v)-nu(v)*N(7,v);
dCdt(v) = (1-d)*eps*N(7,v)+w*(1-a)*N(10,v)-phip*N(8,v)+k(v)*N(3,v)-g(v)*N(8,v)-mu(v)*etap*N(8,v)-nu(v)*N(8,v);
dIdt(v) = d*xi*eps*N(7,v)+k(v)*N(4,v)-g(v)*N(9,v)-mu(v)*etap*N(9,v)-nu(v)*N(9,v);
dTdt(v) = phip*N(8,v)-w*N(10,v)+k(v)*N(5,v)-g(v)*N(10,v)-mu(v)*etap*N(10,v)-nu(v)*N(10,v);

dZdt(v) = psin*Z0-theta(v)*N(11,v)-mu(v)*etaz*N(11,v)-nu(v)*N(11,v);

% Differential equations for the 20-25, 26-29, 30-64 age groups
for v=2:4;
    
dSNdt(v) = w*a*N(5,v)+d*(1-xi)*eps*N(2,v)+g(v)*N(6,v)-k(v)*N(1,v)-mu(v)*death_former*N(1,v)-nu(v)*N(1,v)+nu(v-1)*N(1,v-1);
dANdt(v) = -eps*N(2,v)+g(v)*N(7,v)-k(v)*N(2,v)-mu(v)*death_former*N(2,v)-nu(v)*N(2,v)+nu(v-1)*N(2,v-1);
dCNdt(v) = (1-d)*eps*N(2,v)+w*(1-a)*N(5,v)-phin*N(3,v)+g(v)*N(8,v)-k(v)*N(3,v)-mu(v)*death_former*N(3,v)-nu(v)*N(3,v)+nu(v-1)*N(3,v-1);
dINdt(v) = d*xi*eps*N(2,v)+g(v)*N(9,v)-k(v)*N(4,v)-mu(v)*death_former*N(4,v)-nu(v)*N(4,v)+nu(v-1)*N(4,v-1);
dTNdt(v) = phin*N(3,v)-w*N(5,v)+g(v)*N(10,v)-k(v)*N(5,v)-mu(v)*death_former*N(5,v)-nu(v)*N(5,v)+nu(v-1)*N(5,v-1);

dSdt(v) = theta(v)*N(11,v)-b*N(6,v)*sigma(v)*con(1,v)+w*a*N(10,v)+d*(1-xi)*eps*N(7,v)+k(v)*N(1,v)-g(v)*N(6,v)-mu(v)*etap*N(6,v)-nu(v)*N(6,v)+nu(v-1)*N(6,v-1);
dAdt(v) = b*N(6,v)*sigma(v)*con(1,v)-eps*N(7,v)+k(v)*N(2,v)-g(v)*N(7,v)-mu(v)*etap*N(7,v)-nu(v)*N(7,v)+nu(v-1)*N(7,v-1);
dCdt(v) = (1-d)*eps*N(7,v)+w*(1-a)*N(10,v)-phip*N(8,v)+k(v)*N(3,v)-g(v)*N(8,v)-mu(v)*etap*N(8,v)-nu(v)*N(8,v)+nu(v-1)*N(8,v-1);
dIdt(v) = d*xi*eps*N(7,v)+k(v)*N(4,v)-g(v)*N(9,v)-mu(v)*etap*N(9,v)-nu(v)*N(9,v)+nu(v-1)*N(9,v-1);
dTdt(v) = phip*N(8,v)-w*N(10,v)+k(v)*N(5,v)-g(v)*N(10,v)-mu(v)*etap*N(10,v)-nu(v)*N(10,v)+nu(v-1)*N(10,v-1);

dZdt(v) = nu(v-1)*N(11,v-1)-theta(v)*N(11,v)-mu(v)*etaz*N(11,v)-nu(v)*N(11,v);

end;

% Compile results
hcv = [dSNdt(1); dANdt(1); dCNdt(1); dINdt(1); dTNdt(1); dSdt(1); dAdt(1); dCdt(1); dIdt(1); dTdt(1); r;%psi(1,1);%dZdt(1); 
       dSNdt(2); dANdt(2); dCNdt(2); dINdt(2); dTNdt(2); dSdt(2); dAdt(2); dCdt(2); dIdt(2); dTdt(2); r;%psi(1,2);%dZdt(2); 
       dSNdt(3); dANdt(3); dCNdt(3); dINdt(3); dTNdt(3); dSdt(3); dAdt(3); dCdt(3); dIdt(3); dTdt(3); r;%psi(1,3);%dZdt(3); 
       dSNdt(4); dANdt(4); dCNdt(4); dINdt(4); dTNdt(4); dSdt(4); dAdt(4); dCdt(4); dIdt(4); dTdt(4); r];%psi(1,4)];%dZdt(4); 
   