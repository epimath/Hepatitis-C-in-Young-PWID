%% Ordinary Differential Equations: Gicquelais et al. (2017). Hepatitis C tranmsission model.

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

% The compartments are:
% N[SNi; ANi; CNi; TNi; Si; Ai; Ci; Ti; Zi]                   

function hcv = HCV_DiffEq_04202017(t,N,param,nu,c)

hcv = zeros(27,1);

N=reshape(N,9,3);

% Contact Matrix
con = zeros(3,1);

% con is a 3x1 matrix with each cell as a product of the #
% contacts * sum(infecteds who contribute to transmission 
% [Ci and Ai]) = c*w;

for v=1:3;
    w(v,1)=sum(N(6,v)+N(7,v));
end;

con=c*w;

% Define parameters
b=param(1);
phi=param(2:4);
d=param(5);
k=param(6);
s=param(8);
gn=param(9);
gp=param(10);
a=param(11);
Z0=param(12);
mu=param(16:18);
etap=param(19)*param(21);
etan=param(20);
death_former=param(19)*etan;

% Differential equations for the 15-19 age group
v=1;

dSNdt(v) = s*N(4,v)+2*(1-a)*N(2,v)+d*N(5,v)-k*N(1,v)-nu(v,1)*N(1,v)-mu(v)*death_former*N(1,v);
dANdt(v) = d*N(6,v)-k*N(2,v)-2*N(2,v)-nu(v,1)*N(2,v)-mu(v)*death_former*N(2,v);
dCNdt(v) = 2*a*N(2,v)+d*N(7,v)-k*N(3,v)-gn*N(3,v)-nu(v,1)*N(3,v)-mu(v)*death_former*N(3,v);
dTNdt(v) = gn*N(3,v)+d*N(8,v)-k*N(4,v)-s*N(4,v)-nu(v,1)*N(4,v)-mu(v)*death_former*N(4,v);
dSdt(v) = phi(v)*N(9,v)+s*N(8,v)+2*(1-a)*N(6,v)+k*N(1,v)-d*N(5,v)-b*N(5,v)*con(v,1)-nu(v,1)*N(5,v)-mu(v)*etap*N(5,v);
dAdt(v) = b*N(5,v)*con(v,1)+k*N(2,v)-d*N(6,v)-2*N(6,v)-nu(v,1)*N(6,v)-mu(v)*etap*N(6,v);
dCdt(v) = 2*a*N(6,v)+phi(v)*N(3,v)-d*N(7,v)-gp*N(7,v)-nu(v,1)*N(7,v)-mu(v)*etap*N(7,v);
dTdt(v) = gp*N(7,v)+k*N(4,v)-d*N(8,v)-s*N(8,v)-nu(v,1)*N(8,v)-mu(v)*etap*N(8,v);
dZdt(v) = Z0-phi(v)*N(9,v)-nu(v,1)*N(9,v)-mu(v)*N(9,v);

% Differential equations for the 20-24 and 25-30 age groups
for v=2:3

dSNdt(v) = s*N(4,v)+2*(1-a)*N(2,v)+d*N(5,v)-k*N(1,v)-nu(v,1)*N(1,v)+nu(v-1,1)*N(1,v-1)-mu(v)*death_former*N(1,v);
dANdt(v) = d*N(6,v)-k*N(2,v)-2*N(2,v)-nu(v,1)*N(2,v)+nu(v-1,1)*N(2,v-1)-mu(v)*death_former*N(2,v);
dCNdt(v) = 2*a*N(2,v)+d*N(7,v)-k*N(3,v)-gn*N(3,v)-nu(v,1)*N(3,v)+nu(v-1,1)*N(3,v-1)-mu(v)*death_former*N(3,v);
dTNdt(v) = gn*N(3,v)+d*N(8,v)-k*N(4,v)-s*N(4,v)-nu(v,1)*N(4,v)+nu(v-1,1)*N(4,v-1)-mu(v)*death_former*N(4,v);
dSdt(v) = phi(v)*N(9,v)+s*N(8,v)+2*(1-a)*N(6,v)+k*N(1,v)-d*N(5,v)-b*N(5,v)*con(v,1)-nu(v,1)*N(5,v)+nu(v-1,1)*N(5,v-1)-mu(v)*etap*N(5,v);
dAdt(v) = b*N(5,v)*con(v,1)+k*N(2,v)-d*N(6,v)-2*N(6,v)-nu(v,1)*N(6,v)+nu(v-1,1)*N(6,v-1)-mu(v)*etap*N(6,v);
dCdt(v) = 2*a*N(6,v)+k*N(3,v)-d*N(7,v)-gp*N(7,v)-nu(v,1)*N(7,v)+nu(v-1,1)*N(7,v-1)-mu(v)*etap*N(7,v);
dTdt(v) = gp*N(7,v)+k*N(4,v)-d*N(8,v)-s*N(8,v)-nu(v,1)*N(8,v)+nu(v-1,1)*N(8,v-1)-mu(v)*etap*N(8,v);
dZdt(v) = nu(v-1,1)*N(9,v-1)-phi(v)*N(9,v)-nu(v,1)*N(9,v)-mu(v)*N(9,v);

end;

% Compile results
hcv = [dSNdt(1); dANdt(1); dCNdt(1); dTNdt(1); dSdt(1); dAdt(1); dCdt(1); dTdt(1); dZdt(1);
       dSNdt(2); dANdt(2); dCNdt(2); dTNdt(2); dSdt(2); dAdt(2); dCdt(2); dTdt(2); dZdt(2);
       dSNdt(3); dANdt(3); dCNdt(3); dTNdt(3); dSdt(3); dAdt(3); dCdt(3); dTdt(3); dZdt(3)];
   
