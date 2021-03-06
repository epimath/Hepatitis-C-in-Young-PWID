%% Least Squares Function 

% Gicquelais RE, Foxman B, Coyle J, Eisenberg MC.(2019).
% Hepatitis C transmission in young people who inject drugs: insights using
% a dynamic model informed by state public health surveillance. Epidemics.
% https://doi.org/10.1016/j.epidem.2019.02.003

% Likehood function file used to fit the estimated parameters
% to data and calculate residual sum of squares. 

function RSS = HCV_LS_11152018(time,N0,param1,paramset,nu,pi,popsize,acute1529idu)

param = abs([paramset(1:5);param1(1); paramset(6:18); param1(2:4); paramset(19:end)]);

resid=zeros(3,1);

try

[t,x] = ode15s(@HCV_DiffEq_11152018,time,N0,[],param,nu,pi);

%Calculate the sum of squared residuals for each age group

%define params r, d, eps
r=param(41);
d=param(11);
eps=param(12);

for i=1:3

% Number of acute cases predicted by the model;
y(:,i)=r*(x(:,11*i-4)+x(:,11*i-9));

% Calculate residuals:

% Unweighted least squares;
    resid(i,1) = (acute1529idu(:,i) - y(:,i))'*(acute1529idu(:,i) - y(:,i));

% Other options for calculating the least squares weights (old code):

% Weighted least squares (weighted by size of data point);
   
    %resid(i,1) = ((acute1529idu(:,i) - y(:,i))./acute1529idu(:,i))'*((acute1529idu(:,i) - y(:,i))./acute1529idu(:,i));
        
% Weighted least squares (weighted by max in each age class);
    %resid(i,1) = ((acute1529idu(:,i) - y(:,i))./max(acute1529idu(:,i)))'*((acute1529idu(:,i) - y(:,i))./max(acute1529idu(:,i)));
    
 %Weighted least squares (weighted by sqrt(data));
    %resid(i,1) = ((acute1529idu(:,i) - y(:,i))./sqrt(acute1529idu(:,i)))'*((acute1529idu(:,i) - y(:,i))./sqrt(acute1529idu(:,i)));


end;


% Sum the squared residuals over all age groups;
RSS=sum(resid);

catch
    RSS=0; %Set RSS=0 if model params sampled results in non-convergence
end