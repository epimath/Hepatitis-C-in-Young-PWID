%% Least Squares Function: Gicquelais et al. (2017). Hepatitis C tranmsission model.

% This is the likehood function file used to fit the estimated parameters
% to data and calculate residual sum of squares. 

function RSS = HCV_LS_04202017(time,N0,param1,paramset,nu,c,chronic1530idu)

param=abs([param1;paramset]);

resid=zeros(3,1);

[t,x] = ode15s(@HCV_DiffEq_04202017,time,N0,[],param,nu,c);

%Calculate the sum of squared residuals for each age group
for i=1:3

% Number of new chronic cases predicted by the model;
y(:,i)=param(13)*(2*param(11)*(x(:,(2*i)+7*(i-1))+x(:,(6*i)+3*(i-1))));

% Weighted least squares (weighted by size of data point);
    resid(i,1) = ((chronic1530idu(:,i) - y(:,i))./chronic1530idu(:,i))'*((chronic1530idu(:,i) - y(:,i))./chronic1530idu(:,i));

% Other options for calculating the least squares weights:

% Unweighted least squares;
    %resid(i,1) = (chronic1530idu(:,i) - y(:,i))'*(chronic1530idu(:,i) - y(:,i));
    
% Weighted least squares (weighted by max in each age class);
    %resid(i,1) = ((chronic1530idu(:,i) - y(:,i))./max(chronic1530idu(:,i)))'*((chronic1530idu(:,i) - y(:,i))./max(chronic1530idu(:,i)));
    
 %Weighted least squares (weighted by sqrt(data));
    %resid(i,1) = ((chronic1530idu(:,i) - y(:,i))./sqrt(chronic1530idu(:,i)))'*((chronic1530idu(:,i) - y(:,i))./sqrt(chronic1530idu(:,i)));


end;

resid;

% Sum the squared residuals over all age groups;
RSS=sum(resid);
