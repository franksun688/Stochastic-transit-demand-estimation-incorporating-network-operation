%stochastic demand estimation
%Author: Ran Sun

%Date: DEC 2018

%Result demand_RC- recontruction demand, demand_mean - mean demand
%demand_table, normally simulated demand


%Traffic Assignment using SUE and MNL
%Generate linkflow, pathflow, link proportion, path proportion
%Given cost, demand
%Use SUE and MNL, and compare result


%MAIN----------------------------------------------------------------------
%--------------------------------------------------------------------------

%set random seed
rng('default') % For reproducibility

%days of observation
k = 5;

%k days of observation, demand: k days of demand (k*n)
num_link = 6;
num_od = 6;
num_path = 11;

% delta0 and gamma0 network topology
%delta is link-path
delta0= [0,0,0,1,0,0,0,0,0,0,0;
         1,0,1,0,1,0,1,0,0,0,0;
         1,0,0,0,1,0,0,1,0,1,0;
         1,1,0,0,0,0,0,1,0,0,1;
         0,1,0,0,0,1,0,0,0,0,0;
         0,0,1,0,0,0,0,0,1,0,0];

%gamma is od-path
gamma0=[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0; 
        0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0; 
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0; 
        0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0; 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0; 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
    
%capacity constraint
ka =[2000,2500,2000,2500,1000,1000]';

%in-veh time, time unit: minutes bacause theta in the objective function
ts = [30,24,10,12,10,20];

%frequency, uased to calculate waiting time
f = [5,5,10,5,5,5];

alpha = [0.5,0.5,0.5,0.5,0.5,0.5];
us = 60 * alpha./f;

%link cost = in-veh cost + waiting cost
linkcost = (ts + us)';

% state the dimention of linkflow/pathflow
linkflow_MNL_table = zeros(num_link,k);
linkflow_SUE_table = zeros(num_link,k);

pathflow_MNL_table = zeros(num_link,k);
pathflow_SUE_table = zeros(num_link,k);

proportion_MNL_table = zeros(num_link*k,num_od*k);
proportion_SUE_table = zeros(num_link*k,num_od*k);

link_ob_MNL = zeros(num_link*k,1);
link_ob_SUE = zeros(num_link*k,1);
%record all the randomly generatated demand
demand_table = zeros(num_od,k);

%obtain linkflow, pathflow, link proportion, path proportion of MNL and SUE
%for k day of obesertaion.
for i = 1:k
    %one day observation of damand
    demand = normrnd(60,5,[num_od,1]);
    
    demand_table(:,i) = demand;
    
    %rule out negative demand
    if all(demand<0)
        continue;
    end
    
    %obtain result from MNL and SUE function
    [linkflow_MNL, pathflow_MNL, proportion_path_MNL, proportion_link_MNL]= MNL (demand,num_path,num_od,delta0,gamma0,linkcost);
    [linkflow_SUE, pathflow_SUE, proportion_path_SUE, proportion_link_SUE] = SUE(num_link,num_od,num_path, demand,linkcost,ka,delta0,gamma0);
    
    linkflow_MNL_table(1:num_link,i) = linkflow_MNL;
    linkflow_SUE_table(1:num_link,i) = linkflow_SUE;
    
    pathflow_MNL_table(1:num_path,i) = pathflow_MNL;
    pathflow_SUE_table(1:num_path,i) = pathflow_SUE;
    
    proportion_MNL_table(1+(i-1)*num_link:i*num_link,1+(i-1)*num_od:i*num_od) = proportion_link_MNL;
    proportion_SUE_table(1+(i-1)*num_link:i*num_link,1+(i-1)*num_od:i*num_od) = proportion_link_SUE;
    
    link_ob_MNL(1+(i-1)*num_link:i*num_link,1) =linkflow_MNL ;
    link_ob_SUE(1+(i-1)*num_link:i*num_link,1) =linkflow_SUE ;
end

%result comparison
linkflow_compare2 = linkflow_MNL_table-linkflow_SUE_table;
pathflow_compare2 = pathflow_MNL_table-pathflow_SUE_table;


%demand estimation

[x_sol,demand_mean,demand_RC] = estimation_optim(proportion_MNL_table, link_ob_MNL, num_od,num_link,k);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%functions-----------------------------------------------------------------
%MNL-----------------------------------------------------------------------

% transit assignment using MNL
function [linkflow,pathflow,proportion_path,proportion_link] = MNL (demand,num_path,num_od,delta0,gamma0,linkcost)

%multimonial parameter
theta_MNL = 2;

%pathcost, R*1
pathcost = delta0'*linkcost;

%diag of pathcost, R*R
diagpathcost = diag(pathcost);

%od - pathcost, same dimension of od-path
OD_path_cost = gamma0*diagpathcost;

%od - pathcost, calculate using multinomial logit
pathprob = gamma0;

%the matrix of e^(-theta*Cr)
ecost = exp(1).^(-theta_MNL*OD_path_cost).*gamma0;

for i = 1:num_od
    s = sum(ecost(i,:));
    for j = 1:num_path
        pathprob(i,j) = ecost(i,j)/s;
    end
end
%repeat demand into R columns,n*R
demandrep = repmat(demand,1,num_path);

%od-pathflow n*R
odpathflow = demandrep.*pathprob;
%path-flow
pathflow = sum(odpathflow);

linkflow = delta0*pathflow';
proportion_path = pathprob';

proportion_link = delta0*proportion_path;
end

%SUE-----------------------------------------------------------------------

%transit assignment using SUE
%given one day observation of demand and cost parameters, return link flow
%and path flow, and corresponding proportion matrix
function [linkflow,pathflow,proportion_path, proportion_link] = SUE (num_link,num_od,num_path, demand,linkcost,ka,delta0,gamma0)
%cost parameter
theta_SUE = 2;

%inequality constraints, AX<=b
A1 = zeros(num_link+num_path,num_link+num_path);
for i= 1:num_link
    A1(i,i) = 1;
end
for i = (num_link+1) : (num_link+num_path)
    A1(i,i) = -1;
end

b1 = ka;
b2 = zeros(num_path,1);
A = A1;
b = [b1;b2];

%equality constraints, AeqX=Beq
Beq = [demand;zeros(num_link,1)];
Aeq = zeros(num_od+num_link,num_link+num_path);
Aeq(1:num_od,(num_link+1:num_link+num_path)) = gamma0;
Aeq((num_od+1:num_od+num_link),(num_link+1:num_link+num_path)) = -delta0;
Aeq((num_od+1:num_od+num_link),(1:num_link)) = eye(num_link);

%objective function
fun = @(vh) theta_SUE*(linkcost'*vh(1:num_link))+ vh(1+num_link:num_link+num_path)'*(log(vh(1+num_link:num_link+num_path))-1);

%initial point
x0 = 10*ones([num_link+num_path,1]);

vh = fmincon(fun,x0,A,b,Aeq,Beq);

% vsol1 = [vh(7),vh(8)+vh(9)+vh(10),vh(8),   vh(9)+vh(10),   vh(8)/2+vh(9),vh(8)/2+vh(10)]
% vsol2 = [vh(1),vh(2)+vh(4),vh(2),vh(5)+vh(6),vh(3)/2+vh(6),vh(3)/2]
% vsol3= [vh(1),vh(2)+vh(4),vh(2)+vh(5)/2,vh(5)/2+vh(6),vh(3)/2+vh(6),vh(3)/2+vh(6)]

linkflow = vh(1:num_link);
pathflow = vh(num_link+1:num_link+num_path);

%generate path proportion
proportion_path_0 = pathflow'.*gamma0;
proportion_path = (proportion_path_0./demand)';

%link proportion
proportion_link = delta0*proportion_path;

end


%Estimation Optimization---------------------------------------------------
%--------------------------------------------------------------------------


%proportaional matrix is G_p. link_ob is observed link counts, 
function[xp,demand_mean,demand_RC]= estimation_optim(proportion_MNL, link_ob, num_od,num_link,k)
%x is demand, (n*k) * 1

%xp is [theta;x]
x0 = zeros(num_od*(k+1),1);
c0 = zeros(num_od*k,num_od);
c1 = zeros(num_od*k,num_od*k);
for i = 1:k
    c0(1+(i-1)*num_od:i*num_od,:) = eye(num_od);
    c1(1+(i-1)*num_od:i*num_od , 1+(i-1)*num_od:i*num_od) = -eye(num_od);
end
c = [c0,c1];


fmin = @(xp) norm(c*xp);

%inequality constraints, AX<=b, here does not exist
A = [];
b = [];

%%equality constraints, AeqX=Beq

Aeq =[zeros(num_link*k,num_od), proportion_MNL];
Beq = link_ob;

xp = fmincon(fmin,x0,A,b,Aeq,Beq);
demand_mean = xp(1:num_od);
demand_RC = xp(num_od+1:num_od+k*num_od);
end
