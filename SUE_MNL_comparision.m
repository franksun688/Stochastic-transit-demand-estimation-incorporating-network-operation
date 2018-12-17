%Traffic Assignment using SUE and MNL
%Author: Ran Sun
%Date: 12/16/2018

%Generate linkflow, pathflow, link proportion, path proportion
%Given cost, demand
%Use SUE and MNL, and compare result

%set random seed
rng('default') % For reproducibility

%days of observation
k = 1;

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
linkflow_NML_table = zeros(num_link,k);
linkflow_SUE_table = zeros(num_link,k);

pathflow_NML_table = zeros(num_link,k);
pathflow_SUE_table = zeros(num_link,k);

%record all the randomly generatated demand
demand_table = zeros(num_od,k);

%obtain linkflow, pathflow, link proportion, path proportion of MNL and SUE
%for k day of obesertaion.
for i = 1:k
    %one day observation of damand
    demand = normrnd(60,15,[num_od,1]);
    
    demand_table(:,i) = demand;
    
    %rule out negative demand
    if all(demand<0)
        continue;
    end
    
    %obtain result from MNL and SUE function
    [linkflow_NML, pathflow_NML, proportion_path_MNL, proportion_link_MNL]= MNL (demand,num_path,num_od,delta0,gamma0,linkcost);
    [linkflow_SUE, pathflow_SUE, proportion_path_SUE, proportion_link_SUE] = SUE(num_link,num_od,num_path, demand,linkcost,ka,delta0,gamma0);
    
    linkflow_NML_table(1:num_link,i) = linkflow_NML;
    linkflow_SUE_table(1:num_link,i) = linkflow_SUE;
    
    pathflow_NML_table(1:num_path,i) = pathflow_NML;
    pathflow_SUE_table(1:num_path,i) = pathflow_SUE;
    
    
    
end

%result comparison
linkflow_compare2 = linkflow_NML_table-linkflow_SUE_table;
pathflow_compare2 = pathflow_NML_table-pathflow_SUE_table;


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