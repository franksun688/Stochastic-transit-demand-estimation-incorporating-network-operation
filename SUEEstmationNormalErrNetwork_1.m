%stochastic demand estimation
%Author: Ran Sun

%Date: Feb 2019

%network: Decea 1993
%Result demand_RC- recontruction demand, demand_mean - mean demand
%demand_table, simulated demand, assume Possion distribution
%

%Traffic Assignment using SUE
%Generate linkflow, pathflow, link proportion, path proportion
%Given cost, demand
%Use SUE, and compare result



%set random seed
rng('default') % For reproducibility

%Input---------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%days of observation
k = 10;

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
ka =[250,200,150,150,100,100]';

%in-veh time, time unit: minutes bacause theta in the objective function
ts = [30,24,10,12,10,20];

%frequency, uased to calculate waiting time
f = [5,5,10,5,5,5];

alpha = [0.5,0.5,0.5,0.5,0.5,0.5];
us = 60 * alpha./f;

%link cost = in-veh cost + waiting cost
linkcost = (ts + us)';

M = 20;
demand_mean_0 = [150,120,100,80,60,50]';
variance = 5*ones(num_od,5);
%------------------------------------------------------------------
%------------------------------------------------------------------
%------------------------------------------------------------------
%Main
% state the dimention of linkflow/pathflow
linkflow_SUE_table = zeros(num_link,k);
pathflow_SUE_table = zeros(num_link,k);
proportion_SUE_table = zeros(num_link*k,num_od*k);

link_ob_SUE = zeros(num_link*k,1);
%record all the randomly generatated demand
demand_table = zeros(num_od*k,1);

%obtain linkflow, pathflow, link proportion, path proportion of MNL and SUE
%for k day of obesertaion.

%multiple demand input with variance = sqrt mean value, total M demand sets

multi_demand_input = zeros(num_od,M);
multi_demand_table = zeros(num_od,M);


for j = 1:M
%assume normal demand
    demand_mean = demand_mean_0 + 1*j; %multiple demand increase by 1
    multi_demand_input(:,j) = demand_mean;

    for i = 1:k
        %one day observation of damand, n*1
        demand = zeros(6,1);
        demand(1,1) =  normrnd(demand_mean(1),variance(1));
        demand(2,1) =  normrnd(demand_mean(2),variance(2));
        demand(3,1) =  normrnd(demand_mean(3),variance(3));
        demand(4,1) =  normrnd(demand_mean(4),variance(4));
        demand(5,1) =  normrnd(demand_mean(5),variance(5));
        demand(6,1) =  normrnd(demand_mean(6),variance(6));

        demand_table(1+(i-1)*num_od:i*num_od,1) = demand;

        %rule out negative demand
        if all(demand<0)
            continue;
        end

        %obtain result from MNL and SUE function
        [linkflow_SUE, pathflow_SUE, proportion_path_SUE, proportion_link_SUE] = SUE(num_link,num_od,num_path, demand,linkcost,ka,delta0,gamma0);

        linkflow_SUE_table(1:num_link,i) = linkflow_SUE;

        pathflow_SUE_table(1:num_path,i) = pathflow_SUE;

        proportion_SUE_table(1+(i-1)*num_link:i*num_link,1+(i-1)*num_od:i*num_od) = proportion_link_SUE;

        %link observation from SUE
        %here add normal noise to link-observation
        error = normrnd(0,5,[num_link,1]);
        link_ob_SUE(1+(i-1)*num_link:i*num_link,1) = linkflow_SUE + error;

    end

    %input of estimation: exogeneous proportion matrix, link observation, k
    [x_solErr,demand_meanErr,demand_RCErr,LinkTrue] = estimation_optimErr(proportion_SUE_table, link_ob_SUE, num_od,num_link,k);
    %compare input and reconstruction output
    demand_result = [demand_table,demand_RCErr];
    %demand_resultErr = [demand_meanErr,demand_RCErr];

    %demand mean estimation output
    multi_demand_table(:,j) = demand_meanErr;

end

demand_comparision = [multi_demand_input;multi_demand_table];
csvwrite('DemandComparisonNormalErrNetwork_1.csv',demand_comparision);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%functions-----------------------------------------------------------------

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
function[xpv,demand_mean,demand_RC,LinkTrue]= estimation_optimErr(proportion, link_ob, num_od,num_link,k)
    %use fmincon to solve ||x-theta||^2_2 + ||z-v||^2_2

    %xp is [theta;x;link_ob]
    c0 = zeros(num_od*k,num_od);
    c1 = zeros(num_od*k,num_od*k);
    c2 = eye(num_od*k,num_link*k); 
    for i = 1:k
        c0(1+(i-1)*num_od:i*num_od,:) = eye(num_od);
        c1(1+(i-1)*num_od:i*num_od , 1+(i-1)*num_od:i*num_od) = -eye(num_od);
    end
    c01 = [c0,c1];
    c = zeros(num_link*k+num_link*k,num_od+num_od*k+num_link*k);
    c( 1:num_od*k , 1:num_od+num_od*k) = c01;
    c(num_od*k+1:num_od*k+num_link*k , num_od+num_od*k+1:num_od+num_od*k+num_link*k) = c2;

    d = [zeros(num_od*k,1);link_ob];


    % %inequality constraints, AX<=b, here does not exist
    A = [];
    b = [];
    
    Aeq = [zeros(num_link*k,num_od), proportion, -1*eye(num_link*k)];
    Beq = zeros(num_link*k,1);
    
    lb = zeros(num_od+num_od*k+num_link*k,1);
    ub = Inf*ones(num_od+num_od*k+num_link*k,1);
    
    xpv = lsqlin(c,d,A,b,Aeq,Beq,lb,ub);
    
    demand_mean = xpv(1:num_od);
    demand_RC = xpv(num_od+1:num_od+k*num_od);
    LinkTrue = xpv(num_od+k*num_od+1:num_od+k*num_od+k*num_link);
end