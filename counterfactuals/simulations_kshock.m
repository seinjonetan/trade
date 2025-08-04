% THIS NEEDS TO BE UPDATED!

clear;

% Input necessary data from 1990 US Census.
%pick = table2array(readtable("C:\Users\tsjac\Desktop\research\correlated_location\matlab_input_pi_ckt.csv"));
%yck = ;

% Initiate key parameter inputs.
K = width(T); 
N = height(T);
theta = 2; 
eta = 1.85; 
alpha = 7; 

% Set optimization tolerance levels;
tolfun = 1e-10; 
tolx = 1e-10;

% Set up vector of rho values, results, and shock values.
rho_cf = [0.00, 0.25, 0.50, 0.75];
results = zeros(N,length(rho_cf));
shockvalue = 0.1;
iter = 1;

% Which occupation are we shocking
k = 1;

for i = 1:length(rho_cf)
    disp "Rho"
    rho = rho_cf(i);

    % Calculate initial wage levels, given technology matrix and vector of
    % parameters.
    [w_init] = equilibrium(w,T,alpha,theta,eta,rho,tolfun,tolx);

    % Calculate change in choice shares due to technology shock.
    lambda_init = sum(w_init.^(theta/(1-rho)));
    lambda_sum_init = sum(lambda_init.^(1-rho)); 
    pick_init = (w_init.^(theta/(1-rho))).*(lambda_init.^(-rho))/lambda_sum_init; 
    pic_init = sum(pick_init,2);

    T_cf = T;
    T_cf(:,k) = T_cf(:,k)*(1-shockvalue); 

    % Calculate new technology matrix and solve for new wage matrix.
    [w_cf] = equilibrium(w_init,T_cf,alpha,theta,eta,rho,tolfun,tolx);

    lambda_cf = sum(w_cf.^(theta/(1-rho)));
    lambda_sum_cf = sum(lambda_cf.^(1-rho)); 
    pick_cf = (w_cf.^(theta/(1-rho))).*(lambda_cf.^(-rho))/lambda_sum_cf; 

    % Calculate vector of city-level changes in choice shares due to technology
    % shock.
    
    pic_cf = sum(pick_cf,2);
    dpic = pic_cf ./ pic_init;
    results(:,iter) = dpic;
    iter = iter + 1;
end

writematrix(results,"C:\Users\tsjac\Desktop\research\correlated_location\simulation_results_occupation_shock.csv")

%%
clear;

% Input all data from 1990 US Census.
T = table2array(readtable("C:\Users\tsjac\Desktop\research\correlated_location\matlab_input_T_ckt.csv"));
pick = table2array(readtable("C:\Users\tsjac\Desktop\research\correlated_location\matlab_input_pi_ckt.csv"));
w = table2array(readtable("C:\Users\tsjac\Desktop\research\correlated_location\matlab_input_w_ckt.csv"));

% Initiate key parameter inputs.
K = width(T); 
N = height(T);
theta = 2; 
eta = 1.85; 
alpha = 7; 

options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);

% Generate ranks
pi_c = sum(pick,2);
[~, sortedIndices] = sort(pi_c, 'descend');
ranks = 1:length(pi_c);
ranks(sortedIndices) = ranks;
ranks = ranks';

rho_cf = [0.00, 0.25, 0.50, 0.75];
results_pi = zeros(K,length(rho_cf));
results_wb = zeros(K,length(rho_cf));
shockvalue = 0.1;

for i = 1:length(rho_cf)
    disp "Rho"
    
    rho = rho_cf(i);
    theta_cf = theta;
    
    % Turn this on if you want to run the "false theta" version.
    %theta_cf = theta / (1-rho);
    %rho = 0.00;
    
    % Calculate initial wage levels, given technology matrix and vector of
    % parameters.
    f = @(w)solver(w,T,alpha,theta_cf,eta,rho);
    [w_init] = fsolve(f,w,options);

    % Calculate change in choice shares due to technology shock.
    lambda_init = sum(w_init.^(theta_cf/(1-rho)));
    lambda_sum_init = sum(lambda_init.^(1-rho)); 
    pick_init = (w_init.^(theta_cf/(1-rho))).*(lambda_init.^(-rho))/lambda_sum_init; 
    pik_init = sum(pick_init)';
    zck_init = gamma((theta_cf-1)/theta_cf)*(pick_init.^((1+theta_cf)/theta_cf));
    wb_init = w_init.*zck_init;
    wbk_init = sum(wb_init)';

    for k = 1:K

    T_cf = T;
    T_cf(:,k) = T_cf(:,k)*(1-shockvalue); 

    % Calculate new technology matrix and solve for new wage matrix.
    fcf = @(w)solver(w,T_cf,alpha,theta_cf,eta,rho);
    [w_cf] = fsolve(fcf,w_init,options);
    lambda_cf = sum(w_cf.^(theta_cf/(1-rho)));
    lambda_sum_cf = sum(lambda_cf.^(1-rho)); 
    pick_cf = (w_cf.^(theta_cf/(1-rho))).*(lambda_cf.^(-rho))/lambda_sum_cf; 

    % Calculate vector of city-level changes in choice shares due to technology
    % shock.
    pik_cf = sum(pick_cf)';
    dpik = pik_cf(k) ./ pik_init(k);
    results_pi(k,i) = dpik;
    
    zck_cf = gamma((theta_cf-1)/theta_cf)*(pick_cf.^((1+theta_cf)/theta_cf));
    wb_cf = w_cf.*zck_cf;
    wbk_cf = sum(wb_cf)';
    results_wb(k,i) = wbk_cf(k) / wbk_init(k);
    
    end
end

%%
writematrix(results_pi,"C:\Users\tsjac\Desktop\research\correlated_location\simulation_results_occupation_shock_pi.csv")
writematrix(results_wb,"C:\Users\tsjac\Desktop\research\correlated_location\simulation_results_occupation_shock_wagebill.csv")