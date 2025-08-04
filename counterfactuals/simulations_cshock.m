%THIS NEEDS TO BE UPDATED!

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

%%

options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);

% Generate the ranks
pi_c = sum(pick,2);
[~, sortedIndices] = sort(pi_c, 'descend');
ranks = 1:length(pi_c);
ranks(sortedIndices) = ranks;
ranks = ranks';

rho_cf = [0.00, 0.25, 0.50, 0.75];
results = zeros(N,length(rho_cf));
shockvalue = 0.1;

for i = 1:length(rho_cf)
    disp "Rho"
    rho = rho_cf(i);

    % Calculate initial wage levels, given technology matrix and vector of
    % parameters.
    f = @(w)solver(w,T,alpha,theta,eta,rho);
    [w_init] = fsolve(f,w,options);

    % Calculate change in choice shares due to technology shock.
    lambda_init = sum(w_init.^(theta/(1-rho)));
    lambda_sum_init = sum(lambda_init.^(1-rho)); 
    pick_init = (w_init.^(theta/(1-rho))).*(lambda_init.^(-rho))/lambda_sum_init; 
    pic_init = sum(pick_init,2);

for c=1:N 

    if ranks(c) <= 100
    T_cf = T;
    T_cf(c,:) = T_cf(c,:)*(1-shockvalue); 

    % Calculate new technology matrix and solve for new wage matrix.
    fcf = @(w)solver(w,T_cf,alpha,theta,eta,rho);
    [w_cf] = fsolve(fcf,w_init,options);
    lambda_cf = sum(w_cf.^(theta/(1-rho)));
    lambda_sum_cf = sum(lambda_cf.^(1-rho)); 
    pick_cf = (w_cf.^(theta/(1-rho))).*(lambda_cf.^(-rho))/lambda_sum_cf; 

    % Calculate vector of city-level changes in choice shares due to technology
    % shock.
    pic_cf = sum(pick_cf,2);
    dpic = pic_cf ./ pic_init;
    results(c,i) = dpic(c);
    end
    if ranks(c) > 100 
        results(c,i) = 1;
    end

end

end


writematrix(results,"C:\Users\tsjac\Desktop\research\correlated_location\simulation_results_city_shock.csv")

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

% Generate the ranks
pi_c = sum(pick,2);
[~, sortedIndices] = sort(pi_c, 'descend');
ranks = 1:length(pi_c);
ranks(sortedIndices) = ranks;
ranks = ranks';

rho_cf = [0.00, 0.25, 0.50, 0.75];
results = zeros(N,length(rho_cf));
shockvalue = 0.1;

for i = 1:length(rho_cf)
    disp "Rho"
    
    rho = 0;
    theta = 2 / (1-rho_cf(i));

    % Calculate initial wage levels, given technology matrix and vector of
    % parameters.
    f = @(w)solver(w,T,alpha,theta,eta,rho);
    [w_init] = fsolve(f,w,options);

    % Calculate change in choice shares due to technology shock.
    lambda_init = sum(w_init.^(theta/(1-rho)));
    lambda_sum_init = sum(lambda_init.^(1-rho)); 
    pick_init = (w_init.^(theta/(1-rho))).*(lambda_init.^(-rho))/lambda_sum_init; 
    pic_init = sum(pick_init,2);

for c=1:N 

    if ranks(c) <= 100
    T_cf = T;
    T_cf(c,:) = T_cf(c,:)*(1-shockvalue); 

    % Calculate new technology matrix and solve for new wage matrix.
    fcf = @(w)solver(w,T_cf,alpha,theta,eta,rho);
    [w_cf] = fsolve(fcf,w_init,options);
    lambda_cf = sum(w_cf.^(theta/(1-rho)));
    lambda_sum_cf = sum(lambda_cf.^(1-rho)); 
    pick_cf = (w_cf.^(theta/(1-rho))).*(lambda_cf.^(-rho))/lambda_sum_cf; 

    % Calculate vector of city-level changes in choice shares due to technology
    % shock.
    pic_cf = sum(pick_cf,2);
    dpic = pic_cf ./ pic_init;
    results(c,i) = dpic(c);
    end
    if ranks(c) > 100 
        results(c,i) = 1;
    end

end

end


writematrix(results,"C:\Users\tsjac\Desktop\research\correlated_location\simulation_results_city_shock_falsetheta.csv")